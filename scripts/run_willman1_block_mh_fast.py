#!/usr/bin/env python
"""Fit Willman 1 with block Metropolis-Hastings and make diagnostic plots.

This testing script uses:
- experimental R-z Jeans moment grid + cubic-spline dP_z/dR,
- rho0 linear scaling, while still sampling log10_rho0,
- explicit sampling of systemic_velocity_kms.

The default input is Willman 1. Alternate galaxy/member and center tables can be
passed with --galaxy-csv and --centers-csv.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / "outputs" / ".matplotlib"))
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
try:
    from tqdm.auto import tqdm
except ImportError:  # pragma: no cover - optional progress dependency
    tqdm = None

from experiments.test_halo_force_grid_interpolation import RZMomentGridWithForceGrid
from experiments.test_rz_moment_grid_interpolation import RZMomentGridProjector
from hayashi_jeans.data import GalaxyData, load_galaxy_data
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.mge import MGEPhysicalityConfig, mge_sigma_los2_unit, mge_sigma_los2_unit_checked
from hayashi_jeans.model import build_hernquist_projector
from hayashi_jeans.params import HayashiParameters


@dataclass(frozen=True)
class SlowParams:
    q_halo: float
    log10_b_halo_pc: float
    minus_log10_one_minus_beta_z: float
    alpha: float
    beta: float
    gamma: float
    inclination_deg: float


@dataclass(frozen=True)
class ChainState:
    slow: SlowParams
    log10_rho0_msun_pc3: float
    systemic_velocity_kms: float
    sigma_unit: np.ndarray
    log_probability: float


LIKELIHOOD_MODES = ("fast", "mge", "validation-halo", "validation-strict", "validation-convergence")
_EMCEE_LOG_PROB_CONTEXT: dict[str, object] | None = None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("time", "fit", "validate", "plot", "corner", "all"), default="all")
    parser.add_argument("--sampler", choices=("block-mh", "emcee"), default="block-mh")
    parser.add_argument("--likelihood-mode", choices=LIKELIHOOD_MODES, default="fast")
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/27_Willman_1.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--chain-output", default=str(PROJECT_ROOT / "outputs/willman1_block_mh_fast_chain.csv"))
    parser.add_argument("--validation-output", default=str(PROJECT_ROOT / "outputs/diagnostics/willman1_validation_logprob.csv"))
    parser.add_argument("--figure-output", default=str(PROJECT_ROOT / "outputs/figures/willman1_block_mh_fast_density_profile.png"))
    parser.add_argument("--corner-output", default=str(PROJECT_ROOT / "outputs/figures/willman1_block_mh_fast_corner.png"))
    parser.add_argument("--burn", type=int, default=0, help="Number of initial chain samples to discard before plotting.")
    parser.add_argument("--corner-burn-in", type=int, default=None, help="Deprecated alias for --burn.")
    parser.add_argument("--show-systemic-in-corner", action="store_true")
    parser.add_argument("--n-steps", type=int, default=200)
    parser.add_argument("--n-walkers", type=int, default=24)
    parser.add_argument("--n-chains", type=int, default=1)
    parser.add_argument("--n-processes", type=int, default=1)
    parser.add_argument("--n-fast-updates", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260520)
    parser.add_argument("--n-r", type=int, default=96)
    parser.add_argument("--n-z", type=int, default=192)
    parser.add_argument("--n-los", type=int, default=160)
    parser.add_argument("--n-force-r", type=int, default=48)
    parser.add_argument("--n-force-z", type=int, default=96)
    parser.add_argument("--validation-samples", type=int, default=5)
    parser.add_argument("--strict-epsrel", type=float, default=1.5e-2)
    parser.add_argument("--proposal-scale-slow", default="0.04,0.08,0.05,0.12,0.18,0.06,1.5")
    parser.add_argument("--proposal-scale-log-rho", type=float, default=0.08)
    parser.add_argument("--proposal-scale-systemic", type=float, default=0.7)
    parser.add_argument("--systemic-prior", choices=("adaptive", "legacy", "manual"), default="adaptive")
    parser.add_argument("--systemic-prior-padding", type=float, default=20.0)
    parser.add_argument("--systemic-prior-min-half-width", type=float, default=20.0)
    parser.add_argument(
        "--systemic-prior-bounds",
        default=None,
        help="Manual systemic velocity prior bounds as 'low,high' in km/s; requires --systemic-prior manual.",
    )
    parser.add_argument(
        "--no-mge-physicality-check",
        action="store_true",
        help="Disable the raw MGE/JAM local v_phi^2 physicality guard.",
    )
    parser.add_argument("--no-progress", action="store_true", help="Disable the sampling progress bar.")
    args = parser.parse_args()
    if args.n_chains < 1:
        raise ValueError("--n-chains must be >= 1")
    if args.n_processes < 1:
        raise ValueError("--n-processes must be >= 1")
    burn = args.burn if args.corner_burn_in is None else args.corner_burn_in
    if burn < 0:
        raise ValueError("--burn must be non-negative")

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    systemic_bounds = resolve_systemic_velocity_prior_bounds(
        galaxy,
        mode=args.systemic_prior,
        manual_bounds=args.systemic_prior_bounds,
        padding=args.systemic_prior_padding,
        min_half_width=args.systemic_prior_min_half_width,
    )
    print(f"systemic_velocity_prior=({systemic_bounds[0]:.6g}, {systemic_bounds[1]:.6g})")

    if args.mode in ("time", "all"):
        dt, logp = time_single_likelihood(
            galaxy,
            initial_slow_params(),
            log10_rho0_msun_pc3=-1.961219,
            systemic_velocity_kms=-13.298645,
            likelihood_mode=args.likelihood_mode,
            n_r=args.n_r,
            n_z=args.n_z,
            n_los=args.n_los,
            n_force_r=args.n_force_r,
            n_force_z=args.n_force_z,
            strict_epsrel=args.strict_epsrel,
            systemic_velocity_bounds=systemic_bounds,
            use_mge_physicality_check=not args.no_mge_physicality_check,
        )
        print(f"single_likelihood_seconds={dt:.6f}")
        print(f"single_likelihood_log_probability={logp:.6f}")
        print(f"likelihood_mode={args.likelihood_mode}")

    if args.mode in ("fit", "all"):
        if args.sampler == "block-mh":
            chain = run_parallel_block_mh(
                galaxy,
                galaxy_csv=args.galaxy_csv,
                centers_csv=args.centers_csv,
                n_steps=args.n_steps,
                n_chains=args.n_chains,
                n_processes=args.n_processes,
                n_fast_updates=args.n_fast_updates,
                seed=args.seed,
                likelihood_mode=args.likelihood_mode,
                n_r=args.n_r,
                n_z=args.n_z,
                n_los=args.n_los,
                n_force_r=args.n_force_r,
                n_force_z=args.n_force_z,
                strict_epsrel=args.strict_epsrel,
                proposal_scale_slow=np.array([float(x) for x in args.proposal_scale_slow.split(",")]),
                proposal_scale_log_rho=args.proposal_scale_log_rho,
                proposal_scale_systemic=args.proposal_scale_systemic,
                systemic_velocity_bounds=systemic_bounds,
                use_mge_physicality_check=not args.no_mge_physicality_check,
                show_progress=not args.no_progress,
            )
        else:
            chain = run_emcee_unlayered(
                galaxy,
                n_steps=args.n_steps,
                n_walkers=args.n_walkers,
                n_processes=args.n_processes,
                seed=args.seed,
                likelihood_mode=args.likelihood_mode,
                n_r=args.n_r,
                n_z=args.n_z,
                n_los=args.n_los,
                n_force_r=args.n_force_r,
                n_force_z=args.n_force_z,
                strict_epsrel=args.strict_epsrel,
                systemic_velocity_bounds=systemic_bounds,
                use_mge_physicality_check=not args.no_mge_physicality_check,
                show_progress=not args.no_progress,
            )
        chain_output = Path(args.chain_output)
        chain_output.parent.mkdir(parents=True, exist_ok=True)
        chain.to_csv(chain_output, index=False)
        print(f"wrote {chain_output}")
        print_summary(chain)

    if args.mode in ("validate",):
        chain_path = Path(args.chain_output)
        if not chain_path.exists():
            raise FileNotFoundError(f"{chain_path} does not exist; run --mode fit first")
        chain = pd.read_csv(chain_path)
        validation = validate_chain_log_probability(
            galaxy,
            chain,
            likelihood_mode=args.likelihood_mode,
            n_samples=args.validation_samples,
            seed=args.seed,
            n_r=args.n_r,
            n_z=args.n_z,
            n_los=args.n_los,
            n_force_r=args.n_force_r,
            n_force_z=args.n_force_z,
            strict_epsrel=args.strict_epsrel,
            systemic_velocity_bounds=systemic_bounds,
            use_mge_physicality_check=not args.no_mge_physicality_check,
            show_progress=not args.no_progress,
        )
        validation_output = Path(args.validation_output)
        validation_output.parent.mkdir(parents=True, exist_ok=True)
        validation.to_csv(validation_output, index=False)
        print(f"wrote {validation_output}")
        print(validation.to_string(index=False))

    if args.mode in ("plot", "all"):
        chain_path = Path(args.chain_output)
        if not chain_path.exists():
            raise FileNotFoundError(f"{chain_path} does not exist; run --mode fit first")
        chain = pd.read_csv(chain_path)
        plot_density_profile(galaxy, chain, Path(args.figure_output), burn=burn)

    if args.mode in ("corner", "all"):
        chain_path = Path(args.chain_output)
        if not chain_path.exists():
            raise FileNotFoundError(f"{chain_path} does not exist; run --mode fit first")
        chain = pd.read_csv(chain_path)
        generate_corner_map(
            chain,
            Path(args.corner_output),
            show_systemic_velocity=args.show_systemic_in_corner,
            burn=burn,
        )


def time_single_likelihood(
    galaxy: GalaxyData,
    slow: SlowParams,
    *,
    log10_rho0_msun_pc3: float,
    systemic_velocity_kms: float,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    systemic_velocity_bounds: tuple[float, float],
    use_mge_physicality_check: bool,
) -> tuple[float, float]:
    t0 = time.perf_counter()
    logp = log_prior_slow(slow, galaxy) + log_prior_fast(
        log10_rho0_msun_pc3,
        systemic_velocity_kms,
        systemic_velocity_bounds=systemic_velocity_bounds,
    )
    sigma_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode=likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
        strict_epsrel=strict_epsrel,
        use_mge_physicality_check=use_mge_physicality_check,
    )
    if np.isfinite(logp):
        logp += log_likelihood_from_unit_sigma(
            galaxy,
            sigma_unit,
            log10_rho0_msun_pc3=log10_rho0_msun_pc3,
            systemic_velocity_kms=systemic_velocity_kms,
        )
    return time.perf_counter() - t0, float(logp)


def validation_grid_settings(
    likelihood_mode: str,
    *,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
) -> dict[str, int | None]:
    if likelihood_mode == "fast":
        return {"n_r": n_r, "n_z": n_z, "n_los": n_los, "n_force_r": n_force_r, "n_force_z": n_force_z}
    if likelihood_mode == "validation-halo":
        return {"n_r": 128, "n_z": 256, "n_los": 200, "n_force_r": None, "n_force_z": None}
    if likelihood_mode == "validation-convergence":
        return {"n_r": 128, "n_z": 256, "n_los": 200, "n_force_r": 64, "n_force_z": 128}
    if likelihood_mode == "validation-strict":
        return {"n_r": None, "n_z": None, "n_los": None, "n_force_r": None, "n_force_z": None}
    raise ValueError(f"unknown likelihood_mode={likelihood_mode!r}")


def run_parallel_block_mh(
    galaxy: GalaxyData,
    *,
    galaxy_csv: str,
    centers_csv: str,
    n_steps: int,
    n_chains: int,
    n_processes: int,
    n_fast_updates: int,
    seed: int,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    proposal_scale_slow: np.ndarray,
    proposal_scale_log_rho: float,
    proposal_scale_systemic: float,
    systemic_velocity_bounds: tuple[float, float],
    use_mge_physicality_check: bool,
    show_progress: bool,
) -> pd.DataFrame:
    if n_chains == 1:
        return run_block_mh(
            galaxy,
            n_steps=n_steps,
            n_fast_updates=n_fast_updates,
            seed=seed,
            likelihood_mode=likelihood_mode,
            n_r=n_r,
            n_z=n_z,
            n_los=n_los,
            n_force_r=n_force_r,
            n_force_z=n_force_z,
            strict_epsrel=strict_epsrel,
            proposal_scale_slow=proposal_scale_slow,
            proposal_scale_log_rho=proposal_scale_log_rho,
            proposal_scale_systemic=proposal_scale_systemic,
            systemic_velocity_bounds=systemic_velocity_bounds,
            use_mge_physicality_check=use_mge_physicality_check,
            show_progress=show_progress,
            print_progress=True,
            chain_id=0,
            progress_desc="MH sampling",
        )

    payloads = [
        {
            "chain_id": chain_id,
            "galaxy_csv": galaxy_csv,
            "centers_csv": centers_csv,
            "n_steps": n_steps,
            "n_fast_updates": n_fast_updates,
            "seed": seed + 1_000_003 * chain_id,
            "likelihood_mode": likelihood_mode,
            "n_r": n_r,
            "n_z": n_z,
            "n_los": n_los,
            "n_force_r": n_force_r,
            "n_force_z": n_force_z,
            "strict_epsrel": strict_epsrel,
            "proposal_scale_slow": proposal_scale_slow.tolist(),
            "proposal_scale_log_rho": proposal_scale_log_rho,
            "proposal_scale_systemic": proposal_scale_systemic,
            "systemic_velocity_bounds": list(systemic_velocity_bounds),
            "use_mge_physicality_check": use_mge_physicality_check,
        }
        for chain_id in range(n_chains)
    ]

    if n_processes == 1:
        chains = []
        iterator = payloads
        if show_progress and tqdm:
            iterator = tqdm(payloads, desc="chains", unit="chain", dynamic_ncols=True)
        for payload in iterator:
            chains.append(_run_block_mh_worker(payload))
    else:
        chains = []
        max_workers = min(n_processes, n_chains)
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_run_block_mh_worker, payload) for payload in payloads]
            iterator = as_completed(futures)
            if show_progress and tqdm:
                iterator = tqdm(iterator, total=n_chains, desc="chains", unit="chain", dynamic_ncols=True)
            for future in iterator:
                chains.append(future.result())

    return pd.concat(chains, ignore_index=True).sort_values(["chain_id", "step"]).reset_index(drop=True)


def _run_block_mh_worker(payload: dict[str, object]) -> pd.DataFrame:
    galaxy = load_galaxy_data(
        str(payload["galaxy_csv"]),
        structural_centers_csv=str(payload["centers_csv"]),
    )
    return run_block_mh(
        galaxy,
        n_steps=int(payload["n_steps"]),
        n_fast_updates=int(payload["n_fast_updates"]),
        seed=int(payload["seed"]),
        likelihood_mode=str(payload["likelihood_mode"]),
        n_r=int(payload["n_r"]),
        n_z=int(payload["n_z"]),
        n_los=int(payload["n_los"]),
        n_force_r=int(payload["n_force_r"]),
        n_force_z=int(payload["n_force_z"]),
        strict_epsrel=float(payload["strict_epsrel"]),
        proposal_scale_slow=np.array(payload["proposal_scale_slow"], dtype=float),
        proposal_scale_log_rho=float(payload["proposal_scale_log_rho"]),
        proposal_scale_systemic=float(payload["proposal_scale_systemic"]),
        systemic_velocity_bounds=tuple(float(x) for x in payload["systemic_velocity_bounds"]),
        use_mge_physicality_check=bool(payload["use_mge_physicality_check"]),
        show_progress=False,
        print_progress=False,
        chain_id=int(payload["chain_id"]),
        progress_desc=f"chain {payload['chain_id']}",
    )


def run_block_mh(
    galaxy: GalaxyData,
    *,
    n_steps: int,
    n_fast_updates: int,
    seed: int,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    proposal_scale_slow: np.ndarray,
    proposal_scale_log_rho: float,
    proposal_scale_systemic: float,
    systemic_velocity_bounds: tuple[float, float],
    use_mge_physicality_check: bool,
    show_progress: bool,
    print_progress: bool = True,
    chain_id: int = 0,
    progress_desc: str = "MH sampling",
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    slow = initial_slow_params()
    log_rho = -1.961219
    systemic = float(np.clip(np.median(galaxy.velocity_kms), systemic_velocity_bounds[0], systemic_velocity_bounds[1]))
    sigma_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode=likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
        strict_epsrel=strict_epsrel,
        use_mge_physicality_check=use_mge_physicality_check,
    )
    logp = full_log_probability(
        galaxy,
        slow,
        log_rho,
        systemic,
        sigma_unit,
        systemic_velocity_bounds=systemic_velocity_bounds,
    )
    state = ChainState(slow, log_rho, systemic, sigma_unit, logp)

    rows = []
    accept_slow = 0
    accept_rho = 0
    accept_systemic = 0
    reject_mge_slow = 0
    total_rho = 0
    total_systemic = 0
    progress = tqdm(total=n_steps, desc=progress_desc, unit="step", dynamic_ncols=True) if show_progress and tqdm else None

    for step in range(1, n_steps + 1):
        for _ in range(n_fast_updates):
            total_rho += 1
            proposal_log_rho = state.log10_rho0_msun_pc3 + rng.normal(0.0, proposal_scale_log_rho)
            proposal_logp = full_log_probability(
                galaxy,
                state.slow,
                proposal_log_rho,
                state.systemic_velocity_kms,
                state.sigma_unit,
                systemic_velocity_bounds=systemic_velocity_bounds,
            )
            if accept_move(rng, proposal_logp - state.log_probability):
                state = ChainState(state.slow, proposal_log_rho, state.systemic_velocity_kms, state.sigma_unit, proposal_logp)
                accept_rho += 1

            total_systemic += 1
            proposal_systemic = state.systemic_velocity_kms + rng.normal(0.0, proposal_scale_systemic)
            proposal_logp = full_log_probability(
                galaxy,
                state.slow,
                state.log10_rho0_msun_pc3,
                proposal_systemic,
                state.sigma_unit,
                systemic_velocity_bounds=systemic_velocity_bounds,
            )
            if accept_move(rng, proposal_logp - state.log_probability):
                state = ChainState(state.slow, state.log10_rho0_msun_pc3, proposal_systemic, state.sigma_unit, proposal_logp)
                accept_systemic += 1

        proposal_slow = propose_slow(rng, state.slow, proposal_scale_slow)
        proposal_prior = log_prior_slow(proposal_slow, galaxy)
        if np.isfinite(proposal_prior):
            proposal_sigma_unit = compute_sigma_unit(
                galaxy,
                proposal_slow,
                likelihood_mode=likelihood_mode,
                n_r=n_r,
                n_z=n_z,
                n_los=n_los,
                n_force_r=n_force_r,
                n_force_z=n_force_z,
                strict_epsrel=strict_epsrel,
                use_mge_physicality_check=use_mge_physicality_check,
            )
            proposal_logp = full_log_probability(
                galaxy,
                proposal_slow,
                state.log10_rho0_msun_pc3,
                state.systemic_velocity_kms,
                proposal_sigma_unit,
                systemic_velocity_bounds=systemic_velocity_bounds,
            )
            if likelihood_mode == "mge" and use_mge_physicality_check and not np.isfinite(proposal_logp):
                reject_mge_slow += 1
            if accept_move(rng, proposal_logp - state.log_probability):
                state = ChainState(proposal_slow, state.log10_rho0_msun_pc3, state.systemic_velocity_kms, proposal_sigma_unit, proposal_logp)
                accept_slow += 1

        rows.append(
            {
                "chain_id": chain_id,
                "step": step,
                "q_halo": state.slow.q_halo,
                "log10_b_halo_pc": state.slow.log10_b_halo_pc,
                "log10_rho0_msun_pc3": state.log10_rho0_msun_pc3,
                "minus_log10_one_minus_beta_z": state.slow.minus_log10_one_minus_beta_z,
                "beta_z": beta_z_from_q(state.slow.minus_log10_one_minus_beta_z),
                "alpha": state.slow.alpha,
                "beta": state.slow.beta,
                "gamma": state.slow.gamma,
                "i_deg": state.slow.inclination_deg,
                "systemic_velocity_kms": state.systemic_velocity_kms,
                "log_probability": state.log_probability,
                "likelihood_mode": likelihood_mode,
                "sampler": "block-mh",
                "accept_slow_rate": accept_slow / step,
                "accept_rho_rate": accept_rho / max(total_rho, 1),
                "accept_systemic_rate": accept_systemic / max(total_systemic, 1),
                "mge_physicality_check": bool(use_mge_physicality_check),
                "mge_reject_slow": reject_mge_slow,
                "mge_reject_slow_rate": reject_mge_slow / step,
            }
        )
        progress_text = (
            f"logp={state.log_probability:.3f}, "
            f"logrho={state.log10_rho0_msun_pc3:.3f}, "
            f"u={state.systemic_velocity_kms:.3f}, "
            f"acc_slow={accept_slow/step:.2f}"
        )
        if progress is not None:
            progress.set_postfix_str(progress_text)
            progress.update(1)
        elif print_progress and (step == 1 or step % 10 == 0):
            print(
                f"step={step}/{n_steps} {progress_text}",
                flush=True,
            )

    if progress is not None:
        progress.close()

    return pd.DataFrame(rows)


def compute_sigma_unit(
    galaxy: GalaxyData,
    slow: SlowParams,
    *,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    use_mge_physicality_check: bool = True,
) -> np.ndarray:
    params = HayashiParameters(
        q_halo=slow.q_halo,
        b_halo_pc=10.0**slow.log10_b_halo_pc,
        rho0_msun_pc3=1.0,
        beta_z=beta_z_from_q(slow.minus_log10_one_minus_beta_z),
        alpha=slow.alpha,
        beta=slow.beta,
        gamma=slow.gamma,
        inclination_rad=np.deg2rad(slow.inclination_deg),
        systemic_velocity_kms=0.0,
    )
    if likelihood_mode == "mge":
        mge_func = mge_sigma_los2_unit_checked if use_mge_physicality_check else mge_sigma_los2_unit
        return mge_func(
            galaxy,
            q_halo=params.q_halo,
            b_halo_pc=params.b_halo_pc,
            alpha=params.alpha,
            beta=params.beta,
            gamma=params.gamma,
            beta_z=params.beta_z,
            inclination_rad=params.inclination_rad,
            config=MGEPhysicalityConfig(),
        )
    settings = validation_grid_settings(
        likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
    )
    if likelihood_mode == "validation-strict":
        projector = build_hernquist_projector(
            galaxy,
            params,
            zmax_factor=20.0,
            los_factor=20.0,
            epsrel=strict_epsrel,
        )
        return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    if likelihood_mode == "validation-halo":
        projector = RZMomentGridProjector(
            galaxy=galaxy,
            params=params,
            n_r=int(settings["n_r"]),
            n_z=int(settings["n_z"]),
            n_los=int(settings["n_los"]),
            zmax_factor=20.0,
            los_factor=20.0,
            r_min_pc=1e-3,
            z_min_pc=1e-3,
            grid_padding=1.08,
            radial_derivative="cubic_spline",
        )
        return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)

    projector = RZMomentGridWithForceGrid(
        galaxy=galaxy,
        params=params,
        n_r=int(settings["n_r"]),
        n_z=int(settings["n_z"]),
        n_los=int(settings["n_los"]),
        n_force_r=int(settings["n_force_r"]),
        n_force_z=int(settings["n_force_z"]),
        zmax_factor=20.0,
        los_factor=20.0,
        r_min_pc=1e-3,
        z_min_pc=1e-3,
        grid_padding=1.08,
    )
    return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)


def log_likelihood_from_unit_sigma(
    galaxy: GalaxyData,
    sigma_unit: np.ndarray,
    *,
    log10_rho0_msun_pc3: float,
    systemic_velocity_kms: float,
) -> float:
    sigma_los2 = sigma_unit * 10.0**log10_rho0_msun_pc3
    return GaussianVelocityLikelihood(galaxy).log_likelihood(sigma_los2, systemic_velocity_kms=systemic_velocity_kms)


def full_log_probability(
    galaxy: GalaxyData,
    slow: SlowParams,
    log10_rho0_msun_pc3: float,
    systemic_velocity_kms: float,
    sigma_unit: np.ndarray,
    *,
    systemic_velocity_bounds: tuple[float, float] = (-80.0, 50.0),
) -> float:
    prior = log_prior_slow(slow, galaxy) + log_prior_fast(
        log10_rho0_msun_pc3,
        systemic_velocity_kms,
        systemic_velocity_bounds=systemic_velocity_bounds,
    )
    if not np.isfinite(prior):
        return -np.inf
    return prior + log_likelihood_from_unit_sigma(
        galaxy,
        sigma_unit,
        log10_rho0_msun_pc3=log10_rho0_msun_pc3,
        systemic_velocity_kms=systemic_velocity_kms,
    )


def initial_slow_params() -> SlowParams:
    return SlowParams(
        q_halo=1.0,
        log10_b_halo_pc=3.2,
        minus_log10_one_minus_beta_z=0.3,
        alpha=1.8,
        beta=6.4,
        gamma=1.2,
        inclination_deg=75.0,
    )


def propose_slow(rng: np.random.Generator, slow: SlowParams, scale: np.ndarray) -> SlowParams:
    values = np.array([
        slow.q_halo,
        slow.log10_b_halo_pc,
        slow.minus_log10_one_minus_beta_z,
        slow.alpha,
        slow.beta,
        slow.gamma,
        slow.inclination_deg,
    ])
    proposal = values + rng.normal(0.0, scale)
    return SlowParams(*proposal)


def log_prior_slow(slow: SlowParams, galaxy: GalaxyData) -> float:
    min_i = np.degrees(np.arccos(galaxy.observables.qprime))
    if not (0.1 <= slow.q_halo <= 2.0):
        return -np.inf
    if not (0.0 <= slow.log10_b_halo_pc <= 5.0):
        return -np.inf
    # Hayashi et al. (2023) prior: -1 <= -log10(1 - beta_z) < 1.
    if not (-1.0 <= slow.minus_log10_one_minus_beta_z < 1.0):
        return -np.inf
    if not (0.5 <= slow.alpha <= 3.0):
        return -np.inf
    if not (3.0 <= slow.beta <= 10.0):
        return -np.inf
    if not (0.0 <= slow.gamma <= 2.0):
        return -np.inf
    if not (min_i < slow.inclination_deg <= 90.0):
        return -np.inf
    return 0.0


def beta_z_from_q(minus_log10_one_minus_beta_z: float) -> float:
    return 1.0 - 10.0 ** (-float(minus_log10_one_minus_beta_z))


def q_from_beta_z(beta_z: float) -> float:
    return -float(np.log10(1.0 - float(beta_z)))


def log_prior_fast(
    log10_rho0_msun_pc3: float,
    systemic_velocity_kms: float,
    *,
    systemic_velocity_bounds: tuple[float, float] = (-80.0, 50.0),
) -> float:
    if not (-5.0 <= log10_rho0_msun_pc3 <= 5.0):
        return -np.inf
    low, high = systemic_velocity_bounds
    if not (low <= systemic_velocity_kms <= high):
        return -np.inf
    return 0.0


def systemic_velocity_prior_bounds(
    galaxy: GalaxyData,
    *,
    padding: float = 20.0,
    min_half_width: float = 20.0,
) -> tuple[float, float]:
    velocity = np.asarray(galaxy.velocity_kms, dtype=float)
    velocity = velocity[np.isfinite(velocity)]
    if velocity.size == 0:
        return (-80.0, 50.0)
    center = float(np.median(velocity))
    mad = float(np.median(np.abs(velocity - center)))
    robust_sigma = 1.4826 * mad
    minmax_half_width = 0.5 * float(np.max(velocity) - np.min(velocity))
    half_width = max(5.0 * robust_sigma, minmax_half_width + padding, min_half_width)
    return center - half_width, center + half_width


def parse_systemic_velocity_bounds(value: str) -> tuple[float, float]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) != 2:
        raise ValueError("--systemic-prior-bounds must have the form 'low,high'")
    low, high = (float(parts[0]), float(parts[1]))
    if not np.isfinite(low) or not np.isfinite(high) or not low < high:
        raise ValueError("--systemic-prior-bounds must be finite and satisfy low < high")
    return low, high


def resolve_systemic_velocity_prior_bounds(
    galaxy: GalaxyData,
    *,
    mode: str,
    manual_bounds: str | None,
    padding: float,
    min_half_width: float,
) -> tuple[float, float]:
    if mode == "legacy":
        return (-80.0, 50.0)
    if mode == "manual":
        if manual_bounds is None:
            raise ValueError("--systemic-prior manual requires --systemic-prior-bounds")
        return parse_systemic_velocity_bounds(manual_bounds)
    if mode == "adaptive":
        return systemic_velocity_prior_bounds(galaxy, padding=padding, min_half_width=min_half_width)
    raise ValueError(f"unknown systemic prior mode: {mode!r}")


def accept_move(rng: np.random.Generator, delta_logp: float) -> bool:
    if not np.isfinite(delta_logp):
        return False
    return delta_logp >= 0.0 or np.log(rng.random()) < delta_logp


def run_emcee_unlayered(
    galaxy: GalaxyData,
    *,
    n_steps: int,
    n_walkers: int,
    n_processes: int,
    seed: int,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    systemic_velocity_bounds: tuple[float, float],
    use_mge_physicality_check: bool,
    show_progress: bool,
) -> pd.DataFrame:
    try:
        import emcee
    except ImportError as exc:  # pragma: no cover - optional sampler dependency
        raise ImportError("emcee is required for --sampler emcee; install emcee or use --sampler block-mh") from exc

    ndim = 9
    if n_walkers < 2 * ndim:
        raise ValueError(f"--n-walkers must be at least {2 * ndim} for emcee's default ensemble move")

    rng = np.random.default_rng(seed)
    initial = initial_vector(galaxy)
    scale = np.array([0.03, 0.06, 0.06, 0.04, 0.10, 0.14, 0.05, 1.0, 0.5])
    positions = []
    while len(positions) < n_walkers:
        proposal = initial + rng.normal(0.0, scale)
        if np.isfinite(log_prior_vector(galaxy, proposal, systemic_velocity_bounds=systemic_velocity_bounds)):
            positions.append(proposal)
    positions = np.asarray(positions)

    context = make_emcee_log_prob_context(
        galaxy,
        likelihood_mode=likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
        strict_epsrel=strict_epsrel,
        systemic_velocity_bounds=systemic_velocity_bounds,
        use_mge_physicality_check=use_mge_physicality_check,
    )

    if n_processes == 1:
        sampler = emcee.EnsembleSampler(n_walkers, ndim, lambda vector: emcee_log_prob_from_context(vector, context))
        sampler.run_mcmc(positions, n_steps, progress=show_progress)
    else:
        with Pool(
            processes=n_processes,
            initializer=init_emcee_log_prob_worker,
            initargs=(context,),
        ) as pool:
            sampler = emcee.EnsembleSampler(n_walkers, ndim, emcee_log_prob_worker, pool=pool)
            sampler.run_mcmc(positions, n_steps, progress=show_progress)
    chain = sampler.get_chain(flat=False)
    log_prob_chain = sampler.get_log_prob(flat=False)

    rows = []
    for step_index in range(n_steps):
        for walker in range(n_walkers):
            rows.append(vector_to_row(
                chain[step_index, walker, :],
                log_probability=float(log_prob_chain[step_index, walker]),
                step=step_index + 1,
                chain_id=walker,
                sampler_name="emcee",
                likelihood_mode=likelihood_mode,
            ))
    return pd.DataFrame(rows)


def make_emcee_log_prob_context(
    galaxy: GalaxyData,
    *,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    systemic_velocity_bounds: tuple[float, float],
    use_mge_physicality_check: bool,
) -> dict[str, object]:
    return {
        "galaxy": galaxy,
        "likelihood_mode": likelihood_mode,
        "n_r": n_r,
        "n_z": n_z,
        "n_los": n_los,
        "n_force_r": n_force_r,
        "n_force_z": n_force_z,
        "strict_epsrel": strict_epsrel,
        "systemic_velocity_bounds": tuple(systemic_velocity_bounds),
        "use_mge_physicality_check": bool(use_mge_physicality_check),
    }


def init_emcee_log_prob_worker(context: dict[str, object]) -> None:
    global _EMCEE_LOG_PROB_CONTEXT
    _EMCEE_LOG_PROB_CONTEXT = context


def emcee_log_prob_worker(vector: np.ndarray) -> float:
    if _EMCEE_LOG_PROB_CONTEXT is None:
        raise RuntimeError("emcee worker context has not been initialized")
    return emcee_log_prob_from_context(vector, _EMCEE_LOG_PROB_CONTEXT)


def emcee_log_prob_from_context(vector: np.ndarray, context: dict[str, object]) -> float:
    return full_log_probability_vector(
        context["galaxy"],
        vector,
        likelihood_mode=str(context["likelihood_mode"]),
        n_r=int(context["n_r"]),
        n_z=int(context["n_z"]),
        n_los=int(context["n_los"]),
        n_force_r=int(context["n_force_r"]),
        n_force_z=int(context["n_force_z"]),
        strict_epsrel=float(context["strict_epsrel"]),
        systemic_velocity_bounds=tuple(float(x) for x in context["systemic_velocity_bounds"]),
        use_mge_physicality_check=bool(context["use_mge_physicality_check"]),
    )


def validate_chain_log_probability(
    galaxy: GalaxyData,
    chain: pd.DataFrame,
    *,
    likelihood_mode: str,
    n_samples: int,
    seed: int,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    systemic_velocity_bounds: tuple[float, float],
    use_mge_physicality_check: bool,
    show_progress: bool,
) -> pd.DataFrame:
    if n_samples < 1:
        raise ValueError("--validation-samples must be >= 1")
    required = [
        "q_halo",
        "log10_b_halo_pc",
        "log10_rho0_msun_pc3",
        "minus_log10_one_minus_beta_z",
        "alpha",
        "beta",
        "gamma",
        "i_deg",
        "systemic_velocity_kms",
    ]
    chain = chain_with_anisotropy_columns(chain)
    missing = [column for column in required if column not in chain.columns]
    if missing:
        raise ValueError(f"chain is missing columns required for validation: {missing}")

    finite_chain = chain.replace([np.inf, -np.inf], np.nan).dropna(subset=required)
    if finite_chain.empty:
        raise ValueError("chain has no finite rows to validate")
    rng = np.random.default_rng(seed)
    n_take = min(n_samples, len(finite_chain))
    row_indices = np.sort(rng.choice(finite_chain.index.to_numpy(), size=n_take, replace=False))
    iterator = row_indices
    if show_progress and tqdm:
        iterator = tqdm(row_indices, desc=f"validate {likelihood_mode}", unit="sample", dynamic_ncols=True)

    rows = []
    for row_index in iterator:
        row = finite_chain.loc[row_index]
        vector = row_to_vector(row)
        t0 = time.perf_counter()
        validated_logp = full_log_probability_vector(
            galaxy,
            vector,
            likelihood_mode=likelihood_mode,
            n_r=n_r,
            n_z=n_z,
            n_los=n_los,
            n_force_r=n_force_r,
            n_force_z=n_force_z,
            strict_epsrel=strict_epsrel,
            systemic_velocity_bounds=systemic_velocity_bounds,
            use_mge_physicality_check=use_mge_physicality_check,
        )
        elapsed = time.perf_counter() - t0
        original_logp = float(row["log_probability"]) if "log_probability" in row and pd.notna(row["log_probability"]) else np.nan
        rows.append(
            {
                "source_row": int(row_index),
                "source_step": int(row["step"]) if "step" in row and pd.notna(row["step"]) else np.nan,
                "source_chain_id": int(row["chain_id"]) if "chain_id" in row and pd.notna(row["chain_id"]) else np.nan,
                "validation_mode": likelihood_mode,
                "validation_seconds": elapsed,
                "source_log_probability": original_logp,
                "validation_log_probability": validated_logp,
                "delta_log_probability": validated_logp - original_logp if np.isfinite(original_logp) else np.nan,
                **vector_to_row(
                    vector,
                    log_probability=validated_logp,
                    step=int(row["step"]) if "step" in row and pd.notna(row["step"]) else -1,
                    chain_id=int(row["chain_id"]) if "chain_id" in row and pd.notna(row["chain_id"]) else -1,
                    sampler_name="validation",
                    likelihood_mode=likelihood_mode,
                ),
            }
        )
    return pd.DataFrame(rows)


def full_log_probability_vector(
    galaxy: GalaxyData,
    vector: np.ndarray,
    *,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    systemic_velocity_bounds: tuple[float, float] = (-80.0, 50.0),
    use_mge_physicality_check: bool = True,
) -> float:
    prior = log_prior_vector(galaxy, vector, systemic_velocity_bounds=systemic_velocity_bounds)
    if not np.isfinite(prior):
        return -np.inf
    slow = SlowParams(
        q_halo=float(vector[0]),
        log10_b_halo_pc=float(vector[1]),
        minus_log10_one_minus_beta_z=float(vector[3]),
        alpha=float(vector[4]),
        beta=float(vector[5]),
        gamma=float(vector[6]),
        inclination_deg=float(vector[7]),
    )
    sigma_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode=likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
        strict_epsrel=strict_epsrel,
        use_mge_physicality_check=use_mge_physicality_check,
    )
    return prior + log_likelihood_from_unit_sigma(
        galaxy,
        sigma_unit,
        log10_rho0_msun_pc3=float(vector[2]),
        systemic_velocity_kms=float(vector[8]),
    )


def log_prior_vector(
    galaxy: GalaxyData,
    vector: np.ndarray,
    *,
    systemic_velocity_bounds: tuple[float, float] = (-80.0, 50.0),
) -> float:
    if len(vector) != 9:
        return -np.inf
    slow = SlowParams(
        q_halo=float(vector[0]),
        log10_b_halo_pc=float(vector[1]),
        minus_log10_one_minus_beta_z=float(vector[3]),
        alpha=float(vector[4]),
        beta=float(vector[5]),
        gamma=float(vector[6]),
        inclination_deg=float(vector[7]),
    )
    return log_prior_slow(slow, galaxy) + log_prior_fast(
        float(vector[2]),
        float(vector[8]),
        systemic_velocity_bounds=systemic_velocity_bounds,
    )


def initial_vector(galaxy: GalaxyData | None = None) -> np.ndarray:
    slow = initial_slow_params()
    systemic = -13.298645 if galaxy is None else float(np.median(galaxy.velocity_kms))
    return np.array([
        slow.q_halo,
        slow.log10_b_halo_pc,
        -1.961219,
        slow.minus_log10_one_minus_beta_z,
        slow.alpha,
        slow.beta,
        slow.gamma,
        slow.inclination_deg,
        systemic,
    ])


def row_to_vector(row: pd.Series) -> np.ndarray:
    q_beta = (
        float(row["minus_log10_one_minus_beta_z"])
        if "minus_log10_one_minus_beta_z" in row and pd.notna(row["minus_log10_one_minus_beta_z"])
        else q_from_beta_z(float(row["beta_z"]))
    )
    return np.array([
        float(row["q_halo"]),
        float(row["log10_b_halo_pc"]),
        float(row["log10_rho0_msun_pc3"]),
        q_beta,
        float(row["alpha"]),
        float(row["beta"]),
        float(row["gamma"]),
        float(row["i_deg"]),
        float(row["systemic_velocity_kms"]),
    ])


def vector_to_row(
    vector: np.ndarray,
    *,
    log_probability: float,
    step: int,
    chain_id: int,
    sampler_name: str,
    likelihood_mode: str,
) -> dict[str, float | int | str]:
    return {
        "chain_id": chain_id,
        "step": step,
        "q_halo": float(vector[0]),
        "log10_b_halo_pc": float(vector[1]),
        "log10_rho0_msun_pc3": float(vector[2]),
        "minus_log10_one_minus_beta_z": float(vector[3]),
        "beta_z": beta_z_from_q(float(vector[3])),
        "alpha": float(vector[4]),
        "beta": float(vector[5]),
        "gamma": float(vector[6]),
        "i_deg": float(vector[7]),
        "systemic_velocity_kms": float(vector[8]),
        "log_probability": float(log_probability),
        "sampler": sampler_name,
        "likelihood_mode": likelihood_mode,
    }


def chain_with_anisotropy_columns(chain: pd.DataFrame) -> pd.DataFrame:
    chain = chain.copy()
    if "minus_log10_one_minus_beta_z" not in chain.columns and "beta_z" in chain.columns:
        beta_z = pd.to_numeric(chain["beta_z"], errors="coerce")
        valid = 1.0 - beta_z > 0.0
        chain["minus_log10_one_minus_beta_z"] = np.nan
        chain.loc[valid, "minus_log10_one_minus_beta_z"] = -np.log10(1.0 - beta_z.loc[valid])
    if "beta_z" not in chain.columns and "minus_log10_one_minus_beta_z" in chain.columns:
        q_beta = pd.to_numeric(chain["minus_log10_one_minus_beta_z"], errors="coerce")
        chain["beta_z"] = 1.0 - 10.0 ** (-q_beta)
    return chain


def chain_after_burn(chain: pd.DataFrame, burn: int) -> pd.DataFrame:
    if burn < 0:
        raise ValueError("--burn must be non-negative")
    samples = chain.iloc[burn:].copy()
    if samples.empty:
        raise ValueError(f"no samples remain after burn={burn}; chain has {len(chain)} rows")
    return samples


def plot_density_profile(galaxy: GalaxyData, chain: pd.DataFrame, output: Path, *, burn: int = 0) -> None:
    samples = chain_after_burn(chain, burn)
    radius_kpc = np.geomspace(0.01, 20.0, 300)
    radius_pc = radius_kpc * 1000.0
    profiles = []
    for row in samples.itertuples(index=False):
        from hayashi_jeans.halos import GeneralizedHernquistHalo

        halo = GeneralizedHernquistHalo(
            q=float(row.q_halo),
            b_pc=10.0 ** float(row.log10_b_halo_pc),
            rho0_msun_pc3=10.0 ** float(row.log10_rho0_msun_pc3),
            alpha=float(row.alpha),
            beta=float(row.beta),
            gamma=float(row.gamma),
        )
        profiles.append([halo.density(float(r), 0.0) * 1.0e9 for r in radius_pc])
    profiles = np.asarray(profiles)
    p16, median, p84 = np.percentile(profiles, [16.0, 50.0, 84.0], axis=0)

    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=180)
    ax.fill_between(radius_kpc, p16, p84, color="#4c78a8", alpha=0.28, lw=0, label="68% interval")
    ax.plot(radius_kpc, median, color="#1f5f8b", lw=2.2, label="median")
    ax.axvline(galaxy.observables.b_star_pc / 1000.0, color="#333333", lw=1.4, ls="--", label=r"$b_*$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 20.0)
    ax.set_ylim(1.0e4, 1.0e10)
    ax.set_box_aspect(1.04)
    ax.set_xlabel("Major Axis [kpc]")
    ax.set_ylabel(r"$\rho_{\rm DM}(r)$ [$M_\odot\,{\rm kpc}^{-3}$]")
    ax.set_title(f"{galaxy.observables.galaxy} dark-matter density profile")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, which="both", alpha=0.22, lw=0.6)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    print(f"wrote {output}")


def generate_corner_map(
    chain: pd.DataFrame,
    output: Path,
    *,
    show_systemic_velocity: bool = False,
    burn: int = 0,
) -> None:
    chain = chain_with_anisotropy_columns(chain)
    columns = [
        "q_halo",
        "log10_b_halo_pc",
        "log10_rho0_msun_pc3",
        "minus_log10_one_minus_beta_z",
        "alpha",
        "beta",
        "gamma",
        "i_deg",
    ]
    labels = [
        r"$q_{\rm DM}$",
        r"$\log_{10} b_{\rm halo}\,[{\rm pc}]$",
        r"$\log_{10}\rho_0\,[M_\odot\,{\rm pc}^{-3}]$",
        r"$-\log_{10}(1-\beta_z)$",
        r"$\alpha$",
        r"$\beta$",
        r"$\gamma$",
        r"$i\,[{\rm deg}]$",
    ]
    if show_systemic_velocity:
        columns.append("systemic_velocity_kms")
        labels.append(r"$\langle u\rangle\,[{\rm km\,s}^{-1}]$")

    missing = [column for column in columns if column not in chain.columns]
    if missing:
        raise ValueError(f"chain is missing columns required for corner map: {missing}")
    samples = chain_after_burn(chain, burn)[columns].replace([np.inf, -np.inf], np.nan).dropna()
    if samples.empty:
        raise ValueError("no valid samples remain after burn and finite-value filtering")

    output.parent.mkdir(parents=True, exist_ok=True)
    if len(samples) <= len(columns):
        _plot_scatter_matrix_corner_fallback(
            samples,
            output,
            labels=labels,
            title="too few samples for corner; using scatter-matrix fallback",
        )
        print(f"wrote {output}")
        return

    try:
        import corner

        figure = corner.corner(
            samples.to_numpy(),
            labels=labels,
            quantiles=[0.16, 0.50, 0.84],
            show_titles=True,
            title_fmt=".3g",
            title_kwargs={"fontsize": 9},
            label_kwargs={"fontsize": 10},
            plot_datapoints=False,
            fill_contours=True,
            smooth=0.8,
            color="#1f5f8b",
            hist_kwargs={"density": True},
        )
        figure.savefig(output, dpi=180, bbox_inches="tight")
        plt.close(figure)
    except ImportError:
        _plot_scatter_matrix_corner_fallback(
            samples,
            output,
            labels=labels,
            title="corner package not installed; using scatter-matrix fallback",
        )
    print(f"wrote {output}")


def _plot_scatter_matrix_corner_fallback(
    samples: pd.DataFrame,
    output: Path,
    *,
    labels: list[str],
    title: str,
) -> None:
    axes = pd.plotting.scatter_matrix(
        samples,
        figsize=(1.55 * samples.shape[1], 1.55 * samples.shape[1]),
        diagonal="hist",
        alpha=0.35,
        color="#1f5f8b",
        hist_kwds={"density": True, "color": "#4c78a8"},
    )
    for ax in axes.ravel():
        ax.tick_params(labelsize=7)
        ax.xaxis.label.set_size(8)
        ax.yaxis.label.set_size(8)
    for column_index, label in enumerate(labels):
        axes[-1, column_index].set_xlabel(label)
        axes[column_index, 0].set_ylabel(label)
    fig = axes[0, 0].figure
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def print_summary(chain: pd.DataFrame) -> None:
    chain = chain_with_anisotropy_columns(chain)
    print("parameter medians and 16/84 percentiles:")
    for col in [
        "q_halo",
        "log10_b_halo_pc",
        "log10_rho0_msun_pc3",
        "minus_log10_one_minus_beta_z",
        "beta_z",
        "alpha",
        "beta",
        "gamma",
        "i_deg",
        "systemic_velocity_kms",
        "log_probability",
    ]:
        q16, q50, q84 = chain[col].quantile([0.16, 0.50, 0.84])
        print(f"{col}: {q50:.6g} -{q50-q16:.6g} +{q84-q50:.6g}")


if __name__ == "__main__":
    main()
