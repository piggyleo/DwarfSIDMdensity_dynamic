#!/usr/bin/env python
"""Run Willman 1 with the optional nautilus nested sampler adapter.

This script intentionally reuses the vector-level probability functions from
``run_willman1_block_mh_fast.py`` so sampler comparisons do not change the
Hayashi Jeans likelihood, data loading, projection, or parameter schema.
"""

from __future__ import annotations

import argparse
import inspect
import os
import sys
import time
from pathlib import Path
from typing import Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / "outputs" / ".matplotlib"))
for path in (PROJECT_ROOT, PROJECT_ROOT / "src", PROJECT_ROOT / "scripts"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import numpy as np
import pandas as pd

from hayashi_jeans.data import GalaxyData, load_galaxy_data
try:
    from scripts.run_willman1_block_mh_fast import (
        LIKELIHOOD_MODES,
        full_log_probability_vector,
        initial_vector,
        log_prior_vector,
        print_summary,
        vector_to_row,
    )
except ModuleNotFoundError:
    from run_willman1_block_mh_fast import (
        LIKELIHOOD_MODES,
        full_log_probability_vector,
        initial_vector,
        log_prior_vector,
        print_summary,
        vector_to_row,
    )

PARAMETER_NAMES = [
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

CHAIN_COLUMNS = [
    "chain_id",
    "step",
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
    "sampler",
    "likelihood_mode",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("point", "smoke", "sample", "schema-check", "all"),
        default="all",
        help="point: one log-probability call; smoke: short nautilus run; sample: longer nautilus run.",
    )
    parser.add_argument("--likelihood-mode", choices=LIKELIHOOD_MODES, default="fast")
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/27_Willman_1.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--chain-output", default=str(PROJECT_ROOT / "outputs/willman1_nautilus_chain.csv"))
    parser.add_argument("--checkpoint-output", default=str(PROJECT_ROOT / "outputs/diagnostics/willman1_nautilus_sampler.h5"))
    parser.add_argument("--seed", type=int, default=20260526)
    parser.add_argument("--n-r", type=int, default=96)
    parser.add_argument("--n-z", type=int, default=192)
    parser.add_argument("--n-los", type=int, default=160)
    parser.add_argument("--n-force-r", type=int, default=48)
    parser.add_argument("--n-force-z", type=int, default=96)
    parser.add_argument("--strict-epsrel", type=float, default=1.5e-2)
    parser.add_argument("--n-live", type=int, default=120)
    parser.add_argument("--n-eff", type=float, default=300.0)
    parser.add_argument("--n-like-max", type=float, default=np.inf)
    parser.add_argument("--f-live", type=float, default=0.05)
    parser.add_argument("--n-shell", type=int, default=1)
    parser.add_argument("--n-batch", type=int, default=None)
    parser.add_argument("--n-like-new-bound", type=int, default=None)
    parser.add_argument(
        "--n-processes",
        "--n-workers",
        dest="n_processes",
        type=int,
        default=1,
        help="Number of worker processes used by nautilus for likelihood evaluations.",
    )
    parser.add_argument("--timeout", type=float, default=np.inf)
    parser.add_argument("--resume", action="store_true", help="Resume from --checkpoint-output if it exists.")
    parser.add_argument("--no-checkpoint", action="store_true", help="Do not write a nautilus .h5 checkpoint.")
    parser.add_argument("--weighted-output", action="store_true", help="Write unequal-weight posterior rows plus log weights.")
    parser.add_argument("--equal-weight-boost", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--progress-newlines",
        action="store_true",
        help="With --verbose, print nautilus status updates on separate lines so Slurm .out files keep progress history.",
    )
    parser.add_argument("--smoke-n-live", type=int, default=30)
    parser.add_argument("--smoke-n-eff", type=float, default=20.0)
    parser.add_argument("--smoke-n-like-max", type=float, default=120.0)
    parser.add_argument("--smoke-timeout", type=float, default=300.0)
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)

    if args.mode in ("point", "all"):
        run_point_check(galaxy, args)

    if args.mode in ("schema-check", "all"):
        run_schema_check(galaxy, args)

    if args.mode in ("smoke", "all"):
        smoke_output = Path(args.chain_output).with_name(Path(args.chain_output).stem + "_smoke.csv")
        run_nautilus(
            galaxy,
            args,
            chain_output=smoke_output,
            n_live=args.smoke_n_live,
            n_eff=args.smoke_n_eff,
            n_like_max=args.smoke_n_like_max,
            timeout=args.smoke_timeout,
            checkpoint_output=None,
            resume=False,
            run_label="smoke",
        )

    if args.mode in ("sample", "all"):
        checkpoint_output = None if args.no_checkpoint else Path(args.checkpoint_output)
        run_nautilus(
            galaxy,
            args,
            chain_output=Path(args.chain_output),
            n_live=args.n_live,
            n_eff=args.n_eff,
            n_like_max=args.n_like_max,
            timeout=args.timeout,
            checkpoint_output=checkpoint_output,
            resume=args.resume,
            run_label="nautilus",
        )


def run_point_check(galaxy: GalaxyData, args: argparse.Namespace) -> None:
    vector = initial_vector()
    t0 = time.perf_counter()
    logp = evaluate_log_probability(galaxy, vector, args)
    elapsed = time.perf_counter() - t0
    prior = log_prior_vector(galaxy, vector)
    print("single_point_check")
    print(f"vector={dict(zip(PARAMETER_NAMES, vector))}")
    print(f"log_prior={prior:.6f}")
    print(f"log_probability={logp:.6f}")
    print(f"seconds={elapsed:.6f}")


def run_schema_check(galaxy: GalaxyData, args: argparse.Namespace) -> None:
    vector = initial_vector()
    logp = evaluate_log_probability(galaxy, vector, args)
    row = vector_to_row(
        vector,
        log_probability=logp,
        step=1,
        chain_id=0,
        sampler_name="nautilus",
        likelihood_mode=args.likelihood_mode,
    )
    chain = pd.DataFrame([row])
    missing = [column for column in CHAIN_COLUMNS if column not in chain.columns]
    if missing:
        raise ValueError(f"schema check failed; missing columns: {missing}")
    print("schema_check_passed")
    print(",".join(chain.columns))


def run_nautilus(
    galaxy: GalaxyData,
    args: argparse.Namespace,
    *,
    chain_output: Path,
    n_live: int,
    n_eff: float,
    n_like_max: float,
    timeout: float,
    checkpoint_output: Path | None,
    resume: bool,
    run_label: str,
) -> pd.DataFrame:
    try:
        from nautilus import Prior, Sampler
    except ImportError as exc:  # pragma: no cover - optional sampler dependency
        raise ImportError("nautilus is required for this mode; install the 'nautilus-sampler' package.") from exc

    if n_live < 2:
        raise ValueError("n_live must be >= 2")
    prior = make_nautilus_prior(Prior, galaxy)
    likelihood = make_log_probability_callable(galaxy, args)
    filepath = str(checkpoint_output) if checkpoint_output is not None else None
    if checkpoint_output is not None:
        checkpoint_output.parent.mkdir(parents=True, exist_ok=True)

    sampler = Sampler(
        prior,
        likelihood,
        n_live=n_live,
        n_batch=args.n_batch,
        n_like_new_bound=args.n_like_new_bound,
        pool=args.n_processes if args.n_processes > 1 else None,
        seed=args.seed,
        filepath=filepath,
        resume=resume,
        pass_dict=False,
    )
    if args.progress_newlines:
        enable_progress_newlines(sampler)
    success = sampler.run(
        f_live=args.f_live,
        n_shell=args.n_shell,
        n_eff=n_eff,
        n_like_max=n_like_max,
        timeout=timeout,
        verbose=args.verbose,
    )
    chain = sampler_to_chain(
        sampler,
        likelihood_mode=args.likelihood_mode,
        weighted_output=args.weighted_output,
        equal_weight_boost=args.equal_weight_boost,
    )
    chain_output.parent.mkdir(parents=True, exist_ok=True)
    chain.to_csv(chain_output, index=False)
    print(f"{run_label}_success={success}")
    print(f"wrote {chain_output}")
    print_nautilus_summary(sampler)
    print_summary(chain)
    return chain


def enable_progress_newlines(sampler) -> None:
    original_print_status = sampler.print_status

    def print_status_with_newline(status: str = "", header: bool = False, end: str = "\n") -> None:
        original_print_status(status=status, header=header, end="\n" if end == "\r" else end)

    sampler.print_status = print_status_with_newline


def make_nautilus_prior(prior_class: type, galaxy: GalaxyData):
    min_i = float(np.degrees(np.arccos(galaxy.observables.qprime)))
    prior = prior_class()
    bounds = [
        (0.1, 2.0),
        (0.0, 5.0),
        (-5.0, 5.0),
        (-1.0, 1.0),
        (0.5, 3.0),
        (3.0, 10.0),
        (0.0, 2.0),
        (min_i, 90.0),
        (-80.0, 50.0),
    ]
    for name, bound in zip(PARAMETER_NAMES, bounds):
        prior.add_parameter(name, dist=bound)
    return prior


def make_log_probability_callable(galaxy: GalaxyData, args: argparse.Namespace):
    def log_probability(vector: np.ndarray) -> float:
        return evaluate_log_probability(galaxy, np.asarray(vector, dtype=float), args)

    return log_probability


def evaluate_log_probability(galaxy: GalaxyData, vector: np.ndarray, args: argparse.Namespace) -> float:
    return full_log_probability_vector(
        galaxy,
        vector,
        likelihood_mode=args.likelihood_mode,
        n_r=args.n_r,
        n_z=args.n_z,
        n_los=args.n_los,
        n_force_r=args.n_force_r,
        n_force_z=args.n_force_z,
        strict_epsrel=args.strict_epsrel,
    )


def sampler_to_chain(
    sampler,
    *,
    likelihood_mode: str,
    weighted_output: bool,
    equal_weight_boost: float,
) -> pd.DataFrame:
    if weighted_output:
        points, log_weights, log_prob = sampler.posterior(return_as_dict=False, equal_weight=False)
        rows = rows_from_points(points, log_prob, likelihood_mode=likelihood_mode)
        chain = pd.DataFrame(rows)
        chain["log_weight"] = np.asarray(log_weights, dtype=float)
        chain["weight"] = np.exp(chain["log_weight"] - np.max(chain["log_weight"]))
        return chain

    posterior_kwargs = {"return_as_dict": False, "equal_weight": True}
    if "equal_weight_boost" in inspect.signature(sampler.posterior).parameters:
        posterior_kwargs["equal_weight_boost"] = equal_weight_boost
    points, _log_weights, log_prob = sampler.posterior(**posterior_kwargs)
    return pd.DataFrame(rows_from_points(points, log_prob, likelihood_mode=likelihood_mode))


def rows_from_points(
    points: np.ndarray | dict[str, Iterable[float]],
    log_prob: np.ndarray,
    *,
    likelihood_mode: str,
) -> list[dict[str, float | int | str]]:
    if isinstance(points, dict):
        points_array = np.column_stack([points[name] for name in PARAMETER_NAMES])
    else:
        points_array = np.asarray(points, dtype=float)
    if points_array.ndim == 1:
        points_array = points_array.reshape(1, -1)

    rows = []
    for index, vector in enumerate(points_array):
        rows.append(
            vector_to_row(
                vector,
                log_probability=float(np.asarray(log_prob)[index]),
                step=index + 1,
                chain_id=0,
                sampler_name="nautilus",
                likelihood_mode=likelihood_mode,
            )
        )
    return rows


def print_nautilus_summary(sampler) -> None:
    for name in ("log_z", "n_eff", "f_live", "n_like", "n_like_iter", "n_batch", "n_live"):
        value = getattr(sampler, name, None)
        if value is not None:
            if isinstance(value, (int, np.integer)):
                print(f"{name}={int(value)}")
            else:
                print(f"{name}={float(value):.6g}")


if __name__ == "__main__":
    main()
