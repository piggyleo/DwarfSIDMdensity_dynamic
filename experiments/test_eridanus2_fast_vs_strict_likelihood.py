#!/usr/bin/env python
"""Compare Eridanus II fast likelihood against validation-strict samples."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo
from scripts.run_willman1_block_mh_fast import (
    SlowParams,
    compute_sigma_unit,
    log_likelihood_from_unit_sigma,
    log_prior_vector,
    resolve_systemic_velocity_prior_bounds,
    row_to_vector,
)


PROFILE_RADII_KPC = np.array([0.01, 0.03, 0.1, 0.196, 0.3, 1.0, 3.0, 10.0, 20.0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain", default=str(PROJECT_ROOT / "outputs/eridanus_ii_nautilus_5w_chain.csv"))
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--n-samples", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260528)
    parser.add_argument("--strict-epsrel", type=float, default=1.5e-2)
    parser.add_argument("--selected-output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_vs_strict_selected_samples.csv"))
    parser.add_argument("--detail-output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_vs_strict_likelihood_detail.csv"))
    parser.add_argument("--summary-output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_vs_strict_likelihood_summary.csv"))
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    systemic_bounds = resolve_systemic_velocity_prior_bounds(
        galaxy,
        mode="adaptive",
        manual_bounds=None,
        padding=20.0,
        min_half_width=20.0,
    )
    chain = pd.read_csv(args.chain).replace([np.inf, -np.inf], np.nan)
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
        "log_probability",
    ]
    finite = chain.dropna(subset=required).copy()
    finite = finite.loc[
        finite.apply(lambda row: np.isfinite(log_prior_vector(galaxy, row_to_vector(row), systemic_velocity_bounds=systemic_bounds)), axis=1)
    ].copy()
    if len(finite) < args.n_samples:
        raise ValueError(f"need at least {args.n_samples} finite prior-valid samples; got {len(finite)}")

    selected = select_diverse_density_profiles(finite, n_samples=args.n_samples, seed=args.seed)
    selected_path = Path(args.selected_output)
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(selected_path, index=False)
    print(f"wrote {selected_path}")
    print(selected[["selection_id", "source_row", "step", "log_probability", *density_column_names()]].to_string(index=False))

    rows = []
    for row in selected.itertuples(index=False):
        print(f"checking {row.selection_id}/{len(selected)} source_row={row.source_row}", flush=True)
        result = compare_one_sample(galaxy, row._asdict(), systemic_bounds, strict_epsrel=args.strict_epsrel)
        rows.append(result)
        print(
            f"  dlogL={result['delta_log_likelihood_fast_minus_strict']:.6g}, "
            f"sigma2_rel_max={result['sigma_los2_rel_err_max']:.6g}, "
            f"fast={result['fast_seconds']:.3f}s strict={result['strict_seconds']:.3f}s",
            flush=True,
        )

    detail = pd.DataFrame(rows)
    detail_path = Path(args.detail_output)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(detail_path, index=False)

    summary = summarize(detail, systemic_bounds=systemic_bounds, strict_epsrel=args.strict_epsrel)
    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    print(f"wrote {detail_path}")
    print(detail.to_string(index=False))
    print(f"wrote {summary_path}")
    print(summary.to_string(index=False))


def select_diverse_density_profiles(chain: pd.DataFrame, *, n_samples: int, seed: int) -> pd.DataFrame:
    features = density_profile_features(chain)
    scaled = (features - np.nanmedian(features, axis=0)) / np.nanstd(features, axis=0)
    scaled = np.nan_to_num(scaled, nan=0.0, posinf=0.0, neginf=0.0)

    source_indices = chain.index.to_numpy()
    logp = chain["log_probability"].to_numpy(float)
    selected_positions = [int(np.nanargmax(logp))]
    min_dist = np.linalg.norm(scaled - scaled[selected_positions[0]], axis=1)
    rng = np.random.default_rng(seed)
    jitter = 1.0e-9 * rng.random(len(chain))
    while len(selected_positions) < n_samples:
        score = min_dist + jitter
        score[selected_positions] = -np.inf
        next_pos = int(np.nanargmax(score))
        selected_positions.append(next_pos)
        min_dist = np.minimum(min_dist, np.linalg.norm(scaled - scaled[next_pos], axis=1))

    selected = chain.iloc[selected_positions].copy()
    selected.insert(0, "selection_id", np.arange(1, len(selected) + 1))
    selected.insert(1, "source_row", source_indices[selected_positions])
    density_features = features[selected_positions]
    for index, column in enumerate(density_column_names()):
        selected[column] = density_features[:, index]
    return selected.reset_index(drop=True)


def density_column_names() -> list[str]:
    return [f"log10_rho_msun_kpc3_at_{radius:g}_kpc" for radius in PROFILE_RADII_KPC]


def density_profile_features(chain: pd.DataFrame) -> np.ndarray:
    radius_pc = PROFILE_RADII_KPC * 1000.0
    features = np.empty((len(chain), len(radius_pc)), dtype=float)
    for i, row in enumerate(chain.itertuples(index=False)):
        halo = GeneralizedHernquistHalo(
            q=float(row.q_halo),
            b_pc=10.0 ** float(row.log10_b_halo_pc),
            rho0_msun_pc3=10.0 ** float(row.log10_rho0_msun_pc3),
            alpha=float(row.alpha),
            beta=float(row.beta),
            gamma=float(row.gamma),
        )
        rho = np.array([halo.density(float(r), 0.0) * 1.0e9 for r in radius_pc])
        features[i] = np.log10(np.maximum(rho, 1.0e-300))
    return features


def compare_one_sample(
    galaxy,
    row: dict[str, float],
    systemic_bounds: tuple[float, float],
    *,
    strict_epsrel: float,
) -> dict[str, float | int]:
    slow = SlowParams(
        q_halo=float(row["q_halo"]),
        log10_b_halo_pc=float(row["log10_b_halo_pc"]),
        minus_log10_one_minus_beta_z=float(row["minus_log10_one_minus_beta_z"]),
        alpha=float(row["alpha"]),
        beta=float(row["beta"]),
        gamma=float(row["gamma"]),
        inclination_deg=float(row["i_deg"]),
    )
    log10_rho0 = float(row["log10_rho0_msun_pc3"])
    systemic = float(row["systemic_velocity_kms"])

    t0 = time.perf_counter()
    fast_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode="fast",
        n_r=96,
        n_z=192,
        n_los=160,
        n_force_r=48,
        n_force_z=96,
        strict_epsrel=strict_epsrel,
    )
    fast_seconds = time.perf_counter() - t0
    fast_log_likelihood = log_likelihood_from_unit_sigma(
        galaxy,
        fast_unit,
        log10_rho0_msun_pc3=log10_rho0,
        systemic_velocity_kms=systemic,
    )

    t0 = time.perf_counter()
    strict_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode="validation-strict",
        n_r=96,
        n_z=192,
        n_los=160,
        n_force_r=48,
        n_force_z=96,
        strict_epsrel=strict_epsrel,
    )
    strict_seconds = time.perf_counter() - t0
    strict_log_likelihood = log_likelihood_from_unit_sigma(
        galaxy,
        strict_unit,
        log10_rho0_msun_pc3=log10_rho0,
        systemic_velocity_kms=systemic,
    )

    scale = 10.0**log10_rho0
    fast_sigma2 = fast_unit * scale
    strict_sigma2 = strict_unit * scale
    sigma2_abs = np.abs(fast_sigma2 - strict_sigma2)
    sigma2_rel = sigma2_abs / np.maximum(np.abs(strict_sigma2), 1.0e-300)
    fast_sigma = np.sqrt(np.maximum(fast_sigma2, 0.0))
    strict_sigma = np.sqrt(np.maximum(strict_sigma2, 0.0))
    sigma_abs = np.abs(fast_sigma - strict_sigma)
    sigma_rel = sigma_abs / np.maximum(np.abs(strict_sigma), 1.0e-300)

    vector = row_to_vector(pd.Series(row))
    prior = log_prior_vector(galaxy, vector, systemic_velocity_bounds=systemic_bounds)
    return {
        "selection_id": int(row["selection_id"]),
        "source_row": int(row["source_row"]),
        "step": int(row["step"]),
        "chain_id": int(row["chain_id"]),
        "log_probability_chain": float(row["log_probability"]),
        "log_prior": float(prior),
        "fast_log_likelihood": float(fast_log_likelihood),
        "strict_log_likelihood": float(strict_log_likelihood),
        "delta_log_likelihood_fast_minus_strict": float(fast_log_likelihood - strict_log_likelihood),
        "abs_delta_log_likelihood": float(abs(fast_log_likelihood - strict_log_likelihood)),
        "fast_seconds": float(fast_seconds),
        "strict_seconds": float(strict_seconds),
        "strict_over_fast_seconds": float(strict_seconds / fast_seconds),
        "sigma_los2_abs_err_max": float(np.nanmax(sigma2_abs)),
        "sigma_los2_abs_err_median": float(np.nanmedian(sigma2_abs)),
        "sigma_los2_rel_err_max": float(np.nanmax(sigma2_rel)),
        "sigma_los2_rel_err_p95": float(np.nanpercentile(sigma2_rel, 95.0)),
        "sigma_los2_rel_err_median": float(np.nanmedian(sigma2_rel)),
        "sigma_los_abs_err_max_kms": float(np.nanmax(sigma_abs)),
        "sigma_los_abs_err_median_kms": float(np.nanmedian(sigma_abs)),
        "sigma_los_rel_err_max": float(np.nanmax(sigma_rel)),
        "sigma_los_rel_err_p95": float(np.nanpercentile(sigma_rel, 95.0)),
        "sigma_los_rel_err_median": float(np.nanmedian(sigma_rel)),
        "q_halo": float(row["q_halo"]),
        "log10_b_halo_pc": float(row["log10_b_halo_pc"]),
        "log10_rho0_msun_pc3": log10_rho0,
        "minus_log10_one_minus_beta_z": float(row["minus_log10_one_minus_beta_z"]),
        "alpha": float(row["alpha"]),
        "beta": float(row["beta"]),
        "gamma": float(row["gamma"]),
        "i_deg": float(row["i_deg"]),
        "systemic_velocity_kms": systemic,
    }


def summarize(detail: pd.DataFrame, *, systemic_bounds: tuple[float, float], strict_epsrel: float) -> pd.DataFrame:
    metrics = [
        "delta_log_likelihood_fast_minus_strict",
        "abs_delta_log_likelihood",
        "fast_seconds",
        "strict_seconds",
        "strict_over_fast_seconds",
        "sigma_los2_rel_err_max",
        "sigma_los2_rel_err_p95",
        "sigma_los2_rel_err_median",
        "sigma_los_rel_err_max",
        "sigma_los_rel_err_p95",
        "sigma_los_rel_err_median",
    ]
    row = {
        "n_samples": len(detail),
        "strict_epsrel": strict_epsrel,
        "systemic_prior_low": systemic_bounds[0],
        "systemic_prior_high": systemic_bounds[1],
    }
    for metric in metrics:
        row[f"{metric}_min"] = float(detail[metric].min())
        row[f"{metric}_median"] = float(detail[metric].median())
        row[f"{metric}_max"] = float(detail[metric].max())
    return pd.DataFrame([row])


if __name__ == "__main__":
    main()
