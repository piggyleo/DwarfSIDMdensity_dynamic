#!/usr/bin/env python
"""Compare MGE likelihood against validation-strict for random physical samples."""

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
from scripts.run_willman1_block_mh_fast import (
    SlowParams,
    compute_sigma_unit,
    full_log_probability_vector,
    log_likelihood_from_unit_sigma,
    log_prior_vector,
    resolve_systemic_velocity_prior_bounds,
    row_to_vector,
    slow_params_mge_physical,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain", default=str(PROJECT_ROOT / "outputs/eridanus_ii_nautilus_5w_chain.csv"))
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--n-samples", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--strict-epsrel", type=float, default=1.5e-2)
    parser.add_argument("--output-detail", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_random_physical_mge_vs_strict_detail.csv"))
    parser.add_argument("--output-summary", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_random_physical_mge_vs_strict_summary.csv"))
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
    ]
    finite = chain.dropna(subset=required).copy()
    finite = finite.loc[
        finite.apply(
            lambda row: np.isfinite(
                log_prior_vector(galaxy, row_to_vector(row), systemic_velocity_bounds=systemic_bounds)
            ),
            axis=1,
        )
    ].copy()

    selected = select_random_physical(galaxy, finite, n_samples=args.n_samples, seed=args.seed)
    rows = []
    for i, (source_row, row) in enumerate(selected, start=1):
        row_dict = row.to_dict()
        print(f"sample {i}/{len(selected)} source_row={source_row}", flush=True)
        rows.append(compare_one(galaxy, int(source_row), row_dict, systemic_bounds, strict_epsrel=args.strict_epsrel))
        print(
            f"  dlogL={rows[-1]['mge_minus_strict_log_likelihood']:.6g}; "
            f"dlogP={rows[-1]['mge_minus_strict_log_probability']:.6g}; "
            f"mge={rows[-1]['mge_likelihood_seconds']:.3f}s strict={rows[-1]['strict_likelihood_seconds']:.3f}s",
            flush=True,
        )

    detail = pd.DataFrame(rows)
    detail_path = Path(args.output_detail)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(detail_path, index=False)

    summary = summarize(detail)
    summary_path = Path(args.output_summary)
    summary.to_csv(summary_path, index=False)

    print(f"wrote {detail_path}")
    print(detail.to_string(index=False))
    print(f"wrote {summary_path}")
    print(summary.to_string(index=False))


def select_random_physical(galaxy, finite: pd.DataFrame, *, n_samples: int, seed: int):
    rng = np.random.default_rng(seed)
    source_indices = finite.index.to_numpy()
    order = rng.permutation(len(finite))
    selected = []
    checked = 0
    rejected = 0
    for pos in order:
        row = finite.iloc[int(pos)]
        source_row = source_indices[int(pos)]
        checked += 1
        if slow_params_mge_physical(galaxy, slow_from_row(row)):
            selected.append((source_row, row))
            if len(selected) >= n_samples:
                break
        else:
            rejected += 1
    if len(selected) < n_samples:
        raise ValueError(f"found only {len(selected)} MGE-physical samples after checking {checked} rows")
    print(f"selected {len(selected)} MGE-physical samples after checking {checked} rows; rejected_by_mge={rejected}")
    return selected


def compare_one(galaxy, source_row: int, row: dict[str, float], systemic_bounds: tuple[float, float], *, strict_epsrel: float):
    slow = slow_from_row(row)
    vector = row_to_vector(pd.Series(row))
    prior = log_prior_vector(galaxy, vector, systemic_velocity_bounds=systemic_bounds)
    log10_rho0 = float(row["log10_rho0_msun_pc3"])
    systemic = float(row["systemic_velocity_kms"])

    t0 = time.perf_counter()
    mge_sigma_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode="mge",
        n_r=96,
        n_z=192,
        n_los=160,
        n_force_r=48,
        n_force_z=96,
        strict_epsrel=strict_epsrel,
    )
    mge_likelihood_seconds = time.perf_counter() - t0
    mge_log_likelihood = log_likelihood_from_unit_sigma(
        galaxy,
        mge_sigma_unit,
        log10_rho0_msun_pc3=log10_rho0,
        systemic_velocity_kms=systemic,
    )
    mge_log_probability = full_log_probability_vector(
        galaxy,
        vector,
        likelihood_mode="mge",
        n_r=96,
        n_z=192,
        n_los=160,
        n_force_r=48,
        n_force_z=96,
        strict_epsrel=strict_epsrel,
        systemic_velocity_bounds=systemic_bounds,
    )

    t0 = time.perf_counter()
    strict_sigma_unit = compute_sigma_unit(
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
    strict_likelihood_seconds = time.perf_counter() - t0
    strict_log_likelihood = log_likelihood_from_unit_sigma(
        galaxy,
        strict_sigma_unit,
        log10_rho0_msun_pc3=log10_rho0,
        systemic_velocity_kms=systemic,
    )
    strict_log_probability = prior + strict_log_likelihood if np.isfinite(prior) else -np.inf

    sigma2_scale = 10.0**log10_rho0
    mge_sigma2 = mge_sigma_unit * sigma2_scale
    strict_sigma2 = strict_sigma_unit * sigma2_scale
    rel = np.abs(mge_sigma2 - strict_sigma2) / np.maximum(np.abs(strict_sigma2), 1.0e-300)
    sigma_rel = (
        np.abs(np.sqrt(np.maximum(mge_sigma2, 0.0)) - np.sqrt(np.maximum(strict_sigma2, 0.0)))
        / np.maximum(np.sqrt(np.maximum(strict_sigma2, 0.0)), 1.0e-300)
    )

    return {
        "source_row": source_row,
        "step": int(row["step"]) if "step" in row and pd.notna(row["step"]) else np.nan,
        "q_halo": float(row["q_halo"]),
        "log10_b_halo_pc": float(row["log10_b_halo_pc"]),
        "log10_rho0_msun_pc3": log10_rho0,
        "minus_log10_one_minus_beta_z": float(row["minus_log10_one_minus_beta_z"]),
        "alpha": float(row["alpha"]),
        "beta": float(row["beta"]),
        "gamma": float(row["gamma"]),
        "i_deg": float(row["i_deg"]),
        "systemic_velocity_kms": systemic,
        "log_prior": float(prior),
        "mge_log_likelihood": float(mge_log_likelihood),
        "strict_log_likelihood": float(strict_log_likelihood),
        "mge_minus_strict_log_likelihood": float(mge_log_likelihood - strict_log_likelihood),
        "abs_mge_minus_strict_log_likelihood": float(abs(mge_log_likelihood - strict_log_likelihood)),
        "mge_log_probability": float(mge_log_probability),
        "strict_log_probability": float(strict_log_probability),
        "mge_minus_strict_log_probability": float(mge_log_probability - strict_log_probability),
        "abs_mge_minus_strict_log_probability": float(abs(mge_log_probability - strict_log_probability)),
        "mge_likelihood_seconds": float(mge_likelihood_seconds),
        "strict_likelihood_seconds": float(strict_likelihood_seconds),
        "strict_over_mge_seconds": float(strict_likelihood_seconds / mge_likelihood_seconds),
        "sigma_los2_rel_err_median": float(np.nanmedian(rel)),
        "sigma_los2_rel_err_p95": float(np.nanpercentile(rel, 95.0)),
        "sigma_los2_rel_err_max": float(np.nanmax(rel)),
        "sigma_los_rel_err_median": float(np.nanmedian(sigma_rel)),
        "sigma_los_rel_err_p95": float(np.nanpercentile(sigma_rel, 95.0)),
        "sigma_los_rel_err_max": float(np.nanmax(sigma_rel)),
    }


def slow_from_row(row) -> SlowParams:
    return SlowParams(
        q_halo=float(row["q_halo"]),
        log10_b_halo_pc=float(row["log10_b_halo_pc"]),
        minus_log10_one_minus_beta_z=float(row["minus_log10_one_minus_beta_z"]),
        alpha=float(row["alpha"]),
        beta=float(row["beta"]),
        gamma=float(row["gamma"]),
        inclination_deg=float(row["i_deg"]),
    )


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "abs_mge_minus_strict_log_likelihood",
        "abs_mge_minus_strict_log_probability",
        "sigma_los2_rel_err_median",
        "sigma_los2_rel_err_p95",
        "sigma_los2_rel_err_max",
        "sigma_los_rel_err_median",
        "mge_likelihood_seconds",
        "strict_likelihood_seconds",
        "strict_over_mge_seconds",
    ]
    rows = []
    for metric in metrics:
        rows.append(
            {
                "metric": metric,
                "min": float(detail[metric].min()),
                "median": float(detail[metric].median()),
                "max": float(detail[metric].max()),
            }
        )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
