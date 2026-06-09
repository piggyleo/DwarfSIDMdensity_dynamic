#!/usr/bin/env python
"""Recheck the previous five Eridanus II MGE-physical samples after MGE changes."""

from __future__ import annotations

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
)


SOURCE_ROWS = [2377, 2762, 1240, 1032, 1596]


def main() -> None:
    chain_path = PROJECT_ROOT / "outputs/eridanus_ii_nautilus_5w_chain.csv"
    galaxy_csv = PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"
    centers_csv = PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"
    detail_path = PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fixed5_analytic_mge_vs_strict_detail.csv"
    summary_path = PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fixed5_analytic_mge_vs_strict_summary.csv"

    galaxy = load_galaxy_data(galaxy_csv, structural_centers_csv=centers_csv)
    systemic_bounds = resolve_systemic_velocity_prior_bounds(
        galaxy,
        mode="adaptive",
        manual_bounds=None,
        padding=20.0,
        min_half_width=20.0,
    )
    chain = pd.read_csv(chain_path).replace([np.inf, -np.inf], np.nan)

    rows = []
    for i, source_row in enumerate(SOURCE_ROWS, start=1):
        row = chain.iloc[source_row]
        physical_now = is_mge_physical(galaxy, slow_from_row(row))
        print(f"sample {i}/{len(SOURCE_ROWS)} source_row={source_row} analytic_mge_physical={physical_now}", flush=True)
        result = compare_one(
            galaxy,
            int(source_row),
            row.to_dict(),
            systemic_bounds,
            strict_epsrel=1.5e-2,
        )
        result["analytic_mge_physical"] = bool(physical_now)
        rows.append(result)
        print(
            f"  dlogL={result['mge_minus_strict_log_likelihood']:.6g}; "
            f"dlogP={result['mge_minus_strict_log_probability']:.6g}; "
            f"mge={result['mge_likelihood_seconds']:.3f}s strict={result['strict_likelihood_seconds']:.3f}s",
            flush=True,
        )

    detail = pd.DataFrame(rows)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(detail_path, index=False)

    finite_detail = detail.replace([np.inf, -np.inf], np.nan)
    finite_for_summary = finite_detail.dropna(subset=["abs_mge_minus_strict_log_likelihood"])
    summary = summarize(finite_for_summary) if len(finite_for_summary) else pd.DataFrame()
    summary.to_csv(summary_path, index=False)

    print(f"wrote {detail_path}")
    print(detail.to_string(index=False))
    print(f"wrote {summary_path}")
    print(summary.to_string(index=False))


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
        use_mge_physicality_check=True,
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
        use_mge_physicality_check=True,
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


def is_mge_physical(galaxy, slow: SlowParams) -> bool:
    sigma_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode="mge",
        n_r=96,
        n_z=192,
        n_los=160,
        n_force_r=48,
        n_force_z=96,
        strict_epsrel=1.5e-2,
        use_mge_physicality_check=True,
    )
    return bool(np.all(np.isfinite(sigma_unit)) and np.all(sigma_unit > 0.0))


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
