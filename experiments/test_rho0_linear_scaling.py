#!/usr/bin/env python
"""Experiment: rho0 linear scaling on top of R-z grid + cubic dPz/dR.

This standalone test does not modify production code. It verifies that, for a
fixed halo shape/geometric parameter set, sigma_los^2 scales linearly with the
density normalization rho0 when using the experimental R-z moment grid.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from experiments.test_rz_moment_grid_interpolation import RZMomentGridProjector
from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.params import HayashiParameters


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default="data/galaxies/27_Willman_1.csv")
    parser.add_argument("--centers-csv", default="data/processed/galaxy_structural_centers.csv")
    parser.add_argument("--out-dir", default="outputs/diagnostics")
    parser.add_argument("--log-path", default="logs/likelihood_acceleration_strategy_log.md")
    parser.add_argument("--n-r", type=int, default=96)
    parser.add_argument("--n-z", type=int, default=192)
    parser.add_argument("--n-los", type=int, default=160)
    parser.add_argument("--log-rho-values", default="-2.5,-1.961219,-1.5,-0.5")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    log_rho_values = [float(value) for value in args.log_rho_values.split(",")]

    unit_params = make_params(log_rho0=0.0)
    t0 = time.perf_counter()
    unit_projector = build_projector(galaxy, unit_params, args.n_r, args.n_z, args.n_los)
    unit_build_time = time.perf_counter() - t0
    t0 = time.perf_counter()
    sigma_unit = unit_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    unit_sigma_time = time.perf_counter() - t0

    likelihood = GaussianVelocityLikelihood(galaxy)
    rows = []
    direct_total_time = 0.0
    scaled_total_time = unit_build_time + unit_sigma_time

    for log_rho0 in log_rho_values:
        direct_params = make_params(log_rho0=log_rho0)
        t0 = time.perf_counter()
        direct_projector = build_projector(galaxy, direct_params, args.n_r, args.n_z, args.n_los)
        direct_build_time = time.perf_counter() - t0
        t0 = time.perf_counter()
        direct_sigma = direct_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
        direct_sigma_time = time.perf_counter() - t0
        direct_total_time += direct_build_time + direct_sigma_time

        t0 = time.perf_counter()
        scaled_sigma = sigma_unit * 10.0**log_rho0
        scaled_eval_time = time.perf_counter() - t0
        scaled_total_time += scaled_eval_time

        direct_ln_l = likelihood.log_likelihood(direct_sigma, systemic_velocity_kms=direct_params.systemic_velocity_kms)
        scaled_ln_l = likelihood.log_likelihood(scaled_sigma, systemic_velocity_kms=direct_params.systemic_velocity_kms)
        rel_error = np.abs(scaled_sigma - direct_sigma) / np.maximum(np.abs(direct_sigma), 1e-30)

        for star_id, x_pc, y_pc, direct_value, scaled_value, error in zip(
            galaxy.stars["star_id"],
            galaxy.x_pc,
            galaxy.y_pc,
            direct_sigma,
            scaled_sigma,
            rel_error,
        ):
            rows.append(
                {
                    "log10_rho0_msun_pc3": log_rho0,
                    "star_id": star_id,
                    "x_pc": x_pc,
                    "y_pc": y_pc,
                    "direct_sigma_los2": direct_value,
                    "scaled_sigma_los2": scaled_value,
                    "rel_error": error,
                    "direct_ln_l": direct_ln_l,
                    "scaled_ln_l": scaled_ln_l,
                    "delta_ln_l": scaled_ln_l - direct_ln_l,
                    "direct_build_time_s": direct_build_time,
                    "direct_sigma_time_s": direct_sigma_time,
                    "scaled_eval_time_s": scaled_eval_time,
                }
            )

    detail = pd.DataFrame(rows)
    detail_path = out_dir / "rho0_scaling_willman1_detail.csv"
    detail.to_csv(detail_path, index=False)

    grouped = detail.groupby("log10_rho0_msun_pc3").agg(
        sigma_rel_error_median=("rel_error", "median"),
        sigma_rel_error_max=("rel_error", "max"),
        sigma_rel_error_p95=("rel_error", lambda s: s.quantile(0.95)),
        delta_ln_l=("delta_ln_l", "first"),
        direct_build_time_s=("direct_build_time_s", "first"),
        direct_sigma_time_s=("direct_sigma_time_s", "first"),
        scaled_eval_time_s=("scaled_eval_time_s", "first"),
    ).reset_index()
    summary = {
        "n_star": len(galaxy.stars),
        "n_r": args.n_r,
        "n_z": args.n_z,
        "n_los": args.n_los,
        "n_rho_values": len(log_rho_values),
        "unit_build_time_s": unit_build_time,
        "unit_sigma_time_s": unit_sigma_time,
        "direct_total_time_s": direct_total_time,
        "scaled_total_time_s": scaled_total_time,
        "speedup_for_rho_sweep": direct_total_time / scaled_total_time,
        "overall_sigma_rel_error_median": detail["rel_error"].median(),
        "overall_sigma_rel_error_max": detail["rel_error"].max(),
        "overall_sigma_rel_error_p95": detail["rel_error"].quantile(0.95),
        "max_abs_delta_ln_l": detail.groupby("log10_rho0_msun_pc3")["delta_ln_l"].first().abs().max(),
    }
    grouped_path = out_dir / "rho0_scaling_willman1_by_rho.csv"
    summary_path = out_dir / "rho0_scaling_willman1_summary.csv"
    grouped.to_csv(grouped_path, index=False)
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    append_log(Path(args.log_path), summary, detail_path, grouped_path, summary_path)

    print("By rho0:")
    print(grouped.to_string(index=False))
    print("\nSummary:")
    print(pd.Series(summary).to_string())


def make_params(*, log_rho0: float) -> HayashiParameters:
    return HayashiParameters(
        q_halo=1.0,
        b_halo_pc=10.0**3.2,
        rho0_msun_pc3=10.0**log_rho0,
        beta_z=0.3,
        alpha=1.8,
        beta=6.4,
        gamma=1.2,
        inclination_rad=np.deg2rad(75.0),
        systemic_velocity_kms=-13.298645,
    )


def build_projector(galaxy, params, n_r: int, n_z: int, n_los: int) -> RZMomentGridProjector:
    return RZMomentGridProjector(
        galaxy=galaxy,
        params=params,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        zmax_factor=20.0,
        los_factor=20.0,
        r_min_pc=1e-3,
        z_min_pc=1e-3,
        grid_padding=1.08,
        radial_derivative="cubic_spline",
    )


def append_log(log_path: Path, summary: dict[str, float], detail_path: Path, grouped_path: Path, summary_path: Path) -> None:
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("\n\n### Strategy 2 Test Run 001\n\n")
        handle.write("Status: completed isolated rho0 linear-scaling test using R-z grid + cubic dP_z/dR; production code unchanged.\n\n")
        handle.write("- Method: build one `rho0=1` R-z grid, scale `sigma_los2_unit` by `10**log10_rho0`, and compare to direct grid rebuilds at each rho0.\n")
        handle.write(f"- detail CSV: `{detail_path}`\n")
        handle.write(f"- by-rho CSV: `{grouped_path}`\n")
        handle.write(f"- summary CSV: `{summary_path}`\n")
        for key, value in summary.items():
            handle.write(f"- {key}: {value}\n")


if __name__ == "__main__":
    main()
