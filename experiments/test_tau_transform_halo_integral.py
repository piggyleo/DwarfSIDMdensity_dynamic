#!/usr/bin/env python
"""Compare infinite-tau and unit-interval halo force integrals.

This benchmark keeps the same Willman 1 likelihood setup and compares only the
force integral implementation inside GeneralizedHernquistHalo:

- old path: tau in [0, infinity)
- new path: tau = u / (1 - u), u in [0, 1]
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import IntegrationWarning

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from experiments.test_halo_force_grid_interpolation import RZMomentGridWithForceGrid
from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.params import HayashiParameters
from scripts.run_willman1_block_mh_fast import beta_z_from_q, initial_slow_params


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/27_Willman_1.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--summary-output", default=str(PROJECT_ROOT / "outputs/diagnostics/tau_transform_halo_integral_summary.csv"))
    parser.add_argument("--force-output", default=str(PROJECT_ROOT / "outputs/diagnostics/tau_transform_halo_integral_force_detail.csv"))
    parser.add_argument("--n-r", type=int, default=96)
    parser.add_argument("--n-z", type=int, default=192)
    parser.add_argument("--n-los", type=int, default=160)
    parser.add_argument("--n-force-r", type=int, default=48)
    parser.add_argument("--n-force-z", type=int, default=96)
    parser.add_argument("--n-force-check", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260526)
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    slow = initial_slow_params()
    params = HayashiParameters(
        q_halo=slow.q_halo,
        b_halo_pc=10.0**slow.log10_b_halo_pc,
        rho0_msun_pc3=1.0,
        beta_z=beta_z_from_q(slow.minus_log10_one_minus_beta_z),
        alpha=slow.alpha,
        beta=slow.beta,
        gamma=slow.gamma,
        inclination_rad=np.deg2rad(slow.inclination_deg),
        systemic_velocity_kms=-13.298645,
    )
    log10_rho0 = -1.961219

    force_detail = compare_forces(galaxy, params, n_check=args.n_force_check, seed=args.seed)
    force_path = Path(args.force_output)
    force_path.parent.mkdir(parents=True, exist_ok=True)
    force_detail.to_csv(force_path, index=False)

    likelihood_rows = []
    sigma_by_method: dict[str, np.ndarray] = {}
    for method in ("infinite", "unit_interval"):
        row, sigma_los2 = run_single_likelihood(
            galaxy,
            params,
            force_integral_method=method,
            log10_rho0=log10_rho0,
            n_r=args.n_r,
            n_z=args.n_z,
            n_los=args.n_los,
            n_force_r=args.n_force_r,
            n_force_z=args.n_force_z,
        )
        likelihood_rows.append(row)
        sigma_by_method[method] = sigma_los2

    sigma_old = sigma_by_method["infinite"]
    sigma_new = sigma_by_method["unit_interval"]
    sigma_abs = np.abs(sigma_new - sigma_old)
    sigma_rel = sigma_abs / np.maximum(np.abs(sigma_old), 1e-300)

    old = next(row for row in likelihood_rows if row["force_integral_method"] == "infinite")
    new = next(row for row in likelihood_rows if row["force_integral_method"] == "unit_interval")

    comparison = {
        "force_integral_method": "unit_interval_minus_infinite",
        "single_likelihood_seconds": new["single_likelihood_seconds"] - old["single_likelihood_seconds"],
        "log_likelihood": new["log_likelihood"] - old["log_likelihood"],
        "n_integration_warnings": new["n_integration_warnings"] - old["n_integration_warnings"],
        "sigma_los2_abs_max": float(np.nanmax(sigma_abs)),
        "sigma_los2_abs_median": float(np.nanmedian(sigma_abs)),
        "sigma_los2_rel_max": float(np.nanmax(sigma_rel)),
        "sigma_los2_rel_median": float(np.nanmedian(sigma_rel)),
    }

    force_metrics = summarize_force_errors(force_detail)
    summary = pd.DataFrame([*likelihood_rows, comparison, force_metrics])
    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)

    print(f"wrote {summary_path}")
    print(f"wrote {force_path}")
    print(summary.to_string(index=False))


def compare_forces(galaxy, params: HayashiParameters, *, n_check: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    r_extent = max(float(np.nanmax(np.abs(galaxy.x_pc))), 1.0)
    z_extent = max(float(np.nanmax(np.abs(galaxy.y_pc))), galaxy.observables.b_star_pc, 1.0)
    r_values = 10.0 ** rng.uniform(-3.0, np.log10(20.0 * r_extent), size=n_check)
    z_values = rng.choice([-1.0, 1.0], size=n_check) * 10.0 ** rng.uniform(-3.0, np.log10(20.0 * z_extent), size=n_check)

    old_halo = GeneralizedHernquistHalo(**params.halo_kwargs(), force_integral_method="infinite")
    new_halo = GeneralizedHernquistHalo(**params.halo_kwargs(), force_integral_method="unit_interval")
    rows = []
    for index, (r_pc, z_pc) in enumerate(zip(r_values, z_values, strict=True)):
        with warnings.catch_warnings(record=True) as old_warnings:
            warnings.simplefilter("always", IntegrationWarning)
            t0 = time.perf_counter()
            old_r, old_z = old_halo.potential_gradients(float(r_pc), float(z_pc))
            old_seconds = time.perf_counter() - t0
        with warnings.catch_warnings(record=True) as new_warnings:
            warnings.simplefilter("always", IntegrationWarning)
            t0 = time.perf_counter()
            new_r, new_z = new_halo.potential_gradients(float(r_pc), float(z_pc))
            new_seconds = time.perf_counter() - t0
        rows.append(
            {
                "index": index,
                "r_pc": float(r_pc),
                "z_pc": float(z_pc),
                "old_dphi_dr": old_r,
                "new_dphi_dr": new_r,
                "old_dphi_dz": old_z,
                "new_dphi_dz": new_z,
                "abs_err_dphi_dr": abs(new_r - old_r),
                "abs_err_dphi_dz": abs(new_z - old_z),
                "rel_err_dphi_dr": abs(new_r - old_r) / max(abs(old_r), 1e-300),
                "rel_err_dphi_dz": abs(new_z - old_z) / max(abs(old_z), 1e-300),
                "old_seconds": old_seconds,
                "new_seconds": new_seconds,
                "old_integration_warnings": count_integration_warnings(old_warnings),
                "new_integration_warnings": count_integration_warnings(new_warnings),
            }
        )
    return pd.DataFrame(rows)


def run_single_likelihood(
    galaxy,
    params: HayashiParameters,
    *,
    force_integral_method: str,
    log10_rho0: float,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
) -> tuple[dict[str, float | str], np.ndarray]:
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always", IntegrationWarning)
        t0 = time.perf_counter()
        projector = RZMomentGridWithForceGrid(
            galaxy=galaxy,
            params=params,
            n_r=n_r,
            n_z=n_z,
            n_los=n_los,
            n_force_r=n_force_r,
            n_force_z=n_force_z,
            zmax_factor=20.0,
            los_factor=20.0,
            r_min_pc=1e-3,
            z_min_pc=1e-3,
            grid_padding=1.08,
            force_integral_method=force_integral_method,
        )
        sigma_unit = projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
        sigma_los2 = sigma_unit * 10.0**log10_rho0
        log_likelihood = GaussianVelocityLikelihood(galaxy).log_likelihood(
            sigma_los2,
            systemic_velocity_kms=params.systemic_velocity_kms,
        )
        elapsed = time.perf_counter() - t0
    return (
        {
            "force_integral_method": force_integral_method,
            "single_likelihood_seconds": elapsed,
            "force_grid_build_seconds": getattr(projector, "force_grid_build_time_s", np.nan),
            "log_likelihood": float(log_likelihood),
            "n_integration_warnings": count_integration_warnings(caught_warnings),
        },
        sigma_los2,
    )


def summarize_force_errors(force_detail: pd.DataFrame) -> dict[str, float | str]:
    return {
        "force_integral_method": "force_error_summary",
        "single_likelihood_seconds": np.nan,
        "force_grid_build_seconds": np.nan,
        "log_likelihood": np.nan,
        "n_integration_warnings": int(
            force_detail["new_integration_warnings"].sum() - force_detail["old_integration_warnings"].sum()
        ),
        "force_rel_err_dphi_dr_max": float(force_detail["rel_err_dphi_dr"].max()),
        "force_rel_err_dphi_dr_median": float(force_detail["rel_err_dphi_dr"].median()),
        "force_rel_err_dphi_dz_max": float(force_detail["rel_err_dphi_dz"].max()),
        "force_rel_err_dphi_dz_median": float(force_detail["rel_err_dphi_dz"].median()),
        "force_old_seconds_median": float(force_detail["old_seconds"].median()),
        "force_new_seconds_median": float(force_detail["new_seconds"].median()),
        "force_old_warning_count": int(force_detail["old_integration_warnings"].sum()),
        "force_new_warning_count": int(force_detail["new_integration_warnings"].sum()),
    }


def count_integration_warnings(caught_warnings: list[warnings.WarningMessage]) -> int:
    return sum(issubclass(item.category, IntegrationWarning) for item in caught_warnings)


if __name__ == "__main__":
    main()
