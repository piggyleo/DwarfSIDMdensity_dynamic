#!/usr/bin/env python
"""Decompose Eridanus II fast-vs-strict likelihood errors by control modes."""

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

from experiments.test_halo_force_grid_interpolation import RZMomentGridWithForceGrid
from experiments.test_rz_moment_grid_interpolation import RZMomentGridProjector
from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.params import HayashiParameters
from hayashi_jeans.projection import AxisymmetricJeansProjector
from scripts.run_willman1_block_mh_fast import beta_z_from_q


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_vs_strict_selected_samples.csv"))
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--strict-epsrel", type=float, default=1.5e-2)
    parser.add_argument("--detail-output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_error_source_detail.csv"))
    parser.add_argument("--summary-output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_error_source_summary.csv"))
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    selected = pd.read_csv(args.selected)

    rows = []
    for row in selected.itertuples(index=False):
        row_dict = row._asdict()
        print(f"sample {int(row_dict['selection_id'])}/{len(selected)} source_row={int(row_dict['source_row'])}", flush=True)
        result = diagnose_sample(galaxy, row_dict, strict_epsrel=args.strict_epsrel)
        rows.append(result)
        print(
            f"  total dlogL={result['fast_vs_strict_delta_log_likelihood']:.6g}; "
            f"RZ+other={result['rz96_vs_strict_delta_log_likelihood']:.6g}; "
            f"force-grid={result['fast_vs_rz96_delta_log_likelihood']:.6g}; "
            f"GL128-256={result['fast_gl128_vs_gl256_delta_log_likelihood']:.6g}",
            flush=True,
        )

    detail = pd.DataFrame(rows)
    detail_path = Path(args.detail_output)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(detail_path, index=False)

    summary = summarize(detail)
    summary_path = Path(args.summary_output)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_path, index=False)
    print(f"wrote {detail_path}")
    print(detail.to_string(index=False))
    print(f"wrote {summary_path}")
    print(summary.to_string(index=False))


def diagnose_sample(galaxy, row: dict[str, float], *, strict_epsrel: float) -> dict[str, float | int]:
    params = params_from_row(row)
    log10_rho0 = float(row["log10_rho0_msun_pc3"])
    systemic = float(row["systemic_velocity_kms"])

    modes = {
        "strict": lambda: sigma_strict(galaxy, params, epsrel=strict_epsrel, force_quadrature_order=128),
        "rz96_direct_force": lambda: sigma_rz_direct(galaxy, params, n_r=96, n_z=192, n_los=160, force_quadrature_order=128),
        "rz96_direct_force_gl256": lambda: sigma_rz_direct(galaxy, params, n_r=96, n_z=192, n_los=160, force_quadrature_order=256),
        "rz128_direct_force": lambda: sigma_rz_direct(galaxy, params, n_r=128, n_z=256, n_los=200, force_quadrature_order=128),
        "fast_force48": lambda: sigma_fast(galaxy, params, n_r=96, n_z=192, n_los=160, n_force_r=48, n_force_z=96, force_quadrature_order=128),
        "fast_force48_gl256": lambda: sigma_fast(galaxy, params, n_r=96, n_z=192, n_los=160, n_force_r=48, n_force_z=96, force_quadrature_order=256),
        "fast_force64_moment128": lambda: sigma_fast(galaxy, params, n_r=128, n_z=256, n_los=200, n_force_r=64, n_force_z=128, force_quadrature_order=128),
    }

    evaluated = {}
    for name, function in modes.items():
        t0 = time.perf_counter()
        sigma_unit = function()
        seconds = time.perf_counter() - t0
        sigma2 = sigma_unit * 10.0**log10_rho0
        log_likelihood = GaussianVelocityLikelihood(galaxy).log_likelihood(sigma2, systemic_velocity_kms=systemic)
        evaluated[name] = {"sigma2": sigma2, "log_likelihood": log_likelihood, "seconds": seconds}
        print(f"    {name}: logL={log_likelihood:.6f}, seconds={seconds:.3f}", flush=True)

    result = {
        "selection_id": int(row["selection_id"]),
        "source_row": int(row["source_row"]),
        "step": int(row["step"]),
        "log_probability_chain": float(row["log_probability"]),
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
    for name, values in evaluated.items():
        result[f"{name}_log_likelihood"] = float(values["log_likelihood"])
        result[f"{name}_seconds"] = float(values["seconds"])

    add_comparison(result, evaluated, "fast_force48", "strict", "fast_vs_strict")
    add_comparison(result, evaluated, "rz96_direct_force", "strict", "rz96_vs_strict")
    add_comparison(result, evaluated, "fast_force48", "rz96_direct_force", "fast_vs_rz96")
    add_comparison(result, evaluated, "fast_force48", "fast_force48_gl256", "fast_gl128_vs_gl256")
    add_comparison(result, evaluated, "rz96_direct_force", "rz96_direct_force_gl256", "rz_gl128_vs_gl256")
    add_comparison(result, evaluated, "rz96_direct_force", "rz128_direct_force", "rz96_vs_rz128")
    add_comparison(result, evaluated, "rz128_direct_force", "strict", "rz128_vs_strict")
    add_comparison(result, evaluated, "fast_force48", "fast_force64_moment128", "fast48_vs_fast64_moment128")
    return result


def params_from_row(row: dict[str, float]) -> HayashiParameters:
    return HayashiParameters(
        q_halo=float(row["q_halo"]),
        b_halo_pc=10.0 ** float(row["log10_b_halo_pc"]),
        rho0_msun_pc3=1.0,
        beta_z=beta_z_from_q(float(row["minus_log10_one_minus_beta_z"])),
        alpha=float(row["alpha"]),
        beta=float(row["beta"]),
        gamma=float(row["gamma"]),
        inclination_rad=np.deg2rad(float(row["i_deg"])),
        systemic_velocity_kms=0.0,
    )


def sigma_strict(galaxy, params: HayashiParameters, *, epsrel: float, force_quadrature_order: int) -> np.ndarray:
    halo = GeneralizedHernquistHalo(
        **params.halo_kwargs(),
        force_integral_method="unit_interval",
        force_quadrature_order=force_quadrature_order,
    )
    projector = AxisymmetricJeansProjector(
        b_star_pc=galaxy.observables.b_star_pc,
        qprime=galaxy.observables.qprime,
        inclination_rad=params.inclination_rad,
        beta_z=params.beta_z,
        halo=halo,
        zmax_factor=20.0,
        los_factor=20.0,
        epsrel=epsrel,
    )
    return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)


def sigma_rz_direct(
    galaxy,
    params: HayashiParameters,
    *,
    n_r: int,
    n_z: int,
    n_los: int,
    force_quadrature_order: int,
) -> np.ndarray:
    projector = RZMomentGridProjector(
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
        force_integral_method="unit_interval",
        force_quadrature_order=force_quadrature_order,
    )
    return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)


def sigma_fast(
    galaxy,
    params: HayashiParameters,
    *,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    force_quadrature_order: int,
) -> np.ndarray:
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
        force_integral_method="unit_interval",
        force_quadrature_order=force_quadrature_order,
    )
    return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)


def add_comparison(result: dict[str, float | int], evaluated: dict[str, dict[str, np.ndarray | float]], left: str, right: str, label: str) -> None:
    left_sigma2 = np.asarray(evaluated[left]["sigma2"], dtype=float)
    right_sigma2 = np.asarray(evaluated[right]["sigma2"], dtype=float)
    sigma2_abs = np.abs(left_sigma2 - right_sigma2)
    sigma2_rel = sigma2_abs / np.maximum(np.abs(right_sigma2), 1.0e-300)
    sigma_abs = np.abs(np.sqrt(np.maximum(left_sigma2, 0.0)) - np.sqrt(np.maximum(right_sigma2, 0.0)))
    sigma_rel = sigma_abs / np.maximum(np.sqrt(np.maximum(right_sigma2, 0.0)), 1.0e-300)
    delta_log_likelihood = float(evaluated[left]["log_likelihood"] - evaluated[right]["log_likelihood"])
    result[f"{label}_delta_log_likelihood"] = delta_log_likelihood
    result[f"{label}_abs_delta_log_likelihood"] = abs(delta_log_likelihood)
    result[f"{label}_sigma_los2_rel_err_max"] = float(np.nanmax(sigma2_rel))
    result[f"{label}_sigma_los2_rel_err_p95"] = float(np.nanpercentile(sigma2_rel, 95.0))
    result[f"{label}_sigma_los2_rel_err_median"] = float(np.nanmedian(sigma2_rel))
    result[f"{label}_sigma_los_rel_err_max"] = float(np.nanmax(sigma_rel))
    result[f"{label}_sigma_los_rel_err_p95"] = float(np.nanpercentile(sigma_rel, 95.0))
    result[f"{label}_sigma_los_rel_err_median"] = float(np.nanmedian(sigma_rel))


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    labels = [
        ("fast_vs_strict", "total fast error"),
        ("rz96_vs_strict", "RZMoment+fixed LOS+pressure derivative vs strict"),
        ("fast_vs_rz96", "ForceGrid contribution at 96x192 moments"),
        ("fast_gl128_vs_gl256", "Gauss-Legendre order convergence in fast ForceGrid pipeline"),
        ("rz_gl128_vs_gl256", "Gauss-Legendre order convergence in RZ direct-force pipeline"),
        ("rz96_vs_rz128", "RZMoment resolution contribution"),
        ("rz128_vs_strict", "higher-resolution RZ residual vs strict"),
        ("fast48_vs_fast64_moment128", "fast-grid convergence contribution"),
    ]
    rows = []
    for label, description in labels:
        row = {"label": label, "description": description, "n_samples": len(detail)}
        for metric in [
            "delta_log_likelihood",
            "abs_delta_log_likelihood",
            "sigma_los2_rel_err_max",
            "sigma_los2_rel_err_p95",
            "sigma_los2_rel_err_median",
            "sigma_los_rel_err_max",
            "sigma_los_rel_err_p95",
            "sigma_los_rel_err_median",
        ]:
            column = f"{label}_{metric}"
            row[f"{metric}_min"] = float(detail[column].min())
            row[f"{metric}_median"] = float(detail[column].median())
            row[f"{metric}_max"] = float(detail[column].max())
        rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
