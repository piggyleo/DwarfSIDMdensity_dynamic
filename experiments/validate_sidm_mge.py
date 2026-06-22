#!/usr/bin/env python
"""Validate SIDM MGE density, mass, velocity dispersion, and likelihood."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.special import erf

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.halos import SIDMPSIDM25Halo, SpheroidallyStratifiedMGEHalo
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.mge import (
    MGEPhysicalityConfig,
    decompose_density_mge,
    fit_galaxy_mges_for_halo,
    mge_density_relative_errors,
    mge_los_second_moment,
)
from hayashi_jeans.model import build_projector
from hayashi_jeans.params import HayashiParameters
from hayashi_jeans.tracer import AxisymmetricMGETracer


def mge_enclosed_mass(radius: np.ndarray, amplitudes: np.ndarray, sigmas: np.ndarray) -> np.ndarray:
    r = np.asarray(radius)[:, None]
    s = np.asarray(sigmas)[None, :]
    term = (
        np.sqrt(np.pi / 2.0) * s**3 * erf(r / (np.sqrt(2.0) * s))
        - s**2 * r * np.exp(-0.5 * (r / s) ** 2)
    )
    return 4.0 * np.pi * np.sum(np.asarray(amplitudes)[None, :] * term, axis=1)


def relative_error_summary(prefix: str, approximate: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    error = np.abs(np.asarray(approximate) - np.asarray(reference)) / np.maximum(
        np.abs(reference),
        1.0e-300,
    )
    return {
        f"{prefix}_median": float(np.median(error)),
        f"{prefix}_p95": float(np.percentile(error, 95.0)),
        f"{prefix}_max": float(np.max(error)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/27_Willman_1.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--taus", default="0,0.1,0.5,0.9,1.05,1.08")
    parser.add_argument("--strict-taus", default="0.5,1.08")
    parser.add_argument("--output", default=str(PROJECT_ROOT / "outputs/diagnostics/sidm_mge_validation.csv"))
    parser.add_argument("--n-gauss", type=int, default=60)
    parser.add_argument("--strict-epsrel", type=float, default=1.0e-2)
    parser.add_argument("--strict-factor", type=float, default=20.0)
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    taus = [float(value) for value in args.taus.split(",")]
    strict_taus = {float(value) for value in args.strict_taus.split(",") if value}
    radius = np.geomspace(0.1, 5.0e4, 800)
    integration_radius = np.concatenate([[0.0], radius])
    rows = []

    for tau in taus:
        halo = SIDMPSIDM25Halo(q=1.0, rs0_pc=1600.0, rho_s0_msun_pc3=0.1, tau=tau)
        config = MGEPhysicalityConfig(
            n_gauss_halo=args.n_gauss,
            n_gauss_tracer=60,
            decomposition_terms=28,
            n_u=96,
            r_min_pc=1.0e-2,
            r_max_pc=2.0e5,
            halo_decomposition_method="auto",
            real_lstsq_rcond=1.0e-12,
        )
        analytic = decompose_density_mge(
            halo.density_at_ellipsoidal_radius,
            n_gauss=args.n_gauss,
            decomposition_terms=config.decomposition_terms,
            r_min_pc=config.r_min_pc,
            r_max_pc=config.r_max_pc,
        )
        analytic_density_error = mge_density_relative_errors(
            halo.density_at_ellipsoidal_radius,
            analytic,
            radius,
        )
        q_star, tracer_mge, selected_mge = fit_galaxy_mges_for_halo(
            galaxy,
            halo=halo,
            inclination_rad=np.deg2rad(75.0),
            config=config,
        )
        selected_density_error = mge_density_relative_errors(
            halo.density_at_ellipsoidal_radius,
            selected_mge,
            radius,
        )
        true_density = np.asarray(halo.density_at_ellipsoidal_radius(integration_radius))
        with np.errstate(invalid="ignore"):
            mass_integrand = 4.0 * np.pi * integration_radius**2 * true_density
        mass_integrand[0] = 0.0
        true_mass = cumulative_trapezoid(
            mass_integrand,
            integration_radius,
            initial=0.0,
        )[1:]
        fitted_mass = mge_enclosed_mass(
            radius,
            selected_mge.amplitudes,
            selected_mge.sigmas_major_pc,
        )
        mass_error = np.abs(fitted_mass - true_mass) / np.maximum(np.abs(true_mass), 1.0e-300)
        row = {
            "tau": tau,
            "analytic_density_rel_error_p95": float(np.percentile(analytic_density_error, 95.0)),
            "analytic_density_rel_error_max": float(np.max(analytic_density_error)),
            "selected_decomposition_method": selected_mge.decomposition_method,
            "density_rel_error_median": float(np.median(selected_density_error)),
            "density_rel_error_p95": float(np.percentile(selected_density_error, 95.0)),
            "density_rel_error_max": float(np.max(selected_density_error)),
            "mass_rel_error_median": float(np.median(mass_error[20:])),
            "mass_rel_error_p95": float(np.percentile(mass_error[20:], 95.0)),
            "mass_rel_error_max": float(np.max(mass_error[20:])),
        }

        if tau in strict_taus:
            beta_z = 1.0 - 10.0**-0.3
            inclination_rad = np.deg2rad(75.0)
            params = HayashiParameters(
                q_halo=1.0,
                b_halo_pc=1.0,
                rho0_msun_pc3=1.0,
                beta_z=beta_z,
                alpha=1.0,
                beta=4.0,
                gamma=1.0,
                inclination_rad=inclination_rad,
                systemic_velocity_kms=float(np.median(galaxy.velocity_kms)),
            )
            original_t0 = time.perf_counter()
            original_projector = build_projector(
                galaxy,
                params,
                halo=halo,
                zmax_factor=args.strict_factor,
                los_factor=args.strict_factor,
                epsrel=args.strict_epsrel,
            )
            original_sigma = original_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
            row["original_profile_strict_seconds"] = time.perf_counter() - original_t0

            mge_halo = SpheroidallyStratifiedMGEHalo(
                q=halo.q,
                amplitudes_msun_pc3=selected_mge.amplitudes,
                sigmas_major_pc=selected_mge.sigmas_major_pc,
            )
            mge_tracer = AxisymmetricMGETracer(
                q_intrinsic=q_star,
                amplitudes=tracer_mge.amplitudes,
                sigmas_major_pc=tracer_mge.sigmas_major_pc,
            )
            mge_profile_t0 = time.perf_counter()
            mge_profile_projector = build_projector(
                galaxy,
                params,
                halo=mge_halo,
                tracer=mge_tracer,
                zmax_factor=args.strict_factor,
                los_factor=args.strict_factor,
                epsrel=args.strict_epsrel,
            )
            mge_profile_sigma = mge_profile_projector.sigma_los2_many(
                galaxy.x_pc,
                galaxy.y_pc,
            )
            row["mge_profile_strict_seconds"] = time.perf_counter() - mge_profile_t0

            mge_jam_t0 = time.perf_counter()
            mge_jam_sigma = mge_los_second_moment(
                x_pc=galaxy.x_pc,
                y_pc=galaxy.y_pc,
                halo_mge=selected_mge,
                tracer_mge=tracer_mge,
                q_halo=halo.q,
                q_star=q_star,
                beta_z=beta_z,
                inclination_rad=inclination_rad,
                n_u=config.n_u,
            )
            row["mge_jam_seconds"] = time.perf_counter() - mge_jam_t0
            row.update(
                relative_error_summary(
                    "mge_profile_vs_original_sigma_los2_rel_error",
                    mge_profile_sigma,
                    original_sigma,
                )
            )
            row.update(
                relative_error_summary(
                    "mge_jam_vs_mge_profile_sigma_los2_rel_error",
                    mge_jam_sigma,
                    mge_profile_sigma,
                )
            )
            row.update(
                relative_error_summary(
                    "mge_jam_vs_original_sigma_los2_rel_error",
                    mge_jam_sigma,
                    original_sigma,
                )
            )
            likelihood = GaussianVelocityLikelihood(galaxy)
            systemic = params.systemic_velocity_kms
            original_log_likelihood = likelihood.log_likelihood(
                original_sigma,
                systemic_velocity_kms=systemic,
            )
            mge_profile_log_likelihood = likelihood.log_likelihood(
                mge_profile_sigma,
                systemic_velocity_kms=systemic,
            )
            mge_jam_log_likelihood = likelihood.log_likelihood(
                mge_jam_sigma,
                systemic_velocity_kms=systemic,
            )
            row["original_profile_strict_log_likelihood"] = original_log_likelihood
            row["mge_profile_strict_log_likelihood"] = mge_profile_log_likelihood
            row["mge_jam_log_likelihood"] = mge_jam_log_likelihood
            row["delta_log_likelihood_mge_profile_minus_original"] = (
                mge_profile_log_likelihood - original_log_likelihood
            )
            row["delta_log_likelihood_mge_jam_minus_mge_profile"] = (
                mge_jam_log_likelihood - mge_profile_log_likelihood
            )
            row["delta_log_likelihood_mge_jam_minus_original"] = (
                mge_jam_log_likelihood - original_log_likelihood
            )

        rows.append(row)
        print(pd.Series(row).to_string(), flush=True)

    result = pd.DataFrame(rows)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False)
    print(f"wrote {output}")

    if result["density_rel_error_p95"].max() >= 0.05:
        raise SystemExit("SIDM MGE density validation failed")
    if result["mass_rel_error_p95"].max() >= 0.05:
        raise SystemExit("SIDM MGE mass validation failed")


if __name__ == "__main__":
    main()
