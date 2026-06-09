#!/usr/bin/env python
"""Benchmark direct Jeans likelihood with all integrals using Gauss-Legendre.

This is an isolated experiment. It does not use HaloForceGrid or RZMomentGrid,
and it does not modify the production projector. Compared with
AxisymmetricJeansProjector, the LOS and vertical-pressure integrals are
evaluated by fixed Gauss-Legendre quadrature and vectorized within each star.
The halo force integral uses the same compact unit-interval Gauss-Legendre
formula as GeneralizedHernquistHalo.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import roots_legendre

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hayashi_jeans.constants import G_PC_MSUN_KMS2
from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.params import HayashiParameters
from hayashi_jeans.tracer import AxisymmetricPlummerTracer, intrinsic_q_from_projected
from scripts.run_willman1_block_mh_fast import beta_z_from_q


@dataclass(frozen=True)
class DirectGLSettings:
    n_los: int
    n_pressure: int
    n_force: int
    zmax_factor: float = 20.0
    los_factor: float = 20.0


class DirectGaussLegendreJeansProjector:
    """Direct, no-grid projector with fixed Gauss-Legendre quadrature."""

    def __init__(self, *, galaxy, params: HayashiParameters, settings: DirectGLSettings) -> None:
        self.galaxy = galaxy
        self.params = params
        self.settings = settings
        self.b_star_pc = float(galaxy.observables.b_star_pc)
        self.qprime = float(galaxy.observables.qprime)
        self.inclination_rad = float(params.inclination_rad)
        self.beta_z = float(params.beta_z)
        self.zmax_pc = settings.zmax_factor * self.b_star_pc
        self.los_max_pc = settings.los_factor * self.b_star_pc
        self.tracer = AxisymmetricPlummerTracer(
            b_star_pc=self.b_star_pc,
            q_intrinsic=intrinsic_q_from_projected(self.qprime, self.inclination_rad),
        )

        self.los_nodes, self.los_weights = scaled_legendre_nodes(-self.los_max_pc, self.los_max_pc, settings.n_los)
        self.unit_nodes, self.unit_weights = scaled_legendre_nodes(0.0, 1.0, settings.n_pressure)
        self.force_nodes, self.force_weights = scaled_legendre_nodes(0.0, 1.0, settings.n_force)

    def sigma_los2_many(self, x_pc: np.ndarray, y_pc: np.ndarray) -> np.ndarray:
        return np.array([self.sigma_los2(float(x), float(y)) for x, y in zip(x_pc, y_pc)], dtype=float)

    def sigma_los2(self, x_pc: float, y_pc: float) -> float:
        i = self.inclination_rad
        ell = self.los_nodes
        r = np.sqrt(np.maximum(x_pc**2 + (y_pc * np.cos(i) + ell * np.sin(i)) ** 2, 0.0))
        z = y_pc * np.sin(i) - ell * np.cos(i)
        nu = self.tracer.density(r, z)
        denom = float(np.sum(self.los_weights * nu))
        if denom <= 0.0:
            return np.nan

        v_r2, v_phi2, v_z2 = self.jeans_moments_many(r, z)
        with np.errstate(divide="ignore", invalid="ignore"):
            x2_over_r2 = np.where(r > 1.0e-10, x_pc**2 / r**2, 0.0)
        v_star2 = np.where(r > 1.0e-10, v_phi2 * x2_over_r2 + v_r2 * (1.0 - x2_over_r2), v_r2)
        v_los2 = v_star2 * np.sin(i) ** 2 + v_z2 * np.cos(i) ** 2
        numer = float(np.sum(self.los_weights * nu * np.maximum(v_los2, 0.0)))
        return max(numer / denom, 1.0e-10)

    def jeans_moments_many(self, r_cyl_pc: np.ndarray, z_pc: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        nu = np.asarray(self.tracer.density(r_cyl_pc, z_pc), dtype=float)
        positive = nu > 0.0
        pressure = self.vertical_pressure_many(r_cyl_pc, z_pc)
        v_z2 = np.zeros_like(nu, dtype=float)
        v_z2[positive] = pressure[positive] / nu[positive]
        v_r2 = v_z2 / (1.0 - self.beta_z)

        d_pressure_dr = self.d_vertical_pressure_dr_many(r_cyl_pc, z_pc)
        dphi_dr, _ = self.halo_gradients(r_cyl_pc, z_pc)
        v_phi2 = np.zeros_like(nu, dtype=float)
        v_phi2[positive] = (
            (v_z2[positive] + r_cyl_pc[positive] / nu[positive] * d_pressure_dr[positive])
            / (1.0 - self.beta_z)
            + r_cyl_pc[positive] * dphi_dr[positive]
        )
        return np.maximum(v_r2, 0.0), np.maximum(v_phi2, 0.0), np.maximum(v_z2, 0.0)

    def vertical_pressure_many(self, r_cyl_pc: np.ndarray, z_pc: np.ndarray) -> np.ndarray:
        r = np.asarray(r_cyl_pc, dtype=float)
        z = np.asarray(z_pc, dtype=float)
        sign = np.where(z >= 0.0, 1.0, -1.0)
        z_abs = np.abs(z)
        width = self.zmax_pc - z_abs
        zp_abs = z_abs[:, None] + width[:, None] * self.unit_nodes[None, :]
        z_eval = sign[:, None] * zp_abs
        r_eval = r[:, None] + np.zeros_like(z_eval)
        _, dphi_dz = self.halo_gradients(r_eval, z_eval)
        integrand = self.tracer.density(r_eval, z_eval) * sign[:, None] * dphi_dz
        return width * np.sum(self.unit_weights[None, :] * integrand, axis=1)

    def d_vertical_pressure_dr_many(self, r_cyl_pc: np.ndarray, z_pc: np.ndarray) -> np.ndarray:
        r = np.asarray(r_cyl_pc, dtype=float)
        z = np.asarray(z_pc, dtype=float)
        step = np.maximum.reduce([
            np.full_like(r, 0.01 * self.b_star_pc, dtype=float),
            0.01 * np.abs(r),
            np.full_like(r, 1.0e-2, dtype=float),
        ])
        central = r > step
        derivative = np.empty_like(r, dtype=float)
        if np.any(central):
            p_plus = self.vertical_pressure_many(r[central] + step[central], z[central])
            p_minus = self.vertical_pressure_many(r[central] - step[central], z[central])
            derivative[central] = (p_plus - p_minus) / (2.0 * step[central])
        if np.any(~central):
            p0 = self.vertical_pressure_many(r[~central], z[~central])
            p_plus = self.vertical_pressure_many(r[~central] + step[~central], z[~central])
            derivative[~central] = (p_plus - p0) / step[~central]
        return derivative

    def halo_gradients(self, r_cyl_pc: np.ndarray, z_pc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        r = np.abs(np.asarray(r_cyl_pc, dtype=float))
        z = np.asarray(z_pc, dtype=float)
        q = self.params.q_halo
        b_pc = self.params.b_halo_pc
        alpha = self.params.alpha
        beta = self.params.beta
        gamma = self.params.gamma
        rho0 = self.params.rho0_msun_pc3

        s = self.force_nodes
        weights = self.force_weights
        s2 = s * s
        a_s = 1.0 + (q**2 - 1.0) * s2

        m2 = s2 * (r[..., None] ** 2 + z[..., None] ** 2 / a_s)
        m = np.sqrt(np.maximum(m2, 0.0))
        x = np.maximum(m / b_pc, 1.0e-12)
        rho = rho0 * x ** (-gamma) * (1.0 + x**alpha) ** (-(beta - gamma) / alpha)

        int_r = np.sum(weights * 2.0 * rho * s2 / np.sqrt(a_s), axis=-1)
        int_z = np.sum(weights * 2.0 * rho * s2 / (a_s**1.5), axis=-1)
        prefactor = 2.0 * np.pi * G_PC_MSUN_KMS2 * q
        return prefactor * r * int_r, prefactor * z * int_z


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_vs_strict_selected_samples.csv"))
    parser.add_argument("--strict-detail", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_error_source_detail.csv"))
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--n-los", type=int, default=64)
    parser.add_argument("--n-pressure", type=int, default=64)
    parser.add_argument("--n-force", type=int, default=128)
    parser.add_argument("--sample-limit", type=int, default=None)
    parser.add_argument("--selection-id", type=int, action="append", default=None)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    selected = pd.read_csv(args.selected)
    if args.selection_id is not None:
        selected = selected[selected["selection_id"].isin(args.selection_id)].copy()
    if args.sample_limit is not None:
        selected = selected.head(args.sample_limit).copy()
    strict = pd.read_csv(args.strict_detail).set_index("selection_id")
    likelihood = GaussianVelocityLikelihood(galaxy)
    settings = DirectGLSettings(n_los=args.n_los, n_pressure=args.n_pressure, n_force=args.n_force)

    rows = []
    for row in selected.itertuples(index=False):
        row_dict = row._asdict()
        params = params_from_row(row_dict)
        projector = DirectGaussLegendreJeansProjector(galaxy=galaxy, params=params, settings=settings)
        t0 = time.perf_counter()
        sigma_unit = projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
        seconds = time.perf_counter() - t0
        sigma2 = sigma_unit * 10.0 ** float(row_dict["log10_rho0_msun_pc3"])
        log_likelihood = likelihood.log_likelihood(
            sigma2,
            systemic_velocity_kms=float(row_dict["systemic_velocity_kms"]),
        )

        strict_row = strict.loc[int(row_dict["selection_id"])]
        strict_log_likelihood = float(strict_row["strict_log_likelihood"])
        strict_seconds = float(strict_row["strict_seconds"])
        rows.append(
            {
                "selection_id": int(row_dict["selection_id"]),
                "source_row": int(row_dict["source_row"]),
                "n_los": args.n_los,
                "n_pressure": args.n_pressure,
                "n_force": args.n_force,
                "all_gl_log_likelihood": float(log_likelihood),
                "strict_log_likelihood": strict_log_likelihood,
                "delta_log_likelihood_all_gl_minus_strict": float(log_likelihood - strict_log_likelihood),
                "abs_delta_log_likelihood": float(abs(log_likelihood - strict_log_likelihood)),
                "all_gl_seconds": float(seconds),
                "strict_seconds": strict_seconds,
                "speedup_vs_strict": float(strict_seconds / seconds),
            }
        )
        print(
            f"sample {int(row_dict['selection_id'])}: "
            f"logL={log_likelihood:.6f}, dlogL={log_likelihood - strict_log_likelihood:.6g}, "
            f"seconds={seconds:.3f}, speedup={strict_seconds / seconds:.2f}x",
            flush=True,
        )

    detail = pd.DataFrame(rows)
    if args.output is None:
        output = PROJECT_ROOT / (
            f"outputs/diagnostics/eridanus_ii_all_gl_direct_likelihood_timing_"
            f"los{args.n_los}_pressure{args.n_pressure}_force{args.n_force}.csv"
        )
    else:
        output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(output, index=False)
    print(f"wrote {output}")
    print(detail.to_string(index=False))
    print("\nsummary")
    summary_cols = ["abs_delta_log_likelihood", "all_gl_seconds", "speedup_vs_strict"]
    print(detail[summary_cols].describe(percentiles=[0.25, 0.5, 0.75]).to_string())


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


def scaled_legendre_nodes(a: float, b: float, order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = roots_legendre(order)
    half_width = 0.5 * (b - a)
    center = 0.5 * (a + b)
    return center + half_width * nodes, half_width * weights


if __name__ == "__main__":
    main()
