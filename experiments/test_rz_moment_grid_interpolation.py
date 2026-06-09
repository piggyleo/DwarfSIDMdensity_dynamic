#!/usr/bin/env python
"""Experiment: R-z Jeans moment grid + interpolation for Willman 1.

This script is intentionally standalone and does not modify the production
`hayashi_jeans` implementation. It compares the current slow projector against
an experimental projector that precomputes Jeans moments on an `(R, |z|)` grid
and uses interpolation during LOS projection.
"""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre

from hayashi_jeans.data import GalaxyData, load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo, HaloModel
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.model import build_hernquist_projector
from hayashi_jeans.params import HayashiParameters
from hayashi_jeans.projection import AxisymmetricJeansProjector, JeansMoments
from hayashi_jeans.tracer import AxisymmetricPlummerTracer, intrinsic_q_from_projected


@dataclass(frozen=True)
class Timing:
    label: str
    seconds: float


class RZMomentGridProjector:
    def __init__(
        self,
        *,
        galaxy: GalaxyData,
        params: HayashiParameters,
        n_r: int,
        n_z: int,
        n_los: int,
        zmax_factor: float,
        los_factor: float,
        r_min_pc: float,
        z_min_pc: float,
        grid_padding: float,
        radial_derivative: str = "gradient",
        force_integral_method: str = "unit_interval",
        force_quadrature_order: int = 128,
        halo: HaloModel | None = None,
    ) -> None:
        self.galaxy = galaxy
        self.params = params
        self.n_los = n_los
        self.zmax_pc = zmax_factor * galaxy.observables.b_star_pc
        self.los_max_pc = los_factor * galaxy.observables.b_star_pc
        self.inclination_rad = params.inclination_rad
        self.beta_z = params.beta_z
        self.halo = halo or GeneralizedHernquistHalo(
            **params.halo_kwargs(),
            force_integral_method=force_integral_method,
            force_quadrature_order=force_quadrature_order,
        )
        self.tracer = AxisymmetricPlummerTracer(
            b_star_pc=galaxy.observables.b_star_pc,
            q_intrinsic=intrinsic_q_from_projected(galaxy.observables.qprime, params.inclination_rad),
        )
        self.radial_derivative = radial_derivative

        self.los_nodes, self.los_weights = _scaled_legendre_nodes(self.los_max_pc, n_los)
        path_r, path_z = self._all_los_rz(galaxy.x_pc, galaxy.y_pc)
        self.r_path_max_pc = float(np.nanmax(path_r))
        self.z_path_max_pc = float(np.nanmax(np.abs(path_z)))
        r_max = grid_padding * max(self.r_path_max_pc, np.nanmax(np.abs(galaxy.x_pc)), 1.0)
        z_max = grid_padding * max(self.z_path_max_pc, self.zmax_pc, 1.0)
        self.r_grid = _zero_plus_log_grid(r_min_pc, r_max, n_r)
        self.z_grid = _zero_plus_log_grid(z_min_pc, z_max, n_z)
        self._build_grids()

    def _all_los_rz(self, x_pc: np.ndarray, y_pc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x = np.asarray(x_pc, dtype=float)[:, None]
        y = np.asarray(y_pc, dtype=float)[:, None]
        ell = self.los_nodes[None, :]
        i = self.inclination_rad
        r = np.sqrt(x**2 + (y * np.cos(i) + ell * np.sin(i)) ** 2)
        z = y * np.sin(i) - ell * np.cos(i)
        return r, z

    def _build_grids(self) -> None:
        r_mesh, z_mesh = np.meshgrid(self.r_grid, self.z_grid, indexing="ij")
        self.nu_grid = self.tracer.density(r_mesh, z_mesh)

        force_r = np.empty_like(r_mesh)
        force_z = np.empty_like(z_mesh)
        for index in np.ndindex(r_mesh.shape):
            force_r[index], force_z[index] = self.halo.potential_gradients(
                float(r_mesh[index]),
                float(z_mesh[index]),
            )
        self.force_r_grid = force_r
        self.force_z_grid = force_z

        integrand = self.nu_grid * self.force_z_grid
        pressure = np.zeros_like(integrand)
        dz = np.diff(self.z_grid)
        for iz in range(len(self.z_grid) - 2, -1, -1):
            pressure[:, iz] = pressure[:, iz + 1] + 0.5 * (
                integrand[:, iz] + integrand[:, iz + 1]
            ) * dz[iz]
        self.pressure_grid = pressure

        v_z2 = np.divide(pressure, self.nu_grid, out=np.zeros_like(pressure), where=self.nu_grid > 0.0)
        v_r2 = v_z2 / (1.0 - self.beta_z)
        d_pressure_dr = self._radial_pressure_derivative(pressure)
        v_phi2 = (v_z2 + r_mesh / self.nu_grid * d_pressure_dr) / (1.0 - self.beta_z) + r_mesh * force_r
        self.v_z2_grid = np.maximum(v_z2, 0.0)
        self.v_r2_grid = np.maximum(v_r2, 0.0)
        self.v_phi2_grid = np.maximum(v_phi2, 0.0)

        self.nu_interp = self._interp(self.nu_grid)
        self.v_r2_interp = self._interp(self.v_r2_grid)
        self.v_phi2_interp = self._interp(self.v_phi2_grid)
        self.v_z2_interp = self._interp(self.v_z2_grid)

    def _interp(self, values: np.ndarray) -> RegularGridInterpolator:
        return RegularGridInterpolator(
            (self.r_grid, self.z_grid),
            values,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

    def _radial_pressure_derivative(self, pressure: np.ndarray) -> np.ndarray:
        if self.radial_derivative == "gradient":
            return np.gradient(pressure, self.r_grid, axis=0, edge_order=1)
        if self.radial_derivative == "cubic_spline":
            d_pressure_dr = np.empty_like(pressure)
            for iz in range(pressure.shape[1]):
                spline = CubicSpline(self.r_grid, pressure[:, iz], bc_type="natural", extrapolate=True)
                d_pressure_dr[:, iz] = spline.derivative()(self.r_grid)
            return d_pressure_dr
        raise ValueError(f"unknown radial_derivative={self.radial_derivative!r}")

    def jeans_moments(self, r_cyl_pc: float, z_pc: float) -> JeansMoments:
        point = np.array([[float(r_cyl_pc), abs(float(z_pc))]])
        return JeansMoments(
            v_r2=float(self.v_r2_interp(point)[0]),
            v_phi2=float(self.v_phi2_interp(point)[0]),
            v_z2=float(self.v_z2_interp(point)[0]),
        )

    def sigma_los2_many(self, x_pc: np.ndarray, y_pc: np.ndarray) -> np.ndarray:
        return np.array([self.sigma_los2(float(x), float(y)) for x, y in zip(x_pc, y_pc)], dtype=float)

    def sigma_los2(self, x_pc: float, y_pc: float) -> float:
        i = self.inclination_rad
        ell = self.los_nodes
        r = np.sqrt(x_pc**2 + (y_pc * np.cos(i) + ell * np.sin(i)) ** 2)
        z = y_pc * np.sin(i) - ell * np.cos(i)
        points = np.column_stack([r, np.abs(z)])

        nu = self.nu_interp(points)
        v_r2 = self.v_r2_interp(points)
        v_phi2 = self.v_phi2_interp(points)
        v_z2 = self.v_z2_interp(points)
        valid = np.isfinite(nu) & np.isfinite(v_r2) & np.isfinite(v_phi2) & np.isfinite(v_z2)
        if not np.all(valid):
            return np.nan

        with np.errstate(divide="ignore", invalid="ignore"):
            x_over_r2 = np.where(r > 1e-10, x_pc**2 / r**2, 0.0)
        v_star2 = v_phi2 * x_over_r2 + v_r2 * (1.0 - x_over_r2)
        v_los2 = v_star2 * np.sin(i) ** 2 + v_z2 * np.cos(i) ** 2
        denom = np.sum(self.los_weights * nu)
        if denom <= 0.0:
            return np.nan
        numer = np.sum(self.los_weights * nu * np.maximum(v_los2, 0.0))
        return max(float(numer / denom), 1e-10)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default="data/galaxies/27_Willman_1.csv")
    parser.add_argument("--centers-csv", default="data/processed/galaxy_structural_centers.csv")
    parser.add_argument("--out-dir", default="outputs/diagnostics")
    parser.add_argument("--log-path", default="logs/likelihood_acceleration_strategy_log.md")
    parser.add_argument("--n-r", type=int, default=32)
    parser.add_argument("--n-z", type=int, default=64)
    parser.add_argument("--n-los", type=int, default=80)
    parser.add_argument("--n-moment-points", type=int, default=24)
    parser.add_argument("--zmax-factor", type=float, default=20.0)
    parser.add_argument("--los-factor", type=float, default=20.0)
    parser.add_argument("--epsrel", type=float, default=1.5e-2)
    parser.add_argument("--r-min-pc", type=float, default=1.0e-3)
    parser.add_argument("--z-min-pc", type=float, default=1.0e-3)
    parser.add_argument("--grid-padding", type=float, default=1.08)
    parser.add_argument("--radial-derivative", choices=("gradient", "cubic_spline"), default="gradient")
    parser.add_argument("--seed", type=int, default=20260519)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    params = HayashiParameters(
        q_halo=1.0,
        b_halo_pc=10.0**3.2,
        rho0_msun_pc3=10.0 ** -1.961219,
        beta_z=0.3,
        alpha=1.8,
        beta=6.4,
        gamma=1.2,
        inclination_rad=np.deg2rad(75.0),
        systemic_velocity_kms=-13.298645,
    )

    slow_projector = build_hernquist_projector(
        galaxy,
        params,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        epsrel=args.epsrel,
    )

    t0 = time.perf_counter()
    fast_projector = RZMomentGridProjector(
        galaxy=galaxy,
        params=params,
        n_r=args.n_r,
        n_z=args.n_z,
        n_los=args.n_los,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        r_min_pc=args.r_min_pc,
        z_min_pc=args.z_min_pc,
        grid_padding=args.grid_padding,
        radial_derivative=args.radial_derivative,
    )
    fast_build_time = time.perf_counter() - t0
    print(f"fast grid build time: {fast_build_time:.3f}s", flush=True)
    print(
        "grid coverage: "
        f"R_max={fast_projector.r_grid[-1]:.3f} pc, "
        f"|z|_max={fast_projector.z_grid[-1]:.3f} pc, "
        f"LOS path R_max={fast_projector.r_path_max_pc:.3f} pc, "
        f"LOS path |z|_max={fast_projector.z_path_max_pc:.3f} pc",
        flush=True,
    )

    moment_df = compare_moments(slow_projector, fast_projector, args.n_moment_points, args.seed)
    moment_path = out_dir / "rz_grid_willman1_moment_comparison.csv"
    moment_df.to_csv(moment_path, index=False)

    t0 = time.perf_counter()
    slow_sigma = slow_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    slow_sigma_time = time.perf_counter() - t0
    print(f"slow sigma_los2_many time: {slow_sigma_time:.3f}s", flush=True)

    t0 = time.perf_counter()
    fast_sigma = fast_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    fast_sigma_time = time.perf_counter() - t0
    print(f"fast sigma_los2_many time: {fast_sigma_time:.3f}s", flush=True)

    sigma_df = pd.DataFrame(
        {
            "star_id": galaxy.stars["star_id"],
            "x_pc": galaxy.x_pc,
            "y_pc": galaxy.y_pc,
            "slow_sigma_los2": slow_sigma,
            "fast_sigma_los2": fast_sigma,
        }
    )
    sigma_df["abs_error"] = np.abs(sigma_df["fast_sigma_los2"] - sigma_df["slow_sigma_los2"])
    sigma_df["rel_error"] = sigma_df["abs_error"] / np.abs(sigma_df["slow_sigma_los2"])
    sigma_path = out_dir / "rz_grid_willman1_sigma_los2_comparison.csv"
    sigma_df.to_csv(sigma_path, index=False)

    likelihood = GaussianVelocityLikelihood(galaxy)
    slow_ln_l = likelihood.log_likelihood(slow_sigma, systemic_velocity_kms=params.systemic_velocity_kms)
    fast_ln_l = likelihood.log_likelihood(fast_sigma, systemic_velocity_kms=params.systemic_velocity_kms)
    delta_ln_l = fast_ln_l - slow_ln_l

    summary = {
        "n_star": len(galaxy.stars),
        "n_r": args.n_r,
        "n_z": args.n_z,
        "n_los": args.n_los,
        "radial_derivative": args.radial_derivative,
        "grid_r_max_pc": fast_projector.r_grid[-1],
        "grid_z_max_pc": fast_projector.z_grid[-1],
        "los_path_r_max_pc": fast_projector.r_path_max_pc,
        "los_path_z_abs_max_pc": fast_projector.z_path_max_pc,
        "fast_build_time_s": fast_build_time,
        "slow_sigma_time_s": slow_sigma_time,
        "fast_sigma_time_s": fast_sigma_time,
        "fast_total_time_s": fast_build_time + fast_sigma_time,
        "speedup_excluding_build": slow_sigma_time / fast_sigma_time,
        "speedup_including_build": slow_sigma_time / (fast_build_time + fast_sigma_time),
        "sigma_rel_error_median": sigma_df["rel_error"].median(),
        "sigma_rel_error_max": sigma_df["rel_error"].max(),
        "sigma_rel_error_p95": sigma_df["rel_error"].quantile(0.95),
        "moment_rel_error_median_all": moment_df.filter(like="rel_error").stack().median(),
        "moment_rel_error_max_all": moment_df.filter(like="rel_error").stack().max(),
        "slow_ln_l": slow_ln_l,
        "fast_ln_l": fast_ln_l,
        "delta_ln_l": delta_ln_l,
    }
    summary_path = out_dir / "rz_grid_willman1_summary.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    append_log(Path(args.log_path), summary, moment_path, sigma_path, summary_path)
    print(pd.Series(summary).to_string(), flush=True)


def compare_moments(
    slow_projector: AxisymmetricJeansProjector,
    fast_projector: RZMomentGridProjector,
    n_points: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_inner = n_points // 2
    n_outer = n_points - n_inner
    r_values = np.concatenate(
        [
            10.0 ** rng.uniform(np.log10(0.05), np.log10(100.0), n_inner),
            10.0 ** rng.uniform(np.log10(100.0), np.log10(fast_projector.r_grid[-1] * 0.95), n_outer),
        ]
    )
    z_values = np.concatenate(
        [
            10.0 ** rng.uniform(np.log10(0.05), np.log10(100.0), n_inner),
            10.0 ** rng.uniform(np.log10(100.0), np.log10(fast_projector.z_grid[-1] * 0.95), n_outer),
        ]
    )
    rows = []
    for r, z in zip(r_values, z_values):
        slow = slow_projector.jeans_moments(float(r), float(z))
        fast = fast_projector.jeans_moments(float(r), float(z))
        row = {"R_pc": r, "z_pc": z}
        for name in ("v_r2", "v_phi2", "v_z2"):
            slow_value = getattr(slow, name)
            fast_value = getattr(fast, name)
            row[f"slow_{name}"] = slow_value
            row[f"fast_{name}"] = fast_value
            row[f"{name}_rel_error"] = abs(fast_value - slow_value) / max(abs(slow_value), 1e-30)
        rows.append(row)
    return pd.DataFrame(rows)


def append_log(log_path: Path, summary: dict[str, float], moment_path: Path, sigma_path: Path, summary_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("\n\n### Strategy 1 Test Run 001\n\n")
        handle.write("Status: completed experimental standalone test; production code unchanged.\n\n")
        handle.write(f"- moment comparison CSV: `{moment_path}`\n")
        handle.write(f"- sigma_los2 comparison CSV: `{sigma_path}`\n")
        handle.write(f"- summary CSV: `{summary_path}`\n")
        for key, value in summary.items():
            handle.write(f"- {key}: {value}\n")


def _zero_plus_log_grid(min_positive: float, max_value: float, n: int) -> np.ndarray:
    if n < 3:
        raise ValueError("grid must contain at least 3 points")
    positive = np.geomspace(min_positive, max_value, n - 1)
    return np.concatenate([[0.0], positive])


def _scaled_legendre_nodes(max_abs: float, n: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = roots_legendre(n)
    return max_abs * nodes, max_abs * weights


if __name__ == "__main__":
    main()
