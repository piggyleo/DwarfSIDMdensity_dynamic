#!/usr/bin/env python
"""Experiment: halo force grid + interpolation on top of the R-z moment grid.

This standalone test does not modify production code. It asks whether halo
force evaluations can be moved to a coarser `(R, |z|)` grid and interpolated
onto the Jeans moment grid without materially changing `sigma_los^2` or the
velocity likelihood.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline, RegularGridInterpolator
from scipy.special import roots_legendre

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from experiments.test_rz_moment_grid_interpolation import RZMomentGridProjector
from hayashi_jeans.data import GalaxyData, load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo, HaloModel
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.params import HayashiParameters
from hayashi_jeans.projection import JeansMoments
from hayashi_jeans.tracer import AxisymmetricPlummerTracer, intrinsic_q_from_projected


@dataclass(frozen=True)
class ForceGridSpec:
    n_r: int
    n_z: int

    @property
    def label(self) -> str:
        return f"force{self.n_r}x{self.n_z}"


class RZMomentGridWithForceGrid:
    def __init__(
        self,
        *,
        galaxy: GalaxyData,
        params: HayashiParameters,
        n_r: int,
        n_z: int,
        n_los: int,
        n_force_r: int,
        n_force_z: int,
        zmax_factor: float,
        los_factor: float,
        r_min_pc: float,
        z_min_pc: float,
        grid_padding: float,
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

        self.los_nodes, self.los_weights = _scaled_legendre_nodes(self.los_max_pc, n_los)
        path_r, path_z = self._all_los_rz(galaxy.x_pc, galaxy.y_pc)
        self.r_path_max_pc = float(np.nanmax(path_r))
        self.z_path_max_pc = float(np.nanmax(np.abs(path_z)))
        r_max = grid_padding * max(self.r_path_max_pc, np.nanmax(np.abs(galaxy.x_pc)), 1.0)
        z_max = grid_padding * max(self.z_path_max_pc, self.zmax_pc, 1.0)

        self.r_grid = _zero_plus_log_grid(r_min_pc, r_max, n_r)
        self.z_grid = _zero_plus_log_grid(z_min_pc, z_max, n_z)
        self.force_r_grid_axis = _zero_plus_log_grid(r_min_pc, r_max, n_force_r)
        self.force_z_grid_axis = _zero_plus_log_grid(z_min_pc, z_max, n_force_z)
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
        force_t0 = time.perf_counter()
        force_r_mesh, force_z_mesh = np.meshgrid(
            self.force_r_grid_axis,
            self.force_z_grid_axis,
            indexing="ij",
        )
        force_r_coarse = np.empty_like(force_r_mesh)
        force_z_coarse = np.empty_like(force_z_mesh)
        for index in np.ndindex(force_r_mesh.shape):
            force_r_coarse[index], force_z_coarse[index] = self.halo.potential_gradients(
                float(force_r_mesh[index]),
                float(force_z_mesh[index]),
            )
        self.force_grid_build_time_s = time.perf_counter() - force_t0

        self.force_r_interp = self._force_interp(force_r_coarse)
        self.force_z_interp = self._force_interp(force_z_coarse)

        moment_t0 = time.perf_counter()
        r_mesh, z_mesh = np.meshgrid(self.r_grid, self.z_grid, indexing="ij")
        self.nu_grid = self.tracer.density(r_mesh, z_mesh)
        moment_points = np.column_stack([r_mesh.ravel(), z_mesh.ravel()])
        self.force_r_grid = self.force_r_interp(moment_points).reshape(r_mesh.shape)
        self.force_z_grid = self.force_z_interp(moment_points).reshape(z_mesh.shape)

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
        v_phi2 = (v_z2 + r_mesh / self.nu_grid * d_pressure_dr) / (1.0 - self.beta_z) + r_mesh * self.force_r_grid
        self.v_z2_grid = np.maximum(v_z2, 0.0)
        self.v_r2_grid = np.maximum(v_r2, 0.0)
        self.v_phi2_grid = np.maximum(v_phi2, 0.0)

        self.nu_interp = self._moment_interp(self.nu_grid)
        self.v_r2_interp = self._moment_interp(self.v_r2_grid)
        self.v_phi2_interp = self._moment_interp(self.v_phi2_grid)
        self.v_z2_interp = self._moment_interp(self.v_z2_grid)
        self.moment_grid_build_time_s = time.perf_counter() - moment_t0

    def _force_interp(self, values: np.ndarray) -> RegularGridInterpolator:
        return RegularGridInterpolator(
            (self.force_r_grid_axis, self.force_z_grid_axis),
            values,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

    def _moment_interp(self, values: np.ndarray) -> RegularGridInterpolator:
        return RegularGridInterpolator(
            (self.r_grid, self.z_grid),
            values,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )

    def _radial_pressure_derivative(self, pressure: np.ndarray) -> np.ndarray:
        d_pressure_dr = np.empty_like(pressure)
        for iz in range(pressure.shape[1]):
            spline = CubicSpline(self.r_grid, pressure[:, iz], bc_type="natural", extrapolate=True)
            d_pressure_dr[:, iz] = spline.derivative()(self.r_grid)
        return d_pressure_dr

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
    parser.add_argument("--n-r", type=int, default=96)
    parser.add_argument("--n-z", type=int, default=192)
    parser.add_argument("--n-los", type=int, default=160)
    parser.add_argument("--force-grids", default="32x64,48x96,64x128,96x192")
    parser.add_argument("--n-force-points", type=int, default=200)
    parser.add_argument("--zmax-factor", type=float, default=20.0)
    parser.add_argument("--los-factor", type=float, default=20.0)
    parser.add_argument("--r-min-pc", type=float, default=1.0e-3)
    parser.add_argument("--z-min-pc", type=float, default=1.0e-3)
    parser.add_argument("--grid-padding", type=float, default=1.08)
    parser.add_argument("--seed", type=int, default=20260520)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    params = make_params()
    force_specs = parse_force_grids(args.force_grids)

    print("Building direct-force R-z reference projector...", flush=True)
    t0 = time.perf_counter()
    reference_projector = RZMomentGridProjector(
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
        radial_derivative="cubic_spline",
    )
    reference_build_time = time.perf_counter() - t0
    t0 = time.perf_counter()
    reference_sigma = reference_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    reference_sigma_time = time.perf_counter() - t0
    likelihood = GaussianVelocityLikelihood(galaxy)
    reference_ln_l = likelihood.log_likelihood(reference_sigma, systemic_velocity_kms=params.systemic_velocity_kms)
    print(
        f"Reference build={reference_build_time:.3f}s, "
        f"sigma={reference_sigma_time:.3f}s, lnL={reference_ln_l:.6f}",
        flush=True,
    )

    detail_frames = []
    force_frames = []
    summary_rows = []
    for spec in force_specs:
        print(f"Testing {spec.label}...", flush=True)
        t0 = time.perf_counter()
        projector = RZMomentGridWithForceGrid(
            galaxy=galaxy,
            params=params,
            n_r=args.n_r,
            n_z=args.n_z,
            n_los=args.n_los,
            n_force_r=spec.n_r,
            n_force_z=spec.n_z,
            zmax_factor=args.zmax_factor,
            los_factor=args.los_factor,
            r_min_pc=args.r_min_pc,
            z_min_pc=args.z_min_pc,
            grid_padding=args.grid_padding,
        )
        force_grid_total_build_time = time.perf_counter() - t0
        t0 = time.perf_counter()
        sigma = projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
        sigma_time = time.perf_counter() - t0
        ln_l = likelihood.log_likelihood(sigma, systemic_velocity_kms=params.systemic_velocity_kms)

        sigma_df = make_sigma_detail(galaxy, spec.label, reference_sigma, sigma)
        detail_frames.append(sigma_df)
        force_df = compare_forces(projector, spec.label, args.n_force_points, args.seed)
        force_frames.append(force_df)

        sigma_rel = sigma_df["sigma_rel_error"]
        force_vec_rel = force_df["force_vector_rel_error"]
        summary = {
            "label": spec.label,
            "n_star": len(galaxy.stars),
            "n_r": args.n_r,
            "n_z": args.n_z,
            "n_los": args.n_los,
            "n_force_r": spec.n_r,
            "n_force_z": spec.n_z,
            "reference_build_time_s": reference_build_time,
            "reference_sigma_time_s": reference_sigma_time,
            "force_grid_build_time_s": projector.force_grid_build_time_s,
            "moment_grid_build_time_s": projector.moment_grid_build_time_s,
            "force_grid_total_build_time_s": force_grid_total_build_time,
            "force_grid_sigma_time_s": sigma_time,
            "build_speedup_vs_direct_force": reference_build_time / force_grid_total_build_time,
            "total_speedup_vs_direct_force": (reference_build_time + reference_sigma_time) / (force_grid_total_build_time + sigma_time),
            "force_vector_rel_error_median": force_vec_rel.median(),
            "force_vector_rel_error_max": force_vec_rel.max(),
            "force_vector_rel_error_p95": force_vec_rel.quantile(0.95),
            "sigma_rel_error_median": sigma_rel.median(),
            "sigma_rel_error_max": sigma_rel.max(),
            "sigma_rel_error_p95": sigma_rel.quantile(0.95),
            "reference_ln_l": reference_ln_l,
            "force_grid_ln_l": ln_l,
            "delta_ln_l": ln_l - reference_ln_l,
        }
        summary_rows.append(summary)
        print(pd.Series(summary).to_string(), flush=True)

    detail = pd.concat(detail_frames, ignore_index=True)
    force_detail = pd.concat(force_frames, ignore_index=True)
    summary = pd.DataFrame(summary_rows)

    detail_path = out_dir / "halo_force_grid_willman1_sigma_detail.csv"
    force_path = out_dir / "halo_force_grid_willman1_force_detail.csv"
    summary_path = out_dir / "halo_force_grid_willman1_summary.csv"
    detail.to_csv(detail_path, index=False)
    force_detail.to_csv(force_path, index=False)
    summary.to_csv(summary_path, index=False)
    append_log(Path(args.log_path), summary, detail_path, force_path, summary_path)
    print("\nSummary:")
    print(summary.to_string(index=False), flush=True)


def make_params() -> HayashiParameters:
    return HayashiParameters(
        q_halo=1.0,
        b_halo_pc=10.0**3.2,
        rho0_msun_pc3=10.0**-1.961219,
        beta_z=0.3,
        alpha=1.8,
        beta=6.4,
        gamma=1.2,
        inclination_rad=np.deg2rad(75.0),
        systemic_velocity_kms=-13.298645,
    )


def parse_force_grids(text: str) -> list[ForceGridSpec]:
    specs = []
    for item in text.split(","):
        left, right = item.lower().split("x", maxsplit=1)
        specs.append(ForceGridSpec(int(left), int(right)))
    return specs


def make_sigma_detail(
    galaxy: GalaxyData,
    label: str,
    reference_sigma: np.ndarray,
    force_grid_sigma: np.ndarray,
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "label": label,
            "star_id": galaxy.stars["star_id"],
            "x_pc": galaxy.x_pc,
            "y_pc": galaxy.y_pc,
            "reference_sigma_los2": reference_sigma,
            "force_grid_sigma_los2": force_grid_sigma,
        }
    )
    frame["sigma_abs_error"] = np.abs(frame["force_grid_sigma_los2"] - frame["reference_sigma_los2"])
    frame["sigma_rel_error"] = frame["sigma_abs_error"] / np.maximum(np.abs(frame["reference_sigma_los2"]), 1e-30)
    return frame


def compare_forces(
    projector: RZMomentGridWithForceGrid,
    label: str,
    n_points: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    n_inner = n_points // 2
    n_outer = n_points - n_inner
    r_values = np.concatenate(
        [
            10.0 ** rng.uniform(np.log10(0.01), np.log10(100.0), n_inner),
            10.0 ** rng.uniform(np.log10(100.0), np.log10(projector.r_grid[-1] * 0.95), n_outer),
        ]
    )
    z_values = np.concatenate(
        [
            10.0 ** rng.uniform(np.log10(0.01), np.log10(100.0), n_inner),
            10.0 ** rng.uniform(np.log10(100.0), np.log10(projector.z_grid[-1] * 0.95), n_outer),
        ]
    )
    rows = []
    points = np.column_stack([r_values, z_values])
    interp_force_r = projector.force_r_interp(points)
    interp_force_z = projector.force_z_interp(points)
    for r, z, fast_r, fast_z in zip(r_values, z_values, interp_force_r, interp_force_z):
        exact_r, exact_z = projector.halo.potential_gradients(float(r), float(z))
        exact_vec = np.hypot(exact_r, exact_z)
        fast_vec = np.hypot(fast_r, fast_z)
        rows.append(
            {
                "label": label,
                "R_pc": r,
                "z_pc": z,
                "exact_force_R": exact_r,
                "interp_force_R": fast_r,
                "force_R_rel_error": abs(fast_r - exact_r) / max(abs(exact_r), 1e-30),
                "exact_force_z": exact_z,
                "interp_force_z": fast_z,
                "force_z_rel_error": abs(fast_z - exact_z) / max(abs(exact_z), 1e-30),
                "exact_force_vector": exact_vec,
                "interp_force_vector": fast_vec,
                "force_vector_rel_error": abs(fast_vec - exact_vec) / max(abs(exact_vec), 1e-30),
            }
        )
    return pd.DataFrame(rows)


def append_log(log_path: Path, summary: pd.DataFrame, detail_path: Path, force_path: Path, summary_path: Path) -> None:
    best = summary.sort_values(["sigma_rel_error_p95", "force_grid_total_build_time_s"]).iloc[0]
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("\n\n### Strategy 3 Test Run 001\n\n")
        handle.write("Status: completed isolated halo force-grid interpolation test using the experimental R-z moment grid + cubic `dP_z/dR`; production code unchanged.\n\n")
        handle.write("- Method: build the accepted R-z moment grid, but replace direct halo force calls on every moment-grid point with a coarser halo force grid and linear interpolation.\n")
        handle.write(f"- sigma_los2 detail CSV: `{detail_path}`\n")
        handle.write(f"- force detail CSV: `{force_path}`\n")
        handle.write(f"- summary CSV: `{summary_path}`\n")
        for _, row in summary.iterrows():
            handle.write(
                f"- {row['label']}: total build `{row['force_grid_total_build_time_s']:.3f} s`, "
                f"build speedup `{row['build_speedup_vs_direct_force']:.3f}x`, "
                f"sigma median/max/p95 rel error "
                f"`{row['sigma_rel_error_median']:.6g}` / `{row['sigma_rel_error_max']:.6g}` / `{row['sigma_rel_error_p95']:.6g}`, "
                f"force-vector median/max/p95 rel error "
                f"`{row['force_vector_rel_error_median']:.6g}` / `{row['force_vector_rel_error_max']:.6g}` / `{row['force_vector_rel_error_p95']:.6g}`, "
                f"`Delta lnL={row['delta_ln_l']:.6g}`.\n"
            )
        handle.write(
            f"- Best p95 sigma-error setting in this run: `{best['label']}` with "
            f"p95 sigma relative error `{best['sigma_rel_error_p95']:.6g}` and "
            f"`Delta lnL={best['delta_ln_l']:.6g}`.\n"
        )


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
