#!/usr/bin/env python
"""Experiment: three-layer parameter cache versus the original slow projector.

Layer 1: halo-shape parameters build a unit-rho0 halo force grid.
Layer 2: beta_z and inclination build/project the R-z Jeans moment grid.
Layer 3: rho0 rescales sigma_los^2 and systemic velocity enters likelihood.

This standalone test does not modify production code. Slow-reference results
are checkpointed so interrupted runs do not need to recompute completed
`src/projection.py` likelihoods.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from experiments.test_halo_force_grid_interpolation import RZMomentGridWithForceGrid
from hayashi_jeans.data import GalaxyData, load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.model import build_hernquist_projector
from hayashi_jeans.params import HayashiParameters


@dataclass(frozen=True)
class Layer1Params:
    label: str
    q_halo: float
    log10_b_halo_pc: float
    alpha: float
    beta: float
    gamma: float


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default="data/galaxies/27_Willman_1.csv")
    parser.add_argument("--centers-csv", default="data/processed/galaxy_structural_centers.csv")
    parser.add_argument("--out-dir", default="outputs/diagnostics/three_layer_cache")
    parser.add_argument("--log-path", default="logs/likelihood_acceleration_strategy_log.md")
    parser.add_argument("--n-r", type=int, default=96)
    parser.add_argument("--n-z", type=int, default=192)
    parser.add_argument("--n-los", type=int, default=160)
    parser.add_argument("--n-force-r", type=int, default=48)
    parser.add_argument("--n-force-z", type=int, default=96)
    parser.add_argument("--n-force-points", type=int, default=120)
    parser.add_argument("--n-moment-points", type=int, default=8)
    parser.add_argument("--zmax-factor", type=float, default=20.0)
    parser.add_argument("--los-factor", type=float, default=20.0)
    parser.add_argument("--slow-epsrel", type=float, default=1.5e-2)
    parser.add_argument("--r-min-pc", type=float, default=1.0e-3)
    parser.add_argument("--z-min-pc", type=float, default=1.0e-3)
    parser.add_argument("--grid-padding", type=float, default=1.08)
    parser.add_argument("--beta-z", type=float, default=0.3)
    parser.add_argument("--inclination-deg", type=float, default=75.0)
    parser.add_argument("--log10-rho0", type=float, default=-1.961219)
    parser.add_argument("--systemic-velocity-kms", type=float, default=-13.298645)
    parser.add_argument("--seed", type=int, default=20260520)
    parser.add_argument("--max-groups", type=int, default=5)
    parser.add_argument("--skip-slow", action="store_true", help="Only run the three-layer path and write no slow-reference errors.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    layer1_values = default_layer1_params()[: args.max_groups]

    slow_sigma_path = out_dir / "three_layer_cache_willman1_slow_sigma_checkpoint.csv"
    slow_moment_path = out_dir / "three_layer_cache_willman1_slow_moment_checkpoint.csv"
    slow_sigma_cache = load_cache(slow_sigma_path)
    slow_moment_cache = load_cache(slow_moment_path)

    likelihood = GaussianVelocityLikelihood(galaxy)
    sigma_rows = []
    moment_rows = []
    force_rows = []
    summary_rows = []

    for idx, layer1 in enumerate(layer1_values):
        print(f"\n=== {layer1.label} ({idx + 1}/{len(layer1_values)}) ===", flush=True)
        params = make_params(layer1, args, rho0=10.0**args.log10_rho0)
        unit_params = make_params(layer1, args, rho0=1.0)

        slow_sigma = None
        slow_ln_l = np.nan
        slow_sigma_time = np.nan
        slow_total_time = np.nan
        slow_moment_time = np.nan

        if not args.skip_slow:
            slow_result = get_or_compute_slow_sigma(
                galaxy,
                params,
                args,
                layer1.label,
                likelihood,
                slow_sigma_cache,
                slow_sigma_path,
            )
            slow_sigma = slow_result["sigma"]
            slow_ln_l = slow_result["ln_l"]
            slow_sigma_time = slow_result["sigma_time_s"]
            slow_total_time = slow_result["total_time_s"]

        t0 = time.perf_counter()
        cached_projector = RZMomentGridWithForceGrid(
            galaxy=galaxy,
            params=unit_params,
            n_r=args.n_r,
            n_z=args.n_z,
            n_los=args.n_los,
            n_force_r=args.n_force_r,
            n_force_z=args.n_force_z,
            zmax_factor=args.zmax_factor,
            los_factor=args.los_factor,
            r_min_pc=args.r_min_pc,
            z_min_pc=args.z_min_pc,
            grid_padding=args.grid_padding,
        )
        cached_build_time = time.perf_counter() - t0

        t0 = time.perf_counter()
        sigma_unit = cached_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
        layer2_projection_time = time.perf_counter() - t0

        t0 = time.perf_counter()
        cached_sigma = sigma_unit * params.rho0_msun_pc3
        cached_ln_l = likelihood.log_likelihood(cached_sigma, systemic_velocity_kms=params.systemic_velocity_kms)
        layer3_time = time.perf_counter() - t0

        force_df = compare_forces(cached_projector, layer1.label, args.n_force_points, args.seed + idx)
        force_rows.append(force_df)

        if not args.skip_slow:
            moment_result = get_or_compute_slow_moments(
                cached_projector,
                galaxy,
                params,
                args,
                layer1.label,
                idx,
                slow_moment_cache,
                slow_moment_path,
            )
            moment_df = moment_result["moments"]
            slow_moment_time = moment_result["time_s"]
            moment_rows.append(moment_df)
        else:
            moment_df = pd.DataFrame()

        sigma_df = make_sigma_detail(galaxy, layer1.label, slow_sigma, cached_sigma)
        sigma_rows.append(sigma_df)

        sigma_rel = sigma_df["sigma_rel_error"].dropna()
        moment_rel = moment_df.filter(like="_rel_error").stack() if not moment_df.empty else pd.Series(dtype=float)
        force_vec_rel = force_df["force_vector_rel_error"]
        cached_total_time = cached_build_time + layer2_projection_time + layer3_time
        summary = {
            "label": layer1.label,
            "q_halo": layer1.q_halo,
            "log10_b_halo_pc": layer1.log10_b_halo_pc,
            "alpha": layer1.alpha,
            "beta": layer1.beta,
            "gamma": layer1.gamma,
            "beta_z": args.beta_z,
            "inclination_deg": args.inclination_deg,
            "log10_rho0": args.log10_rho0,
            "systemic_velocity_kms": args.systemic_velocity_kms,
            "n_star": len(galaxy.stars),
            "n_r": args.n_r,
            "n_z": args.n_z,
            "n_los": args.n_los,
            "n_force_r": args.n_force_r,
            "n_force_z": args.n_force_z,
            "slow_sigma_time_s": slow_sigma_time,
            "slow_moment_time_s": slow_moment_time,
            "slow_total_time_s": slow_total_time,
            "cached_force_grid_build_time_s": cached_projector.force_grid_build_time_s,
            "cached_moment_grid_build_time_s": cached_projector.moment_grid_build_time_s,
            "cached_build_time_s": cached_build_time,
            "cached_projection_time_s": layer2_projection_time,
            "cached_layer3_likelihood_time_s": layer3_time,
            "cached_total_likelihood_time_s": cached_total_time,
            "speedup_vs_slow_total": slow_total_time / cached_total_time if np.isfinite(slow_total_time) else np.nan,
            "force_vector_rel_error_median": force_vec_rel.median(),
            "force_vector_rel_error_p95": force_vec_rel.quantile(0.95),
            "force_vector_rel_error_max": force_vec_rel.max(),
            "moment_rel_error_median": moment_rel.median() if not moment_rel.empty else np.nan,
            "moment_rel_error_p95": moment_rel.quantile(0.95) if not moment_rel.empty else np.nan,
            "moment_rel_error_max": moment_rel.max() if not moment_rel.empty else np.nan,
            "sigma_rel_error_median": sigma_rel.median() if not sigma_rel.empty else np.nan,
            "sigma_rel_error_p95": sigma_rel.quantile(0.95) if not sigma_rel.empty else np.nan,
            "sigma_rel_error_max": sigma_rel.max() if not sigma_rel.empty else np.nan,
            "slow_ln_l": slow_ln_l,
            "cached_ln_l": cached_ln_l,
            "delta_ln_l": cached_ln_l - slow_ln_l if np.isfinite(slow_ln_l) else np.nan,
        }
        summary_rows.append(summary)
        print(pd.Series(summary).to_string(), flush=True)

    sigma_detail = pd.concat(sigma_rows, ignore_index=True)
    force_detail = pd.concat(force_rows, ignore_index=True)
    moment_detail = pd.concat(moment_rows, ignore_index=True) if moment_rows else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)

    sigma_detail_path = out_dir / "three_layer_cache_willman1_sigma_detail.csv"
    force_detail_path = out_dir / "three_layer_cache_willman1_force_detail.csv"
    moment_detail_path = out_dir / "three_layer_cache_willman1_moment_detail.csv"
    summary_path = out_dir / "three_layer_cache_willman1_summary.csv"
    sigma_detail.to_csv(sigma_detail_path, index=False)
    force_detail.to_csv(force_detail_path, index=False)
    moment_detail.to_csv(moment_detail_path, index=False)
    summary.to_csv(summary_path, index=False)
    append_log(Path(args.log_path), summary, sigma_detail_path, force_detail_path, moment_detail_path, summary_path)

    print("\nSummary:")
    print(summary.to_string(index=False), flush=True)


def default_layer1_params() -> list[Layer1Params]:
    return [
        Layer1Params("L1_fiducial", q_halo=1.0, log10_b_halo_pc=3.2, alpha=1.8, beta=6.4, gamma=1.2),
        Layer1Params("L1_oblate_compact", q_halo=0.65, log10_b_halo_pc=2.85, alpha=1.4, beta=5.2, gamma=0.8),
        Layer1Params("L1_prolate_extended", q_halo=1.45, log10_b_halo_pc=3.55, alpha=2.2, beta=7.4, gamma=1.4),
        Layer1Params("L1_cored_steep_outer", q_halo=0.9, log10_b_halo_pc=3.0, alpha=2.6, beta=8.2, gamma=0.35),
        Layer1Params("L1_cuspy_shallow_outer", q_halo=1.2, log10_b_halo_pc=3.35, alpha=1.1, beta=4.6, gamma=1.65),
    ]


def make_params(layer1: Layer1Params, args: argparse.Namespace, *, rho0: float) -> HayashiParameters:
    return HayashiParameters(
        q_halo=layer1.q_halo,
        b_halo_pc=10.0**layer1.log10_b_halo_pc,
        rho0_msun_pc3=rho0,
        beta_z=args.beta_z,
        alpha=layer1.alpha,
        beta=layer1.beta,
        gamma=layer1.gamma,
        inclination_rad=np.deg2rad(args.inclination_deg),
        systemic_velocity_kms=args.systemic_velocity_kms,
    )


def load_cache(path: Path) -> pd.DataFrame:
    if path.exists():
        return pd.read_csv(path)
    return pd.DataFrame()


def get_or_compute_slow_sigma(
    galaxy: GalaxyData,
    params: HayashiParameters,
    args: argparse.Namespace,
    label: str,
    likelihood: GaussianVelocityLikelihood,
    cache: pd.DataFrame,
    cache_path: Path,
) -> dict[str, object]:
    if not cache.empty and label in set(cache["label"]):
        rows = cache.loc[cache["label"] == label].sort_values("star_index")
        sigma = rows["slow_sigma_los2"].to_numpy(float)
        return {
            "sigma": sigma,
            "ln_l": float(rows["slow_ln_l"].iloc[0]),
            "sigma_time_s": float(rows["slow_sigma_time_s"].iloc[0]),
            "total_time_s": float(rows["slow_total_time_s"].iloc[0]),
        }

    print(f"Computing original slow projector for {label}...", flush=True)
    t0 = time.perf_counter()
    slow_projector = build_hernquist_projector(
        galaxy,
        params,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        epsrel=args.slow_epsrel,
    )
    build_time = time.perf_counter() - t0
    t0 = time.perf_counter()
    sigma = slow_projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    sigma_time = time.perf_counter() - t0
    ln_l = likelihood.log_likelihood(sigma, systemic_velocity_kms=params.systemic_velocity_kms)
    total_time = build_time + sigma_time

    rows = pd.DataFrame(
        {
            "label": label,
            "star_index": np.arange(len(galaxy.stars)),
            "star_id": galaxy.stars["star_id"],
            "slow_sigma_los2": sigma,
            "slow_ln_l": ln_l,
            "slow_build_time_s": build_time,
            "slow_sigma_time_s": sigma_time,
            "slow_total_time_s": total_time,
        }
    )
    updated = pd.concat([cache, rows], ignore_index=True)
    updated.to_csv(cache_path, index=False)
    return {"sigma": sigma, "ln_l": ln_l, "sigma_time_s": sigma_time, "total_time_s": total_time}


def get_or_compute_slow_moments(
    cached_projector: RZMomentGridWithForceGrid,
    galaxy: GalaxyData,
    params: HayashiParameters,
    args: argparse.Namespace,
    label: str,
    group_index: int,
    cache: pd.DataFrame,
    cache_path: Path,
) -> dict[str, object]:
    if not cache.empty and label in set(cache["label"]):
        rows = cache.loc[cache["label"] == label].copy()
        return {"moments": rows, "time_s": float(rows["slow_moment_time_s"].iloc[0])}

    points = moment_test_points(cached_projector, args.n_moment_points, args.seed + 17 * group_index)
    slow_projector = build_hernquist_projector(
        galaxy,
        params,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        epsrel=args.slow_epsrel,
    )
    rows = []
    t0 = time.perf_counter()
    for point_index, (r_pc, z_pc) in enumerate(points):
        slow = slow_projector.jeans_moments(float(r_pc), float(z_pc))
        fast_unit = cached_projector.jeans_moments(float(r_pc), float(z_pc))
        row = {"label": label, "point_index": point_index, "R_pc": r_pc, "z_pc": z_pc}
        for name in ("v_r2", "v_phi2", "v_z2"):
            slow_value = getattr(slow, name)
            cached_value = getattr(fast_unit, name) * params.rho0_msun_pc3
            row[f"slow_{name}"] = slow_value
            row[f"cached_{name}"] = cached_value
            row[f"{name}_rel_error"] = abs(cached_value - slow_value) / max(abs(slow_value), 1e-30)
        rows.append(row)
    elapsed = time.perf_counter() - t0
    frame = pd.DataFrame(rows)
    frame["slow_moment_time_s"] = elapsed
    updated = pd.concat([cache, frame], ignore_index=True)
    updated.to_csv(cache_path, index=False)
    return {"moments": frame, "time_s": elapsed}


def moment_test_points(projector: RZMomentGridWithForceGrid, n_points: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_inner = n_points // 2
    n_outer = n_points - n_inner
    r_values = np.concatenate(
        [
            10.0 ** rng.uniform(np.log10(0.05), np.log10(100.0), n_inner),
            10.0 ** rng.uniform(np.log10(100.0), np.log10(projector.r_grid[-1] * 0.9), n_outer),
        ]
    )
    z_values = np.concatenate(
        [
            10.0 ** rng.uniform(np.log10(0.05), np.log10(100.0), n_inner),
            10.0 ** rng.uniform(np.log10(100.0), np.log10(projector.z_grid[-1] * 0.9), n_outer),
        ]
    )
    return np.column_stack([r_values, z_values])


def compare_forces(projector: RZMomentGridWithForceGrid, label: str, n_points: int, seed: int) -> pd.DataFrame:
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
    unit_halo = GeneralizedHernquistHalo(
        q=projector.params.q_halo,
        b_pc=projector.params.b_halo_pc,
        rho0_msun_pc3=1.0,
        alpha=projector.params.alpha,
        beta=projector.params.beta,
        gamma=projector.params.gamma,
    )
    for r, z, cached_r, cached_z in zip(r_values, z_values, interp_force_r, interp_force_z):
        exact_r, exact_z = unit_halo.potential_gradients(float(r), float(z))
        exact_vec = np.hypot(exact_r, exact_z)
        cached_vec = np.hypot(cached_r, cached_z)
        rows.append(
            {
                "label": label,
                "R_pc": r,
                "z_pc": z,
                "exact_force_R_unit": exact_r,
                "cached_force_R_unit": cached_r,
                "force_R_rel_error": abs(cached_r - exact_r) / max(abs(exact_r), 1e-30),
                "exact_force_z_unit": exact_z,
                "cached_force_z_unit": cached_z,
                "force_z_rel_error": abs(cached_z - exact_z) / max(abs(exact_z), 1e-30),
                "exact_force_vector_unit": exact_vec,
                "cached_force_vector_unit": cached_vec,
                "force_vector_rel_error": abs(cached_vec - exact_vec) / max(abs(exact_vec), 1e-30),
            }
        )
    return pd.DataFrame(rows)


def make_sigma_detail(
    galaxy: GalaxyData,
    label: str,
    slow_sigma: np.ndarray | None,
    cached_sigma: np.ndarray,
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "label": label,
            "star_id": galaxy.stars["star_id"],
            "x_pc": galaxy.x_pc,
            "y_pc": galaxy.y_pc,
            "cached_sigma_los2": cached_sigma,
        }
    )
    if slow_sigma is not None:
        frame["slow_sigma_los2"] = slow_sigma
        frame["sigma_abs_error"] = np.abs(frame["cached_sigma_los2"] - frame["slow_sigma_los2"])
        frame["sigma_rel_error"] = frame["sigma_abs_error"] / np.maximum(np.abs(frame["slow_sigma_los2"]), 1e-30)
    else:
        frame["slow_sigma_los2"] = np.nan
        frame["sigma_abs_error"] = np.nan
        frame["sigma_rel_error"] = np.nan
    return frame


def append_log(
    log_path: Path,
    summary: pd.DataFrame,
    sigma_detail_path: Path,
    force_detail_path: Path,
    moment_detail_path: Path,
    summary_path: Path,
) -> None:
    if summary.empty:
        return
    speedup = summary["speedup_vs_slow_total"].dropna()
    delta = summary["delta_ln_l"].dropna().abs()
    sigma_p95 = summary["sigma_rel_error_p95"].dropna()
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("\n\n### Three-Layer Cache Test Run 001\n\n")
        handle.write("Status: completed standalone test of layer-1 halo force cache, layer-2 R-z moment/projection cache, and layer-3 rho0/systemic likelihood update; production code unchanged.\n\n")
        handle.write("- Compared against the original slow `src/projection.py` projector for five different layer-1 halo-shape parameter vectors.\n")
        handle.write(f"- sigma detail CSV: `{sigma_detail_path}`\n")
        handle.write(f"- force detail CSV: `{force_detail_path}`\n")
        handle.write(f"- moment detail CSV: `{moment_detail_path}`\n")
        handle.write(f"- summary CSV: `{summary_path}`\n")
        if not speedup.empty:
            handle.write(f"- Speedup versus slow total likelihood time: median `{speedup.median():.3g}x`, min `{speedup.min():.3g}x`, max `{speedup.max():.3g}x`.\n")
        if not sigma_p95.empty:
            handle.write(f"- Per-star sigma_los2 p95 relative error across layer-1 groups: median `{sigma_p95.median():.6g}`, max `{sigma_p95.max():.6g}`.\n")
        if not delta.empty:
            handle.write(f"- Absolute Delta lnL across layer-1 groups: median `{delta.median():.6g}`, max `{delta.max():.6g}`.\n")
        for _, row in summary.iterrows():
            handle.write(
                f"- {row['label']}: cached total `{row['cached_total_likelihood_time_s']:.3f} s`, "
                f"slow total `{row['slow_total_time_s']:.3f} s`, speedup `{row['speedup_vs_slow_total']:.3g}x`, "
                f"sigma p95/max rel error `{row['sigma_rel_error_p95']:.6g}` / `{row['sigma_rel_error_max']:.6g}`, "
                f"moment p95/max rel error `{row['moment_rel_error_p95']:.6g}` / `{row['moment_rel_error_max']:.6g}`, "
                f"force-vector p95/max rel error `{row['force_vector_rel_error_p95']:.6g}` / `{row['force_vector_rel_error_max']:.6g}`, "
                f"`Delta lnL={row['delta_ln_l']:.6g}`.\n"
            )


if __name__ == "__main__":
    main()
