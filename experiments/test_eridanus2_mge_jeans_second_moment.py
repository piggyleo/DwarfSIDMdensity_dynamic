#!/usr/bin/env python
"""Compare Hayashi validation and MGE/JAM LOS second-moment formulae.

This experiment is deliberately standalone. It does not modify the sampler or
likelihood code. For each selected Eridanus II posterior sample it compares:

1. The existing sampler validation projector.
2. A direct-intrinsic-3D MGE numerator divided by projected tracer density.
3. The literal RHS of Shajib (2019) Eq. 3.9 after converting the same 3D MGE
   amplitudes to the paper's 2D projected amplitudes.

The MGE decompositions are positive NNLS fits to the intrinsic 1D density
profiles on logarithmic radii. This isolates the normalization/question around
Eq. 3.9 without changing production code.
"""

from __future__ import annotations

import argparse
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.special import roots_legendre

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hayashi_jeans.constants import G_PC_MSUN_KMS2
from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.tracer import intrinsic_q_from_projected
from scripts.run_willman1_block_mh_fast import SlowParams, beta_z_from_q, compute_sigma_unit


@dataclass(frozen=True)
class MGE1D:
    amplitudes: np.ndarray
    sigmas_major_pc: np.ndarray
    rel_err_median: float
    rel_err_p95: float
    rel_err_max: float
    n_nonzero: int


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_vs_strict_selected_samples.csv"))
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--reference-mode", choices=("validation-convergence", "validation-halo", "validation-strict"), default="validation-convergence")
    parser.add_argument("--strict-epsrel", type=float, default=1.5e-2)
    parser.add_argument("--n-gauss-halo", type=int, default=45)
    parser.add_argument("--n-gauss-tracer", type=int, default=45)
    parser.add_argument("--n-fit-radii", type=int, default=800)
    parser.add_argument("--n-u", type=int, default=96)
    parser.add_argument("--r-min-pc", type=float, default=1e-2)
    parser.add_argument("--r-max-pc", type=float, default=2e5)
    parser.add_argument("--sample-limit", type=int, default=None)
    parser.add_argument("--selection-id", type=int, action="append", default=None)
    parser.add_argument("--top-n-error-samples", type=int, default=None)
    parser.add_argument("--error-detail", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_error_source_detail.csv"))
    parser.add_argument("--error-column", default="fast_vs_strict_abs_delta_log_likelihood")
    parser.add_argument("--output-detail", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_mge_jeans_second_moment_detail.csv"))
    parser.add_argument("--output-summary", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_mge_jeans_second_moment_summary.csv"))
    parser.add_argument("--output-stars", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_mge_jeans_second_moment_per_star.csv"))
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    selected = pd.read_csv(args.selected)
    if args.top_n_error_samples is not None:
        error_detail = pd.read_csv(args.error_detail)
        top_ids = (
            error_detail.sort_values(args.error_column, ascending=False)
            .head(args.top_n_error_samples)["selection_id"]
            .astype(int)
            .tolist()
        )
        selected = selected[selected["selection_id"].astype(int).isin(top_ids)].copy()
        selected["_top_order"] = selected["selection_id"].astype(int).map({sid: i for i, sid in enumerate(top_ids)})
        selected = selected.sort_values("_top_order").drop(columns="_top_order")
    if args.selection_id is not None:
        selected = selected[selected["selection_id"].astype(int).isin(args.selection_id)].copy()
    if args.sample_limit is not None:
        selected = selected.head(args.sample_limit).copy()

    rows = []
    star_rows = []
    for row in selected.itertuples(index=False):
        row_dict = row._asdict()
        print(f"sample {int(row_dict['selection_id'])}/{len(selected)} source_row={int(row_dict['source_row'])}", flush=True)
        result, per_star = compare_sample(galaxy, row_dict, args)
        rows.append(result)
        star_rows.append(per_star)
        print(
            f"  reference_like={result['reference_likelihood_seconds']:.3f}s; "
            f"MGE_like={result['mge_total_seconds_including_fit']:.3f}s; "
            f"dlogL={result['mge_vs_reference_delta_log_likelihood']:.6g}; "
            f"median rel sigma2 MGE/ref={result['mge_vs_reference_sigma_los2_rel_err_median']:.3g}; "
            f"literal/I rel={result['shajib_literal_div_I_vs_mge_sigma_los2_rel_err_max']:.3g}; "
            f"literal/mge median ratio={result['shajib_literal_rhs_to_mge_sigma_los2_ratio_median']:.3g}",
            flush=True,
        )

    detail = pd.DataFrame(rows)
    detail_path = Path(args.output_detail)
    detail_path.parent.mkdir(parents=True, exist_ok=True)
    detail.to_csv(detail_path, index=False)

    summary = summarize(detail)
    summary_path = Path(args.output_summary)
    summary.to_csv(summary_path, index=False)
    stars = pd.concat(star_rows, ignore_index=True)
    stars_path = Path(args.output_stars)
    stars.to_csv(stars_path, index=False)
    print(f"wrote {detail_path}")
    print(detail.to_string(index=False))
    print(f"wrote {summary_path}")
    print(summary.to_string(index=False))
    print(f"wrote {stars_path}")


def compare_sample(galaxy, row: dict[str, float], args: argparse.Namespace) -> tuple[dict[str, float | int | str], pd.DataFrame]:
    log10_rho0 = float(row["log10_rho0_msun_pc3"])
    density_scale = 10.0**log10_rho0
    systemic_velocity_kms = float(row["systemic_velocity_kms"])
    likelihood = GaussianVelocityLikelihood(galaxy)
    slow = SlowParams(
        q_halo=float(row["q_halo"]),
        log10_b_halo_pc=float(row["log10_b_halo_pc"]),
        minus_log10_one_minus_beta_z=float(row["minus_log10_one_minus_beta_z"]),
        alpha=float(row["alpha"]),
        beta=float(row["beta"]),
        gamma=float(row["gamma"]),
        inclination_deg=float(row["i_deg"]),
    )

    t0 = time.perf_counter()
    reference_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode=args.reference_mode,
        n_r=96,
        n_z=192,
        n_los=160,
        n_force_r=48,
        n_force_z=96,
        strict_epsrel=args.strict_epsrel,
    )
    reference_sigma2 = reference_unit * density_scale
    reference_log_likelihood = likelihood.log_likelihood(reference_sigma2, systemic_velocity_kms=systemic_velocity_kms)
    reference_likelihood_seconds = time.perf_counter() - t0

    t_fit = time.perf_counter()
    inclination_rad = np.deg2rad(float(row["i_deg"]))
    q_star = intrinsic_q_from_projected(galaxy.observables.qprime, inclination_rad)
    tracer_mge = fit_positive_mge(
        lambda m: (1.0 + (m / galaxy.observables.b_star_pc) ** 2) ** (-2.5),
        n_gauss=args.n_gauss_tracer,
        n_fit_radii=args.n_fit_radii,
        r_min_pc=args.r_min_pc,
        r_max_pc=args.r_max_pc,
    )
    halo_mge = fit_positive_mge(
        lambda m: generalized_halo_density_unit(
            m,
            b_pc=10.0 ** float(row["log10_b_halo_pc"]),
            alpha=float(row["alpha"]),
            beta=float(row["beta"]),
            gamma=float(row["gamma"]),
        ),
        n_gauss=args.n_gauss_halo,
        n_fit_radii=args.n_fit_radii,
        r_min_pc=args.r_min_pc,
        r_max_pc=args.r_max_pc,
    )
    mge_fit_seconds = time.perf_counter() - t_fit

    t0 = time.perf_counter()
    mge_sigma2_unit, projected_tracer, shajib_literal_rhs_unit, shajib_literal_div_i_unit = mge_los_second_moment(
        x_pc=galaxy.x_pc,
        y_pc=galaxy.y_pc,
        halo_mge=halo_mge,
        tracer_mge=tracer_mge,
        q_halo=float(row["q_halo"]),
        q_star=q_star,
        beta_z=beta_z_from_q(float(row["minus_log10_one_minus_beta_z"])),
        inclination_rad=inclination_rad,
        n_u=args.n_u,
    )
    mge_sigma2 = mge_sigma2_unit * density_scale
    shajib_literal_rhs = shajib_literal_rhs_unit * density_scale
    shajib_literal_div_i = shajib_literal_div_i_unit * density_scale
    mge_log_likelihood = likelihood.log_likelihood(mge_sigma2, systemic_velocity_kms=systemic_velocity_kms)
    shajib_div_i_log_likelihood = likelihood.log_likelihood(shajib_literal_div_i, systemic_velocity_kms=systemic_velocity_kms)
    shajib_literal_rhs_log_likelihood = likelihood.log_likelihood(shajib_literal_rhs, systemic_velocity_kms=systemic_velocity_kms)
    mge_moment_seconds = time.perf_counter() - t0
    mge_total_seconds = mge_fit_seconds + mge_moment_seconds

    result: dict[str, float | int | str] = {
        "selection_id": int(row["selection_id"]),
        "source_row": int(row["source_row"]),
        "reference_mode": args.reference_mode,
        "q_halo": float(row["q_halo"]),
        "log10_b_halo_pc": float(row["log10_b_halo_pc"]),
        "log10_rho0_msun_pc3": log10_rho0,
        "minus_log10_one_minus_beta_z": float(row["minus_log10_one_minus_beta_z"]),
        "alpha": float(row["alpha"]),
        "beta": float(row["beta"]),
        "gamma": float(row["gamma"]),
        "i_deg": float(row["i_deg"]),
        "systemic_velocity_kms": systemic_velocity_kms,
        "q_star_intrinsic": float(q_star),
        "reference_likelihood_seconds": float(reference_likelihood_seconds),
        "reference_log_likelihood": float(reference_log_likelihood),
        "mge_fit_seconds": float(mge_fit_seconds),
        "mge_moment_seconds": float(mge_moment_seconds),
        "mge_total_seconds_including_fit": float(mge_total_seconds),
        "mge_log_likelihood": float(mge_log_likelihood),
        "mge_vs_reference_delta_log_likelihood": safe_delta(mge_log_likelihood, reference_log_likelihood),
        "mge_vs_reference_abs_delta_log_likelihood": safe_abs_delta(mge_log_likelihood, reference_log_likelihood),
        "shajib_div_I_log_likelihood": float(shajib_div_i_log_likelihood),
        "shajib_div_I_vs_reference_delta_log_likelihood": safe_delta(shajib_div_i_log_likelihood, reference_log_likelihood),
        "shajib_div_I_vs_reference_abs_delta_log_likelihood": safe_abs_delta(shajib_div_i_log_likelihood, reference_log_likelihood),
        "shajib_literal_rhs_log_likelihood": float(shajib_literal_rhs_log_likelihood),
        "mge_speedup_vs_reference_moment_only": float(reference_likelihood_seconds / mge_moment_seconds),
        "mge_speedup_vs_reference_including_fit": float(reference_likelihood_seconds / mge_total_seconds),
        "mge_sigma_los2_nonpositive_count": int(np.sum(~np.isfinite(mge_sigma2) | (mge_sigma2 <= 0.0))),
        "reference_sigma_los2_nonpositive_count": int(np.sum(~np.isfinite(reference_sigma2) | (reference_sigma2 <= 0.0))),
        "halo_mge_rel_err_median": halo_mge.rel_err_median,
        "halo_mge_rel_err_p95": halo_mge.rel_err_p95,
        "halo_mge_rel_err_max": halo_mge.rel_err_max,
        "halo_mge_n_nonzero": halo_mge.n_nonzero,
        "tracer_mge_rel_err_median": tracer_mge.rel_err_median,
        "tracer_mge_rel_err_p95": tracer_mge.rel_err_p95,
        "tracer_mge_rel_err_max": tracer_mge.rel_err_max,
        "tracer_mge_n_nonzero": tracer_mge.n_nonzero,
        "projected_tracer_min": float(np.nanmin(projected_tracer)),
        "projected_tracer_median": float(np.nanmedian(projected_tracer)),
    }
    add_sigma_comparison(result, mge_sigma2, reference_sigma2, "mge_vs_reference")
    add_sigma_comparison(result, shajib_literal_div_i, mge_sigma2, "shajib_literal_div_I_vs_mge")

    literal_ratio = shajib_literal_rhs / np.maximum(mge_sigma2, 1.0e-300)
    result["shajib_literal_rhs_to_mge_sigma_los2_ratio_min"] = float(np.nanmin(literal_ratio))
    result["shajib_literal_rhs_to_mge_sigma_los2_ratio_median"] = float(np.nanmedian(literal_ratio))
    result["shajib_literal_rhs_to_mge_sigma_los2_ratio_max"] = float(np.nanmax(literal_ratio))
    result["shajib_literal_rhs_min"] = float(np.nanmin(shajib_literal_rhs))
    result["shajib_literal_rhs_median"] = float(np.nanmedian(shajib_literal_rhs))
    result["mge_sigma_los2_median"] = float(np.nanmedian(mge_sigma2))
    result["reference_sigma_los2_median"] = float(np.nanmedian(reference_sigma2))
    per_star = pd.DataFrame(
        {
            "selection_id": int(row["selection_id"]),
            "source_row": int(row["source_row"]),
            "star_index": np.arange(len(reference_sigma2), dtype=int),
            "x_pc": galaxy.x_pc,
            "y_pc": galaxy.y_pc,
            "reference_sigma_los2_km2_s2": reference_sigma2,
            "mge_derived_sigma_los2_km2_s2": mge_sigma2,
            "shajib_eq39_literal_rhs": shajib_literal_rhs,
            "shajib_eq39_rhs_div_projected_tracer_km2_s2": shajib_literal_div_i,
            "projected_tracer": projected_tracer,
            "mge_vs_reference_sigma_los2_rel_err": np.abs(mge_sigma2 - reference_sigma2)
            / np.maximum(np.abs(reference_sigma2), 1.0e-300),
            "shajib_div_I_vs_mge_sigma_los2_rel_err": np.abs(shajib_literal_div_i - mge_sigma2)
            / np.maximum(np.abs(mge_sigma2), 1.0e-300),
            "shajib_literal_rhs_to_mge_sigma_los2_ratio": shajib_literal_rhs / np.maximum(mge_sigma2, 1.0e-300),
        }
    )
    return result, per_star


def safe_delta(left: float, right: float) -> float:
    if np.isfinite(left) and np.isfinite(right):
        return float(left - right)
    return float("nan")


def safe_abs_delta(left: float, right: float) -> float:
    delta = safe_delta(left, right)
    return float(abs(delta)) if np.isfinite(delta) else float("nan")


def fit_positive_mge(
    density_func,
    *,
    n_gauss: int,
    n_fit_radii: int,
    r_min_pc: float,
    r_max_pc: float,
) -> MGE1D:
    sigmas = np.exp(np.linspace(np.log(r_min_pc / 3.0), np.log(r_max_pc * 3.0), n_gauss))
    radii = np.exp(np.linspace(np.log(r_min_pc), np.log(r_max_pc), n_fit_radii))
    target = np.asarray(density_func(radii), dtype=float)
    basis = np.exp(-0.5 * (radii[:, None] / sigmas[None, :]) ** 2)
    weights = 1.0 / np.maximum(np.abs(target), 1.0e-300)
    amplitudes, _ = nnls(basis * weights[:, None], target * weights, maxiter=10000)

    approx = basis @ amplitudes
    rel_err = np.abs(approx - target) / np.maximum(np.abs(target), 1.0e-300)
    return MGE1D(
        amplitudes=amplitudes,
        sigmas_major_pc=sigmas,
        rel_err_median=float(np.nanmedian(rel_err)),
        rel_err_p95=float(np.nanpercentile(rel_err, 95.0)),
        rel_err_max=float(np.nanmax(rel_err)),
        n_nonzero=int(np.sum(amplitudes > 0.0)),
    )


def generalized_halo_density_unit(m: np.ndarray, *, b_pc: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    x = np.maximum(np.asarray(m, dtype=float) / b_pc, 1.0e-300)
    return x ** (-gamma) * (1.0 + x**alpha) ** (-(beta - gamma) / alpha)


def mge_los_second_moment(
    *,
    x_pc: np.ndarray,
    y_pc: np.ndarray,
    halo_mge: MGE1D,
    tracer_mge: MGE1D,
    q_halo: float,
    q_star: float,
    beta_z: float,
    inclination_rad: float,
    n_u: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nodes, weights = roots_legendre(n_u)
    u = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    x = np.asarray(x_pc, dtype=float)
    y = np.asarray(y_pc, dtype=float)
    cos_i = np.cos(inclination_rad)
    sin_i = np.sin(inclination_rad)
    b_aniso = 1.0 / (1.0 - beta_z)

    qpj2 = cos_i**2 + q_halo**2 * sin_i**2
    qpk2 = cos_i**2 + q_star**2 * sin_i**2
    qpj = np.sqrt(qpj2)
    qpk = np.sqrt(qpk2)

    sigma_h = q_halo * halo_mge.sigmas_major_pc
    rho0 = halo_mge.amplitudes
    tracer_sigma_major = tracer_mge.sigmas_major_pc
    sigma_t = q_star * tracer_sigma_major
    nu0 = tracer_mge.amplitudes

    projected_tracer = np.zeros_like(x, dtype=float)
    for amp_k, sig_k_major in zip(nu0, tracer_sigma_major):
        if amp_k <= 0.0:
            continue
        projected_tracer += (
            np.sqrt(2.0 * np.pi)
            * q_star
            * sig_k_major
            / qpk
            * amp_k
            * np.exp(-(x**2 + y**2 / qpk2) / (2.0 * sig_k_major**2))
        )

    numerator = np.zeros_like(x, dtype=float)
    shajib_literal = np.zeros_like(x, dtype=float)
    u2 = u * u
    x2 = x[None, :] ** 2
    y2 = y[None, :] ** 2

    for amp_j, sig_j in zip(rho0, sigma_h):
        if amp_j <= 0.0:
            continue
        one_minus_qj2 = 1.0 - q_halo**2
        halo_shape = 1.0 - one_minus_qj2 * u2
        for amp_k, sig_k in zip(nu0, sigma_t):
            if amp_k <= 0.0:
                continue
            a = 0.5 * (u2 * q_halo**2 / sig_j**2 + q_star**2 / sig_k**2)
            b = 0.5 * (
                (1.0 - q_star**2) / sig_k**2
                + q_halo**2 * (1.0 - q_halo**2) * u2**2 / (sig_j**2 * halo_shape)
            )
            c = 1.0 - q_halo**2 - q_halo**2 * sig_k**2 / sig_j**2
            d = 1.0 - b_aniso * q_star**2 - ((1.0 - b_aniso) * c + (1.0 - q_halo**2) * b_aniso) * u2
            denom = (1.0 - c * u2) * np.sqrt((a + b * cos_i**2) * halo_shape)
            velocity_factor = sig_k**2 * (cos_i**2 + b_aniso * sin_i**2)
            exponent = -a[:, None] * (x2 + ((a + b) / (a + b * cos_i**2))[:, None] * y2)
            kernel = (
                u2[:, None]
                * (velocity_factor + d[:, None] * x2 * sin_i**2)
                / denom[:, None]
                * np.exp(exponent)
            )
            integral = np.sum(w[:, None] * kernel, axis=0)
            numerator += 4.0 * np.pi ** 1.5 * G_PC_MSUN_KMS2 * q_halo * amp_j * amp_k * integral

            sigma0_j = np.sqrt(2.0 * np.pi) * sig_j * qpj / (q_halo**2) * amp_j
            i0_k = np.sqrt(2.0 * np.pi) * sig_k * qpk / (q_star**2) * amp_k
            shajib_prefactor = (
                2.0
                * np.sqrt(np.pi)
                * G_PC_MSUN_KMS2
                * q_halo**3
                * q_star**2
                * sigma0_j
                * i0_k
                / (sig_j * qpj * sig_k * qpk)
            )
            shajib_literal += shajib_prefactor * integral

    mge_sigma2 = numerator / np.maximum(projected_tracer, 1.0e-300)
    shajib_div_i = shajib_literal / np.maximum(projected_tracer, 1.0e-300)
    return mge_sigma2, projected_tracer, shajib_literal, shajib_div_i


def add_sigma_comparison(result: dict[str, float | int | str], left_sigma2: np.ndarray, right_sigma2: np.ndarray, label: str) -> None:
    rel = np.abs(left_sigma2 - right_sigma2) / np.maximum(np.abs(right_sigma2), 1.0e-300)
    abs_err = np.abs(left_sigma2 - right_sigma2)
    sigma_rel = np.abs(np.sqrt(np.maximum(left_sigma2, 0.0)) - np.sqrt(np.maximum(right_sigma2, 0.0))) / np.maximum(
        np.sqrt(np.maximum(right_sigma2, 0.0)), 1.0e-300
    )
    result[f"{label}_sigma_los2_abs_err_median"] = float(np.nanmedian(abs_err))
    result[f"{label}_sigma_los2_abs_err_max"] = float(np.nanmax(abs_err))
    result[f"{label}_sigma_los2_rel_err_median"] = float(np.nanmedian(rel))
    result[f"{label}_sigma_los2_rel_err_p95"] = float(np.nanpercentile(rel, 95.0))
    result[f"{label}_sigma_los2_rel_err_max"] = float(np.nanmax(rel))
    result[f"{label}_sigma_los_rel_err_median"] = float(np.nanmedian(sigma_rel))
    result[f"{label}_sigma_los_rel_err_p95"] = float(np.nanpercentile(sigma_rel, 95.0))
    result[f"{label}_sigma_los_rel_err_max"] = float(np.nanmax(sigma_rel))


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    rows = []
    rows.append(
        {
            "label": "mge_likelihood_vs_reference",
            "n_samples": len(detail),
            "sigma_los2_rel_err_median_median_across_samples": float(detail["mge_vs_reference_abs_delta_log_likelihood"].median(skipna=True)),
            "sigma_los2_rel_err_median_max_across_samples": float(detail["mge_vs_reference_abs_delta_log_likelihood"].max(skipna=True)),
            "sigma_los2_rel_err_p95_median_across_samples": float(detail["mge_total_seconds_including_fit"].median()),
            "sigma_los2_rel_err_p95_max_across_samples": float(detail["reference_likelihood_seconds"].median()),
            "sigma_los2_rel_err_max_median_across_samples": float(detail["mge_speedup_vs_reference_including_fit"].median()),
            "sigma_los2_rel_err_max_max_across_samples": float(detail["mge_speedup_vs_reference_including_fit"].max()),
            "sigma_los_rel_err_median_median_across_samples": np.nan,
            "sigma_los_rel_err_median_max_across_samples": np.nan,
            "sigma_los_rel_err_p95_median_across_samples": np.nan,
            "sigma_los_rel_err_p95_max_across_samples": np.nan,
            "sigma_los_rel_err_max_median_across_samples": np.nan,
            "sigma_los_rel_err_max_max_across_samples": np.nan,
        }
    )
    for label in ["mge_vs_reference", "shajib_literal_div_I_vs_mge"]:
        row = {"label": label, "n_samples": len(detail)}
        for metric in [
            "sigma_los2_rel_err_median",
            "sigma_los2_rel_err_p95",
            "sigma_los2_rel_err_max",
            "sigma_los_rel_err_median",
            "sigma_los_rel_err_p95",
            "sigma_los_rel_err_max",
        ]:
            col = f"{label}_{metric}"
            row[f"{metric}_median_across_samples"] = float(detail[col].median())
            row[f"{metric}_max_across_samples"] = float(detail[col].max())
        rows.append(row)
    rows.append(
        {
            "label": "shajib_literal_rhs_to_mge_sigma_los2_ratio",
            "n_samples": len(detail),
            "sigma_los2_rel_err_median_median_across_samples": float(detail["shajib_literal_rhs_to_mge_sigma_los2_ratio_median"].median()),
            "sigma_los2_rel_err_median_max_across_samples": float(detail["shajib_literal_rhs_to_mge_sigma_los2_ratio_median"].max()),
            "sigma_los2_rel_err_p95_median_across_samples": np.nan,
            "sigma_los2_rel_err_p95_max_across_samples": np.nan,
            "sigma_los2_rel_err_max_median_across_samples": np.nan,
            "sigma_los2_rel_err_max_max_across_samples": np.nan,
            "sigma_los_rel_err_median_median_across_samples": np.nan,
            "sigma_los_rel_err_median_max_across_samples": np.nan,
            "sigma_los_rel_err_p95_median_across_samples": np.nan,
            "sigma_los_rel_err_p95_max_across_samples": np.nan,
            "sigma_los_rel_err_max_median_across_samples": np.nan,
            "sigma_los_rel_err_max_max_across_samples": np.nan,
        }
    )
    return pd.DataFrame(rows)


if __name__ == "__main__":
    main()
