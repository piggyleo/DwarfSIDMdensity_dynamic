#!/usr/bin/env python
"""Plot Eridanus II density profile with digitized Hayashi Figure 1 curves."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / "outputs" / ".matplotlib"))
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from skimage import measure

from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo
from scripts.plot_willman1_density_profile_with_hayashi_overlay import (
    DEFAULT_PDF,
    DEFAULT_PAGE_RENDER,
    FIG1_XLIM_KPC,
    FIG1_YLIM_MSUN_KPC3,
    parse_crop_box,
    pixel_x_to_kpc,
    pixel_y_to_density,
    render_pdf_page,
)


DEFAULT_CROP_BOX = "511,453,811,765"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--chain-output", default=str(PROJECT_ROOT / "outputs/eridanus_ii_nautilus_chain.csv"))
    parser.add_argument("--figure-output", default=str(PROJECT_ROOT / "outputs/figures/eridanus_ii_nautilus_density_profile_with_hayashi_overlay.png"))
    parser.add_argument("--comparison-output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_nautilus_density_profile_hayashi_comparison.csv"))
    parser.add_argument("--summary-output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_nautilus_density_profile_hayashi_summary.csv"))
    parser.add_argument("--burn", type=int, default=0)
    parser.add_argument("--hayashi-page-render", default=str(DEFAULT_PAGE_RENDER))
    parser.add_argument("--hayashi-pdf", default=DEFAULT_PDF)
    parser.add_argument("--pdf-page-index", type=int, default=4)
    parser.add_argument("--crop-box", default=DEFAULT_CROP_BOX)
    parser.add_argument("--digitized-output", default=str(PROJECT_ROOT / "outputs/diagnostics/hayashi_fig1_eridanus_ii_digitized_density.csv"))
    parser.add_argument("--crop-output", default=str(PROJECT_ROOT / "outputs/diagnostics/hayashi_fig1_eridanus_ii_crop.png"))
    parser.add_argument("--fill-mask-output", default=str(PROJECT_ROOT / "outputs/diagnostics/hayashi_fig1_eridanus_ii_fill_mask.png"))
    parser.add_argument("--line-mask-output", default=str(PROJECT_ROOT / "outputs/diagnostics/hayashi_fig1_eridanus_ii_line_mask.png"))
    args = parser.parse_args()

    if args.burn < 0:
        raise ValueError("--burn must be non-negative")

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    chain = pd.read_csv(args.chain_output)
    samples = chain.iloc[args.burn:].copy()
    if samples.empty:
        raise ValueError(f"no samples remain after burn={args.burn}; chain has {len(chain)} rows")

    hayashi_reference = digitize_eridanus_reference(
        page_render=Path(args.hayashi_page_render),
        pdf_path=Path(args.hayashi_pdf),
        pdf_page_index=args.pdf_page_index,
        crop_box=parse_crop_box(args.crop_box),
        crop_output=Path(args.crop_output),
        digitized_output=Path(args.digitized_output),
        fill_mask_output=Path(args.fill_mask_output),
        line_mask_output=Path(args.line_mask_output),
    )
    profile_summary = plot_density_profile_with_hayashi_reference(
        galaxy,
        samples,
        hayashi_reference,
        Path(args.figure_output),
        burn=args.burn,
    )
    comparison, summary = compare_profiles(galaxy, profile_summary, hayashi_reference, n_samples=len(samples))
    comparison_path = Path(args.comparison_output)
    summary_path = Path(args.summary_output)
    comparison_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(comparison_path, index=False)
    summary.to_csv(summary_path, index=False)
    print(f"wrote {comparison_path}")
    print(comparison.to_string(index=False))
    print(f"wrote {summary_path}")
    print(summary.to_string(index=False))


def digitize_eridanus_reference(
    *,
    page_render: Path,
    pdf_path: Path,
    pdf_page_index: int,
    crop_box: tuple[int, int, int, int],
    crop_output: Path,
    digitized_output: Path,
    fill_mask_output: Path,
    line_mask_output: Path,
) -> pd.DataFrame:
    if not page_render.exists():
        render_pdf_page(pdf_path, pdf_page_index, page_render)

    page = Image.open(page_render).convert("RGB")
    crop = page.crop(crop_box)
    crop_output.parent.mkdir(parents=True, exist_ok=True)
    crop.save(crop_output)

    rgb = np.asarray(crop, dtype=np.uint8)
    r = rgb[:, :, 0].astype(float)
    g = rgb[:, :, 1].astype(float)
    b = rgb[:, :, 2].astype(float)
    height, width = r.shape

    fill_mask = (
        (r > 145.0)
        & (r < 240.0)
        & (g > 105.0)
        & (g < 230.0)
        & (b > 95.0)
        & (b < 225.0)
        & ((r - g) > 8.0)
        & (np.abs(g - b) < 28.0)
    )
    line_mask = (
        (r > 70.0)
        & (r < 185.0)
        & (g > 40.0)
        & (g < 135.0)
        & (b > 35.0)
        & (b < 130.0)
        & ((r - g) > 20.0)
        & ((g - b) > -5.0)
        & ((g - b) < 35.0)
    )

    fill_mask = keep_components(fill_mask, min_area=200)
    line_mask = keep_components(line_mask, min_area=100)
    save_mask(fill_mask, fill_mask_output)
    save_mask(line_mask, line_mask_output)

    shade_pixels = np.linspace(0.0, float(width), width + 1)
    shade_x_kpc = pixel_x_to_kpc(shade_pixels, width)
    upper_y_raw = envelope_from_mask(fill_mask, choose="upper")
    lower_y_raw = envelope_from_mask(fill_mask, choose="lower")
    upper_y = sample_column_envelope(upper_y_raw, shade_pixels, width)
    lower_y = sample_column_envelope(lower_y_raw, shade_pixels, width)
    median_pixel, median_y = trace_thick_line(line_mask)
    median_x_kpc = pixel_x_to_kpc(median_pixel, width)

    frame = pd.DataFrame(
        {
            "x_kpc": shade_x_kpc,
            "hayashi_rho_upper_msun_kpc3": pixel_y_to_density(upper_y, height),
            "hayashi_rho_lower_msun_kpc3": pixel_y_to_density(lower_y, height),
            "median_x_kpc": np.nan,
            "hayashi_rho_median_msun_kpc3": np.nan,
        }
    ).replace([np.inf, -np.inf], np.nan)
    frame.loc[: len(median_x_kpc) - 1, "median_x_kpc"] = median_x_kpc
    frame.loc[: len(median_y) - 1, "hayashi_rho_median_msun_kpc3"] = pixel_y_to_density(median_y, height)

    required = ["hayashi_rho_upper_msun_kpc3", "hayashi_rho_lower_msun_kpc3"]
    if frame[required].isna().any().any():
        raise RuntimeError("digitized shaded region has missing upper/lower bounds")
    if not np.all(frame["hayashi_rho_upper_msun_kpc3"] >= frame["hayashi_rho_lower_msun_kpc3"]):
        raise RuntimeError("digitized shaded region has upper bound below lower bound")

    digitized_output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(digitized_output, index=False)
    print(f"wrote {crop_output}")
    print(f"wrote {fill_mask_output}")
    print(f"wrote {line_mask_output}")
    print(f"wrote {digitized_output}")
    print(f"digitized shaded bounds for {len(frame)} x-samples including both x-axis limits; missing upper/lower=0")
    print(f"digitized thick-line columns from left endpoint to stop: {int(frame['hayashi_rho_median_msun_kpc3'].notna().sum())}")
    return frame


def envelope_from_mask(mask: np.ndarray, *, choose: str) -> np.ndarray:
    height, width = mask.shape
    values = np.full(width, np.nan)
    for col in range(width):
        rows = np.flatnonzero(mask[:, col])
        if len(rows) == 0:
            continue
        if choose == "upper":
            values[col] = rows.min() + 0.5
        elif choose == "lower":
            values[col] = rows.max() + 0.5
        elif choose == "center":
            values[col] = np.median(rows) + 0.5
        else:
            raise ValueError(f"unknown choose={choose!r}")
    return values


def sample_column_envelope(values: np.ndarray, sample_pixels: np.ndarray, width: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    columns = np.arange(width, dtype=float)
    valid = np.isfinite(values)
    if valid.sum() < 2:
        raise RuntimeError("not enough shaded-region columns to sample boundaries")
    sampled = np.interp(np.clip(sample_pixels, 0.0, width - 1.0), columns[valid], values[valid])
    left = sample_pixels < columns[valid][0]
    right = sample_pixels > columns[valid][-1]
    if left.any():
        sampled[left] = extrapolate_endpoint(columns[valid], values[valid], sample_pixels[left], side="left")
    if right.any():
        sampled[right] = extrapolate_endpoint(columns[valid], values[valid], sample_pixels[right], side="right")
    return sampled


def extrapolate_endpoint(x: np.ndarray, y: np.ndarray, xp: np.ndarray, *, side: str) -> np.ndarray:
    n = min(20, len(x))
    if side == "left":
        xs = x[:n]
        ys = y[:n]
    elif side == "right":
        xs = x[-n:]
        ys = y[-n:]
    else:
        raise ValueError(f"unknown side={side!r}")
    slope, intercept = np.polyfit(xs, ys, deg=1)
    return slope * xp + intercept


def trace_thick_line(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    height, width = mask.shape
    centers = envelope_from_mask(mask, choose="center")
    finite_columns = np.flatnonzero(np.isfinite(centers))
    if len(finite_columns) < 2:
        raise RuntimeError("not enough thick-line pixels")
    start = int(finite_columns[0])
    tracked_columns = []
    tracked_y = []
    last_y = centers[start]
    for column in range(start, width):
        rows = np.flatnonzero(mask[:, column])
        if len(rows) == 0:
            tracked_columns.append(column)
            tracked_y.append(np.nan)
            continue
        y = float(rows[np.argmin(np.abs(rows - last_y))]) + 0.5
        tracked_columns.append(column)
        tracked_y.append(y)
        last_y = y
        if y >= height - 1.5:
            break
    columns = np.asarray(tracked_columns, dtype=float)
    y_values = np.asarray(tracked_y, dtype=float)
    missing = ~np.isfinite(y_values)
    if missing.any():
        valid = np.isfinite(y_values)
        y_values[missing] = np.interp(columns[missing], columns[valid], y_values[valid])
    valid_stop = np.flatnonzero(np.isfinite(y_values))
    if len(valid_stop) == 0:
        raise RuntimeError("thick-line tracking failed")
    return columns, y_values


def keep_components(mask: np.ndarray, *, min_area: int) -> np.ndarray:
    labels = measure.label(mask, connectivity=2)
    kept = np.zeros_like(mask, dtype=bool)
    for region in measure.regionprops(labels):
        if region.area >= min_area:
            kept[labels == region.label] = True
    return kept


def save_mask(mask: np.ndarray, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray((mask * 255).astype(np.uint8)).save(output)


def plot_density_profile_with_hayashi_reference(galaxy, chain: pd.DataFrame, hayashi_reference: pd.DataFrame, output: Path, *, burn: int) -> dict[str, np.ndarray]:
    radius_kpc = np.geomspace(*FIG1_XLIM_KPC, 300)
    radius_pc = radius_kpc * 1000.0
    required = ["q_halo", "log10_b_halo_pc", "log10_rho0_msun_pc3", "alpha", "beta", "gamma"]
    samples = chain.replace([np.inf, -np.inf], np.nan).dropna(subset=required)
    if samples.empty:
        raise ValueError("chain has no finite halo samples")

    profiles = []
    for row in samples.itertuples(index=False):
        halo = GeneralizedHernquistHalo(
            q=float(row.q_halo),
            b_pc=10.0 ** float(row.log10_b_halo_pc),
            rho0_msun_pc3=10.0 ** float(row.log10_rho0_msun_pc3),
            alpha=float(row.alpha),
            beta=float(row.beta),
            gamma=float(row.gamma),
        )
        profiles.append([halo.density(float(r), 0.0) * 1.0e9 for r in radius_pc])
    profiles = np.asarray(profiles)
    p16, median, p84 = np.percentile(profiles, [16.0, 50.0, 84.0], axis=0)

    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=180)
    ax.fill_between(
        hayashi_reference["x_kpc"],
        hayashi_reference["hayashi_rho_lower_msun_kpc3"],
        hayashi_reference["hayashi_rho_upper_msun_kpc3"],
        color="#b98278",
        alpha=0.38,
        lw=0,
        label="Hayashi+2023 68% interval",
        zorder=1,
    )
    ax.plot(
        hayashi_reference["median_x_kpc"],
        hayashi_reference["hayashi_rho_median_msun_kpc3"],
        color="#8b554c",
        lw=2.4,
        label="Hayashi+2023 median",
        zorder=2,
    )
    ax.plot(
        hayashi_reference["x_kpc"],
        hayashi_reference["hayashi_rho_upper_msun_kpc3"],
        color="#8b554c",
        lw=0.9,
        alpha=0.75,
        zorder=2,
    )
    lower_boundary = hayashi_reference["hayashi_rho_lower_msun_kpc3"].to_numpy(float).copy()
    lower_boundary[lower_boundary <= 1.05 * FIG1_YLIM_MSUN_KPC3[0]] = np.nan
    ax.plot(
        hayashi_reference["x_kpc"],
        lower_boundary,
        color="#8b554c",
        lw=0.9,
        alpha=0.75,
        zorder=2,
    )
    ax.fill_between(radius_kpc, p16, p84, color="#4c78a8", alpha=0.28, lw=0, label="This work 68% interval", zorder=2)
    ax.plot(radius_kpc, median, color="#1f5f8b", lw=2.2, label="This work median", zorder=3)
    ax.axvline(galaxy.observables.b_star_pc / 1000.0, color="#333333", lw=1.4, ls="--", label=r"$b_*$", zorder=4)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*FIG1_XLIM_KPC)
    ax.set_ylim(*FIG1_YLIM_MSUN_KPC3)
    ax.set_box_aspect(1.04)
    ax.set_xlabel("Major Axis [kpc]")
    ax.set_ylabel(r"$\rho_{\rm DM}(r)$ [$M_\odot\,{\rm kpc}^{-3}$]")
    ax.set_title(f"{galaxy.observables.galaxy} dark-matter density profile")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(False)
    fig.text(0.12, 0.01, f"Posterior samples after burn={burn}; Hayashi curves digitized from Figure 1.", fontsize=7)
    fig.tight_layout(rect=(0.0, 0.03, 1.0, 1.0))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    print(f"wrote {output}")
    print(f"used {len(samples)} finite samples after burn={burn}")
    return {"radius_kpc": radius_kpc, "p16": p16, "median": median, "p84": p84}


def compare_profiles(galaxy, profile_summary: dict[str, np.ndarray], hayashi: pd.DataFrame, *, n_samples: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    radius_kpc = profile_summary["radius_kpc"]
    p16 = profile_summary["p16"]
    median = profile_summary["median"]
    p84 = profile_summary["p84"]
    hayashi = hayashi.sort_values("x_kpc")

    compare_r = np.array([0.01, galaxy.observables.b_star_pc / 1000.0, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0])
    rows = []
    for r in compare_r:
        this16 = interp_log(radius_kpc, p16, r)
        this50 = interp_log(radius_kpc, median, r)
        this84 = interp_log(radius_kpc, p84, r)
        h16 = interp_log(hayashi["x_kpc"], hayashi["hayashi_rho_lower_msun_kpc3"], r)
        h50 = interp_log(hayashi["median_x_kpc"], hayashi["hayashi_rho_median_msun_kpc3"], r)
        h84 = interp_log(hayashi["x_kpc"], hayashi["hayashi_rho_upper_msun_kpc3"], r)
        overlap_values = np.array([this16, this84, h16, h84], dtype=float)
        intervals_overlap = bool(np.all(np.isfinite(overlap_values)) and max(this16, h16) <= min(this84, h84))
        rows.append(
            {
                "radius_kpc": r,
                "this_p16": this16,
                "this_median": this50,
                "this_p84": this84,
                "hayashi_lower": h16,
                "hayashi_median": h50,
                "hayashi_upper": h84,
                "median_ratio_this_over_hayashi": this50 / h50,
                "intervals_overlap": intervals_overlap,
            }
        )
    comparison = pd.DataFrame(rows)

    h16 = interp_log_many(hayashi["x_kpc"], hayashi["hayashi_rho_lower_msun_kpc3"], radius_kpc)
    h50 = interp_log_many(hayashi["median_x_kpc"], hayashi["hayashi_rho_median_msun_kpc3"], radius_kpc)
    h84 = interp_log_many(hayashi["x_kpc"], hayashi["hayashi_rho_upper_msun_kpc3"], radius_kpc)
    finite = np.isfinite(h16) & np.isfinite(h50) & np.isfinite(h84)
    overlap = (np.maximum(p16, h16) <= np.minimum(p84, h84)) & finite
    ratio = median / h50
    summary = pd.DataFrame(
        [
            {
                "n_equal_weight_samples": n_samples,
                "b_star_kpc": galaxy.observables.b_star_pc / 1000.0,
                "common_grid_points": int(finite.sum()),
                "overlap_fraction_on_common_grid": float(overlap.sum() / finite.sum()),
                "median_ratio_min_on_common_grid": float(np.nanmin(ratio[finite])),
                "median_ratio_median_on_common_grid": float(np.nanmedian(ratio[finite])),
                "median_ratio_max_on_common_grid": float(np.nanmax(ratio[finite])),
                "hayashi_shade_missing_bounds": int(hayashi[["hayashi_rho_lower_msun_kpc3", "hayashi_rho_upper_msun_kpc3"]].isna().sum().sum()),
            }
        ]
    )
    return comparison, summary


def interp_log(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray, xp: float) -> float:
    return float(interp_log_many(x, y, np.array([xp]))[0])


def interp_log_many(x: pd.Series | np.ndarray, y: pd.Series | np.ndarray, xp: np.ndarray) -> np.ndarray:
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    good = np.isfinite(x_arr) & np.isfinite(y_arr) & (x_arr > 0.0) & (y_arr > 0.0)
    if good.sum() < 2:
        return np.full_like(xp, np.nan, dtype=float)
    order = np.argsort(x_arr[good])
    x_good = x_arr[good][order]
    y_good = y_arr[good][order]
    return 10.0 ** np.interp(
        np.log10(xp),
        np.log10(x_good),
        np.log10(y_good),
        left=np.nan,
        right=np.nan,
    )


if __name__ == "__main__":
    main()
