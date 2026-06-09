#!/usr/bin/env python
"""Plot Willman 1 density profile with digitized Hayashi Figure 1 curves.

This script is intentionally separate from `run_willman1_block_mh_fast.py`.
It copies the density-profile plotting logic and digitizes the gray Willman 1
median line and shaded-region envelope from Hayashi et al. (2023) Figure 1.
"""

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


DEFAULT_PDF = "/Users/wangkaihao/Zotero/storage/M8AKR2BP/Hayashi 等 - 2023 - Dark Matter Halo Properties of the Galactic Dwarf Satellites Implication for Chemo-dynamical Evolut.pdf"
DEFAULT_PAGE_RENDER = PROJECT_ROOT / "outputs/figures/hayashi_fig1_page_render.png"
# Interior data rectangle of the Willman 1 panel in `hayashi_fig1_page_render.png`.
# The detected panel frame is approximately x=1108..1415, y=449..768.
DEFAULT_CROP_BOX = "1111,453,1411,765"
FIG1_XLIM_KPC = (0.01, 20.0)
FIG1_YLIM_MSUN_KPC3 = (1.0e4, 1.0e10)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/27_Willman_1.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--chain-output", default=str(PROJECT_ROOT / "outputs/willman1_block_mh_fast_chain.csv"))
    parser.add_argument("--figure-output", default=str(PROJECT_ROOT / "outputs/figures/willman1_density_profile_with_hayashi_overlay.png"))
    parser.add_argument("--burn", type=int, default=0)
    parser.add_argument("--hayashi-page-render", default=str(DEFAULT_PAGE_RENDER))
    parser.add_argument("--hayashi-pdf", default=DEFAULT_PDF)
    parser.add_argument("--pdf-page-index", type=int, default=4, help="0-based PDF page index used if --hayashi-page-render is absent.")
    parser.add_argument("--crop-box", default=DEFAULT_CROP_BOX, help="Pixel crop box x0,y0,x1,y1 for the Willman 1 panel.")
    parser.add_argument("--digitized-output", default=str(PROJECT_ROOT / "outputs/diagnostics/hayashi_fig1_willman1_digitized_density.csv"))
    parser.add_argument("--crop-output", default=str(PROJECT_ROOT / "outputs/diagnostics/hayashi_fig1_willman1_crop.png"))
    args = parser.parse_args()

    if args.burn < 0:
        raise ValueError("--burn must be non-negative")

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    chain = pd.read_csv(args.chain_output)
    samples = chain_after_burn(chain, args.burn)
    hayashi_reference = digitize_hayashi_reference(
        page_render=Path(args.hayashi_page_render),
        pdf_path=Path(args.hayashi_pdf),
        pdf_page_index=args.pdf_page_index,
        crop_box=parse_crop_box(args.crop_box),
        crop_output=Path(args.crop_output),
        digitized_output=Path(args.digitized_output),
    )
    plot_density_profile_with_hayashi_reference(galaxy, samples, hayashi_reference, Path(args.figure_output), args.burn)


def chain_after_burn(chain: pd.DataFrame, burn: int) -> pd.DataFrame:
    samples = chain.iloc[burn:].copy()
    if samples.empty:
        raise ValueError(f"no samples remain after burn={burn}; chain has {len(chain)} rows")
    return samples


def digitize_hayashi_reference(
    *,
    page_render: Path,
    pdf_path: Path,
    pdf_page_index: int,
    crop_box: tuple[int, int, int, int],
    crop_output: Path,
    digitized_output: Path,
) -> pd.DataFrame:
    if not page_render.exists():
        render_pdf_page(pdf_path, pdf_page_index, page_render)

    page = Image.open(page_render).convert("RGB")
    crop = page.crop(crop_box)
    crop_output.parent.mkdir(parents=True, exist_ok=True)
    crop.save(crop_output)

    rgb = np.asarray(crop, dtype=float)
    mean = rgb.mean(axis=2)
    std = rgb.std(axis=2)
    height, width = mean.shape

    fill_mask = (std < 4.0) & (mean > 185.0) & (mean < 216.0)
    fill_mask = remove_plot_artifacts(
        fill_mask,
        preserve_left_boundary=True,
        preserve_top_boundary=True,
        preserve_bottom_boundary=True,
    )
    fill_mask = keep_large_components(fill_mask, min_area=250)

    line_mask = (std < 7.0) & (mean > 90.0) & (mean < 172.0)
    line_mask = remove_plot_artifacts(
        line_mask,
        preserve_left_boundary=True,
        preserve_top_boundary=True,
        preserve_bottom_boundary=True,
    )
    line_mask = keep_large_components(line_mask, min_area=200)

    columns = np.arange(width)
    x_kpc = pixel_x_to_kpc(columns + 0.5, width)
    upper = envelope_from_mask(fill_mask, choose="upper")
    lower = envelope_from_mask(fill_mask, choose="lower")
    median = envelope_from_mask(line_mask, choose="center")

    frame = pd.DataFrame(
        {
            "x_kpc": x_kpc,
            "hayashi_rho_upper_msun_kpc3": pixel_y_to_density(upper, height),
            "hayashi_rho_lower_msun_kpc3": pixel_y_to_density(lower, height),
            "hayashi_rho_median_msun_kpc3": pixel_y_to_density(median, height),
        }
    )
    frame = frame.replace([np.inf, -np.inf], np.nan)
    for column in [
        "hayashi_rho_upper_msun_kpc3",
        "hayashi_rho_lower_msun_kpc3",
        "hayashi_rho_median_msun_kpc3",
    ]:
        frame[column] = frame[column].interpolate(limit_area="inside")

    digitized_output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(digitized_output, index=False)
    return frame


def remove_plot_artifacts(
    mask: np.ndarray,
    *,
    preserve_left_boundary: bool,
    preserve_top_boundary: bool,
    preserve_bottom_boundary: bool,
) -> np.ndarray:
    cleaned = mask.copy()
    height, width = cleaned.shape
    if not preserve_top_boundary:
        cleaned[:18, :] = False
    if not preserve_bottom_boundary:
        cleaned[-18:, :] = False
    if preserve_left_boundary:
        cleaned[:, :1] = False
    else:
        cleaned[:, :20] = False
    cleaned[:, -20:] = False
    cleaned[:58, int(0.54 * width):] = False
    bstar_x = int(round(np.log10(0.028 / FIG1_XLIM_KPC[0]) / np.log10(FIG1_XLIM_KPC[1] / FIG1_XLIM_KPC[0]) * width))
    cleaned[:, max(0, bstar_x - 6): min(width, bstar_x + 7)] = False
    return cleaned


def keep_large_components(mask: np.ndarray, *, min_area: int) -> np.ndarray:
    labels = measure.label(mask, connectivity=2)
    kept = np.zeros_like(mask, dtype=bool)
    for region in measure.regionprops(labels):
        if region.area >= min_area:
            kept[labels == region.label] = True
    return kept


def keep_largest_component(mask: np.ndarray) -> np.ndarray:
    labels = measure.label(mask, connectivity=2)
    regions = measure.regionprops(labels)
    if not regions:
        return np.zeros_like(mask, dtype=bool)
    largest = max(regions, key=lambda region: region.area)
    return labels == largest.label


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


def pixel_x_to_kpc(x_pixel: np.ndarray, width: int) -> np.ndarray:
    log_min = np.log10(FIG1_XLIM_KPC[0])
    log_max = np.log10(FIG1_XLIM_KPC[1])
    return 10.0 ** (log_min + x_pixel / width * (log_max - log_min))


def pixel_y_to_density(y_pixel: np.ndarray, height: int) -> np.ndarray:
    log_min = np.log10(FIG1_YLIM_MSUN_KPC3[0])
    log_max = np.log10(FIG1_YLIM_MSUN_KPC3[1])
    return 10.0 ** (log_max - y_pixel / height * (log_max - log_min))


def render_pdf_page(pdf_path: Path, page_index: int, output: Path) -> None:
    import fitz

    doc = fitz.open(pdf_path)
    page = doc[page_index]
    pixmap = page.get_pixmap(matrix=fitz.Matrix(2.0, 2.0), alpha=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    pixmap.save(output)


def parse_crop_box(text: str) -> tuple[int, int, int, int]:
    values = [int(value.strip()) for value in text.split(",")]
    if len(values) != 4:
        raise ValueError("--crop-box must contain four comma-separated integers: x0,y0,x1,y1")
    x0, y0, x1, y1 = values
    if not (x0 < x1 and y0 < y1):
        raise ValueError("--crop-box must satisfy x0 < x1 and y0 < y1")
    return x0, y0, x1, y1


def plot_density_profile_with_hayashi_reference(galaxy, chain: pd.DataFrame, hayashi_reference: pd.DataFrame, output: Path, burn: int) -> None:
    radius_kpc = np.geomspace(*FIG1_XLIM_KPC, 300)
    radius_pc = radius_kpc * 1000.0
    profiles = []
    for row in chain.itertuples(index=False):
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
        color="0.72",
        alpha=0.55,
        lw=0,
        label="Hayashi+2023 68% interval",
        zorder=1,
    )
    ax.plot(
        hayashi_reference["x_kpc"],
        hayashi_reference["hayashi_rho_median_msun_kpc3"],
        color="0.42",
        lw=2.4,
        label="Hayashi+2023 median",
        zorder=2,
    )
    ax.plot(
        hayashi_reference["x_kpc"],
        hayashi_reference["hayashi_rho_upper_msun_kpc3"],
        color="0.60",
        lw=1.0,
        alpha=0.9,
        zorder=2,
    )
    ax.plot(
        hayashi_reference["x_kpc"],
        hayashi_reference["hayashi_rho_lower_msun_kpc3"],
        color="0.60",
        lw=1.0,
        alpha=0.9,
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
    print(f"used {len(chain)} samples after burn={burn}")


if __name__ == "__main__":
    main()
