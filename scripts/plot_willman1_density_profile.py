#!/usr/bin/env python
"""Plot a Hayashi Figure 1 style dark-matter density profile for Willman 1.

Figure 1 in Hayashi et al. is not a plot of a single Table A1 parameter row.
It is a marginalized density profile: evaluate rho_DM(r) for posterior samples,
then plot the median and 68% interval at each radius. This script implements
that semantics. For a real reproduction, pass an MCMC chain with columns:

q_halo, log10_b_halo_pc, log10_rho0_msun_pc3, alpha, beta, gamma
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo

HALO_SAMPLE_COLUMNS = (
    "q_halo",
    "log10_b_halo_pc",
    "log10_rho0_msun_pc3",
    "alpha",
    "beta",
    "gamma",
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default="data/galaxies/27_Willman_1.csv")
    parser.add_argument("--centers-csv", default="data/processed/galaxy_structural_centers.csv")
    parser.add_argument("--posterior-samples", required=True)
    parser.add_argument("--output", default="outputs/figures/willman1_density_profile.png")
    parser.add_argument("--n-radius", type=int, default=300)
    parser.add_argument("--r-min-kpc", type=float, default=0.01)
    parser.add_argument("--r-max-kpc", type=float, default=20.0)
    parser.add_argument("--rho-min-msun-kpc3", type=float, default=1.0e4)
    parser.add_argument("--rho-max-msun-kpc3", type=float, default=1.0e10)
    parser.add_argument("--axis-box-aspect", type=float, default=1.04)
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    samples = load_halo_samples(args.posterior_samples)

    radius_kpc = np.geomspace(args.r_min_kpc, args.r_max_kpc, args.n_radius)
    radius_pc = radius_kpc * 1000.0
    profiles = density_profiles_from_samples(samples, radius_pc)
    p16, median, p84 = np.percentile(profiles, [16.0, 50.0, 84.0], axis=0)

    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=180)
    ax.fill_between(radius_kpc, p16, p84, color="#4c78a8", alpha=0.28, lw=0, label="68% interval")
    ax.plot(radius_kpc, median, color="#1f5f8b", lw=2.2, label="median")
    ax.axvline(galaxy.observables.b_star_pc / 1000.0, color="#333333", lw=1.4, ls="--", label=r"$b_*$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(args.r_min_kpc, args.r_max_kpc)
    ax.set_ylim(args.rho_min_msun_kpc3, args.rho_max_msun_kpc3)
    ax.set_box_aspect(args.axis_box_aspect)
    ax.set_xlabel("Major Axis [kpc]")
    ax.set_ylabel(r"$\rho_{\rm DM}(r)$ [$M_\odot\,{\rm kpc}^{-3}$]")
    ax.set_title("Willman 1 dark-matter density profile")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, which="both", alpha=0.22, lw=0.6)

    subtitle = (
        f"center=({galaxy.observables.center_ra_deg:.4f}, "
        f"{galaxy.observables.center_dec_deg:.4f}) deg; "
        f"n={len(galaxy.stars)}; posterior samples={len(samples)}"
    )
    fig.text(0.12, 0.01, subtitle, fontsize=8, color="#444444")
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    print(f"wrote {output}")
    print(f"galaxy={galaxy.observables.galaxy}")
    print(f"center_ra_deg={galaxy.observables.center_ra_deg:.6f}")
    print(f"center_dec_deg={galaxy.observables.center_dec_deg:.6f}")
    print(f"n_selected={len(galaxy.stars)}")
    print(f"posterior_samples={len(samples)}")
    rho150 = np.interp(150.0, radius_pc, median)
    print(f"median_rho150_log10_msun_kpc3={np.log10(rho150):.6f}")


def load_halo_samples(path: str | Path) -> pd.DataFrame:
    samples = pd.read_csv(path, comment="#")
    missing = [column for column in HALO_SAMPLE_COLUMNS if column not in samples.columns]
    if missing:
        raise ValueError(f"{path} is missing required posterior sample columns: {missing}")
    samples = samples.loc[:, HALO_SAMPLE_COLUMNS].apply(pd.to_numeric, errors="coerce").dropna()
    if samples.empty:
        raise ValueError(f"{path} contains no valid posterior samples")
    return samples


def density_profiles_from_samples(samples: pd.DataFrame, radius_pc: np.ndarray) -> np.ndarray:
    profiles = np.empty((len(samples), len(radius_pc)), dtype=float)
    for row_index, row in enumerate(samples.itertuples(index=False)):
        halo = GeneralizedHernquistHalo(
            q=float(row.q_halo),
            b_pc=10.0 ** float(row.log10_b_halo_pc),
            rho0_msun_pc3=10.0 ** float(row.log10_rho0_msun_pc3),
            alpha=float(row.alpha),
            beta=float(row.beta),
            gamma=float(row.gamma),
        )
        profiles[row_index] = [halo.density(float(radius), 0.0) * 1.0e9 for radius in radius_pc]
    return profiles


if __name__ == "__main__":
    main()
