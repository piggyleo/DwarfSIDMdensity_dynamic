#!/usr/bin/env python
"""Plot raw Jeans physicality diagnostics for selected Eridanus II samples."""

from __future__ import annotations

import sys
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.halos import GeneralizedHernquistHalo
from hayashi_jeans.projection import AxisymmetricJeansProjector
from scripts.run_willman1_block_mh_fast import beta_z_from_q


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection-id", type=int, default=4)
    parser.add_argument("--output-suffix", default=None)
    args = parser.parse_args()

    selected = pd.read_csv(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_vs_strict_selected_samples.csv")
    row = selected.loc[selected["selection_id"].astype(int) == args.selection_id].iloc[0]
    galaxy = load_galaxy_data(
        PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv",
        structural_centers_csv=PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv",
    )

    halo = GeneralizedHernquistHalo(
        q=float(row["q_halo"]),
        b_pc=10.0 ** float(row["log10_b_halo_pc"]),
        rho0_msun_pc3=10.0 ** float(row["log10_rho0_msun_pc3"]),
        alpha=float(row["alpha"]),
        beta=float(row["beta"]),
        gamma=float(row["gamma"]),
    )
    projector = AxisymmetricJeansProjector(
        b_star_pc=galaxy.observables.b_star_pc,
        qprime=galaxy.observables.qprime,
        inclination_rad=np.deg2rad(float(row["i_deg"])),
        beta_z=beta_z_from_q(float(row["minus_log10_one_minus_beta_z"])),
        halo=halo,
        epsrel=3.0e-4,
    )

    radii = np.geomspace(1.0, 2000.0, 180)
    z_slices = [0.0, 50.0, 100.0]
    records: list[dict[str, float]] = []

    for z_pc in z_slices:
        pressure = np.array([projector._vertical_pressure(float(r), z_pc) for r in radii])
        pressure_spline = CubicSpline(np.log(radii), pressure)
        d_pressure_dr = pressure_spline(np.log(radii), 1) / radii
        density = np.array([projector.tracer.density(float(r), z_pc) for r in radii])
        dphi_dr = np.array([halo.potential_gradients(float(r), z_pc)[0] for r in radii])
        delta_phi = cumulative_trapezoid(dphi_dr, radii, initial=0.0)

        vz2 = pressure / density
        b_aniso = 1.0 / (1.0 - projector.beta_z)
        pressure_gradient_term = b_aniso * radii / density * d_pressure_dr
        radial_random_term = b_aniso * vz2
        gravity_term = radii * dphi_dr
        vphi2 = radial_random_term + pressure_gradient_term + gravity_term

        for values in zip(
            radii,
            np.full_like(radii, z_pc),
            pressure,
            d_pressure_dr,
            density,
            vz2,
            radial_random_term,
            pressure_gradient_term,
            gravity_term,
            vphi2,
            dphi_dr,
            delta_phi,
        ):
            records.append(
                dict(
                    R_pc=float(values[0]),
                    z_pc=float(values[1]),
                    Pz=float(values[2]),
                    dPz_dR=float(values[3]),
                    nu=float(values[4]),
                    vz2=float(values[5]),
                    b_vz2=float(values[6]),
                    b_R_over_nu_dPz_dR=float(values[7]),
                    R_dPhi_dR=float(values[8]),
                    vphi2_raw=float(values[9]),
                    dPhi_dR=float(values[10]),
                    delta_Phi=float(values[11]),
                )
            )

    diagnostics = pd.DataFrame.from_records(records)
    suffix = args.output_suffix or f"sample{args.selection_id}"
    output_csv = PROJECT_ROOT / f"outputs/diagnostics/eridanus_ii_{suffix}_jeans_radial_diagnostics.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    diagnostics.to_csv(output_csv, index=False)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), sharex=True)
    colors = {0.0: "#1f77b4", 50.0: "#2ca02c", 100.0: "#d62728"}

    for z_pc in z_slices:
        part = diagnostics[diagnostics["z_pc"] == z_pc]
        label = f"z = {z_pc:.0f} pc"
        color = colors[z_pc]
        axes[0, 0].plot(part["R_pc"], part["Pz"], color=color, label=label)
        axes[0, 1].plot(part["R_pc"], part["vphi2_raw"], color=color, label=label)
        axes[1, 0].plot(part["R_pc"], part["delta_Phi"], color=color, label=label)
        axes[1, 1].plot(part["R_pc"], part["b_vz2"], color=color, linestyle="-", label=f"{label}: b vz2")
        axes[1, 1].plot(part["R_pc"], part["b_R_over_nu_dPz_dR"], color=color, linestyle="--", label=f"{label}: b R/nu dPz/dR")
        axes[1, 1].plot(part["R_pc"], part["R_dPhi_dR"], color=color, linestyle=":", label=f"{label}: R dPhi/dR")

    for ax in axes.ravel():
        ax.axvline(float(row["q_halo"]) * 0.0 + 10.0 ** float(row["log10_b_halo_pc"]), color="0.45", lw=1.0, alpha=0.7)
        ax.axvline(galaxy.observables.b_star_pc, color="0.1", lw=1.0, alpha=0.7)
        ax.set_xscale("log")
        ax.grid(alpha=0.25)

    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylabel(r"$P_z=\nu\,\overline{v_z^2}$")
    axes[0, 0].set_title("Vertical pressure declines with R")

    axes[0, 1].axhline(0.0, color="black", lw=1.0)
    axes[0, 1].set_ylabel(r"raw $\overline{v_\phi^2}$ [(km/s)$^2$]")
    axes[0, 1].set_title("Azimuthal second moment")

    axes[1, 0].set_ylabel(r"$\Phi(R,z)-\Phi(1\,{\rm pc},z)$ [(km/s)$^2$]")
    axes[1, 0].set_xlabel("R [pc]")
    axes[1, 0].set_title("Radial potential rise flattens outside the compact halo")

    axes[1, 1].axhline(0.0, color="black", lw=1.0)
    axes[1, 1].set_ylabel(r"terms in $\overline{v_\phi^2}$ [(km/s)$^2$]")
    axes[1, 1].set_xlabel("R [pc]")
    axes[1, 1].set_title("Term budget")

    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[1, 1].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle(
        f"Eridanus II selection_id={args.selection_id} raw Jeans diagnostics "
        + rf"($\beta_z={projector.beta_z:.3f}$, $q_h={float(row['q_halo']):.2f}$)",
        y=0.99,
    )
    fig.text(
        0.5,
        0.01,
        f"Vertical markers: halo scale b_h={10.0 ** float(row['log10_b_halo_pc']):.1f} pc; stellar scale b_*={galaxy.observables.b_star_pc:.1f} pc.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.965))

    output_png = PROJECT_ROOT / f"outputs/figures/eridanus_ii_{suffix}_jeans_physicality_diagnostics.png"
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=180)
    print(output_png)
    print(output_csv)


if __name__ == "__main__":
    main()
