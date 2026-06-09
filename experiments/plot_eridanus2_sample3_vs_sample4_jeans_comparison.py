#!/usr/bin/env python
"""Overlay Jeans radial diagnostics for physical and unphysical Eri II samples."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    sample3 = pd.read_csv(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_sample3_physical_jeans_radial_diagnostics.csv")
    sample4 = pd.read_csv(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_sample4_unphysical_jeans_radial_diagnostics.csv")
    z_pc = 50.0
    s3 = sample3[sample3["z_pc"] == z_pc]
    s4 = sample4[sample4["z_pc"] == z_pc]

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.2), sharex=True)
    physical = dict(color="#2ca02c", label="sample 3: physical, beta_z=-0.092")
    unphysical = dict(color="#d62728", label="sample 4: unphysical, beta_z=0.872")

    axes[0, 0].plot(s3["R_pc"], s3["Pz"], **physical)
    axes[0, 0].plot(s4["R_pc"], s4["Pz"], **unphysical)
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylabel(r"$P_z=\nu\,\overline{v_z^2}$")
    axes[0, 0].set_title(r"Vertical pressure at $z=50$ pc")

    axes[0, 1].plot(s3["R_pc"], s3["vphi2_raw"], **physical)
    axes[0, 1].plot(s4["R_pc"], s4["vphi2_raw"], **unphysical)
    axes[0, 1].axhline(0.0, color="black", lw=1.0)
    axes[0, 1].set_ylabel(r"raw $\overline{v_\phi^2}$ [(km/s)$^2$]")
    axes[0, 1].set_title("Azimuthal second moment")

    axes[1, 0].plot(s3["R_pc"], s3["delta_Phi"], **physical)
    axes[1, 0].plot(s4["R_pc"], s4["delta_Phi"], **unphysical)
    axes[1, 0].set_ylabel(r"$\Phi(R,z)-\Phi(1\,{\rm pc},z)$ [(km/s)$^2$]")
    axes[1, 0].set_xlabel("R [pc]")
    axes[1, 0].set_title("Radial potential rise")

    for data, color, prefix in [(s3, physical["color"], "sample 3"), (s4, unphysical["color"], "sample 4")]:
        axes[1, 1].plot(data["R_pc"], data["b_vz2"], color=color, linestyle="-", label=f"{prefix}: b vz2")
        axes[1, 1].plot(data["R_pc"], data["b_R_over_nu_dPz_dR"], color=color, linestyle="--", label=f"{prefix}: b R/nu dPz/dR")
        axes[1, 1].plot(data["R_pc"], data["R_dPhi_dR"], color=color, linestyle=":", label=f"{prefix}: R dPhi/dR")
    axes[1, 1].axhline(0.0, color="black", lw=1.0)
    axes[1, 1].set_ylabel(r"terms in $\overline{v_\phi^2}$ [(km/s)$^2$]")
    axes[1, 1].set_xlabel("R [pc]")
    axes[1, 1].set_title("Term budget")

    for ax in axes.ravel():
        ax.set_xscale("log")
        ax.grid(alpha=0.25)

    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[1, 1].legend(frameon=False, fontsize=7, ncol=2)
    fig.suptitle("Eridanus II raw Jeans comparison at z = 50 pc", y=0.99)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))

    output = PROJECT_ROOT / "outputs/figures/eridanus_ii_sample3_vs_sample4_jeans_comparison_z50.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    print(output)


if __name__ == "__main__":
    main()
