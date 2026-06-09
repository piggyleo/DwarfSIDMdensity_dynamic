#!/usr/bin/env python
"""Diagnose corrected MGE/JAM errors for Eridanus II source_row=1032."""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import roots_legendre

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from experiments.test_eridanus2_mge_jeans_second_moment import (
    MGE1D,
    fit_positive_mge,
    generalized_halo_density_unit,
    mge_los_second_moment,
)
from hayashi_jeans.constants import G_PC_MSUN_KMS2
from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.tracer import intrinsic_q_from_projected
from scripts.run_willman1_block_mh_fast import beta_z_from_q


@dataclass(frozen=True)
class MGETracer:
    q: float
    amplitudes: np.ndarray
    sigmas_major_pc: np.ndarray

    def density(self, r_pc: float | np.ndarray, z_pc: float | np.ndarray) -> float | np.ndarray:
        r = np.asarray(r_pc)
        z = np.asarray(z_pc)
        m2 = r**2 + (z / self.q) ** 2
        return np.sum(
            self.amplitudes * np.exp(-0.5 * m2[..., None] / self.sigmas_major_pc[None, :] ** 2),
            axis=-1,
        )


@dataclass(frozen=True)
class MGEHalo:
    q: float
    amplitudes: np.ndarray
    sigmas_major_pc: np.ndarray
    order: int = 160

    def __post_init__(self) -> None:
        nodes, weights = roots_legendre(self.order)
        object.__setattr__(self, "s_nodes", 0.5 * (nodes + 1.0))
        object.__setattr__(self, "s_weights", 0.5 * weights)

    def density_from_m2(self, m2: np.ndarray) -> np.ndarray:
        return np.sum(
            self.amplitudes * np.exp(-0.5 * m2[..., None] / self.sigmas_major_pc[None, :] ** 2),
            axis=-1,
        )

    def potential_gradients(self, r_pc: float, z_pc: float) -> tuple[float, float]:
        r = abs(float(r_pc))
        z = float(z_pc)
        if r == 0.0 and z == 0.0:
            return 0.0, 0.0
        s = self.s_nodes
        s2 = s * s
        a_s = 1.0 + (self.q**2 - 1.0) * s2
        m2 = s2 * (r**2 + z**2 / a_s)
        rho = self.density_from_m2(m2)
        prefactor = 2.0 * np.pi * G_PC_MSUN_KMS2 * self.q
        int_r = float(np.sum(self.s_weights * 2.0 * rho * s2 / np.sqrt(a_s)))
        int_z = float(np.sum(self.s_weights * 2.0 * rho * s2 / (a_s**1.5)))
        return prefactor * r * int_r, prefactor * z * int_z


class DirectMGEProjector:
    def __init__(self, *, tracer: MGETracer, halo: MGEHalo, beta_z: float, inclination_rad: float, b_star_pc: float) -> None:
        self.tracer = tracer
        self.halo = halo
        self.beta_z = beta_z
        self.inclination_rad = inclination_rad
        self.zmax_pc = 20.0 * b_star_pc
        self.los_max_pc = 20.0 * b_star_pc

    def vertical_pressure(self, r_pc: float, z_pc: float) -> float:
        sign = 1.0 if z_pc >= 0.0 else -1.0
        z_abs = abs(z_pc)

        def integrand(zp_abs: float) -> float:
            _, dphi_dz = self.halo.potential_gradients(r_pc, sign * zp_abs)
            return self.tracer.density(r_pc, sign * zp_abs) * sign * dphi_dz

        return quad(integrand, z_abs, self.zmax_pc, epsrel=1.5e-2, limit=100)[0]

    def d_pressure_dr(self, r_pc: float, z_pc: float) -> float:
        step = max(0.01 * 196.0, 0.01 * abs(r_pc), 1.0e-2)
        if r_pc > step:
            return (self.vertical_pressure(r_pc + step, z_pc) - self.vertical_pressure(r_pc - step, z_pc)) / (2.0 * step)
        return (self.vertical_pressure(r_pc + step, z_pc) - self.vertical_pressure(r_pc, z_pc)) / step

    def moments(self, r_pc: float, z_pc: float) -> tuple[float, float, float]:
        nu = self.tracer.density(r_pc, z_pc)
        vz2 = self.vertical_pressure(r_pc, z_pc) / nu
        vr2 = vz2 / (1.0 - self.beta_z)
        dphi_dr, _ = self.halo.potential_gradients(r_pc, z_pc)
        vphi2 = (vz2 + r_pc / nu * self.d_pressure_dr(r_pc, z_pc)) / (1.0 - self.beta_z) + r_pc * dphi_dr
        return float(vr2), float(vphi2), float(vz2)

    def sigma_los2(self, x_pc: float, y_pc: float) -> tuple[float, float, float]:
        i = self.inclination_rad

        def geometry(ell_pc: float) -> tuple[float, float]:
            r = np.sqrt(x_pc**2 + (y_pc * np.cos(i) + ell_pc * np.sin(i)) ** 2)
            z = y_pc * np.sin(i) - ell_pc * np.cos(i)
            return float(r), float(z)

        def surface_integrand(ell_pc: float) -> float:
            r, z = geometry(ell_pc)
            return float(self.tracer.density(r, z))

        def moment_integrand(ell_pc: float) -> float:
            r, z = geometry(ell_pc)
            nu = float(self.tracer.density(r, z))
            vr2, vphi2, vz2 = self.moments(r, z)
            if r > 1.0e-10:
                x_over_r2 = x_pc**2 / r**2
                v_star2 = vphi2 * x_over_r2 + vr2 * (1.0 - x_over_r2)
            else:
                v_star2 = vr2
            v_los2 = v_star2 * np.sin(i) ** 2 + vz2 * np.cos(i) ** 2
            return nu * v_los2

        denom = quad(surface_integrand, -self.los_max_pc, self.los_max_pc, epsrel=1.5e-2, limit=100)[0]
        numer = quad(moment_integrand, -self.los_max_pc, self.los_max_pc, epsrel=1.5e-2, limit=100)[0]
        return float(numer / denom), float(numer), float(denom)


def main() -> None:
    galaxy = load_galaxy_data(
        PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv",
        structural_centers_csv=PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv",
    )
    chain = pd.read_csv(PROJECT_ROOT / "outputs/eridanus_ii_nautilus_5w_chain.csv")
    row = chain.iloc[1032]
    inclination_rad = np.deg2rad(float(row["i_deg"]))
    q_star = intrinsic_q_from_projected(galaxy.observables.qprime, inclination_rad)
    beta_z = beta_z_from_q(float(row["minus_log10_one_minus_beta_z"]))
    density_scale = 10.0 ** float(row["log10_rho0_msun_pc3"])

    tracer_mge = fit_positive_mge(
        lambda m: (1.0 + (m / galaxy.observables.b_star_pc) ** 2) ** (-2.5),
        n_gauss=45,
        n_fit_radii=800,
        r_min_pc=1.0e-2,
        r_max_pc=2.0e5,
    )
    halo_mge = fit_positive_mge(
        lambda m: generalized_halo_density_unit(
            m,
            b_pc=10.0 ** float(row["log10_b_halo_pc"]),
            alpha=float(row["alpha"]),
            beta=float(row["beta"]),
            gamma=float(row["gamma"]),
        ),
        n_gauss=45,
        n_fit_radii=800,
        r_min_pc=1.0e-2,
        r_max_pc=2.0e5,
    )

    stars = pd.read_csv(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_source1032_1596_mge_error_top_stars.csv")
    stars = stars[stars["source_row"] == 1032]
    target = stars[stars["star_index"].isin([72, 91])].copy()

    closed_sigma_unit, projected_tracer, _, _ = mge_los_second_moment(
        x_pc=target["x_pc"].to_numpy(float),
        y_pc=target["y_pc"].to_numpy(float),
        halo_mge=halo_mge,
        tracer_mge=tracer_mge,
        q_halo=float(row["q_halo"]),
        q_star=q_star,
        beta_z=beta_z,
        inclination_rad=inclination_rad,
        n_u=256,
    )

    direct = DirectMGEProjector(
        tracer=MGETracer(q_star, tracer_mge.amplitudes, tracer_mge.sigmas_major_pc),
        halo=MGEHalo(float(row["q_halo"]), halo_mge.amplitudes, halo_mge.sigmas_major_pc),
        beta_z=beta_z,
        inclination_rad=inclination_rad,
        b_star_pc=galaxy.observables.b_star_pc,
    )

    rows = []
    t0 = time.perf_counter()
    for index, star in enumerate(target.itertuples(index=False)):
        direct_unit, direct_num, direct_den = direct.sigma_los2(float(star.x_pc), float(star.y_pc))
        closed_num = closed_sigma_unit[index] * projected_tracer[index]
        rows.append(
            {
                "source_row": 1032,
                "star_index": int(star.star_index),
                "x_pc": float(star.x_pc),
                "y_pc": float(star.y_pc),
                "strict_sigma2": float(star.strict_unit_sigma2 * density_scale),
                "closed_mge_sigma2": float(closed_sigma_unit[index] * density_scale),
                "direct_mge_sigma2": float(direct_unit * density_scale),
                "closed_vs_direct_rel": float(abs(closed_sigma_unit[index] - direct_unit) / abs(direct_unit)),
                "closed_projected_tracer": float(projected_tracer[index]),
                "direct_projected_tracer": float(direct_den),
                "projected_tracer_rel_delta": float((projected_tracer[index] - direct_den) / direct_den),
                "closed_numerator": float(closed_num),
                "direct_numerator": float(direct_num),
                "numerator_rel_delta": float((closed_num - direct_num) / direct_num),
            }
        )
    print(f"direct seconds={time.perf_counter() - t0:.3f}")

    output = PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_source1032_stars72_91_closed_vs_direct_mge.csv"
    pd.DataFrame(rows).to_csv(output, index=False)
    print(output)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
