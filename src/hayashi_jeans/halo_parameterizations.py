"""Parameterizations that resolve to runtime dark-matter halo objects."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.cosmology import Planck15
from astropy import units as u

from .halos import SIDMPSIDM25Halo


@dataclass(frozen=True)
class HaloRunContext:
    redshift: float | None = None
    cosmology_name: str = "planck15"

    def require_redshift(self) -> float:
        if self.redshift is None or not np.isfinite(self.redshift) or self.redshift < 0.0:
            raise ValueError("this halo parameterization requires a fixed non-negative halo redshift")
        if self.cosmology_name.lower() != "planck15":
            raise ValueError("only the planck15 cosmology is currently supported")
        return float(self.redshift)


@dataclass(frozen=True)
class SIDMPhysicalParameters:
    q_halo: float
    rs0_pc: float
    rho_s0_msun_pc3: float
    tau: float
    m200_msun: float | None = None
    c200: float | None = None
    redshift: float | None = None
    concentration_relation: str | None = None
    concentration_scatter_sigma: float | None = None

    def build_halo(self, *, density_scale: float | None = None, **halo_kwargs) -> SIDMPSIDM25Halo:
        return SIDMPSIDM25Halo(
            q=self.q_halo,
            rs0_pc=self.rs0_pc,
            rho_s0_msun_pc3=(
                self.rho_s0_msun_pc3 if density_scale is None else float(density_scale)
            ),
            tau=self.tau,
            **halo_kwargs,
        )


def sidm_from_scale_parameters(
    *,
    q_halo: float,
    rs0_pc: float,
    rho_s0_msun_pc3: float,
    tau: float,
) -> SIDMPhysicalParameters:
    return SIDMPhysicalParameters(
        q_halo=float(q_halo),
        rs0_pc=float(rs0_pc),
        rho_s0_msun_pc3=float(rho_s0_msun_pc3),
        tau=float(tau),
    )


def nfw_scale_parameters_from_m200_c200(
    *,
    m200_msun: float,
    c200: float,
    context: HaloRunContext,
) -> tuple[float, float]:
    redshift = context.require_redshift()
    mass = float(m200_msun)
    concentration = float(c200)
    if mass <= 0.0 or concentration <= 0.0:
        raise ValueError("m200_msun and c200 must be positive")
    rho_crit_msun_pc3 = Planck15.critical_density(redshift).to(u.Msun / u.pc**3).value
    r200_pc = (3.0 * mass / (4.0 * np.pi * 200.0 * rho_crit_msun_pc3)) ** (1.0 / 3.0)
    rs0_pc = r200_pc / concentration
    nfw_mass_factor = np.log1p(concentration) - concentration / (1.0 + concentration)
    rho_s0 = 200.0 / 3.0 * rho_crit_msun_pc3 * concentration**3 / nfw_mass_factor
    return float(rs0_pc), float(rho_s0)


def sidm_from_m200_c200(
    *,
    q_halo: float,
    m200_msun: float,
    c200: float,
    tau: float,
    context: HaloRunContext,
) -> SIDMPhysicalParameters:
    rs0_pc, rho_s0 = nfw_scale_parameters_from_m200_c200(
        m200_msun=m200_msun,
        c200=c200,
        context=context,
    )
    return SIDMPhysicalParameters(
        q_halo=float(q_halo),
        rs0_pc=rs0_pc,
        rho_s0_msun_pc3=rho_s0,
        tau=float(tau),
        m200_msun=float(m200_msun),
        c200=float(c200),
        redshift=context.require_redshift(),
    )


def ludlow16_concentration(
    *,
    m200_msun: float,
    context: HaloRunContext,
    scatter_sigma: float = 0.0,
) -> float:
    redshift = context.require_redshift()
    try:
        from colossus.cosmology import cosmology as col_cosmology
        from colossus.halo.concentration import concentration as col_concentration
    except ImportError as exc:
        raise ImportError("colossus is required for the Ludlow16 parameterization") from exc

    col_cosmo = col_cosmology.setCosmology("planck15")
    concentration = float(
        col_concentration(float(m200_msun) * col_cosmo.h, "200c", redshift, model="ludlow16")
    )
    if not np.isfinite(concentration) or concentration <= 0.0:
        raise ValueError("Ludlow16 returned an invalid concentration")
    return float(10.0 ** (np.log10(concentration) + float(scatter_sigma) * 0.15))


def sidm_from_m200_ludlow(
    *,
    q_halo: float,
    m200_msun: float,
    tau: float,
    context: HaloRunContext,
    scatter_sigma: float = 0.0,
) -> SIDMPhysicalParameters:
    c200 = ludlow16_concentration(
        m200_msun=m200_msun,
        context=context,
        scatter_sigma=scatter_sigma,
    )
    resolved = sidm_from_m200_c200(
        q_halo=q_halo,
        m200_msun=m200_msun,
        c200=c200,
        tau=tau,
        context=context,
    )
    return SIDMPhysicalParameters(
        **{
            **resolved.__dict__,
            "concentration_relation": "ludlow16",
            "concentration_scatter_sigma": float(scatter_sigma),
        }
    )
