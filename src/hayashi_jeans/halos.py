"""Pluggable spheroidally stratified dark-matter halo models."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from .constants import G_PC_MSUN_KMS2


@lru_cache(maxsize=8)
def _unit_interval_legendre_nodes(order: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = roots_legendre(order)
    return 0.5 * (nodes + 1.0), 0.5 * weights


class HaloModel(ABC):
    """Interface expected by the Jeans solver."""

    @abstractmethod
    def density(self, r_cyl_pc: float, z_pc: float) -> float:
        """Return density in Msun / pc^3."""

    @abstractmethod
    def potential_gradients(self, r_cyl_pc: float, z_pc: float) -> tuple[float, float]:
        """Return (dPhi/dR, dPhi/dz) in (km/s)^2 / pc."""


class SpheroidallyStratifiedHalo(HaloModel):
    """Base class for densities stratified on similar spheroids.

    Subclasses provide ``rho(m)`` for ``m^2 = R^2 + z^2 / q^2``.  The common
    force implementation is the exact homeoidal-density integral already used
    by the generalized Hernquist model.
    """

    def __init__(
        self,
        *,
        q: float,
        integration_eps: float = 1e-5,
        force_integral_method: str = "unit_interval",
        force_quadrature_order: int = 128,
    ) -> None:
        self.q = float(q)
        if not np.isfinite(self.q) or self.q <= 0.0:
            raise ValueError("q must be positive and finite")
        self.integration_eps = float(integration_eps)
        if force_integral_method not in {"infinite", "unit_interval"}:
            raise ValueError("force_integral_method must be 'infinite' or 'unit_interval'")
        self.force_integral_method = force_integral_method
        self.force_quadrature_order = int(force_quadrature_order)
        if self.force_quadrature_order < 16:
            raise ValueError("force_quadrature_order must be >= 16")

    @abstractmethod
    def density_at_ellipsoidal_radius(self, m_pc):
        """Return density at ellipsoidal radius ``m_pc`` in Msun / pc^3."""

    @abstractmethod
    def unit_density_profile(self, m_pc):
        """Return the profile with the model's linear density scale set to one."""

    def density(self, r_cyl_pc: float, z_pc: float) -> float:
        m = np.sqrt(float(r_cyl_pc) ** 2 + (float(z_pc) / self.q) ** 2)
        return float(self.density_at_ellipsoidal_radius(m))

    def potential_gradients(self, r_cyl_pc: float, z_pc: float) -> tuple[float, float]:
        r = abs(float(r_cyl_pc))
        z = float(z_pc)
        if r == 0.0 and z == 0.0:
            return 0.0, 0.0

        prefactor = 2.0 * np.pi * G_PC_MSUN_KMS2 * self.q

        def rho_at_tau(tau: float) -> float:
            m2 = r**2 / (1.0 + tau) + z**2 / (self.q**2 + tau)
            return float(self.density_at_ellipsoidal_radius(np.sqrt(max(m2, 0.0))))

        def integrand_r(tau: float) -> float:
            return rho_at_tau(tau) / ((1.0 + tau) ** 2 * np.sqrt(self.q**2 + tau))

        def integrand_z(tau: float) -> float:
            return rho_at_tau(tau) / ((1.0 + tau) * (self.q**2 + tau) ** 1.5)

        if self.force_integral_method == "infinite":
            int_r = quad(integrand_r, 0.0, np.inf, epsrel=self.integration_eps, limit=100)[0]
            int_z = quad(integrand_z, 0.0, np.inf, epsrel=self.integration_eps, limit=100)[0]
            return prefactor * r * int_r, prefactor * z * int_z

        s, weights = _unit_interval_legendre_nodes(self.force_quadrature_order)
        s2 = s * s
        a_s = 1.0 + (self.q**2 - 1.0) * s2
        m = s * np.sqrt(r**2 + z**2 / a_s)
        rho = np.asarray(self.density_at_ellipsoidal_radius(m), dtype=float)
        int_r = float(np.sum(weights * 2.0 * rho * s2 / np.sqrt(a_s)))
        int_z = float(np.sum(weights * 2.0 * rho * s2 / (a_s**1.5)))
        return prefactor * r * int_r, prefactor * z * int_z


class GeneralizedHernquistHalo(SpheroidallyStratifiedHalo):
    """Axisymmetric generalized Hernquist/Zhao halo used by Hayashi et al."""

    def __init__(
        self,
        *,
        q: float,
        b_pc: float,
        rho0_msun_pc3: float,
        alpha: float,
        beta: float,
        gamma: float,
        integration_eps: float = 1e-5,
        force_integral_method: str = "unit_interval",
        force_quadrature_order: int = 128,
    ) -> None:
        self.b_pc = float(b_pc)
        self.rho0_msun_pc3 = float(rho0_msun_pc3)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.gamma = float(gamma)
        if self.b_pc <= 0.0 or self.rho0_msun_pc3 <= 0.0:
            raise ValueError("b_pc and rho0_msun_pc3 must be positive")
        super().__init__(
            q=q,
            integration_eps=integration_eps,
            force_integral_method=force_integral_method,
            force_quadrature_order=force_quadrature_order,
        )

    def unit_density_profile(self, m_pc):
        values = np.asarray(m_pc)
        x = values / self.b_pc
        floor = 1.0e-300 if np.iscomplexobj(x) else 1.0e-12
        x = np.where(np.abs(x) > floor, x, floor + 0.0j if np.iscomplexobj(x) else floor)
        return x ** (-self.gamma) * (1.0 + x**self.alpha) ** (
            -(self.beta - self.gamma) / self.alpha
        )

    def density_at_ellipsoidal_radius(self, m_pc):
        return self.rho0_msun_pc3 * self.unit_density_profile(m_pc)

    def _potential_gradients_unit_interval_quad(self, r_cyl_pc: float, z_pc: float) -> tuple[float, float]:
        """Development-only reference for the compact interval using scipy.quad."""

        r = abs(float(r_cyl_pc))
        z = float(z_pc)
        if r == 0.0 and z == 0.0:
            return 0.0, 0.0

        prefactor = 2.0 * np.pi * G_PC_MSUN_KMS2 * self.q

        def rho_from_m2(m2: float) -> float:
            m = np.sqrt(max(m2, 0.0))
            x = max(m / self.b_pc, 1e-12)
            return self.rho0_msun_pc3 * x ** (-self.gamma) * (
                1.0 + x**self.alpha
            ) ** (-(self.beta - self.gamma) / self.alpha)

        def rho_at_compact_s(s: float) -> tuple[float, float]:
            s2_scalar = s * s
            a_s_scalar = 1.0 + (self.q**2 - 1.0) * s2_scalar
            m2_scalar = s2_scalar * (r**2 + z**2 / a_s_scalar)
            return rho_from_m2(m2_scalar), a_s_scalar

        def integrand_r_unit(s: float) -> float:
            rho, a_s = rho_at_compact_s(s)
            return 2.0 * rho * s * s / np.sqrt(a_s)

        def integrand_z_unit(s: float) -> float:
            rho, a_s = rho_at_compact_s(s)
            return 2.0 * rho * s * s / (a_s**1.5)

        int_r = quad(integrand_r_unit, 0.0, 1.0, epsrel=self.integration_eps, limit=100)[0]
        int_z = quad(integrand_z_unit, 0.0, 1.0, epsrel=self.integration_eps, limit=100)[0]
        return prefactor * r * int_r, prefactor * z * int_z


SIDM_TAU_MAX = 1.08


@dataclass(frozen=True)
class SIDMEvolutionCoefficients:
    density_ratio: float
    scale_radius_ratio: float
    core_radius_ratio: float
    transition_alpha: float
    inner_gamma: float


def sidm_evolution_coefficients(tau: float) -> SIDMEvolutionCoefficients:
    """Return the PSIDM-25 dimensionless profile coefficients."""

    x = float(tau)
    if not np.isfinite(x) or not (0.0 <= x <= SIDM_TAU_MAX):
        raise ValueError(f"tau must be finite and lie in [0, {SIDM_TAU_MAX}]")

    gamma_terms = [2.54741454, -2.75585114, 6.11802594, -7.28902641, 3.41868933, 1.33017751]
    a, b, c, d, g, h = gamma_terms
    inner_gamma = (
        a
        + b * x
        + c * x**2
        + d * x**3
        + g * x**4
        + (1.0 - a) / np.log(0.001) * np.log(x**h + 0.001)
        + 0.1 * (np.arctan(3.0 * (x - 0.5)) - np.arctan(-1.5))
    )

    a, b, c, d, g, h, i = (
        -1.54891697e1,
        7.22099942,
        1.77835229e1,
        -7.21790903e2,
        5.64041228e2,
        -4.52714226e-1,
        -1.24201628,
    )
    density_ratio = np.exp(
        a * x**0.3 + b * x**0.2 + c * x**0.7 + d * x * 0.8 + g * x + h * x**6 + i * x**12
    )

    a, b, c, d, g, h, i = 2.79650252, -0.02224725, -2.68712924, 3.42099585, 0.06969549, 2.02747497, 0.03140837
    scale_radius_ratio = 1.0 + a * x**0.6 + c * x + d * x**3 + g * x**10 + h * x**17 + i * x**80 + b * x**83

    a, b, c, d, g, h, i = -0.41711945, 0.29619283, 0.20237198, 0.93450782, -2.10957393, 2.00795074, -0.57618027
    transition_alpha = 1.0 + a * x**0.6 + c * x**2 + d * x**13 + g * x**26 + h * x**39 + i * x**66 + b * x**71

    a, b, c, d, g, h, i = 6.04642694, -6.95669182, 1.48570521, -0.623443139, 0.209573098, -0.033795816, 5.02603829e-5
    core_radius_ratio = a * x**0.6 + b * x**0.8 + c * x**2 + d * x**5 + g * x**12 + h * x**27 + i * x**80

    return SIDMEvolutionCoefficients(
        density_ratio=float(density_ratio),
        scale_radius_ratio=float(scale_radius_ratio),
        core_radius_ratio=float(core_radius_ratio),
        transition_alpha=float(transition_alpha),
        inner_gamma=float(inner_gamma),
    )


def sidm_dimensionless_density(x, tau: float):
    """PSIDM-25 density divided by the initial NFW scale density."""

    coeff = sidm_evolution_coefficients(tau)
    radius = np.asarray(x)
    complex_input = np.iscomplexobj(radius)
    floor = 1.0e-300
    radius = np.where(np.abs(radius) > floor, radius, floor + 0.0j if complex_input else floor)
    rs = coeff.scale_radius_ratio
    rc = coeff.core_radius_ratio
    alpha = coeff.transition_alpha
    gamma = coeff.inner_gamma
    core = ((radius / rs) ** 4 + (rc / rs) ** 4) ** (gamma / 4.0)
    outer = (1.0 + (radius / rs) ** alpha) ** ((3.0 - gamma) / alpha)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        return coeff.density_ratio / (core * outer)


class SIDMPSIDM25Halo(SpheroidallyStratifiedHalo):
    """Axisymmetric extension of the spherical PSIDM-25 density profile."""

    def __init__(
        self,
        *,
        q: float,
        rs0_pc: float,
        rho_s0_msun_pc3: float,
        tau: float,
        integration_eps: float = 1e-5,
        force_integral_method: str = "unit_interval",
        force_quadrature_order: int = 128,
    ) -> None:
        self.rs0_pc = float(rs0_pc)
        self.rho_s0_msun_pc3 = float(rho_s0_msun_pc3)
        self.tau = float(tau)
        if self.rs0_pc <= 0.0 or self.rho_s0_msun_pc3 <= 0.0:
            raise ValueError("rs0_pc and rho_s0_msun_pc3 must be positive")
        self.evolution = sidm_evolution_coefficients(self.tau)
        super().__init__(
            q=q,
            integration_eps=integration_eps,
            force_integral_method=force_integral_method,
            force_quadrature_order=force_quadrature_order,
        )

    def unit_density_profile(self, m_pc):
        return sidm_dimensionless_density(np.asarray(m_pc) / self.rs0_pc, self.tau)

    def density_at_ellipsoidal_radius(self, m_pc):
        return self.rho_s0_msun_pc3 * self.unit_density_profile(m_pc)
