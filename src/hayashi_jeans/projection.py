"""Axisymmetric Jeans solution and LOS projection."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

from .halos import HaloModel
from .tracer import AxisymmetricPlummerTracer, intrinsic_q_from_projected


@dataclass(frozen=True)
class JeansMoments:
    v_r2: float
    v_phi2: float
    v_z2: float


class AxisymmetricJeansProjector:
    """Compute Hayashi-style sigma_los^2 for individual projected stars."""

    def __init__(
        self,
        *,
        b_star_pc: float,
        qprime: float,
        inclination_rad: float,
        beta_z: float,
        halo: HaloModel,
        zmax_factor: float = 80.0,
        los_factor: float = 80.0,
        epsrel: float = 3e-4,
    ) -> None:
        self.b_star_pc = float(b_star_pc)
        self.qprime = float(qprime)
        self.inclination_rad = float(inclination_rad)
        self.beta_z = float(beta_z)
        self.halo = halo
        self.zmax_pc = zmax_factor * self.b_star_pc
        self.los_max_pc = los_factor * self.b_star_pc
        self.epsrel = float(epsrel)
        self.tracer = AxisymmetricPlummerTracer(
            b_star_pc=self.b_star_pc,
            q_intrinsic=intrinsic_q_from_projected(qprime, inclination_rad),
        )

    def sigma_los2_many(self, x_pc: np.ndarray, y_pc: np.ndarray) -> np.ndarray:
        return np.array([self.sigma_los2(float(x), float(y)) for x, y in zip(x_pc, y_pc)], dtype=float)

    def sigma_los2(self, x_pc: float, y_pc: float) -> float:
        i = self.inclination_rad

        def los_geometry(ell_pc: float) -> tuple[float, float]:
            r2 = x_pc**2 + (y_pc * np.cos(i) + ell_pc * np.sin(i)) ** 2
            z = y_pc * np.sin(i) - ell_pc * np.cos(i)
            return np.sqrt(max(r2, 0.0)), z

        def surface_integrand(ell_pc: float) -> float:
            r, z = los_geometry(ell_pc)
            return self.tracer.density(r, z)

        def moment_integrand(ell_pc: float) -> float:
            r, z = los_geometry(ell_pc)
            nu = self.tracer.density(r, z)
            moments = self.jeans_moments(r, z)
            if r > 1e-10:
                v_star2 = moments.v_phi2 * x_pc**2 / r**2 + moments.v_r2 * (1.0 - x_pc**2 / r**2)
            else:
                v_star2 = moments.v_r2
            v_los2 = v_star2 * np.sin(i) ** 2 + moments.v_z2 * np.cos(i) ** 2
            return nu * max(v_los2, 0.0)

        denom = quad(surface_integrand, -self.los_max_pc, self.los_max_pc, epsrel=self.epsrel, limit=100)[0]
        if denom <= 0.0:
            return np.nan
        numer = quad(moment_integrand, -self.los_max_pc, self.los_max_pc, epsrel=self.epsrel, limit=100)[0]
        return max(float(numer / denom), 1e-10)

    def jeans_moments(self, r_cyl_pc: float, z_pc: float) -> JeansMoments:
        nu = self.tracer.density(r_cyl_pc, z_pc)
        if nu <= 0.0:
            return JeansMoments(0.0, 0.0, 0.0)

        v_z2 = self._vertical_pressure(r_cyl_pc, z_pc) / nu
        v_r2 = v_z2 / (1.0 - self.beta_z)
        d_pressure_dr = self._d_vertical_pressure_dr(r_cyl_pc, z_pc)
        dphi_dr, _ = self.halo.potential_gradients(r_cyl_pc, z_pc)
        v_phi2 = (
            (v_z2 + r_cyl_pc / nu * d_pressure_dr) / (1.0 - self.beta_z)
            + r_cyl_pc * dphi_dr
        )
        return JeansMoments(max(v_r2, 0.0), max(v_phi2, 0.0), max(v_z2, 0.0))

    def _vertical_pressure(self, r_cyl_pc: float, z_pc: float) -> float:
        sign = 1.0 if z_pc >= 0.0 else -1.0
        z_abs = abs(z_pc)

        def integrand(zp_abs: float) -> float:
            _, dphi_dz = self.halo.potential_gradients(r_cyl_pc, sign * zp_abs)
            return self.tracer.density(r_cyl_pc, sign * zp_abs) * sign * dphi_dz

        return quad(integrand, z_abs, self.zmax_pc, epsrel=self.epsrel, limit=100)[0]

    def _d_vertical_pressure_dr(self, r_cyl_pc: float, z_pc: float) -> float:
        step = max(0.01 * self.b_star_pc, 0.01 * abs(r_cyl_pc), 1e-2)
        if r_cyl_pc > step:
            p_plus = self._vertical_pressure(r_cyl_pc + step, z_pc)
            p_minus = self._vertical_pressure(r_cyl_pc - step, z_pc)
            return (p_plus - p_minus) / (2.0 * step)
        p0 = self._vertical_pressure(r_cyl_pc, z_pc)
        p_plus = self._vertical_pressure(r_cyl_pc + step, z_pc)
        return (p_plus - p0) / step

