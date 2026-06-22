"""Axisymmetric Plummer tracer model."""

from __future__ import annotations

import numpy as np


def intrinsic_q_from_projected(qprime: float, inclination_rad: float, *, min_q: float = 1e-3) -> float:
    sin_i = np.sin(inclination_rad)
    if abs(sin_i) < 1e-8:
        raise ValueError("inclination is too close to face-on to infer intrinsic q")
    q2 = (qprime**2 - np.cos(inclination_rad) ** 2) / sin_i**2
    if q2 <= 0:
        raise ValueError("inclination is incompatible with observed qprime")
    return max(float(np.sqrt(q2)), min_q)


class AxisymmetricPlummerTracer:
    """Unnormalized flattened Plummer tracer density."""

    def __init__(self, b_star_pc: float, q_intrinsic: float) -> None:
        self.b_star_pc = float(b_star_pc)
        self.q_intrinsic = float(q_intrinsic)

    def density(self, r_cyl_pc: float | np.ndarray, z_pc: float | np.ndarray) -> float | np.ndarray:
        m2 = np.asarray(r_cyl_pc) ** 2 + (np.asarray(z_pc) / self.q_intrinsic) ** 2
        return (1.0 + m2 / self.b_star_pc**2) ** (-2.5)


class AxisymmetricMGETracer:
    """Axisymmetric tracer represented by a signed intrinsic Gaussian sum."""

    def __init__(self, *, q_intrinsic: float, amplitudes, sigmas_major_pc) -> None:
        self.q_intrinsic = float(q_intrinsic)
        self.amplitudes = np.asarray(amplitudes, dtype=float)
        self.sigmas_major_pc = np.asarray(sigmas_major_pc, dtype=float)
        if not np.isfinite(self.q_intrinsic) or self.q_intrinsic <= 0.0:
            raise ValueError("q_intrinsic must be positive and finite")
        if self.amplitudes.ndim != 1 or self.sigmas_major_pc.ndim != 1:
            raise ValueError("MGE amplitudes and sigmas must be one-dimensional")
        if self.amplitudes.shape != self.sigmas_major_pc.shape:
            raise ValueError("MGE amplitudes and sigmas must have matching shapes")
        if self.amplitudes.size == 0:
            raise ValueError("MGE tracer requires at least one Gaussian")
        if not np.all(np.isfinite(self.amplitudes)):
            raise ValueError("MGE amplitudes must be finite")
        if not np.all(np.isfinite(self.sigmas_major_pc)) or np.any(self.sigmas_major_pc <= 0.0):
            raise ValueError("MGE sigmas must be positive and finite")

    def density(self, r_cyl_pc: float | np.ndarray, z_pc: float | np.ndarray) -> float | np.ndarray:
        r = np.asarray(r_cyl_pc)
        z = np.asarray(z_pc)
        m2 = r**2 + (z / self.q_intrinsic) ** 2
        basis = np.exp(-0.5 * m2[..., None] / self.sigmas_major_pc**2)
        return basis @ self.amplitudes
