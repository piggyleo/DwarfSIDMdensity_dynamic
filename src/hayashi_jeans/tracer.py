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

