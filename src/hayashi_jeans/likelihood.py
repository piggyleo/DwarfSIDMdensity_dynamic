"""Unbinned Gaussian LOS velocity likelihood."""

from __future__ import annotations

import numpy as np

from .data import GalaxyData


def profile_systemic_velocity(velocity_kms: np.ndarray, variance_kms2: np.ndarray) -> float:
    weights = 1.0 / variance_kms2
    return float(np.sum(weights * velocity_kms) / np.sum(weights))


class GaussianVelocityLikelihood:
    """Hayashi et al. single-star Gaussian likelihood."""

    def __init__(self, galaxy: GalaxyData) -> None:
        self.galaxy = galaxy

    def log_likelihood(
        self,
        sigma_los2_kms2: np.ndarray,
        *,
        systemic_velocity_kms: float | None = None,
    ) -> float:
        v = self.galaxy.velocity_kms
        verr = self.galaxy.velocity_err_kms
        sigma_los2 = np.asarray(sigma_los2_kms2, dtype=float)
        if sigma_los2.shape != v.shape:
            raise ValueError(f"sigma_los2 shape {sigma_los2.shape} does not match velocities {v.shape}")
        if not np.all(np.isfinite(sigma_los2)) or np.any(sigma_los2 <= 0.0):
            return -np.inf

        total_var = sigma_los2 + verr**2
        if systemic_velocity_kms is None:
            systemic_velocity_kms = profile_systemic_velocity(v, total_var)
        residual = v - systemic_velocity_kms
        return float(-0.5 * np.sum(residual**2 / total_var + np.log(2.0 * np.pi * total_var)))

