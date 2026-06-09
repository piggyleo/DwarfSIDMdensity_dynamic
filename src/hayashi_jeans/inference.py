"""Minimal inference hooks for the pluggable Jeans likelihood."""

from __future__ import annotations

from collections.abc import Callable

import emcee
import numpy as np

from .data import GalaxyData
from .likelihood import GaussianVelocityLikelihood
from .params import HayashiParameters
from .projection import AxisymmetricJeansProjector

ParameterVector = np.ndarray
PriorFn = Callable[[ParameterVector], float]
ProjectorFactory = Callable[[ParameterVector], AxisymmetricJeansProjector]


def run_inference(
    galaxy: GalaxyData,
    initial_walkers: np.ndarray,
    *,
    log_prior: PriorFn,
    projector_factory: ProjectorFactory,
    n_steps: int,
    systemic_velocity_index: int | None = None,
    progress: bool = True,
) -> emcee.EnsembleSampler:
    """Run MCMC over Hayashi-style parameters.

    ``projector_factory`` is the key plug-in point: it can build a projector
    backed by the generalized Hernquist halo or by a new dark-matter model. If
    ``systemic_velocity_index`` is omitted, the systemic velocity is profiled
    analytically at each step using the total variance.
    """

    initial_walkers = np.asarray(initial_walkers, dtype=float)
    n_walkers, n_dim = initial_walkers.shape
    likelihood = GaussianVelocityLikelihood(galaxy)

    def log_probability(vector: np.ndarray) -> float:
        prior = log_prior(vector)
        if not np.isfinite(prior):
            return -np.inf
        try:
            projector = projector_factory(vector)
            sigma_los2 = projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
        except (FloatingPointError, ValueError, OverflowError):
            return -np.inf
        systemic = None if systemic_velocity_index is None else float(vector[systemic_velocity_index])
        return prior + likelihood.log_likelihood(sigma_los2, systemic_velocity_kms=systemic)

    sampler = emcee.EnsembleSampler(n_walkers, n_dim, log_probability)
    sampler.run_mcmc(initial_walkers, n_steps, progress=progress)
    return sampler


def vector_to_hayashi_params(vector: ParameterVector) -> HayashiParameters:
    """Map a nine-element vector to ``HayashiParameters``.

    Vector order:
    ``Q, log10_b_halo_pc, log10_rho0, -log10(1 - beta_z), alpha, beta, gamma, i_deg, u``.
    """

    q, log_b, log_rho, q_beta, alpha, beta, gamma, i_deg, systemic = vector
    return HayashiParameters(
        q_halo=float(q),
        b_halo_pc=float(10.0**log_b),
        rho0_msun_pc3=float(10.0**log_rho),
        beta_z=float(1.0 - 10.0 ** (-q_beta)),
        alpha=float(alpha),
        beta=float(beta),
        gamma=float(gamma),
        inclination_rad=float(np.deg2rad(i_deg)),
        systemic_velocity_kms=float(systemic),
    )
