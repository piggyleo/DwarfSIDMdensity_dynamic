"""High-level helpers that wire halo, Jeans projection, and likelihood."""

from __future__ import annotations

import numpy as np

from .data import GalaxyData
from .halos import GeneralizedHernquistHalo, HaloModel
from .likelihood import GaussianVelocityLikelihood
from .params import HayashiParameters
from .projection import AxisymmetricJeansProjector


def build_projector(
    galaxy: GalaxyData,
    params: HayashiParameters,
    *,
    halo: HaloModel,
    **projector_kwargs: float,
) -> AxisymmetricJeansProjector:
    """Build a Jeans projector around an already resolved halo model."""

    return AxisymmetricJeansProjector(
        b_star_pc=galaxy.observables.b_star_pc,
        qprime=galaxy.observables.qprime,
        inclination_rad=params.inclination_rad,
        beta_z=params.beta_z,
        halo=halo,
        **projector_kwargs,
    )


def build_hernquist_projector(
    galaxy: GalaxyData,
    params: HayashiParameters,
    **projector_kwargs: float,
) -> AxisymmetricJeansProjector:
    halo = GeneralizedHernquistHalo(**params.halo_kwargs())
    return build_projector(
        galaxy,
        params,
        halo=halo,
        **projector_kwargs,
    )


def log_likelihood_hernquist(galaxy: GalaxyData, params: HayashiParameters) -> float:
    projector = build_hernquist_projector(galaxy, params)
    sigma_los2 = projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    return GaussianVelocityLikelihood(galaxy).log_likelihood(
        sigma_los2,
        systemic_velocity_kms=params.systemic_velocity_kms,
    )


def make_reference_hernquist_params(galaxy: GalaxyData) -> HayashiParameters:
    """A conservative starting point for smoke tests, not a published best fit."""

    velocity = galaxy.velocity_kms
    systemic = float(np.nanmean(velocity))
    return HayashiParameters(
        q_halo=1.0,
        b_halo_pc=5.0 * galaxy.observables.b_star_pc,
        rho0_msun_pc3=0.08,
        beta_z=0.0,
        alpha=1.0,
        beta=4.0,
        gamma=1.0,
        inclination_rad=np.deg2rad(75.0),
        systemic_velocity_kms=systemic,
    )
