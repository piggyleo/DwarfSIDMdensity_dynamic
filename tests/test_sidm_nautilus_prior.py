import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from run_galaxy_nautilus import (
    CONCENTRATION_SCATTER_PRIOR_MEAN,
    CONCENTRATION_SCATTER_PRIOR_STD,
    CONCENTRATION_SCATTER_PRIOR_TRUNCATION,
    concentration_scatter_prior_distribution,
    log_prior_vector,
    make_nautilus_prior,
)


class RecordingPrior:
    def __init__(self) -> None:
        self.parameters = {}

    def add_parameter(self, key, dist) -> None:
        self.parameters[key] = dist


def fake_galaxy():
    return SimpleNamespace(observables=SimpleNamespace(qprime=0.7))


def scatter_vector(scatter: float) -> np.ndarray:
    return np.array([1.0, 9.0, scatter, 0.5, 0.3, 75.0, -20.0])


def test_nautilus_uses_default_truncated_sigma_three_gaussian_scatter_prior():
    prior = make_nautilus_prior(
        RecordingPrior,
        fake_galaxy(),
        (-50.0, 10.0),
        halo_model="sidm",
        sidm_parameterization="m200-ludlow-scatter",
    )
    distribution = prior.parameters["concentration_scatter_sigma"]
    assert distribution.kwds["loc"] == CONCENTRATION_SCATTER_PRIOR_MEAN
    assert distribution.kwds["scale"] == CONCENTRATION_SCATTER_PRIOR_STD
    np.testing.assert_allclose(
        distribution.support(),
        (-CONCENTRATION_SCATTER_PRIOR_TRUNCATION, CONCENTRATION_SCATTER_PRIOR_TRUNCATION),
    )


def test_explicit_scatter_log_prior_is_truncated_gaussian():
    galaxy = fake_galaxy()
    kwargs = {
        "systemic_velocity_bounds": (-50.0, 10.0),
        "halo_model": "sidm",
        "sidm_parameterization": "m200-ludlow-scatter",
    }
    logp_zero = log_prior_vector(galaxy, scatter_vector(0.0), **kwargs)
    logp_one_std = log_prior_vector(galaxy, scatter_vector(3.0), **kwargs)
    logp_lower_edge = log_prior_vector(galaxy, scatter_vector(-4.0), **kwargs)
    logp_upper_edge = log_prior_vector(galaxy, scatter_vector(4.0), **kwargs)
    logp_below = log_prior_vector(galaxy, scatter_vector(-4.01), **kwargs)
    logp_above = log_prior_vector(galaxy, scatter_vector(4.01), **kwargs)

    assert np.isfinite(logp_lower_edge)
    assert np.isfinite(logp_upper_edge)
    assert logp_below == -np.inf
    assert logp_above == -np.inf
    np.testing.assert_allclose(logp_one_std - logp_zero, -0.5)


def test_custom_scatter_truncation_changes_support():
    distribution = concentration_scatter_prior_distribution(2.5)
    np.testing.assert_allclose(distribution.support(), (-2.5, 2.5))
