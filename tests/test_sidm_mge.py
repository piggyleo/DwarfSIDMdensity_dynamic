import numpy as np
from scipy.special import erf
from scipy.integrate import cumulative_trapezoid

from hayashi_jeans.halos import SIDMPSIDM25Halo, SpheroidallyStratifiedMGEHalo
from hayashi_jeans.mge import (
    decompose_density_mge_auto,
    mge_density_relative_errors,
)
from hayashi_jeans.tracer import AxisymmetricMGETracer


def mge_enclosed_mass(radius: np.ndarray, amplitudes: np.ndarray, sigmas: np.ndarray) -> np.ndarray:
    r = np.asarray(radius)[:, None]
    s = np.asarray(sigmas)[None, :]
    term = (
        np.sqrt(np.pi / 2.0) * s**3 * erf(r / (np.sqrt(2.0) * s))
        - s**2 * r * np.exp(-0.5 * (r / s) ** 2)
    )
    return 4.0 * np.pi * np.sum(np.asarray(amplitudes)[None, :] * term, axis=1)


def test_validated_mge_represents_sidm_density_and_mass():
    for tau in (0.0, 0.1, 0.5, 0.9, 1.05, 1.08):
        halo = SIDMPSIDM25Halo(
            q=1.0,
            rs0_pc=1600.0,
            rho_s0_msun_pc3=0.1,
            tau=tau,
        )
        mge = decompose_density_mge_auto(
            halo.density_at_ellipsoidal_radius,
            n_gauss=60,
            decomposition_terms=28,
            r_min_pc=1.0e-2,
            r_max_pc=2.0e5,
            validation_p95_max=0.05,
            n_fit_radii=1600,
            lstsq_rcond=1.0e-12,
        )
        radius = np.geomspace(0.1, 5.0e4, 800)
        density_errors = mge_density_relative_errors(
            halo.density_at_ellipsoidal_radius,
            mge,
            radius,
        )
        assert np.percentile(density_errors, 95.0) < 0.005

        integration_radius = np.concatenate([[0.0], radius])
        density = np.asarray(halo.density_at_ellipsoidal_radius(integration_radius))
        with np.errstate(invalid="ignore"):
            mass_integrand = 4.0 * np.pi * integration_radius**2 * density
        mass_integrand[0] = 0.0
        true_mass = cumulative_trapezoid(
            mass_integrand,
            integration_radius,
            initial=0.0,
        )[1:]
        fitted_mass = mge_enclosed_mass(radius, mge.amplitudes, mge.sigmas_major_pc)
        mass_errors = np.abs(fitted_mass - true_mass) / np.maximum(np.abs(true_mass), 1.0e-300)
        assert np.percentile(mass_errors[20:], 95.0) < 0.05


def test_signed_mge_descriptors_preserve_amplitude_signs():
    amplitudes = np.array([2.0, -0.25])
    sigmas = np.array([10.0, 30.0])
    radius = np.array([0.0, 5.0, 20.0])
    expected = np.exp(-0.5 * (radius[:, None] / sigmas) ** 2) @ amplitudes

    halo = SpheroidallyStratifiedMGEHalo(
        q=0.8,
        amplitudes_msun_pc3=amplitudes,
        sigmas_major_pc=sigmas,
    )
    np.testing.assert_allclose(halo.density_at_ellipsoidal_radius(radius), expected)

    tracer = AxisymmetricMGETracer(
        q_intrinsic=0.8,
        amplitudes=amplitudes,
        sigmas_major_pc=sigmas,
    )
    np.testing.assert_allclose(tracer.density(radius, np.zeros_like(radius)), expected)
