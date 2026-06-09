import numpy as np

from hayashi_jeans.halos import (
    SIDMPSIDM25Halo,
    sidm_dimensionless_density,
    sidm_evolution_coefficients,
)


def test_sidm_tau_zero_is_nfw():
    radius = np.geomspace(1.0e-5, 1.0e3, 400)
    expected = 1.0 / (radius * (1.0 + radius) ** 2)
    actual = sidm_dimensionless_density(radius, 0.0)
    np.testing.assert_allclose(actual, expected, rtol=2.0e-13, atol=0.0)


def test_sidm_density_scales_linearly_with_rho_s0():
    halo_unit = SIDMPSIDM25Halo(q=0.8, rs0_pc=800.0, rho_s0_msun_pc3=1.0, tau=0.5)
    halo_scaled = SIDMPSIDM25Halo(q=0.8, rs0_pc=800.0, rho_s0_msun_pc3=3.7, tau=0.5)
    for r, z in ((0.0, 20.0), (100.0, 40.0), (1200.0, 300.0)):
        np.testing.assert_allclose(halo_scaled.density(r, z), 3.7 * halo_unit.density(r, z), rtol=1e-13)


def test_sidm_fitted_coefficients_are_finite_over_calibrated_range():
    for tau in np.linspace(0.0, 1.08, 55):
        coeff = sidm_evolution_coefficients(float(tau))
        values = np.array(
            [
                coeff.density_ratio,
                coeff.scale_radius_ratio,
                coeff.core_radius_ratio,
                coeff.transition_alpha,
                coeff.inner_gamma,
            ]
        )
        assert np.all(np.isfinite(values))
        assert coeff.density_ratio > 0.0
        assert coeff.scale_radius_ratio > 0.0
        assert coeff.transition_alpha > 0.0


def test_sidm_halo_exposes_positive_gravity_gradient():
    halo = SIDMPSIDM25Halo(q=1.0, rs0_pc=1000.0, rho_s0_msun_pc3=0.1, tau=0.5)
    dphi_dr, dphi_dz = halo.potential_gradients(100.0, 50.0)
    assert dphi_dr > 0.0
    assert dphi_dz > 0.0
