import numpy as np

from hayashi_jeans.halo_parameterizations import (
    HaloRunContext,
    ludlow16_concentration,
    sidm_from_m200_c200,
    sidm_from_m200_ludlow,
    sidm_from_scale_parameters,
)


def test_scale_and_m200_c200_resolve_to_identical_halos():
    context = HaloRunContext(redshift=0.2)
    mass_resolved = sidm_from_m200_c200(
        q_halo=0.9,
        m200_msun=1.0e9,
        c200=15.0,
        tau=0.35,
        context=context,
    )
    scale_resolved = sidm_from_scale_parameters(
        q_halo=mass_resolved.q_halo,
        rs0_pc=mass_resolved.rs0_pc,
        rho_s0_msun_pc3=mass_resolved.rho_s0_msun_pc3,
        tau=mass_resolved.tau,
    )
    mass_halo = mass_resolved.build_halo()
    scale_halo = scale_resolved.build_halo()
    radii = np.geomspace(0.01, 1.0e5, 100)
    np.testing.assert_allclose(
        mass_halo.density_at_ellipsoidal_radius(radii),
        scale_halo.density_at_ellipsoidal_radius(radii),
        rtol=1.0e-13,
    )


def test_ludlow_resolver_matches_explicit_concentration():
    context = HaloRunContext(redshift=0.5)
    c200 = ludlow16_concentration(m200_msun=1.0e10, context=context, scatter_sigma=0.4)
    ludlow = sidm_from_m200_ludlow(
        q_halo=1.0,
        m200_msun=1.0e10,
        tau=0.8,
        context=context,
        scatter_sigma=0.4,
    )
    explicit = sidm_from_m200_c200(
        q_halo=1.0,
        m200_msun=1.0e10,
        c200=c200,
        tau=0.8,
        context=context,
    )
    assert ludlow.concentration_relation == "ludlow16"
    np.testing.assert_allclose(ludlow.rs0_pc, explicit.rs0_pc, rtol=1.0e-13)
    np.testing.assert_allclose(ludlow.rho_s0_msun_pc3, explicit.rho_s0_msun_pc3, rtol=1.0e-13)


def test_mass_parameterization_requires_fixed_redshift():
    context = HaloRunContext()
    try:
        sidm_from_m200_c200(
            q_halo=1.0,
            m200_msun=1.0e9,
            c200=10.0,
            tau=0.5,
            context=context,
        )
    except ValueError as exc:
        assert "redshift" in str(exc)
    else:
        raise AssertionError("missing fixed redshift should be rejected")
