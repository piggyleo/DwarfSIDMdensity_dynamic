from pathlib import Path

import numpy as np

from hayashi_jeans.data import LOW_RISK_GALAXIES, load_galaxy_data, load_low_risk_galaxies
from hayashi_jeans.halos import GeneralizedHernquistHalo
from hayashi_jeans.likelihood import GaussianVelocityLikelihood


ROOT = Path(__file__).resolve().parents[1]


def test_low_risk_counts_match_hayashi_table1():
    galaxies = load_low_risk_galaxies(ROOT / "data" / "galaxies")
    assert set(galaxies) == set(LOW_RISK_GALAXIES)
    for galaxy in galaxies.values():
        assert len(galaxy.stars) == galaxy.observables.n_sample_table1
        assert galaxy.stars["v_los_kms"].notna().all()
        assert galaxy.stars["v_los_err_kms"].notna().all()


def test_loader_derives_projected_coordinates():
    galaxy = load_galaxy_data(ROOT / "data" / "galaxies" / "12_Horologium_I.csv")
    assert len(galaxy.x_pc) == 5
    assert np.isfinite(galaxy.x_pc).all()
    assert np.isfinite(galaxy.y_pc).all()


def test_loader_uses_structural_center_catalog():
    galaxy = load_galaxy_data(ROOT / "data" / "galaxies" / "27_Willman_1.csv")
    assert galaxy.observables.center_ra_deg == 162.3436
    assert galaxy.observables.center_dec_deg == 51.0501


def test_velocity_likelihood_profiles_systemic_velocity():
    galaxy = load_galaxy_data(ROOT / "data" / "galaxies" / "12_Horologium_I.csv")
    sigma_los2 = np.full(len(galaxy.stars), 25.0)
    loglike = GaussianVelocityLikelihood(galaxy).log_likelihood(sigma_los2)
    assert np.isfinite(loglike)


def test_hernquist_halo_exposes_positive_gravity_gradient():
    halo = GeneralizedHernquistHalo(q=1.0, b_pc=100.0, rho0_msun_pc3=0.1, alpha=1.0, beta=4.0, gamma=1.0)
    dphi_dr, dphi_dz = halo.potential_gradients(50.0, 25.0)
    assert dphi_dr > 0.0
    assert dphi_dz > 0.0
