"""Analytic MGE/JAM utilities for Jeans likelihoods and physicality checks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import roots_legendre

from .constants import G_PC_MSUN_KMS2
from .data import GalaxyData
from .halos import GeneralizedHernquistHalo, SIDMPSIDM25Halo, SpheroidallyStratifiedHalo
from .tracer import intrinsic_q_from_projected


@dataclass(frozen=True)
class MGE1D:
    amplitudes: np.ndarray
    sigmas_major_pc: np.ndarray
    decomposition_method: str = "unknown"
    validation_p95_relative_error: float = np.nan


@dataclass(frozen=True)
class MGEPhysicalityConfig:
    n_gauss_halo: int = 45
    n_gauss_tracer: int = 45
    # Kept for compatibility with older CLI/scripts; analytic decomposition does not use it.
    n_fit_radii: int = 800
    decomposition_terms: int = 28
    n_u: int = 96
    r_min_pc: float = 1.0e-2
    r_max_pc: float = 2.0e5
    # Physicality guard settings. n_los_guard is kept as the public/CLI name
    # from the earlier LOS guard; in the current guard it is the number of
    # central radial samples.
    n_los_guard: int = 10
    los_factor: float = 20.0
    n_guard_angles: int = 3
    guard_theta_max_rad: float = np.pi / 3.0
    vphi2_abs_tol: float = 1.0e-8
    vphi2_rel_tol: float = 1.0e-5
    halo_decomposition_method: str = "analytic"
    decomposition_validation_p95_max: float = 0.05
    real_fit_radii: int = 1600
    real_lstsq_rcond: float = 1.0e-12


def mge_physicality_check(
    galaxy: GalaxyData,
    *,
    q_halo: float,
    b_halo_pc: float,
    alpha: float,
    beta: float,
    gamma: float,
    beta_z: float,
    inclination_rad: float,
    config: MGEPhysicalityConfig = MGEPhysicalityConfig(),
) -> bool:
    """Return whether the raw MGE/JAM Jeans solution is locally physical.

    This is a hard physicality guard. It deliberately uses unit halo density:
    multiplying by rho0 cannot change the sign of any second moment.
    """

    sigma_los2_unit = mge_sigma_los2_unit_checked(
        galaxy,
        q_halo=q_halo,
        b_halo_pc=b_halo_pc,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        beta_z=beta_z,
        inclination_rad=inclination_rad,
        config=config,
    )
    return bool(np.all(np.isfinite(sigma_los2_unit)) and np.all(sigma_los2_unit > 0.0))


def mge_sigma_los2_unit(
    galaxy: GalaxyData,
    *,
    q_halo: float,
    b_halo_pc: float,
    alpha: float,
    beta: float,
    gamma: float,
    beta_z: float,
    inclination_rad: float,
    config: MGEPhysicalityConfig = MGEPhysicalityConfig(),
) -> np.ndarray:
    """Compute raw unit-density MGE/JAM sigma_los^2 at observed stars."""

    halo = GeneralizedHernquistHalo(
        q=q_halo,
        b_pc=b_halo_pc,
        rho0_msun_pc3=1.0,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )
    return mge_sigma_los2_for_halo(
        galaxy,
        halo=halo,
        beta_z=beta_z,
        inclination_rad=inclination_rad,
        config=config,
        check_physicality=False,
    )


def mge_sigma_los2_unit_checked(
    galaxy: GalaxyData,
    *,
    q_halo: float,
    b_halo_pc: float,
    alpha: float,
    beta: float,
    gamma: float,
    beta_z: float,
    inclination_rad: float,
    config: MGEPhysicalityConfig = MGEPhysicalityConfig(),
) -> np.ndarray:
    """Compute MGE/JAM sigma_los^2 using one MGE fit, rejecting local v_phi^2 failures."""

    halo = GeneralizedHernquistHalo(
        q=q_halo,
        b_pc=b_halo_pc,
        rho0_msun_pc3=1.0,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )
    return mge_sigma_los2_for_halo(
        galaxy,
        halo=halo,
        beta_z=beta_z,
        inclination_rad=inclination_rad,
        config=config,
        check_physicality=True,
    )


def mge_sigma_los2_for_halo(
    galaxy: GalaxyData,
    *,
    halo: SpheroidallyStratifiedHalo,
    beta_z: float,
    inclination_rad: float,
    config: MGEPhysicalityConfig = MGEPhysicalityConfig(),
    check_physicality: bool = True,
) -> np.ndarray:
    """Compute analytic MGE/JAM moments for a resolved spheroidal halo."""

    failed = np.full_like(np.asarray(galaxy.x_pc, dtype=float), np.nan, dtype=float)
    try:
        q_star, tracer_mge, halo_mge = fit_galaxy_mges_for_halo(
            galaxy,
            halo=halo,
            inclination_rad=inclination_rad,
            config=config,
        )
        if check_physicality and not mge_vphi2_physicality_check(
            galaxy,
            q_halo=halo.q,
            q_star=q_star,
            beta_z=beta_z,
            inclination_rad=inclination_rad,
            tracer_mge=tracer_mge,
            halo_mge=halo_mge,
            config=config,
        ):
            return failed
        return mge_los_second_moment(
            x_pc=galaxy.x_pc,
            y_pc=galaxy.y_pc,
            halo_mge=halo_mge,
            tracer_mge=tracer_mge,
            q_halo=halo.q,
            q_star=q_star,
            beta_z=beta_z,
            inclination_rad=inclination_rad,
            n_u=config.n_u,
        )
    except (FloatingPointError, ValueError, ZeroDivisionError, RuntimeError):
        return failed


def _fit_galaxy_mges(
    galaxy: GalaxyData,
    *,
    q_halo: float,
    b_halo_pc: float,
    alpha: float,
    beta: float,
    gamma: float,
    inclination_rad: float,
    config: MGEPhysicalityConfig,
) -> tuple[float, MGE1D, MGE1D]:
    halo = GeneralizedHernquistHalo(
        q=q_halo,
        b_pc=b_halo_pc,
        rho0_msun_pc3=1.0,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )
    return fit_galaxy_mges_for_halo(
        galaxy,
        halo=halo,
        inclination_rad=inclination_rad,
        config=config,
    )


def fit_galaxy_mges_for_halo(
    galaxy: GalaxyData,
    *,
    halo: SpheroidallyStratifiedHalo,
    inclination_rad: float,
    config: MGEPhysicalityConfig,
) -> tuple[float, MGE1D, MGE1D]:
    """Decompose the tracer and an arbitrary spheroidal halo analytically."""

    q_star = intrinsic_q_from_projected(galaxy.observables.qprime, inclination_rad)
    tracer_mge = decompose_density_mge(
        lambda m: (1.0 + (m / galaxy.observables.b_star_pc) ** 2) ** (-2.5),
        n_gauss=config.n_gauss_tracer,
        decomposition_terms=config.decomposition_terms,
        r_min_pc=config.r_min_pc,
        r_max_pc=config.r_max_pc,
    )
    halo_method = config.halo_decomposition_method
    if halo_method == "auto" or isinstance(halo, SIDMPSIDM25Halo) and halo_method == "analytic":
        halo_mge = decompose_density_mge_auto(
            halo.density_at_ellipsoidal_radius,
            n_gauss=config.n_gauss_halo,
            decomposition_terms=config.decomposition_terms,
            r_min_pc=config.r_min_pc,
            r_max_pc=config.r_max_pc,
            validation_p95_max=config.decomposition_validation_p95_max,
            n_fit_radii=config.real_fit_radii,
            lstsq_rcond=config.real_lstsq_rcond,
        )
    elif halo_method == "real-lstsq":
        halo_mge = decompose_density_mge_real_lstsq(
            halo.density_at_ellipsoidal_radius,
            n_gauss=config.n_gauss_halo,
            n_fit_radii=config.real_fit_radii,
            r_min_pc=config.r_min_pc,
            r_max_pc=config.r_max_pc,
            lstsq_rcond=config.real_lstsq_rcond,
        )
    elif halo_method == "analytic":
        halo_mge = decompose_density_mge(
            halo.density_at_ellipsoidal_radius,
            n_gauss=config.n_gauss_halo,
            decomposition_terms=config.decomposition_terms,
            r_min_pc=config.r_min_pc,
            r_max_pc=config.r_max_pc,
        )
    else:
        raise ValueError(f"unknown halo_decomposition_method={halo_method!r}")
    return q_star, tracer_mge, halo_mge


def mge_vphi2_physicality_check(
    galaxy: GalaxyData,
    *,
    q_halo: float,
    q_star: float,
    beta_z: float,
    inclination_rad: float,
    tracer_mge: MGE1D,
    halo_mge: MGE1D,
    config: MGEPhysicalityConfig = MGEPhysicalityConfig(),
) -> bool:
    """Check that local raw v_phi^2 is not significantly negative on a central grid.

    Under the Hayashi/JAM assumptions, v_z^2=P_z/nu and
    v_R^2=v_z^2/(1-beta_z) are positive when the fitted tracer density,
    vertical pressure, and beta_z<1 are well behaved. The only second moment
    that can become substantially negative is therefore v_phi^2. Checking it
    on a small galaxy-centered meridional grid catches the radial pressure
    gradient failure mode without repeating the expensive calculation along
    every observed stellar LOS.
    """

    if not (beta_z < 1.0):
        return False
    del inclination_rad

    n_radii = max(int(config.n_los_guard), 2)
    n_angles = max(int(config.n_guard_angles), 1)
    b_star_pc = float(galaxy.observables.b_star_pc)
    r_min_pc = max(float(config.r_min_pc), 1.0e-2 * b_star_pc)
    r_max_pc = max(float(config.los_factor) * b_star_pc, r_min_pc * (1.0 + 1.0e-6))
    meridional_radius = np.geomspace(r_min_pc, r_max_pc, n_radii)
    theta = np.linspace(0.0, float(config.guard_theta_max_rad), n_angles)

    r = meridional_radius[:, None] * np.cos(theta)[None, :]
    z = q_star * meridional_radius[:, None] * np.sin(theta)[None, :]
    nu = _mge_tracer_density(r, z, q_star=q_star, tracer_mge=tracer_mge)
    if not np.all(np.isfinite(nu)):
        return False

    pressure, d_pressure_dr = _mge_vertical_pressure_and_dr(
        r,
        z,
        q_halo=q_halo,
        q_star=q_star,
        tracer_mge=tracer_mge,
        halo_mge=halo_mge,
        n_u=config.n_u,
    )
    dphi_dr = _mge_halo_dphi_dr(r, z, q_halo=q_halo, halo_mge=halo_mge, n_u=config.n_u)

    positive = nu > 0.0
    v_z2 = np.full_like(nu, np.nan, dtype=float)
    v_z2[positive] = pressure[positive] / nu[positive]
    v_r2 = v_z2 / (1.0 - beta_z)
    pressure_gradient_term = np.full_like(nu, np.nan, dtype=float)
    pressure_gradient_term[positive] = r[positive] / nu[positive] * d_pressure_dr[positive] / (1.0 - beta_z)
    r_dphi_dr = r * dphi_dr
    v_phi2 = v_r2 + pressure_gradient_term + r_dphi_dr

    scale = np.maximum.reduce(
        [
            np.ones_like(v_phi2),
            np.abs(v_z2),
            np.abs(v_r2),
            np.abs(pressure_gradient_term),
            np.abs(r_dphi_dr),
        ]
    )
    tol = np.maximum(config.vphi2_abs_tol, config.vphi2_rel_tol * scale)
    checked = np.isfinite(v_phi2) & np.isfinite(tol)
    if not np.all(checked):
        return False
    return bool(not np.any(v_phi2[checked] < -tol[checked]))


def _shajib_kesi(func_terms: int) -> np.ndarray:
    """Return the complex sampling abscissae used by the MGE decomposition."""

    n = np.arange(0, 2 * func_terms + 1)
    return np.sqrt(2.0 * func_terms * np.log(10.0) / 3.0 + 2.0j * np.pi * n)


def _shajib_eta(func_terms: int) -> np.ndarray:
    """Return the complex quadrature weights used by the MGE decomposition."""

    i = np.arange(1, func_terms)
    last = 1.0 / 2.0**func_terms
    if i.size:
        middle = last + np.cumsum(np.cumprod((func_terms + 1 - i) / i) * last)
        coeffs = np.hstack([np.array([0.5]), np.ones(func_terms), middle[::-1], np.array([last])])
    else:
        coeffs = np.array([0.5, 1.0, last])
    signs = (-1.0) ** np.arange(0, 2 * func_terms + 1)
    prefactor = 2.0 * np.sqrt(2.0 * np.pi) * 10.0 ** (func_terms / 3.0)
    return signs * prefactor * coeffs


def decompose_density_mge(
    density_func,
    *,
    n_gauss: int,
    decomposition_terms: int,
    r_min_pc: float,
    r_max_pc: float,
) -> MGE1D:
    """Analytically decompose an intrinsic 3D density into Gaussian terms.

    This follows the Shajib/AutoGalaxy MGE kernel for a 3D density directly,
    instead of solving a positive least-squares fit. The resulting amplitudes
    are signed and must be propagated through the moment equations as such.
    """

    if n_gauss < 2:
        raise ValueError("MGE decomposition requires at least two Gaussian terms")
    if decomposition_terms < 1:
        raise ValueError("MGE decomposition requires decomposition_terms >= 1")
    sigmas = np.exp(np.linspace(np.log(r_min_pc), np.log(r_max_pc), n_gauss))
    d_log_sigma = (np.log(r_max_pc) - np.log(r_min_pc)) / (n_gauss - 1)
    kesi = _shajib_kesi(decomposition_terms)
    eta = _shajib_eta(decomposition_terms)
    sampled = density_func(sigmas[:, None] * kesi[None, :])
    f_sigma = np.sum(eta[None, :] * sampled, axis=1)
    amplitudes = np.real(f_sigma) * d_log_sigma / np.sqrt(2.0 * np.pi)
    amplitudes[0] *= 0.5
    amplitudes[-1] *= 0.5
    if not (np.all(np.isfinite(amplitudes)) and np.all(np.isfinite(sigmas))):
        raise ValueError("MGE decomposition returned non-finite terms")
    return MGE1D(
        amplitudes=amplitudes,
        sigmas_major_pc=sigmas,
        decomposition_method="analytic",
    )


def reconstruct_mge_density(radius_pc: np.ndarray, mge: MGE1D) -> np.ndarray:
    radius = np.asarray(radius_pc, dtype=float)
    basis = np.exp(
        -0.5
        * (radius[:, None] / np.asarray(mge.sigmas_major_pc, dtype=float)[None, :]) ** 2
    )
    return basis @ np.asarray(mge.amplitudes, dtype=float)


def mge_density_relative_errors(density_func, mge: MGE1D, radius_pc: np.ndarray) -> np.ndarray:
    true_density = np.asarray(density_func(radius_pc), dtype=float)
    fitted_density = reconstruct_mge_density(radius_pc, mge)
    return np.abs(fitted_density - true_density) / np.maximum(np.abs(true_density), 1.0e-300)


def decompose_density_mge_real_lstsq(
    density_func,
    *,
    n_gauss: int,
    n_fit_radii: int,
    r_min_pc: float,
    r_max_pc: float,
    lstsq_rcond: float = 1.0e-12,
) -> MGE1D:
    """Fit signed Gaussian amplitudes on the real axis with relative weighting."""

    if n_gauss < 2 or n_fit_radii < n_gauss:
        raise ValueError("real-axis MGE requires n_fit_radii >= n_gauss >= 2")
    fit_radius = np.geomspace(r_min_pc, r_max_pc, n_fit_radii)
    density = np.asarray(density_func(fit_radius), dtype=float)
    if not np.all(np.isfinite(density)) or np.any(density <= 0.0):
        raise ValueError("real-axis MGE requires a finite positive density profile")
    sigmas = np.geomspace(r_min_pc / 3.0, r_max_pc * 3.0, n_gauss)
    basis = np.exp(-0.5 * (fit_radius[:, None] / sigmas[None, :]) ** 2)
    weighted_basis = basis / density[:, None]
    amplitudes = np.linalg.lstsq(
        weighted_basis,
        np.ones_like(density),
        rcond=lstsq_rcond,
    )[0]
    mge = MGE1D(
        amplitudes=amplitudes,
        sigmas_major_pc=sigmas,
        decomposition_method="real-lstsq",
    )
    validation_radius = np.geomspace(r_min_pc * 10.0, r_max_pc / 4.0, 600)
    p95 = float(np.percentile(mge_density_relative_errors(density_func, mge, validation_radius), 95.0))
    return MGE1D(
        amplitudes=amplitudes,
        sigmas_major_pc=sigmas,
        decomposition_method="real-lstsq",
        validation_p95_relative_error=p95,
    )


def decompose_density_mge_auto(
    density_func,
    *,
    n_gauss: int,
    decomposition_terms: int,
    r_min_pc: float,
    r_max_pc: float,
    validation_p95_max: float,
    n_fit_radii: int,
    lstsq_rcond: float,
) -> MGE1D:
    """Prefer the analytic transform, falling back after an explicit fit check."""

    analytic = decompose_density_mge(
        density_func,
        n_gauss=n_gauss,
        decomposition_terms=decomposition_terms,
        r_min_pc=r_min_pc,
        r_max_pc=r_max_pc,
    )
    validation_radius = np.geomspace(r_min_pc * 10.0, r_max_pc / 4.0, 600)
    analytic_errors = mge_density_relative_errors(density_func, analytic, validation_radius)
    analytic_fit = reconstruct_mge_density(validation_radius, analytic)
    analytic_p95 = float(np.percentile(analytic_errors, 95.0))
    if (
        np.all(np.isfinite(analytic_fit))
        and np.all(analytic_fit > 0.0)
        and analytic_p95 <= validation_p95_max
    ):
        return MGE1D(
            amplitudes=analytic.amplitudes,
            sigmas_major_pc=analytic.sigmas_major_pc,
            decomposition_method="analytic",
            validation_p95_relative_error=analytic_p95,
        )
    return decompose_density_mge_real_lstsq(
        density_func,
        n_gauss=n_gauss,
        n_fit_radii=n_fit_radii,
        r_min_pc=r_min_pc,
        r_max_pc=r_max_pc,
        lstsq_rcond=lstsq_rcond,
    )


def fit_positive_mge(
    density_func,
    *,
    n_gauss: int,
    n_fit_radii: int,
    r_min_pc: float,
    r_max_pc: float,
) -> MGE1D:
    """Backward-compatible wrapper for the analytic MGE decomposition."""

    del n_fit_radii
    return decompose_density_mge(
        density_func,
        n_gauss=n_gauss,
        decomposition_terms=MGEPhysicalityConfig.decomposition_terms,
        r_min_pc=r_min_pc,
        r_max_pc=r_max_pc,
    )


def generalized_halo_density_unit(
    m: np.ndarray,
    *,
    b_pc: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> np.ndarray:
    x = np.asarray(m) / b_pc
    x = np.where(np.abs(x) > 1.0e-300, x, 1.0e-300 + 0.0j)
    return x ** (-gamma) * (1.0 + x**alpha) ** (-(beta - gamma) / alpha)


def mge_los_second_moment(
    *,
    x_pc: np.ndarray,
    y_pc: np.ndarray,
    halo_mge: MGE1D,
    tracer_mge: MGE1D,
    q_halo: float,
    q_star: float,
    beta_z: float,
    inclination_rad: float,
    n_u: int,
) -> np.ndarray:
    """Return raw MGE/JAM LOS second moment divided by projected tracer."""

    nodes, weights = roots_legendre(n_u)
    u = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    x = np.asarray(x_pc, dtype=float)
    y = np.asarray(y_pc, dtype=float)
    cos_i = np.cos(inclination_rad)
    sin_i = np.sin(inclination_rad)
    b_aniso = 1.0 / (1.0 - beta_z)

    qpk2 = cos_i**2 + q_star**2 * sin_i**2
    qpk = np.sqrt(qpk2)

    sigma_h = q_halo * halo_mge.sigmas_major_pc
    rho0 = halo_mge.amplitudes
    tracer_sigma_major = tracer_mge.sigmas_major_pc
    sigma_t = q_star * tracer_sigma_major
    nu0 = tracer_mge.amplitudes

    projected_tracer = np.zeros_like(x, dtype=float)
    for amp_k, sig_k_major in zip(nu0, tracer_sigma_major):
        if not (np.isfinite(amp_k) and np.isfinite(sig_k_major)) or amp_k == 0.0:
            continue
        projected_tracer += (
            np.sqrt(2.0 * np.pi)
            * q_star
            * sig_k_major
            / qpk
            * amp_k
            * np.exp(-(x**2 + y**2 / qpk2) / (2.0 * sig_k_major**2))
        )

    numerator = np.zeros_like(x, dtype=float)
    u2 = u * u
    x2 = x[None, :] ** 2
    y2 = y[None, :] ** 2

    for amp_j, sig_j in zip(rho0, sigma_h):
        if not (np.isfinite(amp_j) and np.isfinite(sig_j)) or amp_j == 0.0:
            continue
        one_minus_qj2 = 1.0 - q_halo**2
        halo_shape = 1.0 - one_minus_qj2 * u2
        for amp_k, sig_k in zip(nu0, sigma_t):
            if not (np.isfinite(amp_k) and np.isfinite(sig_k)) or amp_k == 0.0:
                continue
            a = 0.5 * (u2 * q_halo**2 / sig_j**2 + q_star**2 / sig_k**2)
            b = 0.5 * (
                (1.0 - q_star**2) / sig_k**2
                + q_halo**2 * (1.0 - q_halo**2) * u2**2 / (sig_j**2 * halo_shape)
            )
            c = 1.0 - q_halo**2 - q_halo**2 * sig_k**2 / sig_j**2
            d = 1.0 - b_aniso * q_star**2 - ((1.0 - b_aniso) * c + (1.0 - q_halo**2) * b_aniso) * u2
            denom = (1.0 - c * u2) * np.sqrt((a + b * cos_i**2) * halo_shape)
            velocity_factor = sig_k**2 * (cos_i**2 + b_aniso * sin_i**2)
            exponent = -a[:, None] * (x2 + ((a + b) / (a + b * cos_i**2))[:, None] * y2)
            kernel = (
                u2[:, None]
                * (velocity_factor + d[:, None] * x2 * sin_i**2)
                / denom[:, None]
                * np.exp(exponent)
            )
            integral = np.sum(w[:, None] * kernel, axis=0)
            numerator += 4.0 * np.pi ** 1.5 * G_PC_MSUN_KMS2 * q_halo * amp_j * amp_k * integral

    sigma_los2 = np.full_like(numerator, np.nan, dtype=float)
    valid_projected = np.isfinite(projected_tracer) & (projected_tracer > 0.0)
    sigma_los2[valid_projected] = numerator[valid_projected] / projected_tracer[valid_projected]
    return sigma_los2


def _mge_tracer_density(
    r_pc: np.ndarray,
    z_pc: np.ndarray,
    *,
    q_star: float,
    tracer_mge: MGE1D,
) -> np.ndarray:
    r = np.asarray(r_pc, dtype=float)
    z = np.asarray(z_pc, dtype=float)
    density = np.zeros_like(r, dtype=float)
    for amp, sigma_major in zip(tracer_mge.amplitudes, tracer_mge.sigmas_major_pc):
        if not (np.isfinite(amp) and np.isfinite(sigma_major)) or amp == 0.0:
            continue
        exponent = -0.5 * (r**2 + (z / q_star) ** 2) / sigma_major**2
        density += amp * np.exp(exponent)
    return density


def _mge_halo_dphi_dr(
    r_pc: np.ndarray,
    z_pc: np.ndarray,
    *,
    q_halo: float,
    halo_mge: MGE1D,
    n_u: int,
) -> np.ndarray:
    r = np.asarray(r_pc, dtype=float)
    z = np.asarray(z_pc, dtype=float)
    nodes, weights = roots_legendre(n_u)
    u = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    u2 = u * u
    shape = 1.0 + (q_halo**2 - 1.0) * u2
    prefactor = 2.0 * np.pi * G_PC_MSUN_KMS2 * q_halo
    force_integral = np.zeros_like(r, dtype=float)
    for amp, sigma_major in zip(halo_mge.amplitudes, halo_mge.sigmas_major_pc):
        if not (np.isfinite(amp) and np.isfinite(sigma_major)) or amp == 0.0:
            continue
        exponent = -0.5 * u2[:, None, None] * (
            r[None, :, :] ** 2 + z[None, :, :] ** 2 / shape[:, None, None]
        ) / sigma_major**2
        kernel = 2.0 * u2[:, None, None] / np.sqrt(shape)[:, None, None] * np.exp(exponent)
        force_integral += amp * np.sum(w[:, None, None] * kernel, axis=0)
    return prefactor * r * force_integral


def _mge_vertical_pressure_and_dr(
    r_pc: np.ndarray,
    z_pc: np.ndarray,
    *,
    q_halo: float,
    q_star: float,
    tracer_mge: MGE1D,
    halo_mge: MGE1D,
    n_u: int,
) -> tuple[np.ndarray, np.ndarray]:
    r = np.asarray(r_pc, dtype=float)
    z_abs = np.abs(np.asarray(z_pc, dtype=float))
    nodes, weights = roots_legendre(n_u)
    u = 0.5 * (nodes + 1.0)
    w = 0.5 * weights
    u2 = u * u
    shape = 1.0 + (q_halo**2 - 1.0) * u2
    force_prefactor = 2.0 * np.pi * G_PC_MSUN_KMS2 * q_halo

    pressure = np.zeros_like(r, dtype=float)
    d_pressure_dr = np.zeros_like(r, dtype=float)
    for halo_amp, halo_sigma_major in zip(halo_mge.amplitudes, halo_mge.sigmas_major_pc):
        if not (np.isfinite(halo_amp) and np.isfinite(halo_sigma_major)) or halo_amp == 0.0:
            continue
        halo_r_coeff = u2 / halo_sigma_major**2
        halo_z_coeff = u2 / (halo_sigma_major**2 * shape)
        force_kernel_coeff = force_prefactor * 2.0 * u2 / (shape**1.5)
        for tracer_amp, tracer_sigma_major in zip(tracer_mge.amplitudes, tracer_mge.sigmas_major_pc):
            if not (np.isfinite(tracer_amp) and np.isfinite(tracer_sigma_major)) or tracer_amp == 0.0:
                continue
            radial_coeff = 1.0 / tracer_sigma_major**2 + halo_r_coeff
            vertical_coeff = 1.0 / (q_star**2 * tracer_sigma_major**2) + halo_z_coeff
            exponent = -0.5 * (
                radial_coeff[:, None, None] * r[None, :, :] ** 2
                + vertical_coeff[:, None, None] * z_abs[None, :, :] ** 2
            )
            term = (
                w[:, None, None]
                * force_kernel_coeff[:, None, None]
                * np.exp(exponent)
                / vertical_coeff[:, None, None]
            )
            summed = halo_amp * tracer_amp * np.sum(term, axis=0)
            pressure += summed
            d_pressure_dr -= halo_amp * tracer_amp * r * np.sum(
                radial_coeff[:, None, None] * term,
                axis=0,
            )
    return pressure, d_pressure_dr
