#!/usr/bin/env python
"""Benchmark numeric vs analytic Plummer surface-density denominators.

This is an experimental comparison script. It leaves the production
AxisymmetricJeansProjector untouched and only replaces the LOS denominator

    I(x, y) = integral nu(R(ell), z(ell)) d ell

with the analytic projected surface density for the same unnormalized
axisymmetric Plummer tracer used by the code.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from hayashi_jeans.data import GalaxyData, load_low_risk_galaxies
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.model import build_hernquist_projector, make_reference_hernquist_params
from hayashi_jeans.projection import AxisymmetricJeansProjector


class AnalyticSurfaceDensityProjector(AxisymmetricJeansProjector):
    """Projector variant using analytic Plummer I(x,y) in the denominator."""

    def surface_density_analytic(self, x_pc: float, y_pc: float) -> float:
        q_intrinsic = self.tracer.q_intrinsic
        qprime = self.qprime
        b = self.b_star_pc
        mprime2 = x_pc**2 + y_pc**2 / qprime**2
        return (4.0 * b / 3.0) * (q_intrinsic / qprime) * (1.0 + mprime2 / b**2) ** (-2.0)

    def sigma_los2(self, x_pc: float, y_pc: float) -> float:
        i = self.inclination_rad

        def los_geometry(ell_pc: float) -> tuple[float, float]:
            r2 = x_pc**2 + (y_pc * np.cos(i) + ell_pc * np.sin(i)) ** 2
            z = y_pc * np.sin(i) - ell_pc * np.cos(i)
            return np.sqrt(max(r2, 0.0)), z

        def moment_integrand(ell_pc: float) -> float:
            r, z = los_geometry(ell_pc)
            nu = self.tracer.density(r, z)
            moments = self.jeans_moments(r, z)
            if r > 1e-10:
                v_star2 = moments.v_phi2 * x_pc**2 / r**2 + moments.v_r2 * (1.0 - x_pc**2 / r**2)
            else:
                v_star2 = moments.v_r2
            v_los2 = v_star2 * np.sin(i) ** 2 + moments.v_z2 * np.cos(i) ** 2
            return nu * max(v_los2, 0.0)

        denom = self.surface_density_analytic(x_pc, y_pc)
        if denom <= 0.0:
            return np.nan
        numer = quad(moment_integrand, -self.los_max_pc, self.los_max_pc, epsrel=self.epsrel, limit=100)[0]
        return max(float(numer / denom), 1e-10)


def numeric_surface_density(projector: AxisymmetricJeansProjector, x_pc: float, y_pc: float) -> float:
    i = projector.inclination_rad

    def los_geometry(ell_pc: float) -> tuple[float, float]:
        r2 = x_pc**2 + (y_pc * np.cos(i) + ell_pc * np.sin(i)) ** 2
        z = y_pc * np.sin(i) - ell_pc * np.cos(i)
        return np.sqrt(max(r2, 0.0)), z

    def surface_integrand(ell_pc: float) -> float:
        r, z = los_geometry(ell_pc)
        return projector.tracer.density(r, z)

    return quad(surface_integrand, -projector.los_max_pc, projector.los_max_pc, epsrel=projector.epsrel, limit=100)[0]


def make_analytic_projector(
    galaxy: GalaxyData,
    *,
    zmax_factor: float,
    los_factor: float,
    epsrel: float,
) -> AnalyticSurfaceDensityProjector:
    params = make_reference_hernquist_params(galaxy)
    numeric_projector = build_hernquist_projector(
        galaxy,
        params,
        zmax_factor=zmax_factor,
        los_factor=los_factor,
        epsrel=epsrel,
    )
    return AnalyticSurfaceDensityProjector(
        b_star_pc=numeric_projector.b_star_pc,
        qprime=numeric_projector.qprime,
        inclination_rad=numeric_projector.inclination_rad,
        beta_z=numeric_projector.beta_z,
        halo=numeric_projector.halo,
        zmax_factor=zmax_factor,
        los_factor=los_factor,
        epsrel=epsrel,
    )


def timed_sigma(projector: AxisymmetricJeansProjector, x_pc: np.ndarray, y_pc: np.ndarray, repeats: int) -> tuple[float, np.ndarray]:
    values = None
    start = time.perf_counter()
    for _ in range(repeats):
        values = projector.sigma_los2_many(x_pc, y_pc)
    elapsed = time.perf_counter() - start
    assert values is not None
    return elapsed / repeats, values


def timed_surface_density(
    numeric_projector: AxisymmetricJeansProjector,
    analytic_projector: AnalyticSurfaceDensityProjector,
    x_pc: np.ndarray,
    y_pc: np.ndarray,
    repeats: int,
) -> tuple[float, np.ndarray, float, np.ndarray]:
    numeric_values = None
    start = time.perf_counter()
    for _ in range(repeats):
        numeric_values = np.array(
            [numeric_surface_density(numeric_projector, float(x), float(y)) for x, y in zip(x_pc, y_pc)],
            dtype=float,
        )
    numeric_elapsed = (time.perf_counter() - start) / repeats

    analytic_values = None
    start = time.perf_counter()
    for _ in range(repeats):
        analytic_values = np.array(
            [analytic_projector.surface_density_analytic(float(x), float(y)) for x, y in zip(x_pc, y_pc)],
            dtype=float,
        )
    analytic_elapsed = (time.perf_counter() - start) / repeats

    assert numeric_values is not None
    assert analytic_values is not None
    return numeric_elapsed, numeric_values, analytic_elapsed, analytic_values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data/galaxies")
    parser.add_argument("--galaxy", default="Horologium I")
    parser.add_argument("--max-stars", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--surface-repeats", type=int, default=20)
    parser.add_argument("--zmax-factor", type=float, default=8.0)
    parser.add_argument("--los-factor", type=float, default=8.0)
    parser.add_argument("--epsrel", type=float, default=5e-2)
    args = parser.parse_args()

    galaxies = load_low_risk_galaxies(Path(args.data_dir))
    galaxy = galaxies[args.galaxy]
    stars = galaxy.stars.head(args.max_stars).copy()
    galaxy = GalaxyData(observables=galaxy.observables, stars=stars)

    params = make_reference_hernquist_params(galaxy)
    numeric_projector = build_hernquist_projector(
        galaxy,
        params,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        epsrel=args.epsrel,
    )
    analytic_projector = make_analytic_projector(
        galaxy,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        epsrel=args.epsrel,
    )

    x_pc = galaxy.x_pc
    y_pc = galaxy.y_pc
    surface_numeric_time, surface_numeric, surface_analytic_time, surface_analytic = timed_surface_density(
        numeric_projector,
        analytic_projector,
        x_pc,
        y_pc,
        args.surface_repeats,
    )
    numeric_time, sigma_numeric = timed_sigma(numeric_projector, x_pc, y_pc, args.repeats)
    analytic_time, sigma_analytic = timed_sigma(analytic_projector, x_pc, y_pc, args.repeats)

    likelihood = GaussianVelocityLikelihood(galaxy)
    ll_numeric = likelihood.log_likelihood(sigma_numeric, systemic_velocity_kms=params.systemic_velocity_kms)
    ll_analytic = likelihood.log_likelihood(sigma_analytic, systemic_velocity_kms=params.systemic_velocity_kms)

    abs_sigma = np.abs(sigma_analytic - sigma_numeric)
    rel_sigma = abs_sigma / np.maximum(np.abs(sigma_numeric), 1e-30)
    abs_surface = np.abs(surface_analytic - surface_numeric)
    rel_surface = abs_surface / np.maximum(np.abs(surface_numeric), 1e-30)
    speedup = numeric_time / analytic_time if analytic_time > 0 else np.inf
    surface_speedup = surface_numeric_time / surface_analytic_time if surface_analytic_time > 0 else np.inf

    print(f"galaxy={galaxy.observables.galaxy}")
    print(f"n_stars={len(galaxy.stars)} repeats={args.repeats}")
    print(f"zmax_factor={args.zmax_factor:g} los_factor={args.los_factor:g} epsrel={args.epsrel:g}")
    print(f"surface_repeats={args.surface_repeats}")
    print(f"surface_numeric_time_s={surface_numeric_time:.9f}")
    print(f"surface_analytic_time_s={surface_analytic_time:.9f}")
    print(f"surface_speedup={surface_speedup:.3f}x")
    print(f"surface_max_abs_delta={abs_surface.max():.9e}")
    print(f"surface_max_rel_delta={rel_surface.max():.9e}")
    print(f"numeric_time_s={numeric_time:.6f}")
    print(f"analytic_time_s={analytic_time:.6f}")
    print(f"speedup={speedup:.3f}x")
    print(f"loglike_numeric={ll_numeric:.9f}")
    print(f"loglike_analytic={ll_analytic:.9f}")
    print(f"loglike_delta={ll_analytic - ll_numeric:.9e}")
    print(f"sigma_los2_numeric={','.join(f'{value:.9g}' for value in sigma_numeric)}")
    print(f"sigma_los2_analytic={','.join(f'{value:.9g}' for value in sigma_analytic)}")
    print(f"sigma_los2_max_abs_delta={abs_sigma.max():.9e}")
    print(f"sigma_los2_max_rel_delta={rel_sigma.max():.9e}")


if __name__ == "__main__":
    main()
