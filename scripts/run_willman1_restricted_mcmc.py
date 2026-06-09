#!/usr/bin/env python
"""Run a restricted Willman 1 Jeans/MCMC chain.

This uses the actual Willman 1 star-by-star LOS velocity likelihood, but fixes
the Hernquist shape/geometric parameters so the current nested-integral Jeans
solver can finish in a reasonable time. The sampled parameters are:

- log10_rho0_msun_pc3
- systemic_velocity_kms

For fixed halo shape, Jeans second moments scale linearly with rho0, so the
expensive sigma_los^2 calculation is done once at rho0=1 Msun/pc^3.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import emcee
import numpy as np
import pandas as pd

from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.model import build_hernquist_projector
from hayashi_jeans.params import HayashiParameters


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy-csv", default="data/galaxies/27_Willman_1.csv")
    parser.add_argument("--centers-csv", default="data/processed/galaxy_structural_centers.csv")
    parser.add_argument("--output", default="outputs/willman1_restricted_posterior_samples.csv")
    parser.add_argument("--chain-output", default="outputs/willman1_restricted_mcmc_chain.csv")
    parser.add_argument("--n-walkers", type=int, default=32)
    parser.add_argument("--n-steps", type=int, default=2500)
    parser.add_argument("--burn-in", type=int, default=600)
    parser.add_argument("--seed", type=int, default=20260515)
    parser.add_argument("--zmax-factor", type=float, default=20.0)
    parser.add_argument("--los-factor", type=float, default=20.0)
    parser.add_argument("--epsrel", type=float, default=1.5e-2)
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    fixed = fixed_willman1_params(galaxy)
    unit_density_params = HayashiParameters(
        q_halo=fixed["q_halo"],
        b_halo_pc=10.0 ** fixed["log10_b_halo_pc"],
        rho0_msun_pc3=1.0,
        beta_z=fixed["beta_z"],
        alpha=fixed["alpha"],
        beta=fixed["beta"],
        gamma=fixed["gamma"],
        inclination_rad=np.deg2rad(fixed["inclination_deg"]),
        systemic_velocity_kms=float(np.mean(galaxy.velocity_kms)),
    )

    print("precomputing sigma_los2 at rho0=1 Msun/pc^3 ...", flush=True)
    projector = build_hernquist_projector(
        galaxy,
        unit_density_params,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        epsrel=args.epsrel,
    )
    sigma_los2_unit = projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    if not np.all(np.isfinite(sigma_los2_unit)) or np.any(sigma_los2_unit <= 0.0):
        raise RuntimeError("invalid sigma_los2_unit values")
    print(
        "sigma_los2_unit range="
        f"{sigma_los2_unit.min():.6g}..{sigma_los2_unit.max():.6g} (km/s)^2",
        flush=True,
    )

    likelihood = GaussianVelocityLikelihood(galaxy)
    rng = np.random.default_rng(args.seed)
    initial = np.column_stack(
        [
            rng.normal(fixed["initial_log10_rho0_msun_pc3"], 0.08, args.n_walkers),
            rng.normal(np.mean(galaxy.velocity_kms), 1.0, args.n_walkers),
        ]
    )

    def log_prior(theta: np.ndarray) -> float:
        log_rho0, systemic = theta
        if -4.5 <= log_rho0 <= 0.5 and -60.0 <= systemic <= 30.0:
            return 0.0
        return -np.inf

    def log_probability(theta: np.ndarray) -> float:
        prior = log_prior(theta)
        if not np.isfinite(prior):
            return -np.inf
        log_rho0, systemic = theta
        sigma_los2 = sigma_los2_unit * 10.0**log_rho0
        return prior + likelihood.log_likelihood(sigma_los2, systemic_velocity_kms=float(systemic))

    print("running emcee ...", flush=True)
    sampler = emcee.EnsembleSampler(args.n_walkers, 2, log_probability)
    sampler.run_mcmc(initial, args.n_steps, progress=False)
    flat = sampler.get_chain(discard=args.burn_in, flat=True)
    flat_log_prob = sampler.get_log_prob(discard=args.burn_in, flat=True)

    chain = pd.DataFrame(
        {
            "log10_rho0_msun_pc3": flat[:, 0],
            "systemic_velocity_kms": flat[:, 1],
            "log_probability": flat_log_prob,
        }
    )
    for key, value in fixed.items():
        chain[f"fixed_{key}"] = value
    chain["n_star"] = len(galaxy.stars)

    posterior_samples = pd.DataFrame(
        {
            "q_halo": fixed["q_halo"],
            "log10_b_halo_pc": fixed["log10_b_halo_pc"],
            "log10_rho0_msun_pc3": flat[:, 0],
            "alpha": fixed["alpha"],
            "beta": fixed["beta"],
            "gamma": fixed["gamma"],
            "systemic_velocity_kms": flat[:, 1],
            "log_probability": flat_log_prob,
        }
    )

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    chain.to_csv(args.chain_output, index=False)
    posterior_samples.to_csv(args.output, index=False)

    print(f"wrote {args.chain_output}")
    print(f"wrote {args.output}")
    print(f"n_selected={len(galaxy.stars)}")
    print(f"posterior_samples={len(posterior_samples)}")
    print(f"acceptance_fraction_mean={np.mean(sampler.acceptance_fraction):.4f}")
    print(f"log10_rho0_median={np.median(flat[:, 0]):.6f}")
    print(f"systemic_velocity_median={np.median(flat[:, 1]):.6f}")


def fixed_willman1_params(galaxy) -> dict[str, float]:
    return {
        "q_halo": 1.0,
        "log10_b_halo_pc": 3.2,
        "initial_log10_rho0_msun_pc3": -1.5,
        "beta_z": 0.3,
        "alpha": 1.8,
        "beta": 6.4,
        "gamma": 1.2,
        "inclination_deg": 75.0,
        "b_star_pc": galaxy.observables.b_star_pc,
        "qprime": galaxy.observables.qprime,
    }


if __name__ == "__main__":
    main()
