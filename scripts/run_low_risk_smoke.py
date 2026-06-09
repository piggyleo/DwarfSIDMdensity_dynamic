#!/usr/bin/env python
"""Run a small Hernquist/Jeans likelihood smoke test on low-risk galaxies."""

from __future__ import annotations

import argparse
from pathlib import Path

from hayashi_jeans.data import load_low_risk_galaxies
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.model import build_hernquist_projector, make_reference_hernquist_params


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data/galaxies")
    parser.add_argument("--galaxy", default="Horologium I")
    parser.add_argument("--max-stars", type=int, default=1)
    parser.add_argument("--zmax-factor", type=float, default=12.0)
    parser.add_argument("--los-factor", type=float, default=12.0)
    parser.add_argument("--epsrel", type=float, default=2e-2)
    args = parser.parse_args()

    galaxies = load_low_risk_galaxies(Path(args.data_dir))
    galaxy = galaxies[args.galaxy]
    stars = galaxy.stars.head(args.max_stars).copy()
    galaxy = type(galaxy)(observables=galaxy.observables, stars=stars)

    params = make_reference_hernquist_params(galaxy)
    projector = build_hernquist_projector(
        galaxy,
        params,
        zmax_factor=args.zmax_factor,
        los_factor=args.los_factor,
        epsrel=args.epsrel,
    )
    sigma_los2 = projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    loglike = GaussianVelocityLikelihood(galaxy).log_likelihood(sigma_los2)

    print(f"galaxy={galaxy.observables.galaxy}")
    print(f"n_stars={len(galaxy.stars)}")
    print(f"table1_n={galaxy.observables.n_sample_table1}")
    print(f"sigma_los2_kms2={','.join(f'{value:.6g}' for value in sigma_los2)}")
    print(f"profiled_log_likelihood={loglike:.6f}")


if __name__ == "__main__":
    main()
