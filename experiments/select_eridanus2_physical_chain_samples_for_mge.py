#!/usr/bin/env python
"""Select prior-valid Eridanus II chain samples with positive corrected MGE moments."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from experiments.test_eridanus2_fast_vs_strict_likelihood import (
    density_column_names,
    density_profile_features,
    select_diverse_density_profiles,
)
from experiments.test_eridanus2_mge_jeans_second_moment import (
    fit_positive_mge,
    generalized_halo_density_unit,
    mge_los_second_moment,
)
from hayashi_jeans.data import load_galaxy_data
from hayashi_jeans.tracer import intrinsic_q_from_projected
from scripts.run_willman1_block_mh_fast import (
    beta_z_from_q,
    log_prior_vector,
    resolve_systemic_velocity_prior_bounds,
    row_to_vector,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain", default=str(PROJECT_ROOT / "outputs/eridanus_ii_nautilus_5w_chain.csv"))
    parser.add_argument("--galaxy-csv", default=str(PROJECT_ROOT / "data/galaxies/08_Eridanus_II.csv"))
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument("--n-samples", type=int, default=5)
    parser.add_argument("--candidate-limit", type=int, default=80)
    parser.add_argument("--seed", type=int, default=20260601)
    parser.add_argument("--exclude-selected", default=None)
    parser.add_argument("--output", default=str(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_chain_physical_mge_selected_samples.csv"))
    args = parser.parse_args()

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    systemic_bounds = resolve_systemic_velocity_prior_bounds(
        galaxy,
        mode="adaptive",
        manual_bounds=None,
        padding=20.0,
        min_half_width=20.0,
    )
    chain = pd.read_csv(args.chain).replace([np.inf, -np.inf], np.nan)
    required = [
        "q_halo",
        "log10_b_halo_pc",
        "log10_rho0_msun_pc3",
        "minus_log10_one_minus_beta_z",
        "alpha",
        "beta",
        "gamma",
        "i_deg",
        "systemic_velocity_kms",
        "log_probability",
    ]
    finite = chain.dropna(subset=required).copy()
    finite = finite.loc[
        finite.apply(
            lambda row: np.isfinite(
                log_prior_vector(galaxy, row_to_vector(row), systemic_velocity_bounds=systemic_bounds)
            ),
            axis=1,
        )
    ].copy()
    if args.exclude_selected is not None:
        excluded = pd.read_csv(args.exclude_selected)
        finite = finite.drop(index=excluded["source_row"].astype(int).tolist(), errors="ignore")

    candidate_count = min(args.candidate_limit, len(finite))
    candidates = select_diverse_density_profiles(finite, n_samples=candidate_count, seed=args.seed)
    tracer_mge = fit_positive_mge(
        lambda m: (1.0 + (m / galaxy.observables.b_star_pc) ** 2) ** (-2.5),
        n_gauss=45,
        n_fit_radii=800,
        r_min_pc=1e-2,
        r_max_pc=2e5,
    )

    selected_rows: list[pd.Series] = []
    physical_flags = []
    for candidate in candidates.itertuples(index=False):
        row = pd.Series(candidate._asdict())
        mge_sigma2_unit = corrected_mge_sigma_unit(galaxy, row, tracer_mge)
        nonpositive = int(np.sum(~np.isfinite(mge_sigma2_unit) | (mge_sigma2_unit <= 0.0)))
        physical_flags.append(
            {
                "candidate_selection_id": int(row["selection_id"]),
                "source_row": int(row["source_row"]),
                "mge_nonpositive_count": nonpositive,
                "mge_sigma_los2_unit_min": float(np.nanmin(mge_sigma2_unit)),
            }
        )
        if nonpositive == 0:
            selected_rows.append(row)
            if len(selected_rows) == args.n_samples:
                break

    if len(selected_rows) < args.n_samples:
        raise RuntimeError(f"found only {len(selected_rows)} physical candidates among {candidate_count}")

    selected = pd.DataFrame(selected_rows).copy()
    selected["original_candidate_selection_id"] = selected["selection_id"].astype(int)
    selected["selection_id"] = np.arange(1, len(selected) + 1)
    feature_columns = density_column_names()
    if not all(column in selected.columns for column in feature_columns):
        features = density_profile_features(selected)
        for index, column in enumerate(feature_columns):
            selected[column] = features[:, index]
    selected.to_csv(args.output, index=False)

    flags_path = Path(args.output).with_name(Path(args.output).stem + "_candidate_flags.csv")
    pd.DataFrame(physical_flags).to_csv(flags_path, index=False)
    print(f"wrote {args.output}")
    print(selected[["selection_id", "source_row", "original_candidate_selection_id", "log_probability"]].to_string(index=False))
    print(f"wrote {flags_path}")


def corrected_mge_sigma_unit(galaxy, row: pd.Series, tracer_mge):
    inclination_rad = np.deg2rad(float(row["i_deg"]))
    q_star = intrinsic_q_from_projected(galaxy.observables.qprime, inclination_rad)
    halo_mge = fit_positive_mge(
        lambda m: generalized_halo_density_unit(
            m,
            b_pc=10.0 ** float(row["log10_b_halo_pc"]),
            alpha=float(row["alpha"]),
            beta=float(row["beta"]),
            gamma=float(row["gamma"]),
        ),
        n_gauss=45,
        n_fit_radii=800,
        r_min_pc=1e-2,
        r_max_pc=2e5,
    )
    mge_sigma2_unit, _, _, _ = mge_los_second_moment(
        x_pc=galaxy.x_pc,
        y_pc=galaxy.y_pc,
        halo_mge=halo_mge,
        tracer_mge=tracer_mge,
        q_halo=float(row["q_halo"]),
        q_star=q_star,
        beta_z=beta_z_from_q(float(row["minus_log10_one_minus_beta_z"])),
        inclination_rad=inclination_rad,
        n_u=96,
    )
    return mge_sigma2_unit


if __name__ == "__main__":
    main()
