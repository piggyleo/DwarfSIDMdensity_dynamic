#!/usr/bin/env python
"""Run one galaxy with the optional nautilus nested sampler adapter.

This script carries its own vector-level probability, prior, and plotting
helpers so it can run independently of ``run_willman1_block_mh_fast.py`` while
keeping the same Hayashi Jeans likelihood, projection code, and chain schema.
"""

from __future__ import annotations

import argparse
import inspect
import os
import re
import sys
import time
from pathlib import Path
from typing import Iterable
from dataclasses import dataclass

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / "outputs" / ".matplotlib"))
for path in (PROJECT_ROOT, PROJECT_ROOT / "src", PROJECT_ROOT / "scripts"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from experiments.test_halo_force_grid_interpolation import RZMomentGridWithForceGrid
from experiments.test_rz_moment_grid_interpolation import RZMomentGridProjector
from hayashi_jeans.data import GalaxyData, load_galaxy_data
from hayashi_jeans.likelihood import GaussianVelocityLikelihood
from hayashi_jeans.mge import (
    MGEPhysicalityConfig,
    mge_sigma_los2_for_halo,
    mge_sigma_los2_unit as shared_mge_sigma_los2_unit,
    mge_sigma_los2_unit_checked as shared_mge_sigma_los2_unit_checked,
)
from hayashi_jeans.halo_parameterizations import (
    HaloRunContext,
    SIDMPhysicalParameters,
    sidm_from_m200_c200,
    sidm_from_m200_ludlow,
    sidm_from_scale_parameters,
)
from hayashi_jeans.halos import HaloModel, SpheroidallyStratifiedHalo
from hayashi_jeans.model import build_hernquist_projector, build_projector
from hayashi_jeans.params import HayashiParameters


@dataclass(frozen=True)
class SlowParams:
    q_halo: float
    log10_b_halo_pc: float
    minus_log10_one_minus_beta_z: float
    alpha: float
    beta: float
    gamma: float
    inclination_deg: float


LIKELIHOOD_MODES = ("fast", "mge", "validation-halo", "validation-strict", "validation-convergence")
HALO_MODELS = ("generalized-hernquist", "sidm")
SIDM_PARAMETERIZATIONS = ("scale", "m200-c200", "m200-ludlow", "m200-ludlow-scatter")

PARAMETER_NAMES = [
    "q_halo",
    "log10_b_halo_pc",
    "log10_rho0_msun_pc3",
    "minus_log10_one_minus_beta_z",
    "alpha",
    "beta",
    "gamma",
    "i_deg",
    "systemic_velocity_kms",
]

CHAIN_COLUMNS = [
    "chain_id",
    "step",
    "q_halo",
    "log10_b_halo_pc",
    "log10_rho0_msun_pc3",
    "minus_log10_one_minus_beta_z",
    "beta_z",
    "alpha",
    "beta",
    "gamma",
    "i_deg",
    "systemic_velocity_kms",
    "log_probability",
    "sampler",
    "likelihood_mode",
]

SIDM_PARAMETER_NAMES = {
    "scale": [
        "q_halo",
        "log10_rs0_pc",
        "log10_rho_s0_msun_pc3",
        "tau",
        "minus_log10_one_minus_beta_z",
        "i_deg",
        "systemic_velocity_kms",
    ],
    "m200-c200": [
        "q_halo",
        "log10_m200_msun",
        "log10_c200",
        "tau",
        "minus_log10_one_minus_beta_z",
        "i_deg",
        "systemic_velocity_kms",
    ],
    "m200-ludlow": [
        "q_halo",
        "log10_m200_msun",
        "tau",
        "minus_log10_one_minus_beta_z",
        "i_deg",
        "systemic_velocity_kms",
    ],
    "m200-ludlow-scatter": [
        "q_halo",
        "log10_m200_msun",
        "concentration_scatter_sigma",
        "tau",
        "minus_log10_one_minus_beta_z",
        "i_deg",
        "systemic_velocity_kms",
    ],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        choices=("point", "smoke", "sample", "schema-check", "plot", "corner", "all"),
        default="all",
        help="point: one log-probability call; smoke: short nautilus run; sample: longer nautilus run.",
    )
    parser.add_argument("--likelihood-mode", choices=LIKELIHOOD_MODES, default="fast")
    parser.add_argument("--halo-model", choices=HALO_MODELS, default="generalized-hernquist")
    parser.add_argument("--sidm-parameterization", choices=SIDM_PARAMETERIZATIONS, default="scale")
    parser.add_argument(
        "--halo-redshift",
        type=float,
        default=None,
        help="Fixed halo redshift required by SIDM M200 parameterizations; it is not sampled.",
    )
    parser.add_argument("--halo-cosmology", choices=("planck15",), default="planck15")
    parser.add_argument(
        "--galaxy",
        default="willman1",
        help="Galaxy to analyze, e.g. willman1, 'Willman 1', coma_berenices, or 'Coma Berenices'.",
    )
    parser.add_argument(
        "--galaxy-csv",
        default=None,
        help="Optional explicit per-galaxy CSV. If omitted, --galaxy is resolved under data/galaxies.",
    )
    parser.add_argument("--centers-csv", default=str(PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv"))
    parser.add_argument(
        "--output-name",
        default=None,
        help=(
            "Optional common output filename stem. Defaults to <galaxy>_nautilus. "
            "Used only for outputs whose explicit --*-output path is not provided."
        ),
    )
    parser.add_argument(
        "--chain-output",
        default=None,
        help="Output CSV path. Defaults to outputs/<output-name>_chain.csv.",
    )
    parser.add_argument(
        "--checkpoint-output",
        default=None,
        help="Nautilus checkpoint path. Defaults to outputs/diagnostics/<output-name>_sampler.h5.",
    )
    parser.add_argument(
        "--figure-output",
        default=None,
        help="Density-profile figure path. Defaults to outputs/figures/<output-name>_density_profile.png.",
    )
    parser.add_argument(
        "--corner-output",
        default=None,
        help="Corner-plot figure path. Defaults to outputs/figures/<output-name>_corner.png.",
    )
    parser.add_argument("--burn", type=int, default=0, help="Rows to discard before plotting. Use 0 for nautilus.")
    systemic_corner = parser.add_mutually_exclusive_group()
    systemic_corner.add_argument(
        "--show-systemic-in-corner",
        dest="show_systemic_in_corner",
        action="store_true",
        help="Include systemic_velocity_kms in corner plots.",
    )
    systemic_corner.add_argument(
        "--hide-systemic-in-corner",
        dest="show_systemic_in_corner",
        action="store_false",
        help="Hide systemic_velocity_kms in corner plots. This is the default.",
    )
    parser.set_defaults(show_systemic_in_corner=False)
    parser.add_argument("--seed", type=int, default=20260526)
    parser.add_argument("--n-r", type=int, default=96)
    parser.add_argument("--n-z", type=int, default=192)
    parser.add_argument("--n-los", type=int, default=160)
    parser.add_argument("--n-force-r", type=int, default=48)
    parser.add_argument("--n-force-z", type=int, default=96)
    parser.add_argument("--strict-epsrel", type=float, default=1.5e-2)
    parser.add_argument(
        "--no-mge-physicality-check",
        action="store_true",
        help="Disable the default central-grid local v_phi^2 MGE/JAM physicality guard.",
    )
    parser.add_argument(
        "--mge-n-gauss-halo",
        type=int,
        default=None,
        help="Halo Gaussian count. Defaults to 60 for SIDM and 45 otherwise.",
    )
    parser.add_argument(
        "--mge-n-gauss-tracer",
        type=int,
        default=None,
        help="Tracer Gaussian count. Defaults to 60 for SIDM and 45 otherwise.",
    )
    parser.add_argument("--mge-decomposition-terms", type=int, default=28)
    parser.add_argument("--mge-n-u", type=int, default=96)
    parser.add_argument("--mge-r-min-pc", type=float, default=1e-2)
    parser.add_argument("--mge-r-max-pc", type=float, default=2e5)
    parser.add_argument("--mge-real-fit-radii", type=int, default=1600)
    parser.add_argument("--mge-real-lstsq-rcond", type=float, default=1e-12)
    parser.add_argument("--systemic-prior", choices=("adaptive", "legacy", "manual"), default="adaptive")
    parser.add_argument("--systemic-prior-padding", type=float, default=20.0)
    parser.add_argument("--systemic-prior-min-half-width", type=float, default=20.0)
    parser.add_argument(
        "--systemic-prior-bounds",
        default=None,
        help="Manual systemic velocity prior bounds as 'low,high' in km/s; requires --systemic-prior manual.",
    )
    parser.add_argument("--n-live", type=int, default=120)
    parser.add_argument("--n-eff", type=float, default=300.0)
    parser.add_argument("--n-like-max", type=float, default=np.inf)
    parser.add_argument("--f-live", type=float, default=0.05)
    parser.add_argument("--n-shell", type=int, default=1)
    parser.add_argument("--n-batch", type=int, default=None)
    parser.add_argument("--n-like-new-bound", type=int, default=None)
    parser.add_argument(
        "--n-processes",
        "--n-workers",
        dest="n_processes",
        type=int,
        default=1,
        help="Number of worker processes used by nautilus for likelihood evaluations.",
    )
    parser.add_argument("--timeout", type=float, default=np.inf)
    parser.add_argument("--resume", action="store_true", help="Resume from --checkpoint-output if it exists.")
    parser.add_argument("--no-checkpoint", action="store_true", help="Do not write a nautilus .h5 checkpoint.")
    parser.add_argument("--weighted-output", action="store_true", help="Write unequal-weight posterior rows plus log weights.")
    parser.add_argument(
        "--use-sample-weights",
        action="store_true",
        help="Use log_weight/weight columns for posterior summaries, density profiles, and corner plots.",
    )
    parser.add_argument("--equal-weight-boost", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--progress-newlines",
        action="store_true",
        help="With --verbose, print nautilus status updates on separate lines so Slurm .out files keep progress history.",
    )
    parser.add_argument(
        "--allow-output-galaxy-mismatch",
        action="store_true",
        help="Allow output paths that do not contain the resolved galaxy slug.",
    )
    parser.add_argument("--smoke-n-live", type=int, default=30)
    parser.add_argument("--smoke-n-eff", type=float, default=20.0)
    parser.add_argument("--smoke-n-like-max", type=float, default=120.0)
    parser.add_argument("--smoke-timeout", type=float, default=300.0)
    args = parser.parse_args()
    args.use_mge_physicality_check = not args.no_mge_physicality_check
    if args.mge_n_gauss_halo is None:
        args.mge_n_gauss_halo = 60 if args.halo_model == "sidm" else 45
    if args.mge_n_gauss_tracer is None:
        args.mge_n_gauss_tracer = 60 if args.halo_model == "sidm" else 45
    args.halo_run_context = HaloRunContext(
        redshift=args.halo_redshift,
        cosmology_name=args.halo_cosmology,
    )
    sampling_modes = {"point", "smoke", "sample", "schema-check", "all"}
    if (
        args.mode in sampling_modes
        and args.halo_model == "sidm"
        and args.sidm_parameterization != "scale"
    ):
        args.halo_run_context.require_redshift()

    galaxy_csv = Path(args.galaxy_csv) if args.galaxy_csv is not None else resolve_galaxy_csv(args.galaxy)
    args.galaxy_csv = str(galaxy_csv)
    args.galaxy_slug = galaxy_slug_from_csv(galaxy_csv)
    output_name = output_slug(args.output_name) if args.output_name is not None else f"{args.galaxy_slug}_nautilus"
    if args.chain_output is None:
        args.chain_output = str(PROJECT_ROOT / "outputs" / f"{output_name}_chain.csv")
    if args.checkpoint_output is None:
        args.checkpoint_output = str(PROJECT_ROOT / "outputs" / "diagnostics" / f"{output_name}_sampler.h5")
    if args.figure_output is None:
        args.figure_output = str(PROJECT_ROOT / "outputs" / "figures" / f"{output_name}_density_profile.png")
    if args.corner_output is None:
        args.corner_output = str(PROJECT_ROOT / "outputs" / "figures" / f"{output_name}_corner.png")
    validate_output_paths(args)

    galaxy = load_galaxy_data(args.galaxy_csv, structural_centers_csv=args.centers_csv)
    args.systemic_velocity_bounds = resolve_systemic_velocity_prior_bounds(
        galaxy,
        mode=args.systemic_prior,
        manual_bounds=args.systemic_prior_bounds,
        padding=args.systemic_prior_padding,
        min_half_width=args.systemic_prior_min_half_width,
    )
    print(f"galaxy={galaxy.observables.galaxy}")
    print(f"galaxy_slug={args.galaxy_slug}")
    print(f"galaxy_csv={args.galaxy_csv}")
    print(f"systemic_velocity_prior=({args.systemic_velocity_bounds[0]:.6g}, {args.systemic_velocity_bounds[1]:.6g})")
    print(f"mge_physicality_check={args.use_mge_physicality_check}")
    print(f"halo_model={args.halo_model}")
    if args.halo_model == "sidm":
        print(f"sidm_parameterization={args.sidm_parameterization}")
        print(f"halo_redshift={args.halo_redshift}")

    if args.mode in ("point", "all"):
        run_point_check(galaxy, args)

    if args.mode in ("schema-check", "all"):
        run_schema_check(galaxy, args)

    if args.mode in ("smoke", "all"):
        smoke_output = Path(args.chain_output).with_name(Path(args.chain_output).stem + "_smoke.csv")
        run_nautilus(
            galaxy,
            args,
            chain_output=smoke_output,
            n_live=args.smoke_n_live,
            n_eff=args.smoke_n_eff,
            n_like_max=args.smoke_n_like_max,
            timeout=args.smoke_timeout,
            checkpoint_output=None,
            resume=False,
            run_label="smoke",
        )

    if args.mode in ("sample", "all"):
        checkpoint_output = None if args.no_checkpoint else Path(args.checkpoint_output)
        run_nautilus(
            galaxy,
            args,
            chain_output=Path(args.chain_output),
            n_live=args.n_live,
            n_eff=args.n_eff,
            n_like_max=args.n_like_max,
            timeout=args.timeout,
            checkpoint_output=checkpoint_output,
            resume=args.resume,
            run_label="nautilus",
        )

    if args.mode in ("plot", "all"):
        chain_path = Path(args.chain_output)
        if not chain_path.exists():
            raise FileNotFoundError(f"{chain_path} does not exist; run --mode sample first or pass --chain-output")
        chain = pd.read_csv(chain_path)
        plot_density_profile(galaxy, chain, Path(args.figure_output), burn=args.burn, use_sample_weights=args.use_sample_weights)

    if args.mode in ("corner", "all"):
        chain_path = Path(args.chain_output)
        if not chain_path.exists():
            raise FileNotFoundError(f"{chain_path} does not exist; run --mode sample first or pass --chain-output")
        chain = pd.read_csv(chain_path)
        generate_corner_map(
            chain,
            Path(args.corner_output),
            show_systemic_velocity=args.show_systemic_in_corner,
            burn=args.burn,
            use_sample_weights=args.use_sample_weights,
        )


def run_point_check(galaxy: GalaxyData, args: argparse.Namespace) -> None:
    vector = initial_vector(
        galaxy,
        halo_model=args.halo_model,
        sidm_parameterization=args.sidm_parameterization,
    )
    t0 = time.perf_counter()
    logp = evaluate_log_probability(galaxy, vector, args)
    elapsed = time.perf_counter() - t0
    prior = log_prior_vector(
        galaxy,
        vector,
        systemic_velocity_bounds=args.systemic_velocity_bounds,
        halo_model=args.halo_model,
        sidm_parameterization=args.sidm_parameterization,
    )
    print("single_point_check")
    print(f"vector={dict(zip(parameter_names_for(args.halo_model, args.sidm_parameterization), vector))}")
    print(f"log_prior={prior:.6f}")
    print(f"log_probability={logp:.6f}")
    print(f"seconds={elapsed:.6f}")


def validate_output_paths(args: argparse.Namespace) -> None:
    if args.allow_output_galaxy_mismatch:
        return
    output_attrs = ("chain_output", "checkpoint_output", "figure_output", "corner_output")
    mismatched = []
    for attr in output_attrs:
        value = getattr(args, attr, None)
        if value is None:
            continue
        if args.galaxy_slug not in Path(value).name:
            mismatched.append((attr.replace("_", "-"), value))
    if mismatched:
        details = "\n".join(f"  --{name}: {value}" for name, value in mismatched)
        raise ValueError(
            f"output path(s) do not include resolved galaxy slug {args.galaxy_slug!r}:\n"
            f"{details}\n"
            "Use the default output paths, rename the files to include the galaxy slug, "
            "or pass --allow-output-galaxy-mismatch if this is intentional."
        )


def normalize_galaxy_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def galaxy_alias_keys(value: str) -> set[str]:
    raw = value.lower()
    aliases = {raw}
    roman_to_number = {
        " viii": " 8",
        " vii": " 7",
        " vi": " 6",
        " v": " 5",
        " iv": " 4",
        " iii": " 3",
        " ii": " 2",
        " i": " 1",
    }
    number_to_roman = {
        " 8": " viii",
        " 7": " vii",
        " 6": " vi",
        " 5": " v",
        " 4": " iv",
        " 3": " iii",
        " 2": " ii",
        " 1": " i",
    }
    spaced = re.sub(r"[_-]+", " ", raw)
    aliases.add(spaced)
    for roman, number in roman_to_number.items():
        if spaced.endswith(roman):
            aliases.add(spaced[: -len(roman)] + number)
    for number, roman in number_to_roman.items():
        if spaced.endswith(number):
            aliases.add(spaced[: -len(number)] + roman)
    return {normalize_galaxy_key(alias) for alias in aliases}


def output_slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    return slug or "galaxy"


def resolve_galaxy_csv(galaxy: str) -> Path:
    data_dir = PROJECT_ROOT / "data" / "galaxies"
    target_keys = galaxy_alias_keys(galaxy)
    candidates: list[tuple[str, Path]] = []
    for csv_path in sorted(data_dir.glob("*.csv")):
        if csv_path.name.startswith("_"):
            continue
        try:
            frame = pd.read_csv(csv_path, comment="#", nrows=1)
        except Exception:
            continue
        if frame.empty or "galaxy" not in frame:
            continue
        galaxy_name = str(frame.loc[0, "galaxy"])
        candidates.append((galaxy_name, csv_path))
        filename_key = re.sub(r"^\d+_", "", csv_path.stem)
        keys = galaxy_alias_keys(galaxy_name) | galaxy_alias_keys(filename_key) | galaxy_alias_keys(csv_path.stem)
        if target_keys & keys:
            return csv_path
    available = ", ".join(name for name, _path in candidates)
    raise ValueError(f"could not resolve --galaxy {galaxy!r}; available galaxies: {available}")


def galaxy_slug_from_csv(galaxy_csv: Path) -> str:
    frame = pd.read_csv(galaxy_csv, comment="#", nrows=1)
    if frame.empty or "galaxy" not in frame:
        return output_slug(re.sub(r"^\d+_", "", galaxy_csv.stem))
    return output_slug(str(frame.loc[0, "galaxy"]))


def validation_grid_settings(
    likelihood_mode: str,
    *,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
) -> dict[str, int | None]:
    if likelihood_mode == "mge":
        return {"n_r": None, "n_z": None, "n_los": None, "n_force_r": None, "n_force_z": None}
    if likelihood_mode == "fast":
        return {"n_r": n_r, "n_z": n_z, "n_los": n_los, "n_force_r": n_force_r, "n_force_z": n_force_z}
    if likelihood_mode == "validation-halo":
        return {"n_r": 128, "n_z": 256, "n_los": 200, "n_force_r": None, "n_force_z": None}
    if likelihood_mode == "validation-convergence":
        return {"n_r": 128, "n_z": 256, "n_los": 200, "n_force_r": 64, "n_force_z": 128}
    if likelihood_mode == "validation-strict":
        return {"n_r": None, "n_z": None, "n_los": None, "n_force_r": None, "n_force_z": None}
    raise ValueError(f"unknown likelihood_mode={likelihood_mode!r}")


def parameter_names_for(halo_model: str, sidm_parameterization: str) -> list[str]:
    if halo_model == "generalized-hernquist":
        return list(PARAMETER_NAMES)
    return list(SIDM_PARAMETER_NAMES[sidm_parameterization])


def sidm_parameters_from_vector(
    vector: np.ndarray,
    *,
    parameterization: str,
    context: HaloRunContext,
) -> tuple[SIDMPhysicalParameters, float, float, float]:
    values = dict(zip(SIDM_PARAMETER_NAMES[parameterization], np.asarray(vector, dtype=float)))
    common = {
        "q_halo": values["q_halo"],
        "tau": values["tau"],
    }
    if parameterization == "scale":
        physical = sidm_from_scale_parameters(
            **common,
            rs0_pc=10.0 ** values["log10_rs0_pc"],
            rho_s0_msun_pc3=10.0 ** values["log10_rho_s0_msun_pc3"],
        )
    elif parameterization == "m200-c200":
        physical = sidm_from_m200_c200(
            **common,
            m200_msun=10.0 ** values["log10_m200_msun"],
            c200=10.0 ** values["log10_c200"],
            context=context,
        )
    else:
        physical = sidm_from_m200_ludlow(
            **common,
            m200_msun=10.0 ** values["log10_m200_msun"],
            context=context,
            scatter_sigma=values.get("concentration_scatter_sigma", 0.0),
        )
    return (
        physical,
        beta_z_from_q(values["minus_log10_one_minus_beta_z"]),
        np.deg2rad(values["i_deg"]),
        values["systemic_velocity_kms"],
    )


def compute_sigma_for_halo(
    galaxy: GalaxyData,
    *,
    halo: SpheroidallyStratifiedHalo,
    beta_z: float,
    inclination_rad: float,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    use_mge_physicality_check: bool,
    mge_config: MGEPhysicalityConfig,
) -> np.ndarray:
    if likelihood_mode == "mge":
        return mge_sigma_los2_for_halo(
            galaxy,
            halo=halo,
            beta_z=beta_z,
            inclination_rad=inclination_rad,
            config=mge_config,
            check_physicality=use_mge_physicality_check,
        )

    params = HayashiParameters(
        q_halo=halo.q,
        b_halo_pc=1.0,
        rho0_msun_pc3=1.0,
        beta_z=beta_z,
        alpha=1.0,
        beta=4.0,
        gamma=1.0,
        inclination_rad=inclination_rad,
        systemic_velocity_kms=0.0,
    )
    settings = validation_grid_settings(
        likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
    )
    if likelihood_mode == "validation-strict":
        projector = build_projector(
            galaxy,
            params,
            halo=halo,
            zmax_factor=20.0,
            los_factor=20.0,
            epsrel=strict_epsrel,
        )
    elif likelihood_mode == "validation-halo":
        projector = RZMomentGridProjector(
            galaxy=galaxy,
            params=params,
            halo=halo,
            n_r=int(settings["n_r"]),
            n_z=int(settings["n_z"]),
            n_los=int(settings["n_los"]),
            zmax_factor=20.0,
            los_factor=20.0,
            r_min_pc=1e-3,
            z_min_pc=1e-3,
            grid_padding=1.08,
            radial_derivative="cubic_spline",
        )
    else:
        projector = RZMomentGridWithForceGrid(
            galaxy=galaxy,
            params=params,
            halo=halo,
            n_r=int(settings["n_r"]),
            n_z=int(settings["n_z"]),
            n_los=int(settings["n_los"]),
            n_force_r=int(settings["n_force_r"]),
            n_force_z=int(settings["n_force_z"]),
            zmax_factor=20.0,
            los_factor=20.0,
            r_min_pc=1e-3,
            z_min_pc=1e-3,
            grid_padding=1.08,
        )
    return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)


def compute_sigma_unit(
    galaxy: GalaxyData,
    slow: SlowParams,
    *,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    use_mge_physicality_check: bool = True,
    mge_config: MGEPhysicalityConfig = MGEPhysicalityConfig(),
) -> np.ndarray:
    params = HayashiParameters(
        q_halo=slow.q_halo,
        b_halo_pc=10.0**slow.log10_b_halo_pc,
        rho0_msun_pc3=1.0,
        beta_z=beta_z_from_q(slow.minus_log10_one_minus_beta_z),
        alpha=slow.alpha,
        beta=slow.beta,
        gamma=slow.gamma,
        inclination_rad=np.deg2rad(slow.inclination_deg),
        systemic_velocity_kms=0.0,
    )
    if likelihood_mode == "mge":
        mge_func = shared_mge_sigma_los2_unit_checked if use_mge_physicality_check else shared_mge_sigma_los2_unit
        return mge_func(
            galaxy,
            q_halo=params.q_halo,
            b_halo_pc=params.b_halo_pc,
            alpha=params.alpha,
            beta=params.beta,
            gamma=params.gamma,
            beta_z=params.beta_z,
            inclination_rad=params.inclination_rad,
            config=mge_config,
        )
    settings = validation_grid_settings(
        likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
    )
    if likelihood_mode == "validation-strict":
        projector = build_hernquist_projector(
            galaxy,
            params,
            zmax_factor=20.0,
            los_factor=20.0,
            epsrel=strict_epsrel,
        )
        return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)
    if likelihood_mode == "validation-halo":
        projector = RZMomentGridProjector(
            galaxy=galaxy,
            params=params,
            n_r=int(settings["n_r"]),
            n_z=int(settings["n_z"]),
            n_los=int(settings["n_los"]),
            zmax_factor=20.0,
            los_factor=20.0,
            r_min_pc=1e-3,
            z_min_pc=1e-3,
            grid_padding=1.08,
            radial_derivative="cubic_spline",
        )
        return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)

    projector = RZMomentGridWithForceGrid(
        galaxy=galaxy,
        params=params,
        n_r=int(settings["n_r"]),
        n_z=int(settings["n_z"]),
        n_los=int(settings["n_los"]),
        n_force_r=int(settings["n_force_r"]),
        n_force_z=int(settings["n_force_z"]),
        zmax_factor=20.0,
        los_factor=20.0,
        r_min_pc=1e-3,
        z_min_pc=1e-3,
        grid_padding=1.08,
    )
    return projector.sigma_los2_many(galaxy.x_pc, galaxy.y_pc)


def log_likelihood_from_unit_sigma(
    galaxy: GalaxyData,
    sigma_unit: np.ndarray,
    *,
    log10_rho0_msun_pc3: float,
    systemic_velocity_kms: float,
) -> float:
    sigma_los2 = sigma_unit * 10.0**log10_rho0_msun_pc3
    return GaussianVelocityLikelihood(galaxy).log_likelihood(sigma_los2, systemic_velocity_kms=systemic_velocity_kms)


def initial_slow_params() -> SlowParams:
    return SlowParams(
        q_halo=1.0,
        log10_b_halo_pc=3.2,
        minus_log10_one_minus_beta_z=0.3,
        alpha=1.8,
        beta=6.4,
        gamma=1.2,
        inclination_deg=75.0,
    )


def log_prior_slow(slow: SlowParams, galaxy: GalaxyData) -> float:
    min_i = np.degrees(np.arccos(galaxy.observables.qprime))
    if not (0.1 <= slow.q_halo <= 2.0):
        return -np.inf
    if not (0.0 <= slow.log10_b_halo_pc <= 5.0):
        return -np.inf
    if not (-1.0 <= slow.minus_log10_one_minus_beta_z < 1.0):
        return -np.inf
    if not (0.5 <= slow.alpha <= 3.0):
        return -np.inf
    if not (3.0 <= slow.beta <= 10.0):
        return -np.inf
    if not (0.0 <= slow.gamma <= 2.0):
        return -np.inf
    if not (min_i < slow.inclination_deg <= 90.0):
        return -np.inf
    return 0.0


def beta_z_from_q(minus_log10_one_minus_beta_z: float) -> float:
    return 1.0 - 10.0 ** (-float(minus_log10_one_minus_beta_z))


def q_from_beta_z(beta_z: float) -> float:
    return -float(np.log10(1.0 - float(beta_z)))


def log_prior_fast(
    log10_rho0_msun_pc3: float,
    systemic_velocity_kms: float,
    *,
    systemic_velocity_bounds: tuple[float, float],
) -> float:
    if not (-5.0 <= log10_rho0_msun_pc3 <= 5.0):
        return -np.inf
    low, high = systemic_velocity_bounds
    if not (low <= systemic_velocity_kms <= high):
        return -np.inf
    return 0.0


def systemic_velocity_prior_bounds(
    galaxy: GalaxyData,
    *,
    padding: float = 20.0,
    min_half_width: float = 20.0,
) -> tuple[float, float]:
    velocity = np.asarray(galaxy.velocity_kms, dtype=float)
    velocity = velocity[np.isfinite(velocity)]
    if velocity.size == 0:
        return (-80.0, 50.0)
    center = float(np.median(velocity))
    mad = float(np.median(np.abs(velocity - center)))
    robust_sigma = 1.4826 * mad
    minmax_half_width = 0.5 * float(np.max(velocity) - np.min(velocity))
    half_width = max(5.0 * robust_sigma, minmax_half_width + padding, min_half_width)
    return center - half_width, center + half_width


def parse_systemic_velocity_bounds(value: str) -> tuple[float, float]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) != 2:
        raise ValueError("--systemic-prior-bounds must have the form 'low,high'")
    low, high = (float(parts[0]), float(parts[1]))
    if not np.isfinite(low) or not np.isfinite(high) or not low < high:
        raise ValueError("--systemic-prior-bounds must be finite and satisfy low < high")
    return low, high


def resolve_systemic_velocity_prior_bounds(
    galaxy: GalaxyData,
    *,
    mode: str,
    manual_bounds: str | None,
    padding: float,
    min_half_width: float,
) -> tuple[float, float]:
    if mode == "legacy":
        return (-80.0, 50.0)
    if mode == "manual":
        if manual_bounds is None:
            raise ValueError("--systemic-prior manual requires --systemic-prior-bounds")
        return parse_systemic_velocity_bounds(manual_bounds)
    if mode == "adaptive":
        return systemic_velocity_prior_bounds(galaxy, padding=padding, min_half_width=min_half_width)
    raise ValueError(f"unknown systemic prior mode: {mode!r}")


def full_log_probability_vector(
    galaxy: GalaxyData,
    vector: np.ndarray,
    *,
    likelihood_mode: str,
    n_r: int,
    n_z: int,
    n_los: int,
    n_force_r: int,
    n_force_z: int,
    strict_epsrel: float,
    systemic_velocity_bounds: tuple[float, float],
    use_mge_physicality_check: bool,
    mge_n_gauss_halo: int,
    mge_n_gauss_tracer: int,
    mge_n_u: int,
    mge_r_min_pc: float,
    mge_r_max_pc: float,
    mge_real_fit_radii: int = 1600,
    mge_real_lstsq_rcond: float = 1.0e-12,
    mge_decomposition_terms: int = 28,
    halo_model: str = "generalized-hernquist",
    sidm_parameterization: str = "scale",
    halo_run_context: HaloRunContext = HaloRunContext(),
) -> float:
    prior = log_prior_vector(
        galaxy,
        vector,
        systemic_velocity_bounds=systemic_velocity_bounds,
        halo_model=halo_model,
        sidm_parameterization=sidm_parameterization,
    )
    if not np.isfinite(prior):
        return -np.inf
    mge_config = MGEPhysicalityConfig(
        n_gauss_halo=mge_n_gauss_halo,
        n_gauss_tracer=mge_n_gauss_tracer,
        decomposition_terms=mge_decomposition_terms,
        n_u=mge_n_u,
        r_min_pc=mge_r_min_pc,
        r_max_pc=mge_r_max_pc,
        real_fit_radii=mge_real_fit_radii,
        real_lstsq_rcond=mge_real_lstsq_rcond,
    )
    if halo_model == "sidm":
        try:
            physical, beta_z, inclination_rad, systemic_velocity = sidm_parameters_from_vector(
                vector,
                parameterization=sidm_parameterization,
                context=halo_run_context,
            )
            sigma_los2 = compute_sigma_for_halo(
                galaxy,
                halo=physical.build_halo(),
                beta_z=beta_z,
                inclination_rad=inclination_rad,
                likelihood_mode=likelihood_mode,
                n_r=n_r,
                n_z=n_z,
                n_los=n_los,
                n_force_r=n_force_r,
                n_force_z=n_force_z,
                strict_epsrel=strict_epsrel,
                use_mge_physicality_check=use_mge_physicality_check,
                mge_config=mge_config,
            )
        except (FloatingPointError, ValueError, ZeroDivisionError, RuntimeError):
            return -np.inf
        return prior + GaussianVelocityLikelihood(galaxy).log_likelihood(
            sigma_los2,
            systemic_velocity_kms=systemic_velocity,
        )

    slow = SlowParams(
        q_halo=float(vector[0]),
        log10_b_halo_pc=float(vector[1]),
        minus_log10_one_minus_beta_z=float(vector[3]),
        alpha=float(vector[4]),
        beta=float(vector[5]),
        gamma=float(vector[6]),
        inclination_deg=float(vector[7]),
    )
    sigma_unit = compute_sigma_unit(
        galaxy,
        slow,
        likelihood_mode=likelihood_mode,
        n_r=n_r,
        n_z=n_z,
        n_los=n_los,
        n_force_r=n_force_r,
        n_force_z=n_force_z,
        strict_epsrel=strict_epsrel,
        use_mge_physicality_check=use_mge_physicality_check,
        mge_config=mge_config,
    )
    return prior + log_likelihood_from_unit_sigma(
        galaxy,
        sigma_unit,
        log10_rho0_msun_pc3=float(vector[2]),
        systemic_velocity_kms=float(vector[8]),
    )


def log_prior_vector(
    galaxy: GalaxyData,
    vector: np.ndarray,
    *,
    systemic_velocity_bounds: tuple[float, float],
    halo_model: str = "generalized-hernquist",
    sidm_parameterization: str = "scale",
) -> float:
    if halo_model == "sidm":
        names = SIDM_PARAMETER_NAMES[sidm_parameterization]
        if len(vector) != len(names):
            return -np.inf
        values = dict(zip(names, np.asarray(vector, dtype=float)))
        min_i = np.degrees(np.arccos(galaxy.observables.qprime))
        if not (0.1 <= values["q_halo"] <= 2.0):
            return -np.inf
        if not (0.0 <= values["tau"] <= 1.08):
            return -np.inf
        if not (-1.0 <= values["minus_log10_one_minus_beta_z"] < 1.0):
            return -np.inf
        if not (min_i < values["i_deg"] <= 90.0):
            return -np.inf
        if not (systemic_velocity_bounds[0] <= values["systemic_velocity_kms"] <= systemic_velocity_bounds[1]):
            return -np.inf
        if sidm_parameterization == "scale":
            if not (0.0 <= values["log10_rs0_pc"] <= 5.0):
                return -np.inf
            if not (-5.0 <= values["log10_rho_s0_msun_pc3"] <= 5.0):
                return -np.inf
        else:
            if not (5.0 <= values["log10_m200_msun"] <= 12.0):
                return -np.inf
            if sidm_parameterization == "m200-c200" and not (0.0 <= values["log10_c200"] <= 2.0):
                return -np.inf
            if "concentration_scatter_sigma" in values and not (-3.0 <= values["concentration_scatter_sigma"] <= 3.0):
                return -np.inf
        return 0.0

    if len(vector) != 9:
        return -np.inf
    slow = SlowParams(
        q_halo=float(vector[0]),
        log10_b_halo_pc=float(vector[1]),
        minus_log10_one_minus_beta_z=float(vector[3]),
        alpha=float(vector[4]),
        beta=float(vector[5]),
        gamma=float(vector[6]),
        inclination_deg=float(vector[7]),
    )
    return log_prior_slow(slow, galaxy) + log_prior_fast(
        float(vector[2]),
        float(vector[8]),
        systemic_velocity_bounds=systemic_velocity_bounds,
    )


def initial_vector(
    galaxy: GalaxyData,
    *,
    halo_model: str = "generalized-hernquist",
    sidm_parameterization: str = "scale",
) -> np.ndarray:
    systemic = float(np.median(galaxy.velocity_kms))
    min_i = float(np.degrees(np.arccos(galaxy.observables.qprime)))
    inclination = min(max(75.0, min_i + 1.0), 90.0)
    if halo_model == "sidm":
        initial = {
            "q_halo": 1.0,
            "log10_rs0_pc": 3.2,
            "log10_rho_s0_msun_pc3": -1.0,
            "log10_m200_msun": 9.0,
            "log10_c200": 1.1,
            "concentration_scatter_sigma": 0.0,
            "tau": 0.5,
            "minus_log10_one_minus_beta_z": 0.3,
            "i_deg": inclination,
            "systemic_velocity_kms": systemic,
        }
        return np.array([initial[name] for name in SIDM_PARAMETER_NAMES[sidm_parameterization]])

    slow = initial_slow_params()
    return np.array([
        slow.q_halo,
        slow.log10_b_halo_pc,
        -1.961219,
        slow.minus_log10_one_minus_beta_z,
        slow.alpha,
        slow.beta,
        slow.gamma,
        slow.inclination_deg,
        systemic,
    ])


def vector_to_row(
    vector: np.ndarray,
    *,
    log_probability: float,
    step: int,
    chain_id: int,
    sampler_name: str,
    likelihood_mode: str,
    halo_model: str = "generalized-hernquist",
    sidm_parameterization: str = "scale",
    halo_run_context: HaloRunContext = HaloRunContext(),
) -> dict[str, float | int | str]:
    if halo_model == "sidm":
        physical, beta_z, _inclination_rad, systemic_velocity = sidm_parameters_from_vector(
            vector,
            parameterization=sidm_parameterization,
            context=halo_run_context,
        )
        values = dict(zip(SIDM_PARAMETER_NAMES[sidm_parameterization], np.asarray(vector, dtype=float)))
        row = {
            "chain_id": chain_id,
            "step": step,
            **values,
            "beta_z": beta_z,
            "systemic_velocity_kms": systemic_velocity,
            "halo_model": "sidm",
            "halo_parameterization": sidm_parameterization,
            "log10_rs0_pc": np.log10(physical.rs0_pc),
            "log10_rho_s0_msun_pc3": np.log10(physical.rho_s0_msun_pc3),
            "c200": physical.c200 if physical.c200 is not None else np.nan,
            "halo_redshift": physical.redshift if physical.redshift is not None else np.nan,
            "concentration_relation": physical.concentration_relation or "",
            "log_probability": float(log_probability),
            "sampler": sampler_name,
            "likelihood_mode": likelihood_mode,
        }
        return row

    return {
        "chain_id": chain_id,
        "step": step,
        "q_halo": float(vector[0]),
        "log10_b_halo_pc": float(vector[1]),
        "log10_rho0_msun_pc3": float(vector[2]),
        "minus_log10_one_minus_beta_z": float(vector[3]),
        "beta_z": beta_z_from_q(float(vector[3])),
        "alpha": float(vector[4]),
        "beta": float(vector[5]),
        "gamma": float(vector[6]),
        "i_deg": float(vector[7]),
        "systemic_velocity_kms": float(vector[8]),
        "log_probability": float(log_probability),
        "sampler": sampler_name,
        "likelihood_mode": likelihood_mode,
        "halo_model": "generalized-hernquist",
        "halo_parameterization": "generalized-hernquist",
    }


def row_to_vector(
    row: pd.Series,
    *,
    halo_model: str | None = None,
    sidm_parameterization: str | None = None,
) -> np.ndarray:
    inferred_model = halo_model or str(row.get("halo_model", "generalized-hernquist"))
    inferred_parameterization = sidm_parameterization or str(
        row.get("halo_parameterization", "scale")
    )
    names = parameter_names_for(inferred_model, inferred_parameterization)
    return np.array([float(row[name]) for name in names], dtype=float)


def run_schema_check(galaxy: GalaxyData, args: argparse.Namespace) -> None:
    vector = initial_vector(
        galaxy,
        halo_model=args.halo_model,
        sidm_parameterization=args.sidm_parameterization,
    )
    logp = evaluate_log_probability(galaxy, vector, args)
    row = vector_to_row(
        vector,
        log_probability=logp,
        step=1,
        chain_id=0,
        sampler_name="nautilus",
        likelihood_mode=args.likelihood_mode,
        halo_model=args.halo_model,
        sidm_parameterization=args.sidm_parameterization,
        halo_run_context=args.halo_run_context,
    )
    chain = pd.DataFrame([row])
    required_columns = CHAIN_COLUMNS if args.halo_model == "generalized-hernquist" else [
        "chain_id",
        "step",
        "halo_model",
        "halo_parameterization",
        "q_halo",
        "tau",
        "log10_rs0_pc",
        "log10_rho_s0_msun_pc3",
        "minus_log10_one_minus_beta_z",
        "beta_z",
        "i_deg",
        "systemic_velocity_kms",
        "log_probability",
        "sampler",
        "likelihood_mode",
    ]
    missing = [column for column in required_columns if column not in chain.columns]
    if missing:
        raise ValueError(f"schema check failed; missing columns: {missing}")
    print("schema_check_passed")
    print(",".join(chain.columns))


def run_nautilus(
    galaxy: GalaxyData,
    args: argparse.Namespace,
    *,
    chain_output: Path,
    n_live: int,
    n_eff: float,
    n_like_max: float,
    timeout: float,
    checkpoint_output: Path | None,
    resume: bool,
    run_label: str,
) -> pd.DataFrame:
    try:
        from nautilus import Prior, Sampler
    except ImportError as exc:  # pragma: no cover - optional sampler dependency
        raise ImportError("nautilus is required for this mode; install the 'nautilus-sampler' package.") from exc

    if n_live < 2:
        raise ValueError("n_live must be >= 2")
    prior = make_nautilus_prior(
        Prior,
        galaxy,
        args.systemic_velocity_bounds,
        halo_model=args.halo_model,
        sidm_parameterization=args.sidm_parameterization,
    )
    likelihood = make_log_probability_callable(galaxy, args)
    filepath = str(checkpoint_output) if checkpoint_output is not None else None
    if checkpoint_output is not None:
        checkpoint_output.parent.mkdir(parents=True, exist_ok=True)

    sampler = Sampler(
        prior,
        likelihood,
        n_live=n_live,
        n_batch=args.n_batch,
        n_like_new_bound=args.n_like_new_bound,
        pool=args.n_processes if args.n_processes > 1 else None,
        seed=args.seed,
        filepath=filepath,
        resume=resume,
        pass_dict=False,
    )
    if args.progress_newlines:
        enable_progress_newlines(sampler)
    success = sampler.run(
        f_live=args.f_live,
        n_shell=args.n_shell,
        n_eff=n_eff,
        n_like_max=n_like_max,
        timeout=timeout,
        verbose=args.verbose,
    )
    chain = sampler_to_chain(
        sampler,
        likelihood_mode=args.likelihood_mode,
        weighted_output=args.weighted_output,
        equal_weight_boost=args.equal_weight_boost,
        halo_model=args.halo_model,
        sidm_parameterization=args.sidm_parameterization,
        halo_run_context=args.halo_run_context,
    )
    chain_output.parent.mkdir(parents=True, exist_ok=True)
    chain.to_csv(chain_output, index=False)
    print(f"{run_label}_success={success}")
    print(f"wrote {chain_output}")
    print_nautilus_summary(sampler)
    print_summary(chain, use_sample_weights=args.use_sample_weights)
    return chain


def enable_progress_newlines(sampler) -> None:
    original_print_status = sampler.print_status

    def print_status_with_newline(status: str = "", header: bool = False, end: str = "\n") -> None:
        original_print_status(status=status, header=header, end="\n" if end == "\r" else end)

    sampler.print_status = print_status_with_newline


def make_nautilus_prior(
    prior_class: type,
    galaxy: GalaxyData,
    systemic_velocity_bounds: tuple[float, float],
    *,
    halo_model: str = "generalized-hernquist",
    sidm_parameterization: str = "scale",
):
    min_i = float(np.degrees(np.arccos(galaxy.observables.qprime)))
    prior = prior_class()
    if halo_model == "sidm":
        bounds_by_name = {
            "q_halo": (0.1, 2.0),
            "log10_rs0_pc": (0.0, 5.0),
            "log10_rho_s0_msun_pc3": (-5.0, 5.0),
            "log10_m200_msun": (5.0, 12.0),
            "log10_c200": (0.0, 2.0),
            "concentration_scatter_sigma": (-3.0, 3.0),
            "tau": (0.0, 1.08),
            "minus_log10_one_minus_beta_z": (-1.0, 1.0),
            "i_deg": (min_i, 90.0),
            "systemic_velocity_kms": systemic_velocity_bounds,
        }
        for name in SIDM_PARAMETER_NAMES[sidm_parameterization]:
            prior.add_parameter(name, dist=bounds_by_name[name])
        return prior

    bounds = [
        (0.1, 2.0),
        (0.0, 5.0),
        (-5.0, 5.0),
        (-1.0, 1.0),
        (0.5, 3.0),
        (3.0, 10.0),
        (0.0, 2.0),
        (min_i, 90.0),
        systemic_velocity_bounds,
    ]
    for name, bound in zip(PARAMETER_NAMES, bounds):
        prior.add_parameter(name, dist=bound)
    return prior


def make_log_probability_callable(galaxy: GalaxyData, args: argparse.Namespace):
    def log_probability(vector: np.ndarray) -> float:
        return evaluate_log_probability(galaxy, np.asarray(vector, dtype=float), args)

    return log_probability


def evaluate_log_probability(galaxy: GalaxyData, vector: np.ndarray, args: argparse.Namespace) -> float:
    return full_log_probability_vector(
        galaxy,
        vector,
        likelihood_mode=args.likelihood_mode,
        n_r=args.n_r,
        n_z=args.n_z,
        n_los=args.n_los,
        n_force_r=args.n_force_r,
        n_force_z=args.n_force_z,
        strict_epsrel=args.strict_epsrel,
        systemic_velocity_bounds=args.systemic_velocity_bounds,
        use_mge_physicality_check=args.use_mge_physicality_check,
        mge_n_gauss_halo=args.mge_n_gauss_halo,
        mge_n_gauss_tracer=args.mge_n_gauss_tracer,
        mge_decomposition_terms=args.mge_decomposition_terms,
        mge_n_u=args.mge_n_u,
        mge_r_min_pc=args.mge_r_min_pc,
        mge_r_max_pc=args.mge_r_max_pc,
        mge_real_fit_radii=getattr(args, "mge_real_fit_radii", 1600),
        mge_real_lstsq_rcond=getattr(args, "mge_real_lstsq_rcond", 1.0e-12),
        halo_model=args.halo_model,
        sidm_parameterization=args.sidm_parameterization,
        halo_run_context=args.halo_run_context,
    )


def sampler_to_chain(
    sampler,
    *,
    likelihood_mode: str,
    weighted_output: bool,
    equal_weight_boost: float,
    halo_model: str = "generalized-hernquist",
    sidm_parameterization: str = "scale",
    halo_run_context: HaloRunContext = HaloRunContext(),
) -> pd.DataFrame:
    if weighted_output:
        points, log_weights, log_prob = sampler.posterior(return_as_dict=False, equal_weight=False)
        rows = rows_from_points(
            points,
            log_prob,
            likelihood_mode=likelihood_mode,
            halo_model=halo_model,
            sidm_parameterization=sidm_parameterization,
            halo_run_context=halo_run_context,
        )
        chain = pd.DataFrame(rows)
        chain["log_weight"] = np.asarray(log_weights, dtype=float)
        chain["weight"] = np.exp(chain["log_weight"] - np.max(chain["log_weight"]))
        return chain

    posterior_kwargs = {"return_as_dict": False, "equal_weight": True}
    if "equal_weight_boost" in inspect.signature(sampler.posterior).parameters:
        posterior_kwargs["equal_weight_boost"] = equal_weight_boost
    points, _log_weights, log_prob = sampler.posterior(**posterior_kwargs)
    return pd.DataFrame(
        rows_from_points(
            points,
            log_prob,
            likelihood_mode=likelihood_mode,
            halo_model=halo_model,
            sidm_parameterization=sidm_parameterization,
            halo_run_context=halo_run_context,
        )
    )


def rows_from_points(
    points: np.ndarray | dict[str, Iterable[float]],
    log_prob: np.ndarray,
    *,
    likelihood_mode: str,
    halo_model: str = "generalized-hernquist",
    sidm_parameterization: str = "scale",
    halo_run_context: HaloRunContext = HaloRunContext(),
) -> list[dict[str, float | int | str]]:
    if isinstance(points, dict):
        points_array = np.column_stack(
            [points[name] for name in parameter_names_for(halo_model, sidm_parameterization)]
        )
    else:
        points_array = np.asarray(points, dtype=float)
    if points_array.ndim == 1:
        points_array = points_array.reshape(1, -1)

    rows = []
    for index, vector in enumerate(points_array):
        rows.append(
            vector_to_row(
                vector,
                log_probability=float(np.asarray(log_prob)[index]),
                step=index + 1,
                chain_id=0,
                sampler_name="nautilus",
                likelihood_mode=likelihood_mode,
                halo_model=halo_model,
                sidm_parameterization=sidm_parameterization,
                halo_run_context=halo_run_context,
            )
        )
    return rows


def chain_with_anisotropy_columns(chain: pd.DataFrame) -> pd.DataFrame:
    chain = chain.copy()
    if "minus_log10_one_minus_beta_z" not in chain.columns and "beta_z" in chain.columns:
        beta_z = pd.to_numeric(chain["beta_z"], errors="coerce")
        valid = 1.0 - beta_z > 0.0
        chain["minus_log10_one_minus_beta_z"] = np.nan
        chain.loc[valid, "minus_log10_one_minus_beta_z"] = -np.log10(1.0 - beta_z.loc[valid])
    if "beta_z" not in chain.columns and "minus_log10_one_minus_beta_z" in chain.columns:
        q_beta = pd.to_numeric(chain["minus_log10_one_minus_beta_z"], errors="coerce")
        chain["beta_z"] = 1.0 - 10.0 ** (-q_beta)
    return chain


def chain_after_burn(chain: pd.DataFrame, burn: int) -> pd.DataFrame:
    if burn < 0:
        raise ValueError("--burn must be non-negative")
    samples = chain.iloc[burn:].copy()
    if samples.empty:
        raise ValueError(f"no samples remain after burn={burn}; chain has {len(chain)} rows")
    return samples


def normalized_sample_weights(samples: pd.DataFrame, *, use_sample_weights: bool) -> np.ndarray | None:
    if not use_sample_weights:
        return None
    if "log_weight" in samples.columns:
        log_weight = pd.to_numeric(samples["log_weight"], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(log_weight)
        if not np.any(finite):
            raise ValueError("--use-sample-weights requested, but log_weight has no finite values")
        weights = np.zeros_like(log_weight, dtype=float)
        weights[finite] = np.exp(log_weight[finite] - np.max(log_weight[finite]))
    elif "weight" in samples.columns:
        weights = pd.to_numeric(samples["weight"], errors="coerce").to_numpy(dtype=float)
    else:
        raise ValueError("--use-sample-weights requested, but chain has no log_weight or weight column")
    weights = np.asarray(weights, dtype=float)
    weights[~np.isfinite(weights)] = 0.0
    if np.any(weights < 0.0):
        raise ValueError("sample weights must be non-negative")
    total = float(np.sum(weights))
    if not (np.isfinite(total) and total > 0.0):
        raise ValueError("sample weights must have a positive finite sum")
    return weights / total


def weighted_quantile_1d(values: np.ndarray, quantiles: Iterable[float], weights: np.ndarray | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    quantile_values = np.asarray(list(quantiles), dtype=float)
    finite = np.isfinite(values)
    if weights is None:
        if not np.any(finite):
            return np.full_like(quantile_values, np.nan, dtype=float)
        return np.nanpercentile(values[finite], 100.0 * quantile_values)
    weights = np.asarray(weights, dtype=float)
    finite &= np.isfinite(weights) & (weights > 0.0)
    if not np.any(finite):
        return np.full_like(quantile_values, np.nan, dtype=float)
    values = values[finite]
    weights = weights[finite]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    cumulative /= cumulative[-1]
    return np.interp(quantile_values, cumulative, values)


def weighted_percentile_axis0(values: np.ndarray, percentiles: Iterable[float], weights: np.ndarray | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if weights is None:
        return np.percentile(values, list(percentiles), axis=0)
    quantiles = [float(percentile) / 100.0 for percentile in percentiles]
    return np.vstack([weighted_quantile_1d(values[:, index], quantiles, weights) for index in range(values.shape[1])]).T


def plot_density_profile(
    galaxy: GalaxyData,
    chain: pd.DataFrame,
    output: Path,
    *,
    burn: int = 0,
    use_sample_weights: bool = False,
) -> None:
    samples = chain_after_burn(chain, burn)
    weights = normalized_sample_weights(samples, use_sample_weights=use_sample_weights)
    radius_kpc = np.geomspace(0.01, 20.0, 300)
    radius_pc = radius_kpc * 1000.0
    profiles = []
    from hayashi_jeans.halos import GeneralizedHernquistHalo, SIDMPSIDM25Halo

    for row in samples.itertuples(index=False):
        if getattr(row, "halo_model", "generalized-hernquist") == "sidm":
            halo = SIDMPSIDM25Halo(
                q=float(row.q_halo),
                rs0_pc=10.0 ** float(row.log10_rs0_pc),
                rho_s0_msun_pc3=10.0 ** float(row.log10_rho_s0_msun_pc3),
                tau=float(row.tau),
            )
        else:
            halo = GeneralizedHernquistHalo(
                q=float(row.q_halo),
                b_pc=10.0 ** float(row.log10_b_halo_pc),
                rho0_msun_pc3=10.0 ** float(row.log10_rho0_msun_pc3),
                alpha=float(row.alpha),
                beta=float(row.beta),
                gamma=float(row.gamma),
            )
        profiles.append([halo.density(float(r), 0.0) * 1.0e9 for r in radius_pc])
    profiles = np.asarray(profiles)
    p16, median, p84 = weighted_percentile_axis0(profiles, [16.0, 50.0, 84.0], weights=weights)

    fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=180)
    ax.fill_between(radius_kpc, p16, p84, color="#4c78a8", alpha=0.28, lw=0, label="68% interval")
    ax.plot(radius_kpc, median, color="#1f5f8b", lw=2.2, label="median")
    ax.axvline(galaxy.observables.b_star_pc / 1000.0, color="#333333", lw=1.4, ls="--", label=r"$b_*$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 20.0)
    ax.set_ylim(1.0e4, 1.0e10)
    ax.set_box_aspect(1.04)
    ax.set_xlabel("Major Axis [kpc]")
    ax.set_ylabel(r"$\rho_{\rm DM}(r)$ [$M_\odot\,{\rm kpc}^{-3}$]")
    ax.set_title(f"{galaxy.observables.galaxy} dark-matter density profile")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, which="both", alpha=0.22, lw=0.6)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)
    print(f"wrote {output}")


def generate_corner_map(
    chain: pd.DataFrame,
    output: Path,
    *,
    show_systemic_velocity: bool = False,
    burn: int = 0,
    use_sample_weights: bool = False,
) -> None:
    chain = chain_with_anisotropy_columns(chain)
    is_sidm = "halo_model" in chain and bool((chain["halo_model"] == "sidm").any())
    if is_sidm:
        parameterization = str(chain.loc[chain["halo_model"] == "sidm", "halo_parameterization"].iloc[0])
        columns = [
            name
            for name in SIDM_PARAMETER_NAMES[parameterization]
            if name != "systemic_velocity_kms"
        ]
        label_map = {
            "q_halo": r"$q_{\rm DM}$",
            "log10_rs0_pc": r"$\log_{10}r_{s,0}\,[{\rm pc}]$",
            "log10_rho_s0_msun_pc3": r"$\log_{10}\rho_{s,0}\,[M_\odot\,{\rm pc}^{-3}]$",
            "log10_m200_msun": r"$\log_{10}M_{200c}\,[M_\odot]$",
            "log10_c200": r"$\log_{10}c_{200c}$",
            "concentration_scatter_sigma": r"$s_c$",
            "tau": r"$\tau$",
            "minus_log10_one_minus_beta_z": r"$-\log_{10}(1-\beta_z)$",
            "i_deg": r"$i\,[{\rm deg}]$",
        }
        labels = [label_map[name] for name in columns]
    else:
        columns = [
            "q_halo",
            "log10_b_halo_pc",
            "log10_rho0_msun_pc3",
            "minus_log10_one_minus_beta_z",
            "alpha",
            "beta",
            "gamma",
            "i_deg",
        ]
        labels = [
            r"$q_{\rm DM}$",
            r"$\log_{10} b_{\rm halo}\,[{\rm pc}]$",
            r"$\log_{10}\rho_0\,[M_\odot\,{\rm pc}^{-3}]$",
            r"$-\log_{10}(1-\beta_z)$",
            r"$\alpha$",
            r"$\beta$",
            r"$\gamma$",
            r"$i\,[{\rm deg}]$",
        ]
    if show_systemic_velocity:
        columns.append("systemic_velocity_kms")
        labels.append(r"$\langle u\rangle\,[{\rm km\,s}^{-1}]$")

    missing = [column for column in columns if column not in chain.columns]
    if missing:
        raise ValueError(f"chain is missing columns required for corner map: {missing}")
    burned = chain_after_burn(chain, burn)
    samples = burned[columns].replace([np.inf, -np.inf], np.nan).dropna()
    if samples.empty:
        raise ValueError("no valid samples remain after burn and finite-value filtering")
    weights = normalized_sample_weights(burned.loc[samples.index], use_sample_weights=use_sample_weights)

    output.parent.mkdir(parents=True, exist_ok=True)
    if len(samples) <= len(columns):
        _plot_scatter_matrix_corner_fallback(
            samples,
            output,
            labels=labels,
            title="too few samples for corner; using scatter-matrix fallback",
            weights=weights,
        )
        print(f"wrote {output}")
        return

    try:
        import corner

        figure = corner.corner(
            samples.to_numpy(),
            labels=labels,
            weights=weights,
            quantiles=[0.16, 0.50, 0.84],
            show_titles=True,
            title_fmt=".3g",
            title_kwargs={"fontsize": 9},
            label_kwargs={"fontsize": 10},
            plot_datapoints=False,
            fill_contours=True,
            smooth=0.8,
            color="#1f5f8b",
            hist_kwargs={"density": True},
        )
        figure.savefig(output, dpi=180, bbox_inches="tight")
        plt.close(figure)
    except ImportError:
        _plot_scatter_matrix_corner_fallback(
            samples,
            output,
            labels=labels,
            title="corner package not installed; using scatter-matrix fallback",
            weights=weights,
        )
    print(f"wrote {output}")


def _plot_scatter_matrix_corner_fallback(
    samples: pd.DataFrame,
    output: Path,
    *,
    labels: list[str],
    title: str,
    weights: np.ndarray | None = None,
) -> None:
    if weights is not None and len(samples) > 1:
        rng = np.random.default_rng(20260526)
        n_resample = min(max(len(samples), 1000), 10000)
        indices = rng.choice(len(samples), size=n_resample, replace=True, p=weights / np.sum(weights))
        samples = samples.iloc[indices].reset_index(drop=True)
    axes = pd.plotting.scatter_matrix(
        samples,
        figsize=(1.55 * samples.shape[1], 1.55 * samples.shape[1]),
        diagonal="hist",
        alpha=0.35,
        color="#1f5f8b",
        hist_kwds={"density": True, "color": "#4c78a8"},
    )
    for ax in axes.ravel():
        ax.tick_params(labelsize=7)
        ax.xaxis.label.set_size(8)
        ax.yaxis.label.set_size(8)
    for column_index, label in enumerate(labels):
        axes[-1, column_index].set_xlabel(label)
        axes[column_index, 0].set_ylabel(label)
    fig = axes[0, 0].figure
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def print_summary(chain: pd.DataFrame, *, use_sample_weights: bool = False) -> None:
    chain = chain_with_anisotropy_columns(chain)
    weights = normalized_sample_weights(chain, use_sample_weights=use_sample_weights)
    print("parameter medians and 16/84 percentiles:")
    preferred_columns = [
        "q_halo",
        "log10_b_halo_pc",
        "log10_rho0_msun_pc3",
        "log10_rs0_pc",
        "log10_rho_s0_msun_pc3",
        "log10_m200_msun",
        "log10_c200",
        "c200",
        "concentration_scatter_sigma",
        "tau",
        "minus_log10_one_minus_beta_z",
        "beta_z",
        "alpha",
        "beta",
        "gamma",
        "i_deg",
        "systemic_velocity_kms",
        "log_probability",
    ]
    for col in [name for name in preferred_columns if name in chain.columns]:
        if not np.any(np.isfinite(pd.to_numeric(chain[col], errors="coerce").to_numpy(dtype=float))):
            continue
        q16, q50, q84 = weighted_quantile_1d(chain[col].to_numpy(dtype=float), [0.16, 0.50, 0.84], weights=weights)
        print(f"{col}: {q50:.6g} -{q50-q16:.6g} +{q84-q50:.6g}")


def print_nautilus_summary(sampler) -> None:
    for name in ("log_z", "n_eff", "f_live", "n_like", "n_like_iter", "n_batch", "n_live"):
        value = getattr(sampler, name, None)
        if value is not None:
            if isinstance(value, (int, np.integer)):
                print(f"{name}={int(value)}")
            else:
                print(f"{name}={float(value):.6g}")


if __name__ == "__main__":
    main()
