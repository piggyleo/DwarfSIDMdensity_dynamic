#!/usr/bin/env python
"""Build source-aware trial membership and unique-star catalogs."""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GALAXY_DIR = PROJECT_ROOT / "data" / "galaxies"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
OUTPUT_DIR = PROCESSED_DIR / "membership"
JI_PATH = PROCESSED_DIR / "kinematics_ref9_ji2021_antlia2_crater2.csv"
SEGUE2_MATCH_AUDIT_PATH = (
    PROJECT_ROOT / "data" / "validation" / "segue2_cross_source_match_audit_2026-06-14.csv"
)
PROBABILITY_THRESHOLD = 0.95
VELOCITY_VARIABILITY_P_THRESHOLD = 0.01
SEGUE2_MATCH_THRESHOLD_ARCSEC = 0.2

SOURCE_VELOCITY_EXCLUSIONS: dict[tuple[str, str], tuple[str, str]] = {
    ("Bootes I", "Boo1_26"): ("known_velocity_variable", "Jenkins et al. 2021 binary candidate."),
    ("Bootes I", "Boo1_61"): ("known_velocity_variable", "Jenkins et al. 2021 binary candidate."),
    ("Bootes I", "Boo1_111"): ("known_velocity_variable", "Jenkins et al. 2021 binary candidate."),
    ("Bootes I", "Boo1_114"): ("known_velocity_variable", "Jenkins et al. 2021 binary candidate."),
    ("Bootes I", "Boo1_32"): ("known_pulsating_variable", "Jenkins et al. 2021 RR Lyrae star."),
    (
        "Grus 1",
        "DESJ225637.05-501024.8",
    ): ("known_velocity_variable", "Chiti et al. 2022 binary candidate."),
    (
        "Grus 1",
        "DESJ225619.67-500913.1",
    ): ("known_velocity_variable", "Chiti et al. 2022 binary candidate."),
    ("Leo IV", "Leo4_1039"): ("known_velocity_variable", "Jenkins et al. 2021 binary candidate."),
    ("Leo IV", "Leo4_1041"): ("known_pulsating_variable", "Jenkins et al. 2021 RR Lyrae star."),
    ("Leo V", "Leo5_1034"): ("known_velocity_variable", "Jenkins et al. 2021 binary candidate."),
    ("Leo V", "Leo5_1038"): ("known_velocity_variable", "Jenkins et al. 2021 binary candidate."),
    ("Leo V", "Leo5_1051"): ("known_pulsating_variable", "Jenkins et al. 2021 RR Lyrae star."),
    (
        "Segue 2",
        "J021900.06+200635.2",
    ): ("known_pulsating_variable", "Kirby et al. 2013 RR Lyrae star."),
    ("Tucana 2", "TucII-078"): ("known_velocity_variable", "Chiti et al. 2023 binary candidate."),
    ("Tucana 2", "TucII-309"): ("known_velocity_variable", "Chiti et al. 2023 binary candidate."),
}

TUCANA2_PREFERRED_INSTRUMENT = {
    "TucII-006": "MIKE",
    "TucII-011": "MIKE",
    "TucII-022": "IMACS",
    "TucII-033": "MIKE",
    "TucII-052": "MIKE",
    "TucII-074": "M2FS",
    "TucII-085": "M2FS",
    "Star12": "IMACS",
    "Star68": "IMACS",
    "TucII-203": "MIKE",
    "TucII-206": "MIKE",
    "TucII-301": "MIKE",
    "TucII-303": "MIKE",
    "TucII-305": "MIKE",
    "TucII-306": "MIKE",
    "TucII-310": "MagE",
    "TucII-320": "MagE",
}


@dataclass(frozen=True)
class GalaxyConfig:
    filename: str
    classification: str


CONFIGS: dict[str, GalaxyConfig] = {
    "Antlia 2": GalaxyConfig("01_Antlia_2.csv", "ji_probability"),
    "Bootes I": GalaxyConfig("02_Bootes_I.csv", "hard_flag"),
    "Crater 2": GalaxyConfig("06_Crater_2.csv", "ji_probability"),
    "Draco 2": GalaxyConfig("07_Draco_2.csv", "hard_flag"),
    "Grus 1": GalaxyConfig("09_Grus_1.csv", "hard_flag_multiepoch"),
    "Grus 2": GalaxyConfig("10_Grus_2.csv", "hard_flag_multiepoch"),
    "Hercules": GalaxyConfig("11_Hercules.csv", "hard_flag"),
    "Leo IV": GalaxyConfig("14_Leo_IV.csv", "hard_flag"),
    "Leo V": GalaxyConfig("15_Leo_V.csv", "hard_flag"),
    "Reticulum II": GalaxyConfig("18_Reticulum_II.csv", "hard_flag"),
    "Segue 2": GalaxyConfig("20_Segue_2.csv", "hard_flag_multiepoch"),
    "Tucana 2": GalaxyConfig("22_Tucana_2.csv", "hard_flag_multiepoch"),
    "Tucana 3": GalaxyConfig("23_Tucana_3.csv", "hard_flag_multiepoch"),
    "Tucana 4": GalaxyConfig("24_Tucana_4.csv", "hard_flag_multiepoch"),
}

OBSERVATION_COLUMNS = [
    "galaxy",
    "canonical_star_id",
    "star_id",
    "source_ref_id",
    "source_catalog",
    "velocity_instrument",
    "ra_deg",
    "dec_deg",
    "v_los_kms",
    "v_los_err_down_kms",
    "v_los_err_up_kms",
    "source_member_label",
    "membership_status",
    "membership_basis",
    "membership_probability",
    "membership_probability_threshold",
    "binary_status",
    "dynamical_velocity_status",
    "classification_notes",
]

UNIQUE_STAR_COLUMNS = [
    "galaxy",
    "canonical_star_id",
    "representative_star_id",
    "ra_deg",
    "dec_deg",
    "n_observations",
    "n_sources",
    "membership_status",
    "membership_basis",
    "membership_probability",
    "binary_status",
    "velocity_variability_p_value",
    "dynamical_velocity_status",
    "v_los_combined_kms",
    "v_los_combined_err_down_kms",
    "v_los_combined_err_up_kms",
    "classification_notes",
]


def _slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def _canonical_star_id(row: pd.Series) -> str:
    ra = pd.to_numeric(row["ra_deg"], errors="coerce")
    dec = pd.to_numeric(row["dec_deg"], errors="coerce")
    if np.isfinite(ra) and np.isfinite(dec):
        return f"coord:{ra:.7f}:{dec:.7f}"
    star_id = str(row.get("star_id", "")).strip()
    if star_id and star_id.lower() != "nan":
        return f"id:{star_id}"
    raise ValueError(f"Cannot construct a canonical identity for row {row.name}")


def _valid_velocity(row: pd.Series) -> bool:
    values = [
        pd.to_numeric(row["v_los_kms"], errors="coerce"),
        pd.to_numeric(row["v_los_err_down_kms"], errors="coerce"),
        pd.to_numeric(row["v_los_err_up_kms"], errors="coerce"),
    ]
    return all(np.isfinite(value) and (i == 0 or value > 0) for i, value in enumerate(values))


def _velocity_instrument(row: pd.Series) -> str:
    notes = str(row.get("notes", ""))
    if "v_los backfilled from ref28" in notes:
        return "backfilled_ref28"
    match = re.search(r"(?:^|[; ])Inst=([^; ]+)", notes)
    if match:
        return match.group(1)
    catalog = str(row.get("source_catalog", ""))
    if catalog == "J/ApJ/892/137/table2":
        return "AAT"
    if catalog == "J/ApJ/892/137/table3":
        return "IMACS"
    return ""


def _ji_binary_lookup() -> dict[tuple[str, float, float], str]:
    table = pd.read_csv(JI_PATH)
    return {
        (str(row.galaxy), round(float(row.ra_deg), 5), round(float(row.dec_deg), 5)): str(row.binary_flag)
        for row in table.itertuples()
    }


def _classify_observations(galaxy: str, config: GalaxyConfig) -> pd.DataFrame:
    table = pd.read_csv(GALAXY_DIR / config.filename, comment="#", dtype={"star_id": str})
    stars = table.loc[table["row_kind"] == "star"].copy()
    for column in [
        "ra_deg",
        "dec_deg",
        "v_los_kms",
        "v_los_err_down_kms",
        "v_los_err_up_kms",
        "memprob",
        "member_flag",
    ]:
        stars[column] = pd.to_numeric(stars[column], errors="coerce")

    ji_binary = _ji_binary_lookup() if config.classification == "ji_probability" else {}
    records: list[dict[str, object]] = []
    for _, row in stars.iterrows():
        probability = row["memprob"] if np.isfinite(row["memprob"]) else np.nan
        hard_flag = int(row["member_flag"]) if np.isfinite(row["member_flag"]) else None
        source_label = ""
        status = "unclassified"
        basis = ""
        binary_status = "unknown"
        notes = ""

        if config.classification == "ji_probability":
            is_ji = str(row["source_catalog"]) in {
                "J/ApJ/921/32/table4",
                "J/ApJ/921/32/table5",
            }
            if is_ji and np.isfinite(probability):
                source_label = "membership_probability_only"
                basis = "source_probability"
                if probability > PROBABILITY_THRESHOLD:
                    status = "confirmed"
                    notes = "Selected by the uniform source-probability rule P_mem > 0.95."
                else:
                    status = "not_selected_probability"
                    notes = "Below the baseline probability threshold; not forced to hard nonmember."
                binary_flag = ji_binary.get(
                    (galaxy, round(float(row["ra_deg"]), 5), round(float(row["dec_deg"]), 5))
                )
                binary_status = {"Y": "known_velocity_variable", "N": "not_flagged_variable"}.get(
                    binary_flag, "unknown"
                )
            else:
                source_label = "no_hard_label"
                basis = "unclassified_source_table"
                notes = "Source row has neither a hard membership label nor the adopted source probability."
        else:
            source_label = {1: "member", 0: "nonmember", 2: "uncertain"}.get(hard_flag, "missing")
            basis = "source_hard_classification"
            status = {1: "confirmed", 0: "nonmember", 2: "unclassified"}.get(hard_flag, "unclassified")
            notes = "Membership follows the source hard classification."

        source_exclusion = SOURCE_VELOCITY_EXCLUSIONS.get((galaxy, str(row["star_id"])))
        if source_exclusion is not None:
            binary_status, exclusion_note = source_exclusion
            notes = f"{notes} {exclusion_note}"

        if not _valid_velocity(row):
            velocity_status = "missing_or_invalid_velocity"
        elif binary_status in {"known_velocity_variable", "known_pulsating_variable"}:
            velocity_status = f"exclude_{binary_status}"
        else:
            velocity_status = "usable_not_known_variable"

        records.append(
            {
                "galaxy": galaxy,
                "canonical_star_id": _canonical_star_id(row),
                "star_id": row["star_id"],
                "source_ref_id": row["source_ref_id"],
                "source_catalog": row["source_catalog"],
                "velocity_instrument": _velocity_instrument(row),
                "ra_deg": row["ra_deg"],
                "dec_deg": row["dec_deg"],
                "v_los_kms": row["v_los_kms"],
                "v_los_err_down_kms": row["v_los_err_down_kms"],
                "v_los_err_up_kms": row["v_los_err_up_kms"],
                "source_member_label": source_label,
                "membership_status": status,
                "membership_basis": basis,
                "membership_probability": probability,
                "membership_probability_threshold": (
                    PROBABILITY_THRESHOLD if basis == "source_probability" else np.nan
                ),
                "binary_status": binary_status,
                "dynamical_velocity_status": velocity_status,
                "classification_notes": notes,
            }
        )
    return pd.DataFrame(records, columns=OBSERVATION_COLUMNS)


def _angular_separation_arcsec(left: pd.DataFrame, right: pd.DataFrame) -> np.ndarray:
    left_ra = np.deg2rad(left["ra_deg"].to_numpy(float))
    left_dec = np.deg2rad(left["dec_deg"].to_numpy(float))
    right_ra = np.deg2rad(right["ra_deg"].to_numpy(float))
    right_dec = np.deg2rad(right["dec_deg"].to_numpy(float))
    delta_ra = (left_ra[:, None] - right_ra[None, :]) * np.cos(
        0.5 * (left_dec[:, None] + right_dec[None, :])
    )
    delta_dec = left_dec[:, None] - right_dec[None, :]
    return np.hypot(delta_ra, delta_dec) * 206264.80624709636


def _apply_segue2_cross_source_identity(
    observations: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    catalogs = sorted(observations["source_catalog"].dropna().unique())
    if catalogs != ["J/ApJ/770/16/table2", "J/MNRAS/397/1748/table2"]:
        raise ValueError(f"Unexpected Segue 2 catalog set: {catalogs}")

    groups = {
        catalog: (
            observations.loc[observations["source_catalog"] == catalog]
            .drop_duplicates("canonical_star_id")
            .reset_index()
        )
        for catalog in catalogs
    }
    left, right = groups[catalogs[0]], groups[catalogs[1]]
    separations = _angular_separation_arcsec(left, right)
    left_nearest = separations.argmin(axis=1)
    right_nearest = separations.argmin(axis=0)

    audit_rows: list[dict[str, object]] = []
    canonical_replacements: dict[str, str] = {}
    for left_position, right_position in enumerate(left_nearest):
        separation = float(separations[left_position, right_position])
        if right_nearest[right_position] != left_position or separation > SEGUE2_MATCH_THRESHOLD_ARCSEC:
            continue
        left_row = left.iloc[left_position]
        right_row = right.iloc[right_position]
        preferred_id = min(left_row["canonical_star_id"], right_row["canonical_star_id"])
        canonical_replacements[left_row["canonical_star_id"]] = preferred_id
        canonical_replacements[right_row["canonical_star_id"]] = preferred_id

        kirby = left_row
        belokurov = right_row
        kirby_has_velocity = (
            np.isfinite(pd.to_numeric(kirby["v_los_kms"], errors="coerce"))
            and kirby["velocity_instrument"] != "backfilled_ref28"
        )
        if kirby["membership_status"] == "confirmed":
            membership_resolution = "kirby_confirmed"
        elif belokurov["membership_status"] == "confirmed" and not kirby_has_velocity:
            membership_resolution = "belokurov_member_kirby_unobserved"
        else:
            membership_resolution = "kirby_classification_priority"
        velocity_resolution = "prefer_kirby_velocity" if kirby_has_velocity else "use_belokurov_velocity"
        audit_rows.append(
            {
                "kirby_star_id": kirby["star_id"],
                "belokurov_star_id": belokurov["star_id"],
                "separation_arcsec": separation,
                "kirby_membership_status": kirby["membership_status"],
                "belokurov_membership_status": belokurov["membership_status"],
                "kirby_v_los_kms": kirby["v_los_kms"],
                "belokurov_v_los_kms": belokurov["v_los_kms"],
                "membership_resolution": membership_resolution,
                "velocity_resolution": velocity_resolution,
                "match_basis": "reciprocal_nearest_neighbour_within_0.2_arcsec",
            }
        )

    result = observations.copy()
    result["canonical_star_id"] = result["canonical_star_id"].replace(canonical_replacements)
    matched = result["canonical_star_id"].isin(set(canonical_replacements.values()))
    result.loc[matched, "classification_notes"] = (
        result.loc[matched, "classification_notes"]
        + " Cross-source identity assigned by reciprocal nearest-neighbour matching."
    )
    audit = pd.DataFrame(audit_rows)
    if len(audit) != 39:
        raise ValueError(f"Expected 39 Segue 2 cross-source matches, found {len(audit)}")
    return result, audit


def _combine_membership_status(group: pd.DataFrame) -> str:
    if group["galaxy"].iloc[0] == "Segue 2" and group["source_catalog"].nunique() > 1:
        kirby = group.loc[group["source_catalog"] == "J/ApJ/770/16/table2"]
        belokurov = group.loc[group["source_catalog"] == "J/MNRAS/397/1748/table2"]
        if (kirby["membership_status"] == "confirmed").any():
            return "confirmed"
        kirby_has_velocity = (
            kirby.apply(_valid_velocity, axis=1)
            & kirby["velocity_instrument"].ne("backfilled_ref28")
        ).any()
        if (belokurov["membership_status"] == "confirmed").any() and not kirby_has_velocity:
            return "confirmed"
        if (kirby["membership_status"] == "nonmember").any():
            return "nonmember"

    unique = set(group["membership_status"])
    if "confirmed" in unique and "nonmember" in unique:
        return "conflicting_source_classification"
    for status in ["confirmed", "unclassified", "not_selected_probability", "nonmember"]:
        if status in unique:
            return status
    return "unclassified"


def _regularized_gamma_q(shape: float, value: float) -> float:
    """Return Q(shape, value) using standard series/continued fractions."""
    if shape <= 0 or value < 0:
        raise ValueError("Gamma arguments must satisfy shape > 0 and value >= 0")
    if value == 0:
        return 1.0
    eps = 3e-14
    max_iter = 500
    log_prefactor = -value + shape * math.log(value) - math.lgamma(shape)

    if value < shape + 1:
        term = 1.0 / shape
        series = term
        denominator = shape
        for _ in range(max_iter):
            denominator += 1.0
            term *= value / denominator
            series += term
            if abs(term) < abs(series) * eps:
                return max(0.0, min(1.0, 1.0 - series * math.exp(log_prefactor)))
        raise RuntimeError("Incomplete-gamma series did not converge")

    tiny = 1e-300
    b = value + 1.0 - shape
    c = 1.0 / tiny
    d = 1.0 / b
    fraction = d
    for iteration in range(1, max_iter + 1):
        coefficient = -iteration * (iteration - shape)
        b += 2.0
        d = coefficient * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + coefficient / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        fraction *= delta
        if abs(delta - 1.0) < eps:
            return max(0.0, min(1.0, math.exp(log_prefactor) * fraction))
    raise RuntimeError("Incomplete-gamma continued fraction did not converge")


def _chi_square_survival(value: float, degrees_of_freedom: int) -> float:
    return _regularized_gamma_q(0.5 * degrees_of_freedom, 0.5 * value)


def _combine_velocity(group: pd.DataFrame) -> tuple[float, float, float, float, str, str]:
    valid = group.loc[group["dynamical_velocity_status"] != "missing_or_invalid_velocity"].copy()
    if valid.empty:
        return np.nan, np.nan, np.nan, np.nan, "unknown", "missing_or_invalid_velocity"
    source_exclusions = valid.loc[
        valid["binary_status"].isin(["known_velocity_variable", "known_pulsating_variable"]),
        "binary_status",
    ]
    if not source_exclusions.empty:
        binary_status = source_exclusions.iloc[0]
        return np.nan, np.nan, np.nan, np.nan, binary_status, f"exclude_{binary_status}"

    galaxy = group["galaxy"].iloc[0]
    representative_ids = set(group["star_id"].astype(str))
    if galaxy == "Tucana 2":
        preferred = {
            TUCANA2_PREFERRED_INSTRUMENT[star_id]
            for star_id in representative_ids
            if star_id in TUCANA2_PREFERRED_INSTRUMENT
        }
        if len(preferred) == 1:
            preferred_instrument = next(iter(preferred))
            valid = valid.loc[valid["velocity_instrument"] == preferred_instrument].copy()
            if valid.empty:
                return (
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    "source_selection_missing",
                    "missing_preferred_source_velocity",
                )

    if galaxy == "Segue 2" and valid["source_catalog"].nunique() > 1:
        kirby = valid.loc[
            (valid["source_catalog"] == "J/ApJ/770/16/table2")
            & valid["velocity_instrument"].ne("backfilled_ref28")
        ]
        if not kirby.empty:
            valid = kirby
        else:
            valid = valid.loc[valid["source_catalog"] == "J/MNRAS/397/1748/table2"]

    if galaxy == "Antlia 2" and valid["source_catalog"].nunique() > 1:
        # Ji et al. re-reduced the Torrealba et al. spectra consistently with
        # the new S5 observations, so its table is the final velocity product.
        valid = valid.loc[valid["source_catalog"] == "J/ApJ/921/32/table4"].copy()

    if galaxy in {"Grus 2", "Tucana 4"} and valid["source_catalog"].nunique() > 1:
        # Simon et al. define membership and derive dynamics from the IMACS
        # table; the AAT table is retained as ancillary provenance.
        valid = valid.loc[valid["source_catalog"] == "J/ApJ/892/137/table3"].copy()

    catalogs = {value for value in valid["source_catalog"].astype(str) if value}
    instruments = {value for value in valid["velocity_instrument"].astype(str) if value}
    if len(catalogs) > 1 or len(instruments) > 1:
        return (
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            "not_evaluated_cross_system",
            "requires_velocity_zero_point_resolution",
        )

    velocity = valid["v_los_kms"].to_numpy(float)
    err_down = valid["v_los_err_down_kms"].to_numpy(float)
    err_up = valid["v_los_err_up_kms"].to_numpy(float)
    symmetric_error = 0.5 * (err_down + err_up)

    if len(valid) == 1:
        return (
            float(velocity[0]),
            float(err_down[0]),
            float(err_up[0]),
            np.nan,
            "untested_single_epoch",
            "usable_not_known_variable",
        )

    weights = 1.0 / np.square(symmetric_error)
    mean_velocity = float(np.sum(weights * velocity) / np.sum(weights))
    mean_error = float(np.sqrt(1.0 / np.sum(weights)))
    chi_square = float(np.sum(np.square((velocity - mean_velocity) / symmetric_error)))
    p_value = _chi_square_survival(chi_square, len(valid) - 1)
    if p_value < VELOCITY_VARIABILITY_P_THRESHOLD:
        return np.nan, np.nan, np.nan, p_value, "velocity_variable", "exclude_velocity_variable"
    return (
        mean_velocity,
        mean_error,
        mean_error,
        p_value,
        "not_flagged_variable",
        "usable_merged_multiepoch",
    )


def _build_unique_stars(observations: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for canonical_id, group in observations.groupby("canonical_star_id", sort=True):
        velocity, err_down, err_up, p_value, binary_status, velocity_status = _combine_velocity(group)
        membership_status = _combine_membership_status(group)
        probabilities = pd.to_numeric(group["membership_probability"], errors="coerce").dropna()
        records.append(
            {
                "galaxy": group["galaxy"].iloc[0],
                "canonical_star_id": canonical_id,
                "representative_star_id": next(
                    (value for value in group["star_id"].astype(str) if value and value.lower() != "nan"),
                    "",
                ),
                "ra_deg": pd.to_numeric(group["ra_deg"], errors="coerce").median(),
                "dec_deg": pd.to_numeric(group["dec_deg"], errors="coerce").median(),
                "n_observations": len(group),
                "n_sources": group["source_catalog"].nunique(),
                "membership_status": membership_status,
                "membership_basis": "|".join(sorted(set(group["membership_basis"].astype(str)))),
                "membership_probability": probabilities.max() if not probabilities.empty else np.nan,
                "binary_status": binary_status,
                "velocity_variability_p_value": p_value,
                "dynamical_velocity_status": velocity_status,
                "v_los_combined_kms": velocity,
                "v_los_combined_err_down_kms": err_down,
                "v_los_combined_err_up_kms": err_up,
                "classification_notes": (
                    "Trial identity uses coordinates rounded to 7 decimal degrees; "
                    "cross-source matches require later manual validation."
                ),
            }
        )
    return pd.DataFrame(records, columns=UNIQUE_STAR_COLUMNS)


def _write_with_notes(frame: pd.DataFrame, path: Path, notes: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    with path.open("a", encoding="utf-8") as handle:
        handle.write("\n# COLUMN_NOTES_BEGIN\n")
        for note in notes:
            handle.write(f"# {note}\n")
        handle.write("# COLUMN_NOTES_END\n")


def _write_summary(summaries: list[dict[str, object]]) -> None:
    frame = pd.DataFrame(summaries)
    _write_with_notes(
        frame,
        OUTPUT_DIR / "trial_membership_summary_2026-06-14.csv",
        [
            "n_confirmed_unique: unique stars classified as confirmed by the stated source-aware rule.",
            "n_known_or_detected_velocity_variable: variables among all stars, including nonmembers.",
            "n_confirmed_velocity_excluded: confirmed members excluded for binary, pulsation, or detected velocity variability.",
            "n_baseline_dynamical: confirmed stars with a usable representative velocity and no known velocity variability.",
            "This is a trial classification and does not overwrite the original observation catalogs.",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--galaxy", action="append", choices=sorted(CONFIGS))
    args = parser.parse_args()
    selected = args.galaxy or list(CONFIGS)

    summaries: list[dict[str, object]] = []
    for galaxy in selected:
        observations = _classify_observations(galaxy, CONFIGS[galaxy])
        if galaxy == "Segue 2":
            observations, match_audit = _apply_segue2_cross_source_identity(observations)
            _write_with_notes(
                match_audit,
                SEGUE2_MATCH_AUDIT_PATH,
                [
                    "Matches are reciprocal nearest neighbours across the two source catalogues within 0.2 arcsec.",
                    "Kirby et al. 2013 is preferred when it contains an effective observation; otherwise a Belokurov et al. 2009 member is retained.",
                    "The audit records identity resolution only and does not overwrite either raw source row.",
                ],
            )
        unique_stars = _build_unique_stars(observations)
        slug = _slug(galaxy)
        _write_with_notes(
            observations,
            OUTPUT_DIR / f"{slug}_trial_classified_observations.csv",
            [
                "One row is one source observation.",
                "membership_status separates source membership from dynamical velocity usability.",
                "Probability-only rows use the strict rule P_mem > 0.95; lower probabilities are not forced to hard nonmember.",
            ],
        )
        _write_with_notes(
            unique_stars,
            OUTPUT_DIR / f"{slug}_trial_unique_stars.csv",
            [
                "One row is one trial unique star.",
                "Multi-epoch velocities use inverse-variance means only when a constant-velocity chi-square test has p >= 0.01.",
                "Single-epoch stars are untested, not proven non-variable.",
            ],
        )
        summaries.append(
            {
                "galaxy": galaxy,
                "n_observation_rows": len(observations),
                "n_unique_stars": len(unique_stars),
                "n_confirmed_unique": int((unique_stars["membership_status"] == "confirmed").sum()),
                "n_known_or_detected_velocity_variable": int(
                    unique_stars["binary_status"].isin(
                        ["known_velocity_variable", "known_pulsating_variable", "velocity_variable"]
                    ).sum()
                ),
                "n_confirmed_velocity_excluded": int(
                    (
                        (unique_stars["membership_status"] == "confirmed")
                        & unique_stars["binary_status"].isin(
                            [
                                "known_velocity_variable",
                                "known_pulsating_variable",
                                "velocity_variable",
                            ]
                        )
                    ).sum()
                ),
                "n_baseline_dynamical": int(
                    (
                        (unique_stars["membership_status"] == "confirmed")
                        & unique_stars["dynamical_velocity_status"].isin(
                            ["usable_not_known_variable", "usable_merged_multiepoch"]
                        )
                    ).sum()
                ),
                "classification_mode": CONFIGS[galaxy].classification,
            }
        )
    _write_summary(summaries)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
