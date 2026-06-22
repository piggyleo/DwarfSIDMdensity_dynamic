#!/usr/bin/env python
"""Write model-ready unique-star rows for the 14 audited galaxies.

The original observation-level galaxy CSVs are expected to be backed up before
running this script. It uses the audited trial unique-star and classified
observation tables to keep one model-ready stellar row per unique star.
"""

from __future__ import annotations

import csv
import io
import math
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GALAXY_DIR = PROJECT_ROOT / "data" / "galaxies"
MEMBERSHIP_DIR = PROJECT_ROOT / "data" / "processed" / "membership"
AUDIT_PATH = PROJECT_ROOT / "data" / "validation" / "model_member_flag_application_2026-06-16.csv"

TRIAL_GALAXIES: dict[str, str] = {
    "Antlia 2": "01_Antlia_2.csv",
    "Bootes I": "02_Bootes_I.csv",
    "Crater 2": "06_Crater_2.csv",
    "Draco 2": "07_Draco_2.csv",
    "Grus 1": "09_Grus_1.csv",
    "Grus 2": "10_Grus_2.csv",
    "Hercules": "11_Hercules.csv",
    "Leo IV": "14_Leo_IV.csv",
    "Leo V": "15_Leo_V.csv",
    "Reticulum II": "18_Reticulum_II.csv",
    "Segue 2": "20_Segue_2.csv",
    "Tucana 2": "22_Tucana_2.csv",
    "Tucana 3": "23_Tucana_3.csv",
    "Tucana 4": "24_Tucana_4.csv",
}

SLUGS: dict[str, str] = {
    galaxy: filename.removesuffix(".csv").split("_", 1)[1].lower()
    for galaxy, filename in TRIAL_GALAXIES.items()
}

SENSITIVITY_STARS: dict[tuple[str, str], str] = {
    ("Grus 2", "J220420.12-462341.3"): "weak_velocity_variable_candidate_p_0.020",
    ("Grus 2", "J220437.47-461955.0"): "weak_velocity_variable_candidate_p_0.030",
    ("Grus 2", "J220403.28-462520.0"): "weak_velocity_variable_candidate_p_0.083",
    ("Hercules", "308_734"): "borderline_source_member_velocity_nai_metallicity_edge",
}

GALAXY_WIDE_SENSITIVITY: dict[str, str] = {
    "Draco 2": "single_epoch_no_source_binary_screen",
    "Tucana 2": "source_text_16_vs_17_dynamical_count_inconsistency",
}


def _read_project_csv(path: Path) -> tuple[list[str], list[dict[str, str]], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    comment_start = next((i for i, line in enumerate(lines) if line.startswith("#")), len(lines))
    reader = csv.DictReader(io.StringIO("\n".join(lines[:comment_start])))
    if reader.fieldnames is None:
        raise ValueError(f"{path} has no CSV header")
    return reader.fieldnames, list(reader), lines[comment_start:]


def _write_project_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]], comments: list[str]) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    text = buffer.getvalue()
    if comments:
        text += "\n".join(_comments_with_model_notes(comments)) + "\n"
    path.write_text(text, encoding="utf-8")


def _comments_with_model_notes(comments: list[str]) -> list[str]:
    notes = [
        "# model_member_flag: 建模样本标记；1 为 baseline clean member，2 为默认纳入但需敏感性测试的成员，3 为不进入默认 Jeans 建模的非成员或速度不可用成员。",
        "# model_selection_reason: model_member_flag 的人工可读原因，记录边缘成员、弱速度变量、双星排除、非成员等口径。",
        "# model_velocity_rule: 写入 v_los_kms 和误差列时采用的速度选择、合并或排除规则。",
    ]
    filtered = [line for line in comments if not any(key in line for key in ("model_member_flag:", "model_selection_reason:", "model_velocity_rule:"))]
    try:
        end = filtered.index("# COLUMN_NOTES_END")
    except ValueError:
        return filtered + ["# COLUMN_NOTES_BEGIN", *notes, "# COLUMN_NOTES_END"]
    return filtered[:end] + notes + filtered[end:]


def _text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return str(value)


def _number(value: object) -> str:
    numeric = pd.to_numeric(value, errors="coerce")
    if pd.isna(numeric):
        return ""
    return f"{float(numeric):.12g}"


def _source_values(values: pd.Series) -> str:
    cleaned: list[str] = []
    for value in values.dropna().astype(str):
        if value.endswith(".0"):
            value = value[:-2]
        if value and value not in cleaned:
            cleaned.append(value)
    return ";".join(cleaned)


def _member_flag_from_status(status: str) -> str:
    if status == "confirmed":
        return "1"
    if status in {"nonmember", "not_selected_probability"}:
        return "0"
    return "2"


def _model_flag_and_reason(galaxy: str, unique_row: pd.Series) -> tuple[str, str]:
    star_id = str(unique_row["representative_star_id"])
    membership = str(unique_row["membership_status"])
    velocity_status = str(unique_row["dynamical_velocity_status"])
    binary_status = str(unique_row["binary_status"])

    if membership != "confirmed":
        return "3", f"excluded_{membership}"
    if not velocity_status.startswith("usable"):
        return "3", f"confirmed_member_velocity_excluded_{velocity_status}"

    reasons: list[str] = []
    if galaxy in GALAXY_WIDE_SENSITIVITY:
        reasons.append(GALAXY_WIDE_SENSITIVITY[galaxy])
    star_reason = SENSITIVITY_STARS.get((galaxy, star_id))
    if star_reason is not None:
        reasons.append(star_reason)
    p_value = pd.to_numeric(unique_row.get("velocity_variability_p_value"), errors="coerce")
    if (
        galaxy in {"Grus 2", "Tucana 4"}
        and pd.notna(p_value)
        and 0.01 <= float(p_value) <= 0.10
        and "weak_velocity_variable" not in " ".join(reasons)
    ):
        reasons.append(f"weak_velocity_variable_candidate_p_{float(p_value):.3f}")

    if reasons:
        return "2", ";".join(reasons)
    if binary_status == "untested_single_epoch":
        return "1", "baseline_member_single_epoch_velocity"
    return "1", "baseline_member"


def _row_notes(galaxy: str, unique_row: pd.Series, observations: pd.DataFrame) -> str:
    parts = [
        "model-ready unique-star row generated from audited trial membership tables",
        f"membership_status={unique_row['membership_status']}",
        f"binary_status={unique_row['binary_status']}",
        f"dynamical_velocity_status={unique_row['dynamical_velocity_status']}",
        f"n_observations={int(unique_row['n_observations'])}",
        f"source_catalogs={_source_values(observations['source_catalog'])}",
    ]
    if galaxy == "Hercules" and unique_row["representative_star_id"] == "309_528":
        parts.append("Kirby metallicity retained in original audit because upstream Simon & Geha [Fe/H] is conflicting")
    return "; ".join(parts)


def _build_rows(galaxy: str, filename: str) -> tuple[list[str], list[dict[str, object]], list[dict[str, object]]]:
    path = GALAXY_DIR / filename
    fieldnames, existing_rows, comments = _read_project_csv(path)
    for col in ("model_member_flag", "model_selection_reason", "model_velocity_rule"):
        if col not in fieldnames:
            fieldnames.append(col)

    global_rows = [row for row in existing_rows if row["row_kind"] == "global"]
    if len(global_rows) != 1:
        raise ValueError(f"{path} should contain exactly one global row")
    global_row = {col: global_rows[0].get(col, "") for col in fieldnames}
    for col in ("model_member_flag", "model_selection_reason", "model_velocity_rule"):
        global_row[col] = ""

    slug = SLUGS[galaxy]
    unique = pd.read_csv(MEMBERSHIP_DIR / f"{slug}_trial_unique_stars.csv", comment="#")
    observations = pd.read_csv(MEMBERSHIP_DIR / f"{slug}_trial_classified_observations.csv", comment="#")

    rows: list[dict[str, object]] = [global_row]
    audit_rows: list[dict[str, object]] = []
    for unique_row in unique.sort_values(["membership_status", "representative_star_id"]).itertuples(index=False):
        unique_series = pd.Series(unique_row._asdict())
        group = observations.loc[observations["canonical_star_id"] == unique_series["canonical_star_id"]]
        if group.empty:
            raise ValueError(f"No observation rows for {galaxy} {unique_series['representative_star_id']}")

        model_flag, model_reason = _model_flag_and_reason(galaxy, unique_series)
        row = {col: global_row.get(col, "") for col in fieldnames}
        row.update(
            {
                "row_kind": "star",
                "source_ref_id": _source_values(group["source_ref_id"]),
                "source_catalog": _source_values(group["source_catalog"]),
                "star_id": unique_series["representative_star_id"],
                "ra_deg": _number(unique_series["ra_deg"]),
                "dec_deg": _number(unique_series["dec_deg"]),
                "v_los_kms": _number(unique_series["v_los_combined_kms"]),
                "v_los_err_down_kms": _number(unique_series["v_los_combined_err_down_kms"]),
                "v_los_err_up_kms": _number(unique_series["v_los_combined_err_up_kms"]),
                "memprob": _number(unique_series["membership_probability"]),
                "memprob_err_down": "",
                "memprob_err_up": "",
                "member_flag": _member_flag_from_status(str(unique_series["membership_status"])),
                "data_status": "model_ready_unique_star_table",
                "notes": _row_notes(galaxy, unique_series, group),
                "model_member_flag": model_flag,
                "model_selection_reason": model_reason,
                "model_velocity_rule": str(unique_series["dynamical_velocity_status"]),
            }
        )
        rows.append(row)
        audit_rows.append(
            {
                "galaxy": galaxy,
                "star_id": unique_series["representative_star_id"],
                "membership_status": unique_series["membership_status"],
                "binary_status": unique_series["binary_status"],
                "dynamical_velocity_status": unique_series["dynamical_velocity_status"],
                "model_member_flag": model_flag,
                "model_selection_reason": model_reason,
                "model_velocity_rule": unique_series["dynamical_velocity_status"],
            }
        )
    return fieldnames, rows, audit_rows, comments


def main() -> None:
    all_audit_rows: list[dict[str, object]] = []
    for galaxy, filename in TRIAL_GALAXIES.items():
        fieldnames, rows, audit_rows, comments = _build_rows(galaxy, filename)
        _write_project_csv(GALAXY_DIR / filename, fieldnames, rows, comments)
        all_audit_rows.extend(audit_rows)

    AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    audit = pd.DataFrame(all_audit_rows)
    audit.to_csv(AUDIT_PATH, index=False)
    summary = (
        audit.groupby(["galaxy", "model_member_flag"])
        .size()
        .unstack(fill_value=0)
        .rename(columns={"1": "n_model_flag_1", "2": "n_model_flag_2", "3": "n_model_flag_3"})
    )
    print(summary.to_string())
    print(f"wrote {AUDIT_PATH}")


if __name__ == "__main__":
    main()
