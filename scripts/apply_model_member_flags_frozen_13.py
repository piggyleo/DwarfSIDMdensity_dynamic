#!/usr/bin/env python
"""Add model_member_flag columns to the 13 frozen Hayashi-equivalent galaxies."""

from __future__ import annotations

import csv
import io
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GALAXY_DIR = PROJECT_ROOT / "data" / "galaxies"
AUDIT_PATH = PROJECT_ROOT / "data" / "validation" / "model_member_flag_frozen_13_2026-06-16.csv"

FROZEN_GALAXIES: dict[str, str] = {
    "Canes Venatici I": "03_Canes_Venatici_I.csv",
    "Canes Venatici II": "04_Canes_Venatici_II.csv",
    "Coma Berenices": "05_Coma_Berenices.csv",
    "Eridanus II": "08_Eridanus_II.csv",
    "Horologium I": "12_Horologium_I.csv",
    "Hydra II": "13_Hydra_II.csv",
    "Leo T": "16_Leo_T.csv",
    "Pisces II": "17_Pisces_II.csv",
    "Segue 1": "19_Segue_1.csv",
    "Triangulum II": "21_Triangulum_II.csv",
    "Ursa Major I": "25_Ursa_Major_I.csv",
    "Ursa Major II": "26_Ursa_Major_II.csv",
    "Willman 1": "27_Willman_1.csv",
}


def _read_project_csv(path: Path) -> tuple[list[str], list[dict[str, str]], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    comment_start = next((i for i, line in enumerate(lines) if line.startswith("#")), len(lines))
    reader = csv.DictReader(io.StringIO("\n".join(lines[:comment_start])))
    if reader.fieldnames is None:
        raise ValueError(f"{path} has no CSV header")
    return reader.fieldnames, list(reader), lines[comment_start:]


def _comments_with_model_notes(comments: list[str]) -> list[str]:
    notes = [
        "# model_member_flag: 建模样本标记；1 为 baseline clean member，2 为默认纳入但需敏感性测试的成员，3 为不进入默认 Jeans 建模的非成员或速度不可用成员。",
        "# model_selection_reason: model_member_flag 的人工可读原因，记录成员、非成员、不确定行或速度缺失等口径。",
        "# model_velocity_rule: 写入 v_los_kms 和误差列时采用的速度选择、合并或排除规则。",
    ]
    filtered = [
        line
        for line in comments
        if not any(key in line for key in ("model_member_flag:", "model_selection_reason:", "model_velocity_rule:"))
    ]
    try:
        end = filtered.index("# COLUMN_NOTES_END")
    except ValueError:
        return filtered + ["# COLUMN_NOTES_BEGIN", *notes, "# COLUMN_NOTES_END"]
    return filtered[:end] + notes + filtered[end:]


def _write_project_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]], comments: list[str]) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    text = buffer.getvalue()
    if comments:
        text += "\n".join(_comments_with_model_notes(comments)) + "\n"
    path.write_text(text, encoding="utf-8")


def _has_valid_velocity(row: dict[str, str]) -> bool:
    for column in ("v_los_kms", "v_los_err_down_kms", "v_los_err_up_kms"):
        try:
            value = float(row.get(column, ""))
        except ValueError:
            return False
        if column != "v_los_kms" and value <= 0:
            return False
    return True


def _model_assignment(row: dict[str, str]) -> tuple[str, str, str]:
    if row.get("row_kind") != "star":
        return "", "", ""

    try:
        member_flag = int(float(row.get("member_flag", "")))
    except ValueError:
        return "3", "excluded_missing_member_flag", "excluded_missing_member_flag"

    valid_velocity = _has_valid_velocity(row)
    if member_flag == 1 and valid_velocity:
        return "1", "baseline_hayashi_equivalent_hard_member", "source_single_row_velocity"
    if member_flag == 1:
        return "3", "confirmed_member_missing_model_velocity", "missing_or_invalid_velocity"
    if member_flag == 2:
        return "3", "excluded_nonhard_or_uncertain_source_member_status", "excluded_from_default_model"
    if member_flag == 0:
        return "3", "excluded_source_nonmember", "excluded_from_default_model"
    return "3", f"excluded_unrecognized_member_flag_{member_flag}", "excluded_from_default_model"


def main() -> None:
    audit_rows: list[dict[str, str]] = []
    for galaxy, filename in FROZEN_GALAXIES.items():
        path = GALAXY_DIR / filename
        fieldnames, rows, comments = _read_project_csv(path)
        for column in ("model_member_flag", "model_selection_reason", "model_velocity_rule"):
            if column not in fieldnames:
                fieldnames.append(column)

        for row in rows:
            model_flag, reason, velocity_rule = _model_assignment(row)
            row["model_member_flag"] = model_flag
            row["model_selection_reason"] = reason
            row["model_velocity_rule"] = velocity_rule
            if row.get("row_kind") == "star":
                audit_rows.append(
                    {
                        "galaxy": galaxy,
                        "star_id": row.get("star_id", ""),
                        "member_flag": row.get("member_flag", ""),
                        "model_member_flag": model_flag,
                        "model_selection_reason": reason,
                        "model_velocity_rule": velocity_rule,
                    }
                )

        _write_project_csv(path, fieldnames, rows, comments)

    audit = pd.DataFrame(audit_rows)
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
