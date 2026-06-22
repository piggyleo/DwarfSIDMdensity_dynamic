#!/usr/bin/env python
"""Rebuild the Hercules project CSV from the complete Simon & Geha table."""

from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GALAXY_PATH = PROJECT_ROOT / "data" / "galaxies" / "11_Hercules.csv"
RAW_PATH = PROJECT_ROOT / "data" / "raw" / "jsimon" / "Herc_feh.dat"
AUDIT_PATH = (
    PROJECT_ROOT
    / "data"
    / "validation"
    / "hercules_full_source_rebuild_2026-06-14.csv"
)


def _read_project_csv(path: Path) -> tuple[list[str], list[dict[str, str]], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    comment_start = next((i for i, line in enumerate(lines) if line.startswith("#")), len(lines))
    reader = csv.DictReader(io.StringIO("\n".join(lines[:comment_start])))
    if reader.fieldnames is None:
        raise ValueError(f"{path} has no CSV header")
    return reader.fieldnames, list(reader), lines[comment_start:]


def _write_project_csv(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    comments: list[str],
) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    text = buffer.getvalue()
    if comments:
        text += "\n".join(comments) + "\n"
    path.write_text(text, encoding="utf-8")


def _optional_number(value: str) -> str:
    try:
        number = float(value)
    except ValueError:
        return ""
    return str(number) if math.isfinite(number) else ""


def _read_raw_table() -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for line in RAW_PATH.read_text(encoding="utf-8").splitlines()[2:]:
        fields = line.split()
        if not fields:
            continue
        if len(fields) != 10:
            raise ValueError(f"Unexpected Hercules raw row with {len(fields)} fields: {line}")
        star_id, ra, dec, r_mag, g_minus_r, velocity, velocity_error, feh, feh_error, member = fields
        records.append(
            {
                "star_id": star_id,
                "ra_deg": str(float(ra)),
                "dec_deg": str(float(dec)),
                "r_mag": str(float(r_mag)),
                "g_minus_r": str(float(g_minus_r)),
                "v_los_kms": str(float(velocity)),
                "v_los_err_kms": str(float(velocity_error)),
                "feh_star_dex": _optional_number(feh),
                "feh_star_err_dex": _optional_number(feh_error),
                "member_flag": member,
            }
        )
    if len(records) != 86:
        raise ValueError(f"Expected 86 Hercules source rows, found {len(records)}")
    if sum(record["member_flag"] == "1" for record in records) != 30:
        raise ValueError("Expected 30 Hercules source members")
    return records


def _build_rows(
    fieldnames: list[str],
    old_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    global_rows = [row for row in old_rows if row["row_kind"] == "global"]
    if len(global_rows) != 1:
        raise ValueError(f"Expected one Hercules global row, found {len(global_rows)}")
    global_row = global_rows[0].copy()
    global_row["data_status"] = "star_data_available"
    global_row["notes"] = (
        "global parameters from Hayashi Table 1; all 86 Simon & Geha (2007) "
        "machine-table targets are retained, including 30 source hard members and "
        "56 nonmembers; source_ref_id=5 on star rows follows the Hayashi Table 1 "
        "Kirby et al. 2013b reference chain, while source_catalog and row notes "
        "preserve the upstream kinematic and membership provenance."
    )
    old_stars = {
        row["star_id"]: row for row in old_rows if row["row_kind"] == "star"
    }
    if len(old_stars) not in {21, 86}:
        raise ValueError(f"Expected 21 legacy or 86 rebuilt Hercules rows, found {len(old_stars)}")
    kirby_overlap = {
        star_id: row
        for star_id, row in old_stars.items()
        if (
            "ref5 metallicity catalog row" in row.get("notes", "")
            or "Kirby et al. 2013b ref5 metallicity overlap" in row.get("notes", "")
        )
    }
    if len(kirby_overlap) != 21:
        raise ValueError(f"Expected 21 Kirby metallicity overlap rows, found {len(kirby_overlap)}")

    rebuilt: list[dict[str, str]] = [global_row]
    audit: list[dict[str, str]] = []
    for source in _read_raw_table():
        old = kirby_overlap.get(source["star_id"])
        row = {column: global_row.get(column, "") for column in fieldnames}
        row.update(
            {
                "row_kind": "star",
                "source_ref_id": "5",
                "source_catalog": "jsimon_data/Herc_feh.dat",
                "star_id": source["star_id"],
                "ra_deg": source["ra_deg"],
                "dec_deg": source["dec_deg"],
                "v_los_kms": source["v_los_kms"],
                "v_los_err_down_kms": source["v_los_err_kms"],
                "v_los_err_up_kms": source["v_los_err_kms"],
                "memprob": "1.0" if source["member_flag"] == "1" else "0.0",
                "memprob_err_down": "",
                "memprob_err_up": "",
                "member_flag": source["member_flag"],
                "data_status": "star_data_available",
            }
        )

        if old is not None:
            row["feh_star_dex"] = old["feh_star_dex"]
            row["feh_star_err_down_dex"] = old["feh_star_err_down_dex"]
            row["feh_star_err_up_dex"] = old["feh_star_err_up_dex"]
            feh_provenance = "Kirby et al. 2013b ref5 metallicity overlap"
            row["notes"] = (
                "complete Simon & Geha (2007) source-table row; source Mem mapped to "
                f"member_flag={source['member_flag']}; velocity and membership traced through "
                "jsimon_data/Herc_feh.dat; [Fe/H] preserved from the Kirby et al. 2013b ref5 "
                "metallicity overlap; source_ref_id=5 follows Hayashi Table 1."
            )
        else:
            row["feh_star_dex"] = source["feh_star_dex"]
            row["feh_star_err_down_dex"] = source["feh_star_err_dex"]
            row["feh_star_err_up_dex"] = source["feh_star_err_dex"]
            feh_provenance = "Simon & Geha 2007 upstream machine table"
            row["notes"] = (
                "complete Simon & Geha (2007) source-table row newly restored; source Mem "
                f"mapped to member_flag={source['member_flag']}; velocity, membership, and "
                "[Fe/H] traced through jsimon_data/Herc_feh.dat; source_ref_id=5 records the "
                "Hayashi Table 1 primary-reference chain, while this note preserves the "
                "upstream machine-table provenance."
            )

        rebuilt.append(row)
        audit.append(
            {
                "star_id": source["star_id"],
                "present_in_previous_21_row_catalog": "1" if old is not None else "0",
                "source_member_flag": source["member_flag"],
                "ra_deg": source["ra_deg"],
                "dec_deg": source["dec_deg"],
                "v_los_kms": source["v_los_kms"],
                "v_los_err_down_kms": source["v_los_err_kms"],
                "v_los_err_up_kms": source["v_los_err_kms"],
                "feh_provenance": feh_provenance,
                "action": "preserved_and_retraced" if old is not None else "restored_from_complete_source",
            }
        )

    return rebuilt, audit


def _write_audit(rows: list[dict[str, str]]) -> None:
    AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with AUDIT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.write("\n# COLUMN_NOTES_BEGIN\n")
        handle.write("# source_member_flag: hard Mem value from jsimon_data/Herc_feh.dat.\n")
        handle.write("# present_in_previous_21_row_catalog: 1 if the star was in the Kirby metallicity overlap.\n")
        handle.write("# feh_provenance: star-level metallicity source retained in the rebuilt project row.\n")
        handle.write("# COLUMN_NOTES_END\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    fieldnames, old_rows, comments = _read_project_csv(GALAXY_PATH)
    rebuilt, audit = _build_rows(fieldnames, old_rows)
    if not args.dry_run:
        _write_project_csv(GALAXY_PATH, fieldnames, rebuilt, comments)
        _write_audit(audit)

    print(
        "Hercules rebuild: "
        f"source_rows={len(audit)} "
        f"members={sum(row['source_member_flag'] == '1' for row in audit)} "
        f"preserved_kirby_overlap={sum(row['present_in_previous_21_row_catalog'] == '1' for row in audit)} "
        f"dry_run={args.dry_run}"
    )


if __name__ == "__main__":
    main()
