#!/usr/bin/env python
"""Repair source membership labels that were flattened incorrectly."""

from __future__ import annotations

import argparse
import csv
import io
from dataclasses import dataclass
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
GALAXY_DIR = PROJECT_ROOT / "data" / "galaxies"
AUDIT_PATH = PROJECT_ROOT / "data" / "validation" / "source_member_semantics_repair_2026-06-14.csv"

SEGUE2_HORIZONTAL_BRANCH_MEMBERS = {
    "J021907.59+201220.8",
    "J021918.49+201021.9",
    "J021921.13+200740.2",
    "J021934.68+201144.3",
}

BOOTES1_CATALOG_PAPER_CONFLICTS = {
    "Boo1_9",
    "Boo1_35",
    "Boo1_106",
    "Boo1_113",
}


@dataclass(frozen=True)
class Repair:
    galaxy: str
    star_id: str
    source_catalog: str
    source_member_label: str
    old_member_flag: str
    new_member_flag: str
    evidence: str


def _read_project_csv(path: Path) -> tuple[list[str], list[dict[str, str]], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    comment_start = next((i for i, line in enumerate(lines) if line.startswith("#")), len(lines))
    data_lines = lines[:comment_start]
    comments = lines[comment_start:]
    reader = csv.DictReader(io.StringIO("\n".join(data_lines)))
    if reader.fieldnames is None:
        raise ValueError(f"{path} has no CSV header")
    return reader.fieldnames, list(reader), comments


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


def _append_note(note: str, message: str) -> str:
    if message in note:
        return note
    return f"{message}; {note}" if note else message


def _remove_obsolete_note_fragments(note: str) -> str:
    fragments = [
        "source member status is unknown/blank for this row",
        "member_flag set to 2",
    ]
    parts = [part.strip() for part in note.split(";")]
    return "; ".join(part for part in parts if part and part not in fragments)


def repair_grus1(*, dry_run: bool) -> list[Repair]:
    path = GALAXY_DIR / "09_Grus_1.csv"
    fieldnames, rows, comments = _read_project_csv(path)
    repairs: list[Repair] = []

    for row in rows:
        if row.get("row_kind") != "star" or "MEM=CM" not in row.get("notes", ""):
            continue
        old_flag = row["member_flag"]
        if old_flag not in {"0", "2"}:
            raise ValueError(f"Unexpected Grus 1 CM member_flag={old_flag!r} for {row['star_id']}")
        row["member_flag"] = "2"
        row["notes"] = _append_note(
            row.get("notes", ""),
            "source_member_label=CM; candidate member with compatible velocity but no metallicity",
        )
        repairs.append(
            Repair(
                galaxy="Grus 1",
                star_id=row["star_id"],
                source_catalog=row["source_catalog"],
                source_member_label="CM",
                old_member_flag=old_flag,
                new_member_flag="2",
                evidence="Chiti et al. 2022: CM denotes a candidate member lacking a metallicity measurement.",
            )
        )

    if len(repairs) != 4:
        raise ValueError(f"Expected 4 Grus 1 CM rows, found {len(repairs)}")
    if not dry_run:
        _write_project_csv(path, fieldnames, rows, comments)
    return repairs


def repair_segue2(*, dry_run: bool) -> list[Repair]:
    path = GALAXY_DIR / "20_Segue_2.csv"
    fieldnames, rows, comments = _read_project_csv(path)
    repairs: list[Repair] = []

    for row in rows:
        if (
            row.get("row_kind") != "star"
            or row.get("source_catalog") != "J/ApJ/770/16/table2"
            or row.get("star_id") not in SEGUE2_HORIZONTAL_BRANCH_MEMBERS
        ):
            continue
        old_flag = row["member_flag"]
        if old_flag not in {"1", "2"}:
            raise ValueError(f"Unexpected Segue 2 B member_flag={old_flag!r} for {row['star_id']}")
        row["member_flag"] = "1"
        row["notes"] = _remove_obsolete_note_fragments(row.get("notes", ""))
        row["notes"] = _append_note(
            row.get("notes", ""),
            "source_member_label=B; confirmed Segue 2 horizontal-branch member",
        )
        repairs.append(
            Repair(
                galaxy="Segue 2",
                star_id=row["star_id"],
                source_catalog=row["source_catalog"],
                source_member_label="B",
                old_member_flag=old_flag,
                new_member_flag="1",
                evidence="Kirby et al. 2013 VizieR ReadMe: B means yes, horizontal-branch star.",
            )
        )

    if len(repairs) != 4:
        raise ValueError(f"Expected 4 Segue 2 B rows, found {len(repairs)}")
    if not dry_run:
        _write_project_csv(path, fieldnames, rows, comments)
    return repairs


def repair_bootes1_catalog_conflicts(*, dry_run: bool) -> list[Repair]:
    path = GALAXY_DIR / "02_Bootes_I.csv"
    fieldnames, rows, comments = _read_project_csv(path)
    repairs: list[Repair] = []

    for row in rows:
        if row.get("row_kind") != "star" or row.get("star_id") not in BOOTES1_CATALOG_PAPER_CONFLICTS:
            continue
        old_flag = row["member_flag"]
        if old_flag not in {"1", "2"}:
            raise ValueError(
                f"Unexpected Bootes I catalog-conflict member_flag={old_flag!r} for {row['star_id']}"
            )
        row["member_flag"] = "2"
        row["notes"] = _append_note(
            row.get("notes", ""),
            (
                "source_member_label=published_catalog_conflict; VizieR Member=1 "
                "but star is absent from Jenkins et al. 2021 Table 2 final M/VCNM list"
            ),
        )
        repairs.append(
            Repair(
                galaxy="Bootes I",
                star_id=row["star_id"],
                source_catalog=row["source_catalog"],
                source_member_label="published_catalog_conflict",
                old_member_flag=old_flag,
                new_member_flag="2",
                evidence=(
                    "Jenkins et al. 2021 Table 2 contains 69 final M stars; these four "
                    "VizieR Member=1 rows are absent from both the final M and VCNM list."
                ),
            )
        )

    if len(repairs) != 4:
        raise ValueError(f"Expected 4 Bootes I catalog conflicts, found {len(repairs)}")
    if not dry_run:
        _write_project_csv(path, fieldnames, rows, comments)
    return repairs


def repair_tucana2_member_compilation(*, dry_run: bool) -> list[Repair]:
    path = GALAXY_DIR / "22_Tucana_2.csv"
    fieldnames, rows, comments = _read_project_csv(path)
    repairs: list[Repair] = []

    for row in rows:
        if row.get("row_kind") != "star" or row.get("source_catalog") != "J/AJ/165/55/table6":
            continue
        old_flag = row["member_flag"]
        if old_flag not in {"1", "2"}:
            raise ValueError(
                f"Unexpected Tucana 2 compiled-member flag={old_flag!r} for {row['star_id']}"
            )
        row["member_flag"] = "1"
        row["notes"] = _remove_obsolete_note_fragments(row.get("notes", ""))
        row["notes"] = _append_note(
            row.get("notes", ""),
            "source_member_label=confirmed_member_velocity_compilation",
        )
        repairs.append(
            Repair(
                galaxy="Tucana 2",
                star_id=row["star_id"],
                source_catalog=row["source_catalog"],
                source_member_label="confirmed_member_velocity_compilation",
                old_member_flag=old_flag,
                new_member_flag="1",
                evidence=(
                    "Chiti et al. 2023 Appendix: Table 6 compiles all literature "
                    "velocity measurements of Tucana II members."
                ),
            )
        )

    if len(repairs) != 60:
        raise ValueError(f"Expected 60 Tucana 2 member-observation rows, found {len(repairs)}")
    if not dry_run:
        _write_project_csv(path, fieldnames, rows, comments)
    return repairs


def write_audit(repairs: list[Repair]) -> None:
    AUDIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    existing: dict[tuple[str, str, str], dict[str, str]] = {}
    if AUDIT_PATH.exists():
        with AUDIT_PATH.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                existing[(row["galaxy"], row["star_id"], row["source_catalog"])] = row

    with AUDIT_PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(Repair.__dataclass_fields__))
        writer.writeheader()
        for repair in repairs:
            row = repair.__dict__.copy()
            key = (repair.galaxy, repair.star_id, repair.source_catalog)
            if key in existing:
                row["old_member_flag"] = existing[key]["old_member_flag"]
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    repairs = (
        repair_grus1(dry_run=args.dry_run)
        + repair_segue2(dry_run=args.dry_run)
        + repair_bootes1_catalog_conflicts(dry_run=args.dry_run)
        + repair_tucana2_member_compilation(dry_run=args.dry_run)
    )
    if not args.dry_run:
        write_audit(repairs)
    changed = sum(repair.old_member_flag != repair.new_member_flag for repair in repairs)
    print(f"validated={len(repairs)} changed={changed} dry_run={args.dry_run}")


if __name__ == "__main__":
    main()
