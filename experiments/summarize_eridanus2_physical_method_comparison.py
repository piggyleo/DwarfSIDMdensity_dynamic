#!/usr/bin/env python
"""Summarize physical-sample timing/error comparison against validation-strict."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    mge = pd.read_csv(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_physical_mge_vs_validation_strict_detail.csv")
    fast = pd.read_csv(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_fast_error_source_detail.csv")
    per_star = pd.read_csv(PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_mge_jeans_second_moment_per_star.csv")

    physical_ids = (
        per_star.groupby("selection_id")["mge_derived_sigma_los2_km2_s2"]
        .apply(lambda values: bool(np.all(np.asarray(values) > 0.0)))
        .loc[lambda s: s]
        .index.astype(int)
        .tolist()
    )

    rows: list[dict[str, float | int | str]] = []
    fast = fast[fast["selection_id"].astype(int).isin(physical_ids)].copy()
    mge = mge[mge["selection_id"].astype(int).isin(physical_ids)].copy()
    strict_seconds_by_id = {
        int(row.selection_id): float(row.reference_likelihood_seconds)
        for row in mge.itertuples(index=False)
    }
    strict_log_likelihood_by_id = {
        int(row.selection_id): float(row.reference_log_likelihood)
        for row in mge.itertuples(index=False)
    }

    for row in fast.itertuples(index=False):
        data = row._asdict()
        selection_id = int(data["selection_id"])
        strict_seconds = strict_seconds_by_id[selection_id]
        method_seconds = float(data["fast_force48_seconds"])
        rows.append(
            {
                "selection_id": selection_id,
                "source_row": int(data["source_row"]),
                "method": "RZMomentGrid+HaloForceGrid",
                "method_config": "n_r=96,n_z=192,n_los=160,n_force_r=48,n_force_z=96",
                "reference": "validation-strict",
                "reference_seconds": strict_seconds,
                "method_seconds": method_seconds,
                "speedup_vs_reference": strict_seconds / method_seconds,
                "method_log_likelihood": float(data["fast_force48_log_likelihood"]),
                "reference_log_likelihood": strict_log_likelihood_by_id[selection_id],
                "delta_log_likelihood": float(data["fast_vs_strict_delta_log_likelihood"]),
                "abs_delta_log_likelihood": float(data["fast_vs_strict_abs_delta_log_likelihood"]),
                "sigma_los2_rel_err_median": float(data["fast_vs_strict_sigma_los2_rel_err_median"]),
                "sigma_los2_rel_err_p95": float(data["fast_vs_strict_sigma_los2_rel_err_p95"]),
                "sigma_los2_rel_err_max": float(data["fast_vs_strict_sigma_los2_rel_err_max"]),
                "sigma_los_rel_err_median": float(data["fast_vs_strict_sigma_los_rel_err_median"]),
                "sigma_los_rel_err_p95": float(data["fast_vs_strict_sigma_los_rel_err_p95"]),
                "sigma_los_rel_err_max": float(data["fast_vs_strict_sigma_los_rel_err_max"]),
                "raw_nonpositive_count": 0,
            }
        )

    for row in mge.itertuples(index=False):
        data = row._asdict()
        strict_seconds = float(data["reference_likelihood_seconds"])
        method_seconds = float(data["mge_total_seconds_including_fit"])
        rows.append(
            {
                "selection_id": int(data["selection_id"]),
                "source_row": int(data["source_row"]),
                "method": "MGE/JAM",
                "method_config": "n_gauss_halo=45,n_gauss_tracer=45,n_u=96,including_fit",
                "reference": "validation-strict",
                "reference_seconds": strict_seconds,
                "method_seconds": method_seconds,
                "speedup_vs_reference": strict_seconds / method_seconds,
                "method_log_likelihood": float(data["mge_log_likelihood"]),
                "reference_log_likelihood": float(data["reference_log_likelihood"]),
                "delta_log_likelihood": float(data["mge_vs_reference_delta_log_likelihood"]),
                "abs_delta_log_likelihood": float(data["mge_vs_reference_abs_delta_log_likelihood"]),
                "sigma_los2_rel_err_median": float(data["mge_vs_reference_sigma_los2_rel_err_median"]),
                "sigma_los2_rel_err_p95": float(data["mge_vs_reference_sigma_los2_rel_err_p95"]),
                "sigma_los2_rel_err_max": float(data["mge_vs_reference_sigma_los2_rel_err_max"]),
                "sigma_los_rel_err_median": float(data["mge_vs_reference_sigma_los_rel_err_median"]),
                "sigma_los_rel_err_p95": float(data["mge_vs_reference_sigma_los_rel_err_p95"]),
                "sigma_los_rel_err_max": float(data["mge_vs_reference_sigma_los_rel_err_max"]),
                "raw_nonpositive_count": int(data["mge_sigma_los2_nonpositive_count"]),
            }
        )

    detail = pd.DataFrame(rows).sort_values(["selection_id", "method"]).reset_index(drop=True)
    detail_path = PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_physical_method_comparison_vs_validation_strict_detail.csv"
    detail.to_csv(detail_path, index=False)

    metric_cols = [
        "reference_seconds",
        "method_seconds",
        "speedup_vs_reference",
        "abs_delta_log_likelihood",
        "sigma_los2_rel_err_median",
        "sigma_los2_rel_err_p95",
        "sigma_los2_rel_err_max",
        "sigma_los_rel_err_median",
        "sigma_los_rel_err_p95",
        "sigma_los_rel_err_max",
    ]
    summary_rows = []
    for method, group in detail.groupby("method", sort=False):
        row: dict[str, float | int | str] = {
            "method": method,
            "n_physical_samples": int(group["selection_id"].nunique()),
            "selection_ids": ",".join(str(x) for x in sorted(group["selection_id"].astype(int).unique())),
        }
        for metric in metric_cols:
            row[f"{metric}_median"] = float(group[metric].median())
            row[f"{metric}_max"] = float(group[metric].max())
            row[f"{metric}_min"] = float(group[metric].min())
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary_path = PROJECT_ROOT / "outputs/diagnostics/eridanus_ii_physical_method_comparison_vs_validation_strict_summary.csv"
    summary.to_csv(summary_path, index=False)

    print(f"physical_ids={physical_ids}")
    print(f"wrote {detail_path}")
    print(detail.to_string(index=False))
    print(f"wrote {summary_path}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
