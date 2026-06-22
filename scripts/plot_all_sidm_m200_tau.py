#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_validated_sidm_m200_tau import (
    GALAXY_CONFIG,
    M200_PARAMETERIZATIONS,
    PROJECT_ROOT,
    SIDM_PARAMETERIZATIONS,
    Y_QUANTITIES,
    galaxy_display_name,
    halo_output_tag,
    summarize_chain,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot weighted median log10(M200/Msun) and SIDM halo constraints for "
            "all galaxies listed in config/sidm_galaxies.tsv."
        )
    )
    parser.add_argument(
        "--halo-model",
        choices=("generalized-hernquist", "sidm"),
        default="sidm",
    )
    parser.add_argument(
        "--sidm-parameterization",
        choices=SIDM_PARAMETERIZATIONS,
        default="m200-ludlow-scatter",
    )
    parser.add_argument(
        "--y-quantity",
        choices=tuple(Y_QUANTITIES),
        default="tau",
        help=(
            "Quantity plotted on the vertical axis: tau, c200, or scatter "
            "(s_c = delta log(c) / sigma_log(c)). Default: tau."
        ),
    )
    parser.add_argument(
        "--chain-dir",
        type=Path,
        default=PROJECT_ROOT / "outputs",
        help="Directory containing the per-galaxy chain CSV files.",
    )
    parser.add_argument(
        "--figure-output",
        type=Path,
        default=None,
        help="Optional output PNG path.",
    )
    parser.add_argument(
        "--label-galaxies",
        action="store_true",
        help="Add galaxy names beside the plotted points. Labels are hidden by default.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Skip missing or unreadable chains instead of failing.",
    )
    return parser.parse_args()


def all_task_ids(config: pd.DataFrame) -> list[int]:
    task_ids = [int(value) for value in config.index]
    if not task_ids:
        raise ValueError(f"No task IDs found in {GALAXY_CONFIG}")
    return sorted(task_ids)


def main() -> None:
    args = parse_args()
    if args.halo_model != "sidm":
        raise ValueError("M200 and tau are SIDM quantities; use --halo-model sidm")
    if args.sidm_parameterization not in M200_PARAMETERIZATIONS:
        raise ValueError(
            "--sidm-parameterization scale does not contain M200; choose "
            "m200-c200, m200-ludlow, or m200-ludlow-scatter"
        )
    if args.y_quantity == "scatter" and args.sidm_parameterization != "m200-ludlow-scatter":
        raise ValueError(
            "--y-quantity scatter requires "
            "--sidm-parameterization m200-ludlow-scatter"
        )

    config = pd.read_csv(GALAXY_CONFIG, sep="\t").set_index("task_id")
    task_ids = all_task_ids(config)
    tag = halo_output_tag(args.halo_model, args.sidm_parameterization)
    y_column = str(Y_QUANTITIES[args.y_quantity]["column"])
    y_label = str(Y_QUANTITIES[args.y_quantity]["label"])
    rows: list[dict[str, float | str]] = []
    errors: list[str] = []

    for task_id in task_ids:
        galaxy_slug = str(config.at[task_id, "galaxy_slug"])
        chain_path = args.chain_dir / f"{galaxy_slug}_{tag}_chain.csv"
        try:
            m200, y_summary = summarize_chain(
                chain_path,
                args.sidm_parameterization,
                y_column,
            )
            rows.append(
                {
                    "galaxy": galaxy_display_name(task_id),
                    "log10_m200_q16": m200[0],
                    "log10_m200_median": m200[1],
                    "log10_m200_q84": m200[2],
                    f"{args.y_quantity}_q16": y_summary[0],
                    f"{args.y_quantity}_median": y_summary[1],
                    f"{args.y_quantity}_q84": y_summary[2],
                }
            )
        except (FileNotFoundError, OSError, ValueError, pd.errors.ParserError) as exc:
            errors.append(f"{galaxy_slug}: {chain_path}: {exc}")

    if errors and not args.allow_missing:
        details = "\n".join(f"  - {error}" for error in errors)
        raise RuntimeError(f"Could not analyze every configured target:\n{details}")
    if errors:
        print("Skipped targets:")
        for error in errors:
            print(f"  - {error}")
    if not rows:
        raise RuntimeError("No usable chains were found")

    summary = pd.DataFrame(rows)
    print(
        summary.to_string(
            index=False,
            columns=[
                "galaxy",
                "log10_m200_median",
                "log10_m200_q16",
                "log10_m200_q84",
                f"{args.y_quantity}_median",
                f"{args.y_quantity}_q16",
                f"{args.y_quantity}_q84",
            ],
            float_format=lambda value: f"{value:.5g}",
        )
    )

    x = summary["log10_m200_median"].to_numpy(dtype=float)
    y = summary[f"{args.y_quantity}_median"].to_numpy(dtype=float)
    xerr = np.vstack(
        [
            x - summary["log10_m200_q16"].to_numpy(dtype=float),
            summary["log10_m200_q84"].to_numpy(dtype=float) - x,
        ]
    )
    yerr = np.vstack(
        [
            y - summary[f"{args.y_quantity}_q16"].to_numpy(dtype=float),
            summary[f"{args.y_quantity}_q84"].to_numpy(dtype=float) - y,
        ]
    )

    fig, ax = plt.subplots(figsize=(9.8, 7.2), dpi=180)
    ax.errorbar(
        x,
        y,
        xerr=xerr,
        yerr=yerr,
        fmt="o",
        ms=5.2,
        color="#176b87",
        ecolor="#6a7880",
        elinewidth=1.05,
        capsize=2.4,
        capthick=0.95,
        alpha=0.9,
        zorder=3,
    )
    if args.label_galaxies:
        annotation_order = np.argsort(x)
        for rank, row_index in enumerate(annotation_order):
            row = summary.iloc[row_index]
            level = 10 + 12 * (rank // 2)
            dy = level if rank % 2 == 0 else -level
            dx = 7 if rank % 4 < 2 else -7
            ax.annotate(
                str(row["galaxy"]),
                (
                    float(row["log10_m200_median"]),
                    float(row[f"{args.y_quantity}_median"]),
                ),
                xytext=(dx, dy),
                textcoords="offset points",
                ha="left" if dx > 0 else "right",
                va="bottom" if dy > 0 else "top",
                fontsize=7.5,
                arrowprops={
                    "arrowstyle": "-",
                    "color": "#8a969c",
                    "linewidth": 0.5,
                    "shrinkA": 2,
                    "shrinkB": 3,
                },
            )

    ax.set_xlabel(r"$\log_{10}(M_{200}/M_\odot)$")
    ax.set_ylabel(y_label)
    ax.set_title(
        f"All configured dwarf galaxies: {args.sidm_parameterization}, "
        f"{args.y_quantity}"
    )
    ax.grid(True, alpha=0.22, linewidth=0.7)
    ax.margins(x=0.12, y=0.15)
    fig.tight_layout()

    figure_output = args.figure_output
    if figure_output is None:
        figure_output = (
            PROJECT_ROOT
            / "outputs/figures"
            / f"all_galaxies_{tag}_m200_{args.y_quantity}.png"
        )
    figure_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_output, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {figure_output}")


if __name__ == "__main__":
    main()
