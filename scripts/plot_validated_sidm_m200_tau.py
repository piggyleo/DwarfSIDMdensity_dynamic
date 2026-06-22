#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SUBMIT_SCRIPT = PROJECT_ROOT / "scripts/submit_validated_sidm.sh"
GALAXY_CONFIG = PROJECT_ROOT / "config/sidm_galaxies.tsv"
SIDM_PARAMETERIZATIONS = ("scale", "m200-c200", "m200-ludlow", "m200-ludlow-scatter")
M200_PARAMETERIZATIONS = ("m200-c200", "m200-ludlow", "m200-ludlow-scatter")
Y_QUANTITIES = {
    "tau": {
        "column": "tau",
        "label": r"$\tau$",
    },
    "c200": {
        "column": "c200",
        "label": r"$c_{200c}$",
    },
    "scatter": {
        "column": "concentration_scatter_sigma",
        "label": r"$\Delta\log c/\sigma_{\log c}$",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot weighted median log10(M200/Msun) and SIDM halo constraints for the "
            "galaxies listed in scripts/submit_validated_sidm.sh."
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
    return parser.parse_args()


def validated_task_ids(path: Path) -> list[int]:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"^TASK_IDS=\(([^)]*)\)", text, flags=re.MULTILINE)
    if match is None:
        raise ValueError(f"Could not find TASK_IDS in {path}")
    task_ids = [int(value) for value in match.group(1).split()]
    if not task_ids:
        raise ValueError(f"TASK_IDS is empty in {path}")
    return task_ids


def halo_output_tag(halo_model: str, sidm_parameterization: str) -> str:
    tag = halo_model.replace("-", "_")
    if halo_model == "sidm" and sidm_parameterization != "m200-ludlow-scatter":
        tag += "_" + sidm_parameterization.replace("-", "_")
    return tag


def galaxy_display_name(task_id: int) -> str:
    matches = sorted((PROJECT_ROOT / "data/galaxies").glob(f"{task_id:02d}_*.csv"))
    if len(matches) != 1:
        raise ValueError(f"Expected one galaxy CSV for task {task_id}, found {len(matches)}")
    return matches[0].stem.split("_", 1)[1].replace("_", " ")


def normalized_weights(chain: pd.DataFrame) -> np.ndarray:
    if "log_weight" in chain:
        log_weights = pd.to_numeric(chain["log_weight"], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(log_weights)
        if not np.any(finite):
            raise ValueError("log_weight contains no finite values")
        weights = np.zeros_like(log_weights)
        weights[finite] = np.exp(log_weights[finite] - np.max(log_weights[finite]))
    elif "weight" in chain:
        weights = pd.to_numeric(chain["weight"], errors="coerce").to_numpy(dtype=float)
    else:
        raise ValueError("chain has neither log_weight nor weight")

    weights = np.asarray(weights, dtype=float)
    weights[~np.isfinite(weights)] = 0.0
    if np.any(weights < 0.0):
        raise ValueError("chain contains negative sample weights")
    total = float(np.sum(weights))
    if not (np.isfinite(total) and total > 0.0):
        raise ValueError("sample weights do not have a positive finite sum")
    return weights / total


def weighted_quantiles(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    finite = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(finite):
        raise ValueError("no finite positive-weight samples remain")
    values = values[finite]
    weights = weights[finite]
    order = np.argsort(values)
    values = values[order]
    cumulative = np.cumsum(weights[order])
    cumulative /= cumulative[-1]
    return np.interp([0.16, 0.50, 0.84], cumulative, values)


def summarize_chain(
    chain_path: Path,
    expected_parameterization: str,
    y_column: str,
) -> tuple[np.ndarray, np.ndarray]:
    chain = pd.read_csv(chain_path)
    required = {"log10_m200_msun", y_column}
    missing = sorted(required.difference(chain.columns))
    if missing:
        raise ValueError(f"missing columns: {', '.join(missing)}")

    if "halo_model" in chain:
        models = set(chain["halo_model"].dropna().astype(str))
        if models != {"sidm"}:
            raise ValueError(f"unexpected halo_model values: {sorted(models)}")
    if "halo_parameterization" in chain:
        parameterizations = set(chain["halo_parameterization"].dropna().astype(str))
        if parameterizations != {expected_parameterization}:
            raise ValueError(
                "unexpected halo_parameterization values: "
                f"{sorted(parameterizations)}; expected {expected_parameterization}"
            )

    values = chain[["log10_m200_msun", y_column]].apply(pd.to_numeric, errors="coerce")
    valid = np.isfinite(values).all(axis=1).to_numpy()
    if not np.any(valid):
        raise ValueError(f"no rows have finite log10_m200_msun and {y_column}")
    values = values.loc[valid]
    weights = normalized_weights(chain)[valid]
    if not np.any(weights > 0.0):
        raise ValueError(
            f"no positive-weight rows have finite log10_m200_msun and {y_column}"
        )
    weights /= np.sum(weights)
    return (
        weighted_quantiles(values["log10_m200_msun"].to_numpy(), weights),
        weighted_quantiles(values[y_column].to_numpy(), weights),
    )


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

    task_ids = validated_task_ids(SUBMIT_SCRIPT)
    config = pd.read_csv(GALAXY_CONFIG, sep="\t").set_index("task_id")
    tag = halo_output_tag(args.halo_model, args.sidm_parameterization)
    y_column = str(Y_QUANTITIES[args.y_quantity]["column"])
    y_label = str(Y_QUANTITIES[args.y_quantity]["label"])
    rows: list[dict[str, float | str]] = []
    errors: list[str] = []

    for task_id in task_ids:
        if task_id not in config.index:
            errors.append(f"task {task_id}: missing from {GALAXY_CONFIG}")
            continue
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

    if errors:
        details = "\n".join(f"  - {error}" for error in errors)
        raise RuntimeError(f"Could not analyze every validated target:\n{details}")

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

    fig, ax = plt.subplots(figsize=(9.0, 6.8), dpi=180)
    ax.errorbar(
        x,
        y,
        xerr=xerr,
        yerr=yerr,
        fmt="o",
        ms=5.5,
        color="#176b87",
        ecolor="#6a7880",
        elinewidth=1.1,
        capsize=2.5,
        capthick=1.0,
        alpha=0.9,
        zorder=3,
    )
    if args.label_galaxies:
        annotation_order = np.argsort(x)
        for rank, row_index in enumerate(annotation_order):
            row = summary.iloc[row_index]
            level = 12 + 14 * (rank // 2)
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
                fontsize=8,
                arrowprops={
                    "arrowstyle": "-",
                    "color": "#8a969c",
                    "linewidth": 0.55,
                    "shrinkA": 2,
                    "shrinkB": 3,
                },
            )

    ax.set_xlabel(r"$\log_{10}(M_{200}/M_\odot)$")
    ax.set_ylabel(y_label)
    ax.set_title(
        f"Validated dwarf galaxies: {args.sidm_parameterization}, "
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
            / f"validated_galaxies_{tag}_m200_{args.y_quantity}.png"
        )
    figure_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_output, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {figure_output}")


if __name__ == "__main__":
    main()
