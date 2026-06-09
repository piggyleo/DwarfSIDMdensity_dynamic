#!/usr/bin/env python
"""MCMC convergence diagnostics and trace plots for the Willman 1 chain."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "outputs" / ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from emcee.autocorr import AutocorrError, integrated_time


PROJECT_ROOT = Path(__file__).resolve().parents[1]

PARAMETERS = [
    ("q_halo", r"$Q$"),
    ("log10_b_halo_pc", r"$\log_{10} b_{\rm halo}$"),
    ("log10_rho0_msun_pc3", r"$\log_{10}\rho_0$"),
    ("minus_log10_one_minus_beta_z", r"$-\log_{10}(1-\beta_z)$"),
    ("alpha", r"$\alpha$"),
    ("beta", r"$\beta$"),
    ("gamma", r"$\gamma$"),
    ("i_deg", r"$i\,[{\rm deg}]$"),
    ("systemic_velocity_kms", r"$\langle u\rangle$"),
    ("log_probability", r"$\log p$"),
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--chain", default=str(PROJECT_ROOT / "outputs/willman1_block_mh_fast_chain.csv"))
    parser.add_argument("--burn", type=int, default=0)
    parser.add_argument("--figure-output", default=str(PROJECT_ROOT / "outputs/figures/willman1_mcmc_trace_diagnostics.png"))
    parser.add_argument("--summary-output", default=str(PROJECT_ROOT / "outputs/diagnostics/willman1_mcmc_convergence_summary.csv"))
    args = parser.parse_args()

    chain = pd.read_csv(args.chain)
    chain = add_anisotropy_columns(chain)
    post_burn = discard_burn(chain, args.burn)
    chain_array, steps, chain_ids, columns, labels = make_chain_array(post_burn)

    acceptance = estimate_acceptance_rates(chain_array, chain_ids)
    tau, ess = estimate_autocorr_and_ess(chain_array)
    rhat = estimate_rank_normalized_split_rhat(chain_array)

    summary = pd.DataFrame(
        {
            "parameter": columns,
            "label": labels,
            "tau_integrated": tau,
            "effective_sample_size": ess,
            "rank_normalized_split_rhat": rhat,
        }
    )
    summary["mean_acceptance_rate"] = acceptance["mean_acceptance_rate"]
    summary["min_walker_acceptance_rate"] = acceptance["min_acceptance_rate"]
    summary["max_walker_acceptance_rate"] = acceptance["max_acceptance_rate"]

    summary_output = Path(args.summary_output)
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(summary_output, index=False)

    plot_trace_diagnostics(
        chain_array,
        steps,
        chain_ids,
        labels,
        tau,
        ess,
        rhat,
        acceptance,
        burn=args.burn,
        output=Path(args.figure_output),
    )

    print(f"wrote {args.figure_output}")
    print(f"wrote {args.summary_output}")
    print(summary.to_string(index=False, float_format=lambda value: f"{value:.4g}"))


def add_anisotropy_columns(chain: pd.DataFrame) -> pd.DataFrame:
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


def discard_burn(chain: pd.DataFrame, burn: int) -> pd.DataFrame:
    if burn < 0:
        raise ValueError("--burn must be non-negative")
    if "step" not in chain.columns or "chain_id" not in chain.columns:
        raise ValueError("chain must contain chain_id and step columns")
    post_burn = chain.loc[chain["step"] > burn].copy()
    if post_burn.empty:
        raise ValueError(f"no samples remain after burn={burn}")
    return post_burn


def make_chain_array(chain: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[str]]:
    columns = [column for column, _label in PARAMETERS if column in chain.columns]
    labels = [label for column, label in PARAMETERS if column in chain.columns]
    missing = [column for column, _label in PARAMETERS[:8] if column not in chain.columns]
    if missing:
        raise ValueError(f"chain is missing required parameter columns: {missing}")

    chain_ids = np.sort(chain["chain_id"].unique())
    common_steps = None
    by_chain: list[np.ndarray] = []
    for chain_id in chain_ids:
        sub = chain.loc[chain["chain_id"] == chain_id].sort_values("step")
        steps = sub["step"].to_numpy()
        if common_steps is None:
            common_steps = steps
        elif not np.array_equal(common_steps, steps):
            raise ValueError("all chain_id traces must share the same step grid")
        by_chain.append(sub[columns].to_numpy(dtype=float))

    if common_steps is None:
        raise ValueError("empty chain")
    array = np.stack(by_chain, axis=0)
    return array, common_steps, chain_ids, columns, labels


def estimate_acceptance_rates(chain_array: np.ndarray, chain_ids: np.ndarray) -> dict[str, float]:
    if chain_array.shape[1] < 2:
        rates = np.zeros(len(chain_ids), dtype=float)
    else:
        parameter_values = chain_array[:, :, :-1]  # exclude log_probability if present
        moved = np.any(np.diff(parameter_values, axis=1) != 0.0, axis=2)
        rates = moved.mean(axis=1)
    return {
        "mean_acceptance_rate": float(np.mean(rates)),
        "min_acceptance_rate": float(np.min(rates)),
        "max_acceptance_rate": float(np.max(rates)),
    }


def estimate_autocorr_and_ess(chain_array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    emcee_array = np.transpose(chain_array, (1, 0, 2))
    n_raw = chain_array.shape[0] * chain_array.shape[1]
    try:
        tau = np.asarray(integrated_time(emcee_array, tol=0, quiet=True), dtype=float)
    except AutocorrError as exc:
        tau = np.asarray(exc.tau, dtype=float)
    tau = np.where(np.isfinite(tau) & (tau > 0.0), tau, np.nan)
    ess = n_raw / tau
    return tau, ess


def estimate_rank_normalized_split_rhat(chain_array: np.ndarray) -> np.ndarray:
    try:
        from scipy.stats import norm, rankdata
    except ImportError:
        return estimate_split_rhat(chain_array)

    n_chain, n_draw, n_param = chain_array.shape
    result = np.full(n_param, np.nan)
    for param_index in range(n_param):
        values = chain_array[:, :, param_index]
        flat = values.reshape(-1)
        finite = np.isfinite(flat)
        if finite.sum() != len(flat):
            continue
        ranks = rankdata(flat, method="average")
        z = norm.ppf((ranks - 0.375) / (len(ranks) + 0.25)).reshape(n_chain, n_draw)
        folded = np.abs(z - np.median(z))
        result[param_index] = np.nanmax([split_rhat(z), split_rhat(folded)])
    return result


def estimate_split_rhat(chain_array: np.ndarray) -> np.ndarray:
    return np.array([split_rhat(chain_array[:, :, i]) for i in range(chain_array.shape[2])])


def split_rhat(values: np.ndarray) -> float:
    n_chain, n_draw = values.shape
    half = n_draw // 2
    if half < 2:
        return np.nan
    split = np.concatenate([values[:, :half], values[:, -half:]], axis=0)
    m, n = split.shape
    chain_means = split.mean(axis=1)
    chain_vars = split.var(axis=1, ddof=1)
    within = chain_vars.mean()
    if within <= 0.0 or not np.isfinite(within):
        return np.nan
    between = n * chain_means.var(ddof=1)
    var_hat = ((n - 1) / n) * within + between / n
    return float(np.sqrt(var_hat / within))


def plot_trace_diagnostics(
    chain_array: np.ndarray,
    steps: np.ndarray,
    chain_ids: np.ndarray,
    labels: list[str],
    tau: np.ndarray,
    ess: np.ndarray,
    rhat: np.ndarray,
    acceptance: dict[str, float],
    *,
    burn: int,
    output: Path,
) -> None:
    n_param = chain_array.shape[2]
    n_cols = 2
    n_rows = int(np.ceil(n_param / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12.0, 1.95 * n_rows), dpi=180, sharex=True)
    axes_flat = np.atleast_1d(axes).ravel()
    colors = plt.cm.tab20(np.linspace(0, 1, max(len(chain_ids), 2)))

    for param_index, ax in enumerate(axes_flat[:n_param]):
        for chain_index, _chain_id in enumerate(chain_ids):
            ax.plot(steps, chain_array[chain_index, :, param_index], color=colors[chain_index % len(colors)], alpha=0.28, lw=0.55)
        median = np.nanmedian(chain_array[:, :, param_index], axis=0)
        ax.plot(steps, median, color="black", lw=1.0, alpha=0.85)
        ax.set_xlim(float(steps[0]), float(steps[-1]))
        ax.set_ylabel(labels[param_index], fontsize=9)
        ax.grid(alpha=0.2, lw=0.5)
        text = (
            rf"$\tau$={format_number(tau[param_index])}   "
            rf"ESS={format_number(ess[param_index])}   "
            rf"$\hat R$={format_number(rhat[param_index])}"
        )
        ax.text(
            0.99,
            0.96,
            text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            bbox={"facecolor": "white", "edgecolor": "0.82", "alpha": 0.86, "pad": 2.5},
        )

    for ax in axes_flat[n_param:]:
        ax.axis("off")
    for ax in axes_flat[-n_cols:]:
        ax.set_xlabel("step")

    fig.suptitle(
        "Willman 1 MCMC trace diagnostics "
        f"(burn={burn}; chains/walkers={len(chain_ids)}; "
        f"acceptance mean={acceptance['mean_acceptance_rate']:.3f}, "
        f"range={acceptance['min_acceptance_rate']:.3f}-{acceptance['max_acceptance_rate']:.3f})",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def format_number(value: float) -> str:
    if not np.isfinite(value):
        return "nan"
    if abs(value) >= 1000.0:
        return f"{value:.2g}"
    if abs(value) >= 10.0:
        return f"{value:.1f}"
    return f"{value:.3g}"


if __name__ == "__main__":
    main()
