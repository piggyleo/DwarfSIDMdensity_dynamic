# Hayashi Jeans/MCMC Analysis Handoff

This document summarizes the current code structure and command-line workflow for the Willman 1 Hayashi-style Jeans analysis. It is intended as a quick-start handoff for a new conversation that will test alternative samplers such as `nautilus`.

## Important Constraint for New Sampler Tests

When testing `nautilus` or another sampler, do not modify the existing `emcee` sampler, block-MH sampler, likelihood, data loader, or Jeans/projection implementation unless there is a necessary bug fix.

Prefer adding a new standalone adapter script, for example:

- `scripts/run_willman1_nautilus.py`
- or `experiments/test_nautilus_willman1.py`

The new sampler should reuse the existing log-probability and parameter-conversion functions from `scripts/run_willman1_block_mh_fast.py` as the baseline. This keeps sampler comparisons attributable to the sampler itself rather than to hidden changes in the likelihood or model implementation.

## Project Context

The project reproduces the Hayashi et al. 2023 dwarf-galaxy dark-matter density-profile workflow with a pluggable halo model. The current working example is `Willman 1`.

Core modeling choices currently in use:

- unbinned per-star line-of-sight velocity likelihood;
- axisymmetric Plummer stellar tracer;
- generalized Hernquist dark halo, following the Hayashi paper baseline;
- projected sky coordinates are derived from each star's RA/Dec using structural centers from `data/processed/galaxy_structural_centers.csv`;
- default member selection is `member_flag == 1` with finite velocity and velocity uncertainty;
- `rho0` is still sampled, but LOS dispersion is accelerated using linear scaling with `rho0`;
- anisotropy is sampled as `q = -log10(1 - beta_z)`, not directly as `beta_z`.

## Key Input Files

- `PROJECT_MEMORY.md`: project state and recovery protocol. New sessions should read sections 2.2, 2.3, 5, and 6 first.
- `hayashi2023_method_notes.md`: detailed Hayashi method notes.
- `hayashi2023_table1_observables.csv`: Hayashi Table 1 observable parameters.
- `data/galaxies/27_Willman_1.csv`: current Willman 1 per-galaxy table.
- `data/processed/galaxy_structural_centers.csv`: structural centers used for RA/Dec projection.
- `data/validation/member_count_discrepancy_log_2026-05-06.md`: data/member-count caveats.

## Current Code Structure

### Core package: `src/hayashi_jeans/`

- `data.py`
  - `load_galaxy_data(...)`
  - `load_low_risk_galaxies(...)`
  - Reads one galaxy CSV with `pd.read_csv(..., comment="#")`.
  - Selects hard members by default: `member_flag == 1`.
  - Reads structural centers from `galaxy_structural_centers.csv`.
  - Converts RA/Dec into tangent-plane projected coordinates `x_pc`, `y_pc`.

- `params.py`
  - `HayashiParameters`
  - Container for halo, anisotropy, geometry, and systemic velocity parameters.

- `halos.py`
  - `HaloModel`
  - `GeneralizedHernquistHalo`
  - Provides the Hayashi-style generalized Hernquist halo and potential-gradient interface.

- `tracer.py`
  - `AxisymmetricPlummerTracer`
  - `intrinsic_q_from_projected(...)`
  - Implements the axisymmetric stellar tracer.

- `projection.py`
  - `AxisymmetricJeansProjector`
  - `JeansMoments`
  - Slow/reference LOS projection and Jeans moment calculation.

- `likelihood.py`
  - `GaussianVelocityLikelihood`
  - `profile_systemic_velocity(...)`
  - Implements the unbinned Gaussian LOS velocity likelihood:

```text
s_i^2 = sigma_los_i^2 + delta_v_i^2
ln L = -1/2 sum_i [ (v_i - <u>)^2 / s_i^2 + ln(2 pi s_i^2) ]
```

- `model.py`
  - `build_hernquist_projector(...)`
  - `log_likelihood_hernquist(...)`
  - `make_reference_hernquist_params(...)`
  - Wires halo, tracer/projection, and likelihood for the baseline Hernquist model.

### Optimization experiments: `experiments/`

These are experimental acceleration layers used by the current fast Willman 1 script.

- `test_rz_moment_grid_interpolation.py`
  - Defines `RZMomentGridProjector`.
  - Builds an R-z grid of Jeans moments and interpolates them for LOS projection.
  - Uses cubic-spline radial derivative for improved `dP_z/dR`.

- `test_halo_force_grid_interpolation.py`
  - Defines `RZMomentGridWithForceGrid`.
  - Adds halo force-grid interpolation under the R-z moment-grid layer.

- `test_rho0_linear_scaling.py`
  - Isolated tests showing that `rho0` can be handled by linear scaling of `sigma_los^2`.

- `test_three_layer_parameter_cache.py`
  - Experimental three-layer cache tests. This idea is not part of the current emcee production-style workflow.

## Main Sampling Script

The current main script is:

```text
scripts/run_willman1_block_mh_fast.py
```

It contains both the current block-MH sampler and an unlayered `emcee` mode.

Important functions to reuse for a new sampler:

- `load_galaxy_data(...)` from `hayashi_jeans.data`
- `initial_vector()`
- `full_log_probability_vector(...)`
- `log_prior_vector(...)`
- `vector_to_row(...)`
- `row_to_vector(...)`
- `compute_sigma_unit(...)`
- `log_likelihood_from_unit_sigma(...)`
- `validation_grid_settings(...)`

The parameter vector order is:

```text
[q_halo,
 log10_b_halo_pc,
 log10_rho0_msun_pc3,
 minus_log10_one_minus_beta_z,
 alpha,
 beta,
 gamma,
 i_deg,
 systemic_velocity_kms]
```

Column names expected in chain outputs:

```text
chain_id
step
q_halo
log10_b_halo_pc
log10_rho0_msun_pc3
minus_log10_one_minus_beta_z
beta_z
alpha
beta
gamma
i_deg
systemic_velocity_kms
log_probability
sampler
likelihood_mode
```

For nautilus, preserving this output schema is strongly recommended because the existing plotting and diagnostics scripts can then be reused directly.

## Likelihood Modes

The current script supports these modes:

- `fast`
  - Uses the caller-specified grid settings.
  - Default fast mode combines `HaloForceGrid + RZMomentGrid`.

- `validation-convergence`
  - Uses `force64x128 + moment128x256`.
  - This was the proposed convergence-validation setting.

- `validation-halo`
  - Uses R-z moment grid without the halo force grid.
  - Intended to isolate halo-force interpolation error.

- `validation-strict`
  - Uses the slow/reference projector from `build_hernquist_projector`.
  - Most expensive; use only for small validation samples.

The default fast grid values in `run_willman1_block_mh_fast.py` are:

```text
n_r = 96
n_z = 192
n_los = 160
n_force_r = 48
n_force_z = 96
```

## Example Calls

Run a single-likelihood timing test:

```bash
python scripts/run_willman1_block_mh_fast.py \
  --mode time \
  --likelihood-mode fast
```

Run the current block-MH sampler:

```bash
python scripts/run_willman1_block_mh_fast.py \
  --mode fit \
  --sampler block-mh \
  --likelihood-mode fast \
  --n-steps 10000 \
  --n-chains 24 \
  --n-processes 8 \
  --chain-output outputs/willman1_block_mh_fast_chain.csv
```

Run the existing unlayered emcee mode:

```bash
python scripts/run_willman1_block_mh_fast.py \
  --mode fit \
  --sampler emcee \
  --likelihood-mode fast \
  --n-steps 1000 \
  --n-walkers 24 \
  --n-processes 8 \
  --chain-output outputs/diagnostics/willman1_emcee_test_chain.csv
```

Plot the density profile from a chain:

```bash
python scripts/run_willman1_block_mh_fast.py \
  --mode plot \
  --chain-output outputs/willman1_block_mh_fast_chain.csv \
  --figure-output outputs/figures/willman1_block_mh_fast_density_profile.png \
  --burn 0
```

Make a corner map:

```bash
python scripts/run_willman1_block_mh_fast.py \
  --mode corner \
  --chain-output outputs/willman1_block_mh_fast_chain.csv \
  --corner-output outputs/figures/willman1_block_mh_fast_corner.png \
  --burn 0 \
  --show-systemic-in-corner
```

Validate randomly selected chain rows against another likelihood mode:

```bash
python scripts/run_willman1_block_mh_fast.py \
  --mode validate \
  --chain-output outputs/willman1_block_mh_fast_chain.csv \
  --likelihood-mode validation-convergence \
  --validation-samples 10 \
  --validation-output outputs/diagnostics/willman1_validation_convergence.csv
```

## Hayashi Figure 1 Overlay Plot

The script below creates the Willman 1 density-profile plot overlaid with digitized Hayashi Figure 1 median and 68 percent region:

```text
scripts/plot_willman1_density_profile_with_hayashi_overlay.py
```

Current call:

```bash
python scripts/plot_willman1_density_profile_with_hayashi_overlay.py \
  --chain outputs/willman1_block_mh_fast_chain.csv \
  --burn 0 \
  --figure-output outputs/figures/willman1_density_profile_with_hayashi_overlay_10000step.png
```

Important outputs used by this overlay workflow:

- `outputs/diagnostics/hayashi_fig1_willman1_crop.png`
- `outputs/diagnostics/hayashi_fig1_willman1_digitized_density.csv`
- `outputs/figures/willman1_density_profile_with_hayashi_overlay_10000step.png`

## Trace and Convergence Diagnostics

Trace diagnostics are generated by:

```text
scripts/diagnose_willman1_mcmc_trace.py
```

Current call:

```bash
python scripts/diagnose_willman1_mcmc_trace.py \
  --chain outputs/willman1_block_mh_fast_chain.csv \
  --burn 0 \
  --figure-output outputs/figures/willman1_mcmc_trace_diagnostics_10000step.png \
  --summary-output outputs/diagnostics/willman1_mcmc_convergence_summary_10000step.csv
```

The diagnostics script computes:

- integrated autocorrelation time `tau`;
- effective sample size `ESS`;
- rank-normalized split R-hat;
- acceptance-rate summaries where applicable.

Current comparison outputs:

- `outputs/diagnostics/willman1_mcmc_convergence_summary.csv`: earlier 2000-step run.
- `outputs/diagnostics/willman1_mcmc_convergence_summary_updated.csv`: earlier 5000-step run.
- `outputs/diagnostics/willman1_mcmc_convergence_summary_10000step.csv`: current 10000-step run.
- `outputs/diagnostics/willman1_mcmc_convergence_comparison_2000_5000_10000.csv`: merged comparison table.

## Current Main Chain and Figures

Current updated chain:

```text
outputs/willman1_block_mh_fast_chain.csv
```

It currently contains:

```text
24 walkers/chains x 10000 steps = 240000 rows
```

Recent figures:

- `outputs/figures/willman1_density_profile_with_hayashi_overlay_10000step.png`
- `outputs/figures/willman1_mcmc_trace_diagnostics_10000step.png`
- `outputs/figures/willman1_block_mh_fast_corner.png`

Current convergence summary:

- R-hat improved substantially from 2000 to 10000 steps.
- Slowest parameters remain `log10_b_halo_pc` and `log10_rho0_msun_pc3`.
- The 10000-step chain is improved but should not yet be treated as strictly converged.

## Suggested Nautilus Integration Strategy

Add a new script that imports the existing log-probability machinery:

```python
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for path in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from hayashi_jeans.data import load_galaxy_data
from scripts.run_willman1_block_mh_fast import (
    full_log_probability_vector,
    initial_vector,
    vector_to_row,
)
```

Recommended adapter responsibilities:

1. Load Willman 1 with:

```python
galaxy = load_galaxy_data(
    PROJECT_ROOT / "data/galaxies/27_Willman_1.csv",
    structural_centers_csv=PROJECT_ROOT / "data/processed/galaxy_structural_centers.csv",
)
```

2. Define the same 9-dimensional parameter vector and the same prior domain.

3. Let the sampler call:

```python
logp = full_log_probability_vector(
    galaxy,
    vector,
    likelihood_mode="fast",
    n_r=96,
    n_z=192,
    n_los=160,
    n_force_r=48,
    n_force_z=96,
    strict_epsrel=1.5e-2,
)
```

4. Convert samples into the standard CSV schema using `vector_to_row(...)`.

5. Write the chain to a new output path, for example:

```text
outputs/willman1_nautilus_chain.csv
```

6. Reuse existing plotting and diagnostics scripts by passing the new chain path.

## Notes for Comparing Samplers

For a fair comparison, keep these fixed unless explicitly testing them:

- input galaxy file;
- structural center file;
- member-selection rule;
- likelihood mode;
- grid sizes;
- prior ranges;
- parameter vector definition;
- output schema;
- density-profile plotting code;
- trace/convergence diagnostic code.

The first nautilus tests should focus on:

- single likelihood-call compatibility;
- sample output schema compatibility;
- posterior density-profile comparison;
- posterior marginal comparison against the current block-MH/emcee baseline;
- likelihood-mode validation on a small random subset of samples.

## Quick Recovery Prompt for a New Conversation

Use this prompt when starting the new conversation:

```text
请读取 /Users/wangkaihao/Documents/New project 4/docs/nautilus_sampler_handoff.md 和 PROJECT_MEMORY.md 的第 2.2、2.3、5、6 节。我要测试 nautilus 等替代采样器。请优先新增独立适配脚本，不要在非必要情况下修改现有 emcee 采样器、block-MH 采样器、likelihood、Jeans/projection 或绘图诊断代码。新采样器应复用 scripts/run_willman1_block_mh_fast.py 中现有的 full_log_probability_vector、prior、parameter vector 和输出 schema，以便结果可直接接入已有 density profile、corner 和 trace diagnostics 流程。
```
