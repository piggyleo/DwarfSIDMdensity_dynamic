# Nautilus Sampler and Pluggable Halo Handoff

Updated: 2026-06-09

This document describes the current production-facing Nautilus entry point,
the generalized Hernquist baseline, the SIDM PSIDM-25 integration, and the
real-axis signed-MGE acceleration/validation workflow.

## Current Branch and Entry Point

Development and validation were performed on:

```text
codex/setup-repo
```

The current core sampling entry point is:

```text
scripts/run_galaxy_nautilus.py
```

The older `scripts/run_willman1_block_mh_fast.py` remains useful as historical
comparison code, but it is not the interface to extend for new halo models.

## Compatibility Contract

The generalized Hernquist behavior remains the default:

```bash
PYTHONPATH=src python scripts/run_galaxy_nautilus.py \
  --mode point \
  --galaxy willman1
```

The model is changed only when `--halo-model sidm` is supplied. Existing
generalized Hernquist parameter names, priors, chain columns, fast-grid modes,
and the 45-component MGE defaults are preserved.

The shared runtime halo contract is defined by
`src/hayashi_jeans/halos.py`:

```python
class HaloModel:
    def density(self, r_cyl_pc, z_pc): ...
    def potential_gradients(self, r_cyl_pc, z_pc): ...
```

`SpheroidallyStratifiedHalo` implements the common homeoidal force integral for
profiles stratified on

```text
m^2 = R^2 + z^2 / q^2.
```

Both `GeneralizedHernquistHalo` and `SIDMPSIDM25Halo` use this interface.
`build_projector(...)` in `src/hayashi_jeans/model.py` is the model-independent
constructor; `build_hernquist_projector(...)` remains the compatibility wrapper.

## SIDM Model

`SIDMPSIDM25Halo` implements the density-only PSIDM-25 profile ported from:

```text
/Users/wangkaihao/sidm/SIDM_Lensing_Model-main/lib/SIDM_Parametric_Model_jax.py
```

No lensing code is included. The spherical profile is extended to the existing
axisymmetric Jeans convention by evaluating the density at ellipsoidal radius
`m`.

The physical profile is

```text
rho(m, tau) = rho_s0 * f(m / rs0, tau),
```

where the fitted evolution functions determine the time-dependent density
normalization, scale radius, core radius, transition sharpness, and inner
slope. The supported domain is:

```text
0 <= tau <= 1.08
```

### SIDM Parameterizations

The sampler supports four parameter vectors.

Direct scale parameters:

```text
q_halo
log10_rs0_pc
log10_rho_s0_msun_pc3
tau
minus_log10_one_minus_beta_z
i_deg
systemic_velocity_kms
```

Free `M200` and concentration:

```text
q_halo
log10_m200_msun
log10_c200
tau
minus_log10_one_minus_beta_z
i_deg
systemic_velocity_kms
```

Ludlow16 concentration:

```text
q_halo
log10_m200_msun
tau
minus_log10_one_minus_beta_z
i_deg
systemic_velocity_kms
```

Ludlow16 with sampled log-concentration scatter coordinate:

```text
q_halo
log10_m200_msun
concentration_scatter_sigma
tau
minus_log10_one_minus_beta_z
i_deg
systemic_velocity_kms
```

`concentration_scatter_sigma` has the truncated Gaussian prior

```text
concentration_scatter_sigma ~ Normal(0, 3^2), truncated to [-4, 4].
```

The symmetric truncation is configured with
`--concentration-scatter-truncation`; its default is `4`. Nautilus applies the
normalized truncated distribution through its `Prior` transform, so samples
outside the interval are never proposed. The likelihood callback does not add
the prior density a second time.

Changing this prior makes an existing Nautilus checkpoint incompatible. Start
with a new checkpoint/output name instead of using `--resume`.

The conversion layer is in
`src/hayashi_jeans/halo_parameterizations.py`. It resolves every
parameterization to one `SIDMPhysicalParameters` descriptor and then to the
same `SIDMPSIDM25Halo`.

For all `M200` parameterizations, redshift is fixed before sampling and is not
a sampled parameter:

```bash
--halo-redshift 0.01
```

Planck15 is currently the only supported cosmology. Ludlow16 evaluation uses
`colossus`.

## Signed-MGE Method

The MGE/JAM implementation lives in `src/hayashi_jeans/mge.py`.

For generalized Hernquist, the existing Shajib complex-plane analytic
decomposition remains available. For SIDM, `auto` validates that decomposition
on the real radial axis and falls back to a real-axis least-squares fit when
needed. This avoids evaluating the SIDM profile near or across its complex
branch structure.

The fallback solves a real linear problem:

```text
rho(r_i) ~= sum_j A_j exp[-r_i^2 / (2 sigma_j^2)].
```

The coefficients `A_j` are signed. Negative individual Gaussian amplitudes are
allowed; physicality applies to the reconstructed total density and Jeans
moments, not to each basis coefficient separately.

SIDM sampling defaults:

```text
n_gauss_halo = 60
n_gauss_tracer = 60
real_fit_radii = 1600
real_lstsq_rcond = 1e-12
n_u = 96
r_min_pc = 1e-2
r_max_pc = 2e5
```

Generalized Hernquist retains `45/45` Gaussian defaults.

The MGE likelihood computes the projected second moment using the
one-dimensional JAM quadrature in `mge_los_second_moment(...)`. It does not run
the slow line-of-sight Jeans projector for each sample.

### Physicality Guard

The MGE physicality check is enabled by default. It evaluates the local raw
`v_phi^2` on a compact central meridional grid and rejects non-finite or
significantly negative solutions. Multiplying the halo density by a positive
normalization cannot change these signs, so this check can use a unit-density
halo where linear scaling applies.

Disable it only for diagnostics:

```bash
--no-mge-physicality-check
```

## Error Attribution

`experiments/validate_sidm_mge.py` uses three distinct calculations at the
same physical parameter point:

1. Original SIDM density plus original Plummer tracer through the slow
   projector. This is the reference.
2. Signed-MGE halo plus signed-MGE tracer through the same slow projector.
   The difference from step 1 isolates density-decomposition error.
3. The same MGE arrays through the fast one-dimensional JAM calculation.
   The difference from step 2 isolates the MGE/JAM moment implementation and
   finite integration-boundary error.

Do not validate the fast MGE result only against an MGE-based reference: that
would hide density-decomposition error.

High-accuracy Willman 1 checks used:

```bash
PYTHONPATH=src python experiments/validate_sidm_mge.py \
  --taus 0.5,1.08 \
  --strict-taus 0.5,1.08 \
  --strict-factor 40 \
  --strict-epsrel 1e-2 \
  --output outputs/diagnostics/sidm_mge_three_layer_factor40.csv
```

Results:

| tau | MGE profile - original logL | MGE/JAM - MGE profile logL | MGE/JAM - original logL |
|---:|---:|---:|---:|
| 0.50 | -0.000005 | +0.000687 | +0.000682 |
| 1.08 | -0.000057 | -0.000939 | -0.000995 |

The signed-MGE profile substitution is therefore much smaller than
`1e-3` in log likelihood at these checks. The remaining error near `1e-3`
comes mainly from the fast JAM calculation relative to the finite-bound slow
reference.

The corresponding 95th-percentile density errors were about `2.5e-4` and
`2.2e-4`. Slow original-profile evaluations took roughly 87-98 seconds,
slow MGE-profile evaluations 119-139 seconds, and the isolated JAM moment
calculation about 0.27 seconds on this machine. These are diagnostic timings,
not full sampler-throughput measurements.

## Commands

SIDM direct-scale point check:

```bash
PYTHONPATH=src python scripts/run_galaxy_nautilus.py \
  --mode point \
  --likelihood-mode mge \
  --halo-model sidm \
  --sidm-parameterization scale \
  --galaxy willman1
```

Free `M200-c200` point check:

```bash
PYTHONPATH=src python scripts/run_galaxy_nautilus.py \
  --mode point \
  --likelihood-mode mge \
  --halo-model sidm \
  --sidm-parameterization m200-c200 \
  --halo-redshift 0.01 \
  --galaxy willman1
```

Ludlow16 forms use `--sidm-parameterization m200-ludlow` or
`m200-ludlow-scatter`, also with a fixed `--halo-redshift`.

Set a different symmetric scatter truncation with:

```bash
--concentration-scatter-truncation 3
```

For the Slurm array worker, the equivalent environment variable is
`CONCENTRATION_SCATTER_TRUNCATION`; it defaults to `4`.

Schema check:

```bash
PYTHONPATH=src python scripts/run_galaxy_nautilus.py \
  --mode schema-check \
  --likelihood-mode mge \
  --halo-model sidm \
  --sidm-parameterization scale \
  --galaxy willman1
```

Short Nautilus smoke:

```bash
PYTHONPATH=src python scripts/run_galaxy_nautilus.py \
  --mode smoke \
  --likelihood-mode mge \
  --halo-model sidm \
  --sidm-parameterization scale \
  --galaxy willman1 \
  --output-name willman_1_sidm_signed_mge_smoke \
  --smoke-n-live 10 \
  --smoke-n-eff 2 \
  --smoke-n-like-max 120 \
  --smoke-timeout 180 \
  --f-live 1 \
  --n-batch 10
```

Longer run:

```bash
PYTHONPATH=src python scripts/run_galaxy_nautilus.py \
  --mode sample \
  --likelihood-mode mge \
  --halo-model sidm \
  --sidm-parameterization scale \
  --galaxy willman1 \
  --output-name willman_1_sidm_signed_mge \
  --n-live 120 \
  --n-eff 300
```

## Smoke-Test Record

The signed-MGE SIDM smoke test completed successfully on 2026-06-09 with the
physicality guard enabled:

```text
smoke_success = True
n_like = 50
n_live = 10
n_eff = 2.27483
log_z = -146.801
```

Output:

```text
outputs/willman_1_sidm_signed_mge_smoke_chain_smoke.csv
```

The chain uses the dynamic SIDM schema and includes resolved metadata such as
`halo_model`, `halo_parameterization`, `c200`, `halo_redshift`, and
`concentration_relation`. Derived values that do not apply to a chosen
parameterization are retained as nullable/NaN metadata and must be summarized
with finite-value checks.

## Tests

Focused package tests:

```bash
PYTHONPATH=src pytest -q tests
```

The current focused suite covers SIDM profile consistency, all parameterization
conversions, MGE density and enclosed-mass accuracy, signed descriptor
semantics, generalized Hernquist defaults, dynamic schema conversion, and
Nautilus adapter behavior.
