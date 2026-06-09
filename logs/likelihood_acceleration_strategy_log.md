# Likelihood Acceleration Strategy Log

Created: 2026-05-19

Purpose: Track planned acceleration strategies for the Hayashi/Willman 1 Jeans likelihood, with separate tests for speedup and numerical error. Do not execute tests until the user gives an explicit instruction.

## Baseline Context

- Target bottleneck: one full likelihood evaluation for Willman 1 currently takes roughly 126 seconds in the full-parameter path.
- Scientific target: preserve the Hayashi et al. 2023 axisymmetric Jeans/LOS likelihood as the reference model.
- Current slow reference implementation:
  - `src/hayashi_jeans/projection.py`
  - nested LOS integration, Jeans vertical integration, finite-difference radial derivative, and halo force integration.
- Primary test galaxy:
  - `Willman 1`
  - galaxy data: `data/galaxies/27_Willman_1.csv`
  - structural center: `data/processed/galaxy_structural_centers.csv`

## Execution Status

Status: Strategies 1-3 testing completed on Willman 1 with standalone experimental code; Strategy 4 remains next in the planned first-priority block. Strategies 5-6 are not yet formally tested in this log.

## Test Metrics To Record

For each strategy, record:

- date/time
- code branch/file(s) changed
- parameter vector(s) tested
- number of stars used
- integration/grid settings
- wall time per likelihood evaluation
- speedup relative to baseline
- per-star `sigma_los2` error versus slow reference:
  - max relative error
  - median relative error
  - worst affected star id/location if available
- log-likelihood error versus slow reference:
  - absolute `Delta lnL`
  - relative/context note
- posterior-facing impact, if tested:
  - median/68% changes for key parameters
  - density profile changes at selected radii, especially 10 pc, `b_*`, 150 pc, 1 kpc
- conclusion:
  - accept/reject/needs refinement

## Planned Strategy Order

### 1. R-z Jeans Moment Grid + Interpolation

Priority: highest.

Idea:
- For a fixed parameter vector, build a 2D grid in intrinsic cylindrical coordinates `(R,z)`.
- Precompute:
  - tracer density `nu(R,z)`
  - Jeans vertical pressure `P_z = nu * v_z2`
  - `v_R2`, `v_phi2`, `v_z2`
- During LOS projection, interpolate moments instead of solving Jeans moments at each LOS point.

Targeted tests:
- Compare fast grid moments against current slow `jeans_moments()` at random `(R,z)` points.
- Compare `sigma_los2_many()` for all 40 Willman 1 stars.
- Test multiple grid resolutions and spacing choices, especially inner-region resolution below 10 pc.

Expected risk:
- Grid/interpolation may smooth the central cusp and bias inner density constraints.

Status: Tested on Willman 1 with standalone experimental code; production code unchanged. Strategy is validated as promising. Cubic-spline radial derivatives for `dP_z/dR` strongly improve accuracy over direct `np.gradient` on the nonuniform R grid.

Results:
- Slow reference for all 40 Willman 1 members took `430.39 s` for `sigma_los2_many`.
- All tested R-z grids covered the full member-star LOS path; highest-resolution test used `R_max=603.85 pc`, `|z|_max=604.80 pc`, while the LOS path reached `R_max=559.12 pc`, `|z|_max=210.83 pc`.
- Test Run 001 (`32 x 64`, `n_los=80`) was extremely fast (`0.94 s`, `458x` including build) but too inaccurate for production use: median `sigma_los2` relative error `9.73%`, max `46.92%`.
- Test Run 002 (`64 x 128`, `n_los=120`) improved accuracy but still had noticeable outliers: total `3.69 s`, median error `2.23%`, max `10.30%`, `Delta lnL=+0.0717`.
- Test Run 003 (`96 x 192`, `n_los=160`, gradient `dP_z/dR`) is a usable exploratory baseline: total `8.44 s`, median error `0.81%`, max `4.37%`, p95 `3.80%`, `Delta lnL=+0.0335`.
- Test Run 004 (`128 x 256`, `n_los=200`, gradient `dP_z/dR`) improves accuracy: total `15.19 s`, median error `0.41%`, max `2.38%`, p95 `2.09%`, `Delta lnL=+0.0179`, speedup `28.3x` including grid build.
- Test Run 005 (`96 x 192`, `n_los=160`, cubic-spline `dP_z/dR`) shows that the radial derivative was the main accuracy limiter: total `8.01 s`, median error `0.078%`, max `0.812%`, p95 `0.551%`, `Delta lnL=+0.00185`.
- Test Run 006 (`128 x 256`, `n_los=200`, cubic-spline `dP_z/dR`) is the best current configuration: total `14.40 s`, median error `0.034%`, max `0.302%`, p95 `0.194%`, `Delta lnL=+0.000115`, speedup about `29.9x` including grid build.
- Conclusion: R-z moment grid + LOS interpolation can reduce one Willman 1 likelihood-like evaluation from several minutes to seconds while keeping `lnL` error negligible when `dP_z/dR` is computed with per-z-slice cubic-spline derivatives. The `96 x 192` cubic-spline grid is likely sufficient for exploratory MCMC; the `128 x 256` cubic-spline grid is preferred for validation/final checks.

Derivative-improvement summary:
- Motivation: the original experimental grid used `np.gradient(P_z, R)` on a nonuniform R grid. Since `dP_z/dR` enters `v_phi2` directly, this derivative was suspected to dominate the remaining interpolation error.
- Tested change: compute `dP_z/dR` by fitting `CubicSpline(P_z(R,z_k))` independently at each `z_k` slice and evaluating the spline derivative on the R grid.
- `96 x 192`, `n_los=160`: median `sigma_los2` error improved from `0.81%` to `0.078%`; max error from `4.37%` to `0.812%`; `Delta lnL` from `+0.0335` to `+0.00185`.
- `128 x 256`, `n_los=200`: median `sigma_los2` error improved from `0.41%` to `0.034%`; max error from `2.38%` to `0.302%`; `Delta lnL` from `+0.0179` to `+0.000115`.
- Conclusion: cubic-spline radial derivatives substantially reduce both per-star and likelihood-level errors at essentially no extra runtime cost. This derivative method should be used in any production implementation of the R-z grid projector, subject to additional convergence checks across more parameter vectors.

### 2. Use rho0 Linear Scaling

Priority: high.

Idea:
- For fixed halo shape/geometric parameters, compute `sigma_los2` once at `rho0=1`.
- For arbitrary density normalization:
  `sigma_los2(rho0) = rho0 * sigma_los2(rho0=1)`.

Targeted tests:
- Numerically verify against slow reference at several `rho0` values.
- Confirm `Delta lnL` is negligible compared with direct recomputation.
- Measure reduction in per-likelihood cost when only `rho0` changes.

Expected risk:
- Very low under Hayashi assumptions: dark-matter-dominated system, linear density normalization, no stellar self-gravity.

Status: Isolated test completed on top of the experimental R-z grid + cubic `dP_z/dR`; production code unchanged. Strategy is numerically validated.

Results:
- Test Run 001 used `96 x 192`, `n_los=160`, cubic-spline `dP_z/dR`, with `log10_rho0_msun_pc3 = -2.5, -1.961219, -1.5, -0.5`.
- Direct comparison rebuilt the R-z grid separately for each rho0 value; scaled comparison built one `rho0=1` grid and used `sigma_los2(rho0)=10**log10_rho0 * sigma_los2_unit`.
- Per-star `sigma_los2` differences were at floating-point precision: overall median relative error `3.51e-16`, max `1.64e-15`, p95 `9.62e-16`.
- Likelihood differences were negligible: maximum absolute `Delta lnL = 2.84e-14`.
- For a four-rho sweep, direct rebuild time was `53.80 s`; unit-grid scaling time was `13.20 s`, giving `4.07x` speedup. For larger rho grids or repeated rho updates, the speedup approaches the number of avoided grid rebuilds.
- Conclusion: rho0 linear scaling is exact to numerical precision under the current Hayashi assumptions and should be used in production. It is especially useful for profiling/marginalizing rho0 or block-updating rho0 without recomputing Jeans moments.

### 3. Halo Force Grid + Interpolation

Priority: high.

Idea:
- Precompute halo potential gradients:
  - `dPhi/dR(R,z)`
  - `dPhi/dz(R,z)`
- Interpolate force values during Jeans moment construction.

Targeted tests:
- Compare force-grid interpolation against current `GeneralizedHernquistHalo.potential_gradients()`.
- Evaluate how force-grid errors propagate into moments, `sigma_los2`, and `lnL`.
- Test grid choices compatible with Strategy 1.

Expected risk:
- Force errors propagate directly into Jeans moments; inner grid resolution is critical.

Status: Isolated test completed on top of the experimental R-z moment grid + cubic `dP_z/dR`; production code unchanged. Strategy is useful as a secondary acceleration layer, but the force-grid resolution must not be too coarse.

Results:
- Formal Test Run 001 used Willman 1, moment grid `96 x 192`, `n_los=160`, and compared against the accepted direct-force R-z grid reference at the same moment resolution. The direct-force R-z reference build took `28.26 s`.
- `force32x64` gave the strongest speedup: total build `3.19 s`, build speedup `8.85x`, but accuracy was marginal for production-facing use: median/p95/max `sigma_los2` relative error `0.60%` / `1.77%` / `3.69%`, `Delta lnL=+0.0318`.
- `force48x96` is the best current tradeoff: total build `7.80 s`, build speedup `3.62x`, median/p95/max `sigma_los2` relative error `0.16%` / `0.37%` / `0.51%`, force-vector p95 relative error `0.29%`, `Delta lnL=+0.00286`.
- `force64x128` further reduces force interpolation error but gives little downstream improvement over `48x96`: total build `12.61 s`, build speedup `2.24x`, median/p95/max `sigma_los2` relative error `0.11%` / `0.34%` / `0.41%`, `Delta lnL=-0.00409`.
- `force96x192`, matching the moment grid resolution, is a consistency check rather than an acceleration: total build `28.82 s`, speedup `0.98x`, and the projected dispersions match the direct-force R-z reference to machine precision because both evaluate the force on the same grid nodes.
- Conclusion: adopt `force48x96` as the exploratory default if this strategy is connected to the MCMC likelihood; use `force64x128` or direct force for validation runs. Reject very coarse force grids unless a later posterior test shows the `lnL` perturbation is harmless.

### 4. Profile or Analytically Marginalize Systemic Velocity <u>

Priority: high.

Idea:
- Avoid explicit sampling of `systemic_velocity_kms`.
- For a given set of model variances:
  `s_i^2 = sigma_los_i^2 + delta_v_i^2`,
  profile:
  `<u> = sum(u_i/s_i^2) / sum(1/s_i^2)`.
- Prefer analytic marginalization over pure profiling for better uncertainty treatment:
  `ln L_marg = ln L(<u>_profiled) + 0.5 * ln(2*pi / sum(1/s_i^2))`
  under a flat prior on `<u>`.

Targeted tests:
- Compare explicit sampling/profile/marginalized likelihood on the restricted model.
- Measure posterior changes for `rho0`.
- Confirm dimension reduction and runtime impact.

Expected risk:
- Profiling may understate uncertainty in low-sample systems.
- Analytic marginalization is preferred for final use.

Status: Not started as an isolated test.

Results: TBD.

### 5. Parallelize Walkers

Priority: medium.

Idea:
- Use multiprocessing or another pool with `emcee.EnsembleSampler(..., pool=pool)`.
- This does not shorten a single likelihood call, but reduces wall-clock time per ensemble step.

Targeted tests:
- Compare serial vs parallel wall time for a fixed number of likelihood calls.
- Check reproducibility/seed behavior.
- Prefer running as a `.py` script, because Jupyter multiprocessing can be fragile.

Expected risk:
- No mathematical bias.
- Operational risk from pickling/import issues and memory overhead.

Status: Not started.

Results: TBD.

### 6. Use AGAMA For Potential/Interpolation Acceleration

Priority: exploratory.

Idea:
- Investigate using AGAMA for potential/force representation or interpolation.
- Keep Hayashi Jeans/projection likelihood if possible; only replace low-level force evaluation.

Targeted tests:
- Verify AGAMA can represent the generalized Hernquist/Zhao axisymmetric halo with required accuracy.
- Compare AGAMA forces against current reference force implementation.
- Measure downstream `sigma_los2` and `lnL` differences.

Expected risk:
- Depends on how closely AGAMA represents the exact Hayashi halo.
- More dependency/setup complexity.

Status: Not started.

Results: TBD.

## Notes

- Strategies 1-4 are the first block to try.
- Strategy 5 is useful after a faster likelihood exists, or for batch validation.
- Strategy 6 should be treated as exploratory and only adopted after force-level validation.



### Strategy 1 Test Run 001

Status: completed experimental standalone test; production code unchanged.

- moment comparison CSV: `outputs/diagnostics/rz_grid_willman1_moment_comparison.csv`
- sigma_los2 comparison CSV: `outputs/diagnostics/rz_grid_willman1_sigma_los2_comparison.csv`
- summary CSV: `outputs/diagnostics/rz_grid_willman1_summary.csv`
- n_star: 40
- n_r: 32
- n_z: 64
- n_los: 80
- grid_r_max_pc: 603.6317267802842
- grid_z_max_pc: 604.8000000000001
- los_path_r_max_pc: 558.9182655373002
- los_path_z_abs_max_pc: 210.7802051648668
- fast_build_time_s: 0.9251152500073658
- slow_sigma_time_s: 430.3948537920078
- fast_sigma_time_s: 0.013730250007938594
- fast_total_time_s: 0.9388455000153044
- speedup_excluding_build: 31346.46882199238
- speedup_including_build: 458.4299054370413
- sigma_rel_error_median: 0.09734869400749621
- sigma_rel_error_max: 0.46920550864580757
- sigma_rel_error_p95: 0.4003698435530401
- moment_rel_error_median_all: 0.07155741653461817
- moment_rel_error_max_all: 0.4868485040642426
- slow_ln_l: -129.83601965984826
- fast_ln_l: -129.75295854113742
- delta_ln_l: 0.08306111871084454


### Strategy 1 Test Run 002

Status: completed fine-grid follow-up using saved slow sigma reference; production code unchanged.

- sigma_los2 comparison CSV: `outputs/diagnostics/rz_grid_willman1_sigma_los2_comparison_grid64x128_los120.csv`
- summary CSV: `outputs/diagnostics/rz_grid_willman1_summary_grid64x128_los120.csv`
- label: grid64x128_los120
- n_star: 40
- n_r: 64
- n_z: 128
- n_los: 120
- grid_r_max_pc: 603.7758971639918
- grid_z_max_pc: 604.8000000000001
- los_path_r_max_pc: 559.0517566333257
- los_path_z_abs_max_pc: 210.81601124861015
- fast_build_time_s: 3.681334583001444
- fast_sigma_time_s: 0.0047843329957686365
- fast_total_time_s: 3.6861189159972128
- sigma_rel_error_median: 0.022272481995873034
- sigma_rel_error_max: 0.10298534610128672
- sigma_rel_error_p95: 0.09009388502018104
- slow_ln_l: -129.83601965984826
- fast_ln_l: -129.76430216273255
- delta_ln_l: 0.07171749711571351


### Strategy 1 Test Run 003

Status: completed finer-grid follow-up using saved slow sigma reference; production code unchanged.

- sigma_los2 comparison CSV: `outputs/diagnostics/rz_grid_willman1_sigma_los2_comparison_grid96x192_los160.csv`
- summary CSV: `outputs/diagnostics/rz_grid_willman1_summary_grid96x192_los160.csv`
- label: grid96x192_los160
- n_star: 40
- n_r: 96
- n_z: 192
- n_los: 160
- grid_r_max_pc: 603.8266028314492
- grid_z_max_pc: 604.8000000000001
- los_path_r_max_pc: 559.098706325416
- los_path_z_abs_max_pc: 210.82860447835674
- fast_build_time_s: 8.433832082999288
- fast_sigma_time_s: 0.007725834002485499
- fast_total_time_s: 8.441557917001774
- sigma_rel_error_median: 0.008075233565221045
- sigma_rel_error_max: 0.04369806975790097
- sigma_rel_error_p95: 0.03800814159173505
- slow_ln_l: -129.83601965984826
- fast_ln_l: -129.80248356139765
- delta_ln_l: 0.03353609845061101


### Strategy 1 Test Run 004

Status: completed higher-resolution grid follow-up using saved slow sigma reference; production code unchanged.

- sigma_los2 comparison CSV: `outputs/diagnostics/rz_grid_willman1_sigma_los2_comparison_grid128x256_los200.csv`
- summary CSV: `outputs/diagnostics/rz_grid_willman1_summary_grid128x256_los200.csv`
- label: grid128x256_los200
- n_star: 40
- n_r: 128
- n_z: 256
- n_los: 200
- grid_r_max_pc: 603.8501323792941
- grid_z_max_pc: 604.8000000000001
- los_path_r_max_pc: 559.1204929437909
- los_path_z_abs_max_pc: 210.8344482622719
- fast_build_time_s: 15.179379957990022
- fast_sigma_time_s: 0.007427375006955117
- fast_total_time_s: 15.186807332996977
- speedup_including_build_vs_slow430s: 28.340048329767892
- sigma_rel_error_median: 0.004120203555954686
- sigma_rel_error_max: 0.023794735021041825
- sigma_rel_error_p95: 0.02092831149809898
- slow_ln_l: -129.83601965984826
- fast_ln_l: -129.81807694270765
- delta_ln_l: 0.01794271714061324


### Strategy 1 Test Run 005

Status: tested cubic-spline radial derivative for dP_z/dR using saved slow sigma reference; production code unchanged.

Improvement idea: replace direct `np.gradient(P_z, R)` on the nonuniform R grid with per-z-slice `CubicSpline(P_z(R))` derivatives.

- sigma_los2 comparison CSV: `outputs/diagnostics/rz_grid_willman1_sigma_los2_comparison_grid96x192_los160_cubic_dPdr.csv`
- summary CSV: `outputs/diagnostics/rz_grid_willman1_summary_grid96x192_los160_cubic_dPdr.csv`
- label: grid96x192_los160_cubic_dPdr
- n_star: 40
- n_r: 96
- n_z: 192
- n_los: 160
- radial_derivative: cubic_spline
- grid_r_max_pc: 603.8266028314492
- grid_z_max_pc: 604.8000000000001
- los_path_r_max_pc: 559.098706325416
- los_path_z_abs_max_pc: 210.82860447835674
- fast_build_time_s: 8.008999750003568
- fast_sigma_time_s: 0.005717084000934847
- fast_total_time_s: 8.014716834004503
- sigma_rel_error_median: 0.0007810965916943743
- sigma_rel_error_max: 0.00812244116607262
- sigma_rel_error_p95: 0.005507411547689685
- slow_ln_l: -129.83601965984826
- fast_ln_l: -129.8341717103089
- delta_ln_l: 0.0018479495393535217


### Strategy 1 Test Run 006

Status: tested high-resolution cubic-spline radial derivative for dP_z/dR using saved slow sigma reference; production code unchanged.

Improvement idea: same as Run 005, with denser R-z grid.

- sigma_los2 comparison CSV: `outputs/diagnostics/rz_grid_willman1_sigma_los2_comparison_grid128x256_los200_cubic_dPdr.csv`
- summary CSV: `outputs/diagnostics/rz_grid_willman1_summary_grid128x256_los200_cubic_dPdr.csv`
- label: grid128x256_los200_cubic_dPdr
- n_star: 40
- n_r: 128
- n_z: 256
- n_los: 200
- radial_derivative: cubic_spline
- grid_r_max_pc: 603.8501323792941
- grid_z_max_pc: 604.8000000000001
- los_path_r_max_pc: 559.1204929437909
- los_path_z_abs_max_pc: 210.8344482622719
- fast_build_time_s: 14.393573584005935
- fast_sigma_time_s: 0.006334165998850949
- fast_total_time_s: 14.399907750004786
- sigma_rel_error_median: 0.0003392181844923543
- sigma_rel_error_max: 0.003021514834402893
- sigma_rel_error_p95: 0.0019407790447803698
- slow_ln_l: -129.83601965984826
- fast_ln_l: -129.83590447532373
- delta_ln_l: 0.00011518452453174177


### Strategy 2 Test Run 001

Status: completed isolated rho0 linear-scaling test using R-z grid + cubic dP_z/dR; production code unchanged.

- Method: build one `rho0=1` R-z grid, scale `sigma_los2_unit` by `10**log10_rho0`, and compare to direct grid rebuilds at each rho0.
- detail CSV: `outputs/diagnostics/rho0_scaling_willman1_detail.csv`
- by-rho CSV: `outputs/diagnostics/rho0_scaling_willman1_by_rho.csv`
- summary CSV: `outputs/diagnostics/rho0_scaling_willman1_summary.csv`
- n_star: 40
- n_r: 96
- n_z: 192
- n_los: 160
- n_rho_values: 4
- unit_build_time_s: 13.195415750000393
- unit_sigma_time_s: 0.008363916000234894
- direct_total_time_s: 53.80429691700556
- scaled_total_time_s: 13.20378870901186
- speedup_for_rho_sweep: 4.07491350420376
- overall_sigma_rel_error_median: 3.5105642924613e-16
- overall_sigma_rel_error_max: 1.6350542349384328e-15
- overall_sigma_rel_error_p95: 9.621133872659828e-16
- max_abs_delta_ln_l: 2.842170943040401e-14


### Strategy 3 Smoke Run 000

Status: completed isolated halo force-grid interpolation test using the experimental R-z moment grid + cubic `dP_z/dR`; production code unchanged.

- Method: build the accepted R-z moment grid, but replace direct halo force calls on every moment-grid point with a coarser halo force grid and linear interpolation.
- sigma_los2 detail CSV: `outputs/diagnostics/halo_force_smoke/halo_force_grid_willman1_sigma_detail.csv`
- force detail CSV: `outputs/diagnostics/halo_force_smoke/halo_force_grid_willman1_force_detail.csv`
- summary CSV: `outputs/diagnostics/halo_force_smoke/halo_force_grid_willman1_summary.csv`
- force16x32: total build `0.847 s`, build speedup `3.889x`, sigma median/max/p95 rel error `0.0269015` / `0.0643738` / `0.0608249`, force-vector median/max/p95 rel error `0.0119838` / `0.0340339` / `0.030101`, `Delta lnL=0.110724`.
- force24x48: total build `1.868 s`, build speedup `1.764x`, sigma median/max/p95 rel error `0.00653799` / `0.026579` / `0.0206009`, force-vector median/max/p95 rel error `0.00343838` / `0.0124627` / `0.0107236`, `Delta lnL=0.0112254`.
- Best p95 sigma-error setting in this run: `force24x48` with p95 sigma relative error `0.0206009` and `Delta lnL=0.0112254`.


### Strategy 3 Test Run 001

Status: completed isolated halo force-grid interpolation test using the experimental R-z moment grid + cubic `dP_z/dR`; production code unchanged.

- Method: build the accepted R-z moment grid, but replace direct halo force calls on every moment-grid point with a coarser halo force grid and linear interpolation.
- sigma_los2 detail CSV: `outputs/diagnostics/halo_force_grid_willman1_sigma_detail.csv`
- force detail CSV: `outputs/diagnostics/halo_force_grid_willman1_force_detail.csv`
- summary CSV: `outputs/diagnostics/halo_force_grid_willman1_summary.csv`
- force32x64: total build `3.192 s`, build speedup `8.854x`, sigma median/max/p95 rel error `0.00598171` / `0.0368673` / `0.0177315`, force-vector median/max/p95 rel error `0.00279153` / `0.00869608` / `0.00674072`, `Delta lnL=0.031753`.
- force48x96: total build `7.798 s`, build speedup `3.624x`, sigma median/max/p95 rel error `0.00157252` / `0.0050552` / `0.00366958`, force-vector median/max/p95 rel error `0.00120172` / `0.00358535` / `0.00288414`, `Delta lnL=0.00285867`.
- force64x128: total build `12.614 s`, build speedup `2.241x`, sigma median/max/p95 rel error `0.00105069` / `0.00410139` / `0.00344813`, force-vector median/max/p95 rel error `0.000597088` / `0.00196065` / `0.00158553`, `Delta lnL=-0.00409013`.
- force96x192: total build `28.817 s`, build speedup `0.981x`, sigma median/max/p95 rel error `0` / `0` / `0`, force-vector median/max/p95 rel error `0.00027373` / `0.000868998` / `0.0007452`, `Delta lnL=0`.
- Best p95 sigma-error setting in this run: `force96x192` with p95 sigma relative error `0` and `Delta lnL=0`.


### Three-Layer Cache Smoke Run 000

Status: completed standalone test of layer-1 halo force cache, layer-2 R-z moment/projection cache, and layer-3 rho0/systemic likelihood update; production code unchanged.

- Compared against the original slow `src/projection.py` projector for five different layer-1 halo-shape parameter vectors.
- sigma detail CSV: `outputs/diagnostics/three_layer_cache_smoke/three_layer_cache_willman1_sigma_detail.csv`
- force detail CSV: `outputs/diagnostics/three_layer_cache_smoke/three_layer_cache_willman1_force_detail.csv`
- moment detail CSV: `outputs/diagnostics/three_layer_cache_smoke/three_layer_cache_willman1_moment_detail.csv`
- summary CSV: `outputs/diagnostics/three_layer_cache_smoke/three_layer_cache_willman1_summary.csv`
- Speedup versus slow total likelihood time: median `1.57e+03x`, min `1.57e+03x`, max `1.57e+03x`.
- Per-star sigma_los2 p95 relative error across layer-1 groups: median `0.0697903`, max `0.0697903`.
- Absolute Delta lnL across layer-1 groups: median `0.0561857`, max `0.0561857`.
- L1_fiducial: cached total `0.494 s`, slow total `776.914 s`, speedup `1.57e+03x`, sigma p95/max rel error `0.0697903` / `0.0846044`, moment p95/max rel error `0.088096` / `0.088096`, force-vector p95/max rel error `0.0239221` / `0.0251924`, `Delta lnL=0.0561857`.


### Three-Layer Cache Test Run 001

Status: completed standalone test of layer-1 halo force cache, layer-2 R-z moment/projection cache, and layer-3 rho0/systemic likelihood update; production code unchanged.

- Compared against the original slow `src/projection.py` projector for five different layer-1 halo-shape parameter vectors.
- sigma detail CSV: `outputs/diagnostics/three_layer_cache/three_layer_cache_willman1_sigma_detail.csv`
- force detail CSV: `outputs/diagnostics/three_layer_cache/three_layer_cache_willman1_force_detail.csv`
- moment detail CSV: `outputs/diagnostics/three_layer_cache/three_layer_cache_willman1_moment_detail.csv`
- summary CSV: `outputs/diagnostics/three_layer_cache/three_layer_cache_willman1_summary.csv`
- Speedup versus slow total likelihood time: median `209x`, min `181x`, max `259x`.
- Per-star sigma_los2 p95 relative error across layer-1 groups: median `0.00520151`, max `0.0111798`.
- Absolute Delta lnL across layer-1 groups: median `0.00946769`, max `0.0472559`.
- L1_fiducial: cached total `3.845 s`, slow total `776.914 s`, speedup `202x`, sigma p95/max rel error `0.00439137` / `0.00809702`, moment p95/max rel error `0.0820546` / `0.0846505`, force-vector p95/max rel error `0.00308095` / `0.0035057`, `Delta lnL=0.00470662`.
- L1_oblate_compact: cached total `4.180 s`, slow total `754.870 s`, speedup `181x`, sigma p95/max rel error `0.00520151` / `0.00758172`, moment p95/max rel error `0.220401` / `0.244618`, force-vector p95/max rel error `0.00291038` / `0.0063167`, `Delta lnL=-0.0163395`.
- L1_prolate_extended: cached total `3.680 s`, slow total `851.549 s`, speedup `231x`, sigma p95/max rel error `0.00865795` / `0.0125903`, moment p95/max rel error `0.240788` / `0.277524`, force-vector p95/max rel error `0.00488734` / `0.00578385`, `Delta lnL=-0.00946769`.
- L1_cored_steep_outer: cached total `2.765 s`, slow total `578.799 s`, speedup `209x`, sigma p95/max rel error `0.00361086` / `0.00691676`, moment p95/max rel error `0.286538` / `0.329343`, force-vector p95/max rel error `0.00798464` / `0.0118146`, `Delta lnL=-0.00245688`.
- L1_cuspy_shallow_outer: cached total `2.136 s`, slow total `552.448 s`, speedup `259x`, sigma p95/max rel error `0.0111798` / `0.014561`, moment p95/max rel error `0.0465208` / `0.0488333`, force-vector p95/max rel error `0.00959524` / `0.0109516`, `Delta lnL=-0.0472559`.
