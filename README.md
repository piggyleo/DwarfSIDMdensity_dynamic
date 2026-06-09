# Hayashi Jeans Framework

This repository now includes a first-pass, pluggable implementation of the
Hayashi et al. (2023) unbinned LOS-velocity likelihood.

The default modeling sample intentionally uses only the low-risk galaxies listed
in `PROJECT_MEMORY.md` and selects `member_flag == 1` stars with finite velocities
and velocity errors. Galaxies with unresolved member-count or repeated-observation
口径 issues are left out of the default loader.

## Main Pieces

- `src/hayashi_jeans/data.py`: per-galaxy CSV loader, low-risk sample selector,
  and RA/Dec to projected pc coordinates.
- `src/hayashi_jeans/halos.py`: plug-in halo interface plus Hayashi generalized
  Hernquist/Zhao halo with spheroidal potential gradients.
- `src/hayashi_jeans/tracer.py`: flattened Plummer tracer.
- `src/hayashi_jeans/projection.py`: axisymmetric Jeans second moments and LOS
  projection for each star.
- `src/hayashi_jeans/likelihood.py`: Hayashi single-star Gaussian likelihood,
  including optional profiled systemic velocity.
- `src/hayashi_jeans/inference.py`: minimal `emcee` runner; alternate dark matter
  profiles plug in through a projector factory.

## Quick Checks

```bash
PYTHONPATH=src pytest -q
PYTHONPATH=src python scripts/run_low_risk_smoke.py --galaxy "Horologium I"
PYTHONPATH=src python scripts/make_willman1_test_posterior.py
PYTHONPATH=src python scripts/plot_willman1_density_profile.py --posterior-samples outputs/willman1_test_posterior_samples.csv
PYTHONPATH=src python scripts/run_willman1_restricted_mcmc.py
PYTHONPATH=src python scripts/plot_willman1_density_profile.py --posterior-samples outputs/willman1_restricted_posterior_samples.csv
```

The smoke script uses deliberately loose integration settings so it can confirm
the chain of data loading -> Hernquist halo -> Jeans projection -> likelihood.
For publication-level reproduction, tighten `--epsrel` and increase the
line-of-sight / vertical integration factors.
