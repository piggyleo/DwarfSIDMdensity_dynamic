#!/bin/bash

python scripts/run_galaxy_nautilus.py \
  --mode corner \
  --galaxy eridanus2 \
  --halo-model sidm \
  --sidm-parameterization m200-ludlow-scatter \
  --halo-redshift 0.0 \
  --chain-output outputs/eridanus_2_sidm_chain.csv \
  --figure-output outputs/figures/eridanus_2_sidm_density_profile.png \
  --corner-output outputs/figures/eridanus_2_sidm_corner.png \
  --use-sample-weights

python scripts/run_galaxy_nautilus.py \
  --mode plot \
  --galaxy eridanus2 \
  --halo-model sidm \
  --sidm-parameterization m200-ludlow-scatter \
  --halo-redshift 0.0 \
  --chain-output outputs/eridanus_2_sidm_chain.csv \
  --figure-output outputs/figures/eridanus_2_sidm_density_profile.png \
  --corner-output outputs/figures/eridanus_2_sidm_corner.png \
  --use-sample-weights
