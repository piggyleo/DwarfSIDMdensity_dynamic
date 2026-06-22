#!/bin/bash

python scripts/run_galaxy_nautilus.py \
  --mode corner \
  --galaxy willman1 \
  --halo-model sidm \
  --sidm-parameterization m200-ludlow-scatter \
  --halo-redshift 0.000 \
  --chain-output outputs/willman_1_weighted_mcrscatter_chain.csv \
  --figure-output outputs/figures/willman_1_weighted_mcrscatter_density_profile.png \
  --corner-output outputs/figures/willman_1_weighted_mcrscatter_corner.png \
  --use-sample-weights

python scripts/run_galaxy_nautilus.py \
  --mode plot \
  --galaxy willman1 \
  --halo-model sidm \
  --sidm-parameterization m200-ludlow-scatter \
  --halo-redshift 0.000 \
  --chain-output outputs/willman_1_weighted_mcrscatter_chain.csv \
  --figure-output outputs/figures/willman_1_weighted_mcrscatter_density_profile.png \
  --corner-output outputs/figures/willman_1_weighted_mcrscatter_corner.png \
  --use-sample-weights
