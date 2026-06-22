#!/bin/bash
#SBATCH --job-name=sidm_plot
#SBATCH --qos=normal
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --array=1-27%2
#SBATCH --output=logs/sidm_plot_%A_%a.out
#SBATCH --error=logs/sidm_plot_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
CONFIG_FILE="$PROJECT_ROOT/config/sidm_galaxies.tsv"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"
HALO_MODEL="${HALO_MODEL:-sidm}"
SIDM_PARAMETERIZATION="${SIDM_PARAMETERIZATION:-m200-ludlow-scatter}"
halo_output_tag="${HALO_MODEL//-/_}"
if [[ "$HALO_MODEL" == "sidm" && "$SIDM_PARAMETERIZATION" != "m200-ludlow-scatter" ]]; then
  halo_output_tag+="_${SIDM_PARAMETERIZATION//-/_}"
fi

cd "$PROJECT_ROOT"

row="$(awk -F $'\t' -v task_id="$TASK_ID" 'NR > 1 && $1 == task_id {print; exit}' "$CONFIG_FILE")"
if [[ -z "$row" ]]; then
  echo "No galaxy configuration found for array task $TASK_ID" >&2
  exit 2
fi

IFS=$'\t' read -r configured_task_id galaxy_input galaxy_slug redshift <<< "$row"
if [[ "$configured_task_id" != "$TASK_ID" ]]; then
  echo "Configuration task mismatch: expected $TASK_ID, got $configured_task_id" >&2
  exit 2
fi

chain_output="outputs/${galaxy_slug}_${halo_output_tag}_chain.csv"
figure_output="outputs/figures/${galaxy_slug}_${halo_output_tag}_density_profile.png"
corner_output="outputs/figures/${galaxy_slug}_${halo_output_tag}_corner.png"

mkdir -p logs outputs/figures

export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

echo "array_job_id=${SLURM_ARRAY_JOB_ID:-unknown}"
echo "array_task_id=$TASK_ID"
echo "galaxy=$galaxy_input"
echo "galaxy_slug=$galaxy_slug"
echo "redshift=$redshift"
echo "halo_model=$HALO_MODEL"
echo "sidm_parameterization=$SIDM_PARAMETERIZATION"
echo "chain_output=$chain_output"
echo "figure_output=$figure_output"
echo "corner_output=$corner_output"

if [[ ! -s "$chain_output" ]]; then
  echo "Weighted chain is missing or empty: $chain_output" >&2
  exit 3
fi

python -u scripts/run_galaxy_nautilus.py \
  --mode corner \
  --galaxy "$galaxy_input" \
  --halo-model "$HALO_MODEL" \
  --sidm-parameterization "$SIDM_PARAMETERIZATION" \
  --halo-redshift "$redshift" \
  --chain-output "$chain_output" \
  --figure-output "$figure_output" \
  --corner-output "$corner_output" \
  --use-sample-weights

python -u scripts/run_galaxy_nautilus.py \
  --mode plot \
  --galaxy "$galaxy_input" \
  --halo-model "$HALO_MODEL" \
  --sidm-parameterization "$SIDM_PARAMETERIZATION" \
  --halo-redshift "$redshift" \
  --chain-output "$chain_output" \
  --figure-output "$figure_output" \
  --corner-output "$corner_output" \
  --use-sample-weights
