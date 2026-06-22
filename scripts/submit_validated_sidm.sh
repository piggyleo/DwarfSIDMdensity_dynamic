#!/bin/bash

set -euo pipefail

# Task IDs from config/sidm_galaxies.tsv for the galaxies that have passed
# the current SIDM fitting checks:
#   3  Canes Venatici I       4  Canes Venatici II
#   5  Coma Berenices         8  Eridanus II
#  12  Horologium I          13  Hydra II
#  16  Leo T                 17  Pisces II
#  19  Segue 1               21  Triangulum II
#  25  Ursa Major I          26  Ursa Major II
#  27  Willman 1
TASK_IDS=(3 4 5 8 12 13 16 17 19 21 25 26 27)
SERVER_MAX_CONCURRENT=2
MAX_CONCURRENT="${MAX_CONCURRENT:-$SERVER_MAX_CONCURRENT}"
SUBMIT_BATCH_SIZE="${SUBMIT_BATCH_SIZE:-$SERVER_MAX_CONCURRENT}"
HALO_MODEL="sidm"
SIDM_PARAMETERIZATION="m200-ludlow-scatter"

usage() {
  cat <<'EOF'
Usage: scripts/submit_validated_sidm.sh [options]

Options:
  --halo-model MODEL
      generalized-hernquist or sidm (default: sidm)
  --sidm-parameterization PARAMETERIZATION
      scale, m200-c200, m200-ludlow, or m200-ludlow-scatter
      (default: m200-ludlow-scatter)
  -h, --help
      Show this help message.
EOF
}

while (($# > 0)); do
  case "$1" in
    --halo-model)
      [[ $# -ge 2 ]] || {
        echo "--halo-model requires a value." >&2
        exit 2
      }
      HALO_MODEL="$2"
      shift 2
      ;;
    --sidm-parameterization)
      [[ $# -ge 2 ]] || {
        echo "--sidm-parameterization requires a value." >&2
        exit 2
      }
      SIDM_PARAMETERIZATION="$2"
      shift 2
      ;;
    -h | --help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

case "$HALO_MODEL" in
  generalized-hernquist | sidm) ;;
  *)
    echo "Unsupported --halo-model: $HALO_MODEL" >&2
    exit 2
    ;;
esac

case "$SIDM_PARAMETERIZATION" in
  scale | m200-c200 | m200-ludlow | m200-ludlow-scatter) ;;
  *)
    echo "Unsupported --sidm-parameterization: $SIDM_PARAMETERIZATION" >&2
    exit 2
    ;;
esac

if [[ ! "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]]; then
  echo "MAX_CONCURRENT must be a positive integer." >&2
  exit 2
fi
if [[ ! "$SUBMIT_BATCH_SIZE" =~ ^[1-9][0-9]*$ ]]; then
  echo "SUBMIT_BATCH_SIZE must be a positive integer." >&2
  exit 2
fi
if ((MAX_CONCURRENT > SERVER_MAX_CONCURRENT)); then
  echo "MAX_CONCURRENT cannot exceed the server limit of $SERVER_MAX_CONCURRENT." >&2
  exit 2
fi
if ((SUBMIT_BATCH_SIZE > SERVER_MAX_CONCURRENT)); then
  echo "SUBMIT_BATCH_SIZE cannot exceed the server limit of $SERVER_MAX_CONCURRENT." >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"
mkdir -p logs outputs outputs/diagnostics outputs/figures

echo "Submitting validated SIDM galaxy tasks: ${TASK_IDS[*]}"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"
echo "Tasks submitted per batch: $SUBMIT_BATCH_SIZE"
echo "Resume existing checkpoints: ${RESUME_EXISTING:-0}"
echo "Halo model: $HALO_MODEL"
echo "SIDM parameterization: $SIDM_PARAMETERIZATION"

for ((offset = 0; offset < ${#TASK_IDS[@]}; offset += SUBMIT_BATCH_SIZE)); do
  batch_ids=("${TASK_IDS[@]:offset:SUBMIT_BATCH_SIZE}")
  task_spec="$(IFS=,; echo "${batch_ids[*]}")"

  echo "Submitting SIDM batch: $task_spec"
  if submission="$(
    sbatch \
      --parsable \
      --wait \
      --array="${task_spec}%${MAX_CONCURRENT}" \
      --export=ALL,RESUME_EXISTING="${RESUME_EXISTING:-0}",HALO_MODEL="$HALO_MODEL",SIDM_PARAMETERIZATION="$SIDM_PARAMETERIZATION" \
      scripts/slurm_fit_all_sidm.sh
  )"; then
    job_id="${submission%%;*}"
    echo "Completed SIDM batch $task_spec: $job_id"
  else
    status=$?
    echo "SIDM batch $task_spec failed with sbatch status $status." >&2
    echo "No later batches were submitted." >&2
    exit "$status"
  fi
done

echo "All validated SIDM galaxy batches completed."
