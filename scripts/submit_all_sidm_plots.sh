#!/bin/bash

set -euo pipefail

HALO_MODEL="sidm"
SIDM_PARAMETERIZATION="m200-ludlow-scatter"
POSITIONAL_ARGS=()

usage() {
  cat <<'EOF'
Usage: scripts/submit_all_sidm_plots.sh [options] [START_TASK [END_TASK]]

Submit plotting-only Slurm array jobs for already completed SIDM chains.

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
    --)
      shift
      while (($# > 0)); do
        POSITIONAL_ARGS+=("$1")
        shift
      done
      ;;
    -*)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done

if [[ ${#POSITIONAL_ARGS[@]} -gt 2 ]]; then
  usage >&2
  exit 2
fi

start_task="${POSITIONAL_ARGS[0]:-1}"
end_task="${POSITIONAL_ARGS[1]:-27}"

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

if [[ ! "$start_task" =~ ^[0-9]+$ || ! "$end_task" =~ ^[0-9]+$ ]]; then
  echo "START_TASK and END_TASK must be integers." >&2
  exit 2
fi
if ((start_task < 1 || start_task > 27 || end_task < 1 || end_task > 27)); then
  echo "START_TASK and END_TASK must be between 1 and 27." >&2
  exit 2
fi
if ((start_task > end_task)); then
  echo "START_TASK must not be greater than END_TASK." >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"
mkdir -p logs outputs/figures

echo "Halo model: $HALO_MODEL"
echo "SIDM parameterization: $SIDM_PARAMETERIZATION"

for batch_start in $(seq "$start_task" 2 "$end_task"); do
  batch_end=$((batch_start + 1))
  if ((batch_end > end_task)); then
    batch_end=$end_task
  fi
  if ((batch_start == batch_end)); then
    task_spec="${batch_start}%2"
  else
    task_spec="${batch_start}-${batch_end}%2"
  fi

  echo "Submitting SIDM plotting batch: $task_spec"
  if submission="$(
    sbatch \
      --parsable \
      --wait \
      --array="$task_spec" \
      --export=ALL,HALO_MODEL="$HALO_MODEL",SIDM_PARAMETERIZATION="$SIDM_PARAMETERIZATION" \
      scripts/slurm_plot_all_sidm.sh
  )"; then
    job_id="${submission%%;*}"
    echo "Completed SIDM plotting batch $batch_start-$batch_end: $job_id"
  else
    status=$?
    echo "SIDM plotting batch $batch_start-$batch_end failed with sbatch status $status." >&2
    echo "No later plotting batches were submitted." >&2
    exit "$status"
  fi
done

echo "All SIDM plotting batches completed for tasks $start_task-$end_task."
