#!/bin/bash

set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 TASK_SPEC [TIME_LIMIT]" >&2
  exit 2
fi

task_spec="$1"
time_limit="${2:-04:00:00}"

if [[ ! "$task_spec" =~ ^[0-9,-]+$ ]]; then
  echo "TASK_SPEC must contain task numbers, commas, and ranges only: $task_spec" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"
mkdir -p logs outputs outputs/diagnostics outputs/figures

fit_submission="$(
  sbatch \
    --parsable \
    --wait \
    --array="${task_spec}%3" \
    --time="$time_limit" \
    --export=ALL,RESUME_EXISTING=1 \
    scripts/slurm_fit_all_sidm.sh
)"
fit_job_id="${fit_submission%%;*}"

echo "Completed resubmitted SIDM fit-and-plot tasks $task_spec: $fit_job_id"
echo "Resume existing checkpoints: enabled"
echo "Fit time limit: $time_limit"
