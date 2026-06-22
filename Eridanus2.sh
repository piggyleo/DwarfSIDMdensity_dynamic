#!/bin/bash
#SBATCH --job-name=img2imgTEST          # 作业名称
#SBATCH --qos=normal             # QoS 等级
#SBATCH --mem=20G                  # 内存
#SBATCH --cpus-per-task=24      # CPU 核心数
#SBATCH --time=2:00:00             # 运行时间
#SBATCH --output=job_%j.out        # 输出文件（%j 是作业ID）
#SBATCH --error=job_%j.err         # 错误文件

export PYTHONUNBUFFERED=1

python -u scripts/run_galaxy_nautilus.py \
  --mode sample \
  --galaxy eridanus2 \
  --likelihood-mode mge \
  --halo-model sidm \
  --sidm-parameterization m200-ludlow-scatter \
  --halo-redshift 0.000 \
  --n-live 1000 \
  --n-eff 10000 \
  --f-live 0.01 \
  --n-processes 24 \
  --n-batch 1000 \
  --n-like-max 200000 \
  --verbose \
  --progress-newlines \
  --seed 20260526 \
  --weighted-output \
  --use-sample-weights
  # --output-name eridanus_ii_weighted_sidm \
  # --no-mge-physicality-check \
