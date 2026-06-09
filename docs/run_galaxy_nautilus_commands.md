# `run_galaxy_nautilus.py` 常用指令

本文记录 `scripts/run_galaxy_nautilus.py` 的常用调用方式。该脚本支持用 `--galaxy` 自动选择星系数据、结构中心、默认输出文件名，并默认使用自适应系统速度先验。

当前 likelihood mode 包括：

```text
fast
mge
validation-halo
validation-strict
validation-convergence
```

其中 `mge` 使用共享模块 `hayashi_jeans.mge` 的 MGE/JAM 计算；`fast` 和 `validation-*` 仍使用原 Jeans/projector 路径，不额外运行 MGE guard。

## 1. 单点 log-probability 检查

用于确认星系名解析、数据读取、结构中心、prior 和 likelihood 都能正常工作。

```bash
python scripts/run_galaxy_nautilus.py \
  --mode point \
  --galaxy willman1 \
  --likelihood-mode fast
```

例如测量 Eridanus II：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode point \
  --galaxy eridanus2 \
  --likelihood-mode fast
```

## 2. Schema 检查

用于确认输出 chain 的列名可以接入现有 density profile 和 corner 诊断流程。

```bash
python scripts/run_galaxy_nautilus.py \
  --mode schema-check \
  --galaxy eridanus2
```

## 3. 快速 smoke test

只用于测试 nautilus 调用、CSV 写出和输出 schema，不用于科学结果。

```bash
python scripts/run_galaxy_nautilus.py \
  --mode smoke \
  --galaxy eridanus2 \
  --likelihood-mode fast \
  --n-r 24 --n-z 48 --n-los 48 \
  --n-force-r 16 --n-force-z 32 \
  --smoke-n-live 30 \
  --smoke-n-eff 20 \
  --smoke-n-like-max 120 \
  --verbose \
  --progress-newlines
```

## 4. 正式采样模板

默认使用 fast 网格：

```text
n_r=96, n_z=192, n_los=160, n_force_r=48, n_force_z=96
```

正式采样示例：

```bash
python -u scripts/run_galaxy_nautilus.py \
  --mode sample \
  --galaxy eridanus2 \
  --likelihood-mode fast \
  --n-live 1000 \
  --n-eff 2000 \
  --f-live 0.01 \
  --n-processes 16 \
  --n-batch 1000 \
  --verbose \
  --progress-newlines \
  --seed 20260526
```

如果要用 MGE/JAM likelihood 采样，把 likelihood mode 改成：

```bash
--likelihood-mode mge
```

例如：

```bash
python -u scripts/run_galaxy_nautilus.py \
  --mode sample \
  --galaxy eridanus2 \
  --likelihood-mode mge \
  --n-live 1000 \
  --n-eff 2000 \
  --f-live 0.01 \
  --n-processes 16 \
  --n-batch 1000 \
  --verbose \
  --progress-newlines \
  --seed 20260526
```

默认输出文件会自动包含星系名：

```text
outputs/eridanus_ii_nautilus_chain.csv
outputs/diagnostics/eridanus_ii_nautilus_sampler.h5
```

如果手动指定输出路径，文件名应包含当前星系 slug，例如：

```bash
--chain-output outputs/eridanus_ii_nautilus_chain.csv \
--checkpoint-output outputs/diagnostics/eridanus_ii_nautilus_sampler.h5
```

## 5. 恢复采样

如果已有 checkpoint，可以继续运行：

```bash
python -u scripts/run_galaxy_nautilus.py \
  --mode sample \
  --galaxy eridanus2 \
  --likelihood-mode fast \
  --n-live 1000 \
  --n-eff 2000 \
  --f-live 0.01 \
  --n-processes 16 \
  --n-batch 1000 \
  --resume \
  --verbose \
  --progress-newlines
```

## 6. 作图

Density profile：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode plot \
  --galaxy eridanus2 \
  --burn 0
```

默认读取：

```text
outputs/eridanus_ii_nautilus_chain.csv
```

默认写出：

```text
outputs/figures/eridanus_ii_nautilus_density_profile.png
```

如果要指定作图时读取的 chain 文件，用 `--chain-output`。在 `--mode plot` / `--mode corner` 中，这个参数表示输入 chain 路径：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode plot \
  --galaxy eridanus2 \
  --chain-output outputs/eridanus_ii_mge_run1_chain.csv \
  --figure-output outputs/figures/eridanus_ii_mge_run1_density_profile.png
```

Corner map：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode corner \
  --galaxy eridanus2 \
  --burn 0 \
  --show-systemic-in-corner
```

默认写出：

```text
outputs/figures/eridanus_ii_nautilus_corner.png
```

指定输入 chain 和输出 corner 文件：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode corner \
  --galaxy eridanus2 \
  --chain-output outputs/eridanus_ii_mge_run1_chain.csv \
  --corner-output outputs/figures/eridanus_ii_mge_run1_corner.png
```

如果 chain 是用 `--weighted-output` 生成的，作图时建议加：

```bash
--use-sample-weights
```

例如：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode corner \
  --galaxy eridanus2 \
  --chain-output outputs/eridanus_ii_mge_run1_chain.csv \
  --corner-output outputs/figures/eridanus_ii_mge_run1_corner.png \
  --use-sample-weights
```

如果 chain 是默认 equal-weight 输出，则不要加 `--use-sample-weights`。

nautilus 输出已经是 posterior 样本，通常不需要 MCMC burn-in，因此使用 `--burn 0`。

## 7. 系统速度 prior

默认：

```bash
--systemic-prior adaptive
```

adaptive prior 会根据当前星系成员星速度自动生成：

```text
center = median(v_los)
half_width = max(5 * robust_sigma, 0.5 * (vmax - vmin) + padding, min_half_width)
prior = [center - half_width, center + half_width]
```

可调参数：

```bash
--systemic-prior-padding 20
--systemic-prior-min-half-width 20
```

复现旧 Willman 1 测试的系统速度 prior：

```bash
--systemic-prior legacy
```

手动指定：

```bash
--systemic-prior manual \
--systemic-prior-bounds "20,130"
```

## 8. MGE likelihood 与 physicality check

`--likelihood-mode mge` 会调用共享模块：

```python
hayashi_jeans.mge
```

默认使用：

```text
mge_n_gauss_halo   = 45
mge_n_gauss_tracer = 45
mge_n_u            = 96
mge_r_min_pc       = 1e-2
mge_r_max_pc       = 2e5
```

可通过命令行调整：

```bash
--mge-n-gauss-halo 45
--mge-n-gauss-tracer 45
--mge-n-u 96
--mge-r-min-pc 1e-2
--mge-r-max-pc 2e5
```

在 `likelihood_mode == "mge"` 时，默认使用：

```python
mge_sigma_los2_unit_checked()
```

它会在一次 MGE 拟合中同时做局部 `v_phi^2` physicality check 和 `sigma_los^2` 计算。若发现非物理局部二阶矩或非正/非有限 `sigma_los^2`，该参数点会在 likelihood 中得到 `-inf`。

可以关闭 MGE physicality check：

```bash
--no-mge-physicality-check
```

此时 `mge` mode 会使用未检查版本：

```python
mge_sigma_los2_unit()
```

注意：`--no-mge-physicality-check` 只影响 `--likelihood-mode mge`。`fast` 和 `validation-*` 模式不会额外运行 MGE guard。

## 9. Slurm 建议

`sbatch` 脚本中建议设置：

```bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONUNBUFFERED=1
```

运行时让 worker 数匹配 CPU 数：

```bash
python -u scripts/run_galaxy_nautilus.py \
  --mode sample \
  --galaxy eridanus2 \
  --likelihood-mode fast \
  --n-live 1000 \
  --n-eff 2000 \
  --f-live 0.01 \
  --n-processes 16 \
  --n-batch 1000 \
  --verbose \
  --progress-newlines
```

`--progress-newlines` 会让 nautilus 中间状态逐行写入 `.out` 文件，而不是用终端刷新覆盖同一行。

如果在 Slurm 上跑 MGE likelihood，只需把上面的模式改成：

```bash
--likelihood-mode mge
```

## 10. 支持的星系名输入示例

脚本会自动匹配罗马数字和阿拉伯数字，例如：

```bash
--galaxy willman1
--galaxy eridanus2
--galaxy coma_berenices
--galaxy ursamajor2
--galaxy reticulum2
```

也可以使用带空格的正式名称：

```bash
--galaxy "Willman 1"
--galaxy "Eridanus II"
--galaxy "Coma Berenices"
```
