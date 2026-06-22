# `run_galaxy_nautilus.py` 常用指令

本文记录 `scripts/run_galaxy_nautilus.py` 的常用调用方式。该脚本支持用 `--galaxy` 自动选择星系数据、结构中心、默认输出文件名，并默认使用自适应系统速度先验。

## 0. 数据读取口径

`run_galaxy_nautilus.py` 通过 `hayashi_jeans.data.load_galaxy_data()` 读取每星系 CSV。

当前默认规则是：

```text
读取 model_member_flag in [1, 2]
```

其中：

```text
model_member_flag = 1  baseline clean member
model_member_flag = 2  默认纳入，但需要敏感性测试的成员
model_member_flag = 3  不进入默认 Jeans 建模
```

所有 27 个星系都已迁移到 `model_member_flag`。其中 14 个统一试分类星系的主 CSV `star` 行已经是一颗唯一恒星一行；多历元速度、跨源去重和速度变量排除已在星表整理阶段完成。冻结备用星系的科学成员样本未改变，只新增了 `model_member_flag` 以满足统一读取接口。缺少 `model_member_flag` 时 loader 会报错，这是预期行为。

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

默认输出 stem 会包含星系名和 halo 模型，不再包含采样器名称。默认
`--halo-model generalized-hernquist` 时：

```text
outputs/eridanus_2_generalized_hernquist_chain.csv
outputs/diagnostics/eridanus_2_generalized_hernquist_sampler.h5
```

使用 `--halo-model sidm` 时，默认改为：

```text
outputs/eridanus_2_sidm_chain.csv
outputs/diagnostics/eridanus_2_sidm_sampler.h5
outputs/figures/eridanus_2_sidm_density_profile.png
outputs/figures/eridanus_2_sidm_corner.png
```

可以继续用 `--output-name` 覆盖共同 stem。

如果手动指定输出路径，文件名应包含当前星系 slug，例如：

```bash
--chain-output outputs/eridanus_2_generalized_hernquist_chain.csv \
--checkpoint-output outputs/diagnostics/eridanus_2_generalized_hernquist_sampler.h5
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
outputs/eridanus_2_generalized_hernquist_chain.csv
```

默认写出：

```text
outputs/figures/eridanus_2_generalized_hernquist_density_profile.png
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
outputs/figures/eridanus_2_generalized_hernquist_corner.png
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

## 11. 批量运行全部 SIDM 星系

批量配置位于：

```text
config/sidm_galaxies.tsv
```

采样和作图由同一个 Slurm array task 顺序执行：

```text
scripts/slurm_fit_all_sidm.sh
```

驱动脚本将 27 个目标拆成每批 3 个：

```bash
scripts/submit_all_sidm.sh
```

每批使用 `sbatch --wait`。当前批的 3 个星系全部结束后，才提交下一批，因此
Slurm 中最多登记和运行 3 个目标任务。每个 task 在采样成功后立即生成 density
profile 和 corner；如果当前批失败或超时，驱动脚本停止，不再提交后续批次。

查看运行状态：

```bash
squeue -u "$USER"

sacct -j <fit_job_id> \
  --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

如果 task `3,8,17` 超时，只恢复并重跑这些 checkpoint，同时重新提交对应作图：

```bash
scripts/resubmit_sidm_tasks.sh 3,8,17 04:00:00
```

第二个参数是新的采样时间上限；省略时默认使用 `04:00:00`。

如果分批提交因某个批次失败而停止，可以在处理失败 task 后，从后续编号继续：

```bash
scripts/submit_all_sidm.sh 10
```

这表示从 task 10 继续到 task 27。也可以指定起止范围：

```bash
scripts/submit_all_sidm.sh 10 18
```

不传参数时仍执行全部 task：

```bash
scripts/submit_all_sidm.sh
```

日志文件按 array job 和 task 编号保存：

```text
logs/sidm_fit_plot_<array-job-id>_<task-id>.out
logs/sidm_fit_plot_<array-job-id>_<task-id>.err
```
