# MGE Physicality Guard Implementation Notes

## 1. 非物理解原因

当前 Hayashi-style Jeans / JAM 问题中的非物理解不是数值归一化错误，而是 Jeans 方程本身在某些参数组合下给出了不可接受的二阶矩解。具体地，垂直 Jeans 方程给出 \(P_z=\nu \overline{v_z^2}\)，径向 Jeans 方程进一步要求

\[
\overline{v_\phi^2}
=
\frac{\overline{v_z^2}+(R/\nu)\partial_R P_z}{1-\beta_z}
+
R\partial_R\Phi .
\]

由于 \(P_z\) 通常随 \(R\) 下降，\(\partial_R P_z\) 可以为负。当 halo/tracer/geometry 使压力下降很陡，同时 \(\beta_z\) 很高时，负的压力梯度项会被 \(1/(1-\beta_z)\) 放大，压过 \(\overline{v_z^2}/(1-\beta_z)\) 和 \(R\partial_R\Phi\)，导致 \(\overline{v_\phi^2}<0\)，进而可能出现 \(\sigma_{\rm los}^2\le0\)。这说明该参数点不是物理 Jeans 解。当前 direct validation projector 中的 `max(..., 0.0)` clipping 会掩盖这个问题，因此采样时应在 likelihood 前加入 raw MGE/JAM physicality check：如果 corrected MGE raw \(\sigma_{\rm los}^2\) 对任一观测星非有限或 \(\le0\)，直接返回 \(-\infty\)。

## 2. 如何修改 `run_willman1_block_mh_fast.py` 以应用 MGE 方法和排除非物理解

在 `scripts/run_willman1_block_mh_fast.py` 中加入一个 MGE physicality helper，复用 `experiments/test_eridanus2_mge_jeans_second_moment.py` 中已经修正 projected tracer denominator 的 `fit_positive_mge`, `generalized_halo_density_unit`, `mge_los_second_moment` 逻辑，或者将这些函数提升到 `src/hayashi_jeans/mge.py` 后从 production 模块导入。helper 接收 `galaxy` 和 `SlowParams`，构造 tracer MGE、halo MGE，计算 unit-density corrected MGE \(\sigma_{\rm los}^2\)，并返回：

```python
np.all(np.isfinite(sigma2_unit)) and np.all(sigma2_unit > 0.0)
```

然后在 slow proposal 接受前、调用 `compute_sigma_unit()` 之前插入检查。也就是在 block-MH slow block 里，当前流程大致是 `proposal_prior = log_prior_slow(...)`，如果 prior finite 就计算 `proposal_sigma_unit = compute_sigma_unit(...)`。改成：

```python
proposal_prior = log_prior_slow(proposal_slow, galaxy)
if np.isfinite(proposal_prior) and mge_physicality_check(galaxy, proposal_slow):
    proposal_sigma_unit = compute_sigma_unit(...)
    proposal_logp = full_log_probability(...)
else:
    proposal_logp = -np.inf
```

初始化 state 时也要做同样检查，避免初始点本身非物理。rho0 和 systemic velocity 的 fast updates 不需要重复 MGE 检查，因为 physicality 与密度归一化和 systemic velocity 无关，只依赖 slow parameters、tracer geometry 和 halo shape。建议给 CLI 加一个开关，例如 `--mge-physicality-check/--no-mge-physicality-check`，默认开启；并在 chain metadata 或日志里记录被 physicality check 拒绝的 slow proposals 数量。

## 3. 如何修改 `run_galaxy_nautilus.py` 以应用 MGE 方法和排除非物理解

在 `scripts/run_galaxy_nautilus.py` 中同样加入 corrected MGE physicality helper，最好与 `run_willman1_block_mh_fast.py` 共用同一个 production helper，避免两边公式漂移。最合适的接入点是 `full_log_probability_vector()`：该函数已经先调用 `log_prior_vector()`，然后构造 `SlowParams`，再调用 `compute_sigma_unit()`。在 `compute_sigma_unit()` 之前插入：

```python
if args.use_mge_physicality_check:
    if not mge_physicality_check(galaxy, slow):
        return -np.inf
```

因为 Nautilus 只看返回的 log-probability，这样非物理点会被当作 hard-prior 外部区域处理，不会进入 posterior。注意这个 check 应使用 unit-density MGE \(\sigma_{\rm los}^2\)，不要乘 `rho0`，因为正负性不依赖密度归一化。建议在 `evaluate_log_probability()` 或 `full_log_probability_vector()` 参数中传入开关，例如 `mge_physicality_check=True`，并在 CLI 加：

```python
parser.add_argument("--no-mge-physicality-check", action="store_true")
```

默认启用检查。对于性能，tracer MGE 只依赖 galaxy 和 inclination 以外的固定 tracer profile，但 `q_star` 依赖 inclination；仍可缓存 tracer radial MGE 与 halo MGE fit，至少用 `functools.lru_cache` 按 rounded slow parameters 缓存。第一版可以先不缓存，因为 corrected MGE 单次约 0.1-0.2 秒，已经比 validation-strict 快很多；之后如果 Nautilus 调用量很大，再把 MGE fit 和 physicality check 做成缓存层。
