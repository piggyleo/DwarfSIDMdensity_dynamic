# 暗物质晕模型替换功能 Handoff

本文档供新对话快速接手当前 Hayashi-style Jeans / Nautilus 管线，并在此基础上加入“可替换暗物质晕模型”的功能。

请新对话优先阅读本文档，再进入代码。当前核心采样功能主要通过：

```text
scripts/run_galaxy_nautilus.py
```

实现。旧的 `scripts/run_willman1_block_mh_fast.py` 仍有参考价值，但当前多星系、Nautilus 采样、输出 schema、plot/corner 主流程应以 `run_galaxy_nautilus.py` 为准。

## 0. 重要约束

新对话的目标是给管线加入暗物质晕模型替换能力。建议保持以下原则：

- 优先做可插拔接口，不要直接大改 likelihood 数学。
- 非必要不要调整现有 Nautilus sampler 调用逻辑、参数先验、数据读取、plotting schema。
- 现有 generalized Hernquist/Hayashi 模型必须作为默认模型继续可用。
- 可插拔接口应是统一入口：`validation-strict`、`fast`、`validation-*`、`mge` 不应各自独立硬接 halo model。理想状态是所有路径从同一个 halo/model factory 或参数对象获得模型。
- `mge` 应作为优先考虑项之一：即使当前 MGE 数值仍需继续验证，新 halo model 接口设计也必须覆盖 MGE 所需的 density/MGE decomposition 信息，避免后续再返工。
- 所有新增模型都应能输出同样的 chain/diagnostics 表，或在文档中明确新增列。

## 1. 当前管线的整体结构

### 1.1 数据输入

主要数据位于：

```text
data/galaxies/
data/processed/galaxy_structural_centers.csv
```

每个星系一个 CSV，例如：

```text
data/galaxies/08_Eridanus_II.csv
data/galaxies/27_Willman_1.csv
```

结构中心来自：

```text
data/processed/galaxy_structural_centers.csv
```

当前代码会用其中的：

```text
ra_center_deg
dec_center_deg
```

作为星系中心，把逐星 RA/Dec 投影为相对平面坐标 `x_pc`, `y_pc`。

核心读取代码：

```text
src/hayashi_jeans/data.py
```

重要对象：

- `GalaxyObservables`
- `GalaxyData`
- `load_galaxy_data(...)`
- `load_low_risk_galaxies(...)`

默认成员选择：

```text
member_flag == 1
```

并要求逐星速度和速度误差有限。

### 1.2 模型参数

当前参数容器：

```text
src/hayashi_jeans/params.py
```

核心类：

```python
HayashiParameters
```

当前字段：

```text
q_halo
b_halo_pc
rho0_msun_pc3
beta_z
alpha
beta
gamma
inclination_rad
systemic_velocity_kms
```

其中 halo 参数默认对应 Hayashi 2023 generalized Hernquist/Zhao profile。

采样向量在 `run_galaxy_nautilus.py` 中定义，顺序为：

```text
[
  q_halo,
  log10_b_halo_pc,
  log10_rho0_msun_pc3,
  minus_log10_one_minus_beta_z,
  alpha,
  beta,
  gamma,
  i_deg,
  systemic_velocity_kms
]
```

链文件列名为：

```text
chain_id
step
q_halo
log10_b_halo_pc
log10_rho0_msun_pc3
minus_log10_one_minus_beta_z
beta_z
alpha
beta
gamma
i_deg
systemic_velocity_kms
log_probability
sampler
likelihood_mode
```

如果加入新 halo model，建议不要破坏这些基础列。可以新增：

```text
halo_model
```

以及模型特有参数列。

### 1.3 暗晕模型接口

当前 halo 接口位于：

```text
src/hayashi_jeans/halos.py
```

已有抽象类：

```python
class HaloModel(ABC):
    def density(self, r_cyl_pc: float, z_pc: float) -> float
    def potential_gradients(self, r_cyl_pc: float, z_pc: float) -> tuple[float, float]
```

当前默认实现：

```python
GeneralizedHernquistHalo
```

它提供：

- `density(R,z)`：返回 \(\rho(R,z)\)，单位 `Msun / pc^3`
- `potential_gradients(R,z)`：返回 \((\partial\Phi/\partial R,\partial\Phi/\partial z)\)，单位 `(km/s)^2 / pc`

这是加入可替换 halo model 的最自然入口。任何新 halo model 如果能实现这两个方法，就能被 slow/reference Jeans projector 使用。

### 1.4 Stellar tracer 与 Jeans projection

stellar tracer：

```text
src/hayashi_jeans/tracer.py
```

核心：

- `AxisymmetricPlummerTracer`
- `intrinsic_q_from_projected(...)`

Jeans/LOS 投影：

```text
src/hayashi_jeans/projection.py
```

核心：

- `AxisymmetricJeansProjector`
- `JeansMoments`

`AxisymmetricJeansProjector` 接收：

```python
halo: HaloModel
```

因此它本身已经具备一定的 halo 可插拔能力。真正还没有完全可插拔的是 sampler 参数 schema、prior、fast grid 和 plotting 对 halo 参数的硬编码。

### 1.5 Likelihood

速度 likelihood 位于：

```text
src/hayashi_jeans/likelihood.py
```

核心类：

```python
GaussianVelocityLikelihood
```

对每颗星：

\[
s_i^2 = \sigma_{{\rm los},i}^2 + \delta v_i^2
\]

\[
\ln L =
-\frac{1}{2}
\sum_i
\left[
\frac{(v_i-u_{\rm sys})^2}{s_i^2}
+
\ln(2\pi s_i^2)
\right].
\]

暗晕替换一般不应修改这里。新 halo model 应通过改变 \(\sigma_{\rm los}^2\) 进入 likelihood。

### 1.6 高层模型 wiring

当前高层 helper：

```text
src/hayashi_jeans/model.py
```

核心函数：

```python
build_hernquist_projector(...)
log_likelihood_hernquist(...)
make_reference_hernquist_params(...)
```

这里目前直接构造：

```python
GeneralizedHernquistHalo(**params.halo_kwargs())
```

如果要让 strict/reference 路径支持多 halo model，建议把这里扩展为：

```python
build_projector(galaxy, params, halo_model="generalized_hernquist", ...)
```

或者新增一个 halo factory，例如：

```python
build_halo_model(model_name, params)
```

并让 `build_hernquist_projector(...)` 保留为向后兼容 wrapper。

## 2. 当前核心采样入口：run_galaxy_nautilus.py

文件：

```text
scripts/run_galaxy_nautilus.py
```

这是当前最重要的管线入口。它承担：

- CLI 参数解析；
- 星系 CSV 自动解析；
- systemic velocity prior；
- parameter vector schema；
- prior；
- likelihood mode dispatch；
- Nautilus sampler adapter；
- chain CSV 写出；
- density profile plot；
- corner plot；
- smoke/schema check。

### 2.1 CLI 模式

`--mode` 支持：

```text
point
smoke
sample
schema-check
plot
corner
all
```

常用采样调用示例：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode sample \
  --galaxy eridanus2 \
  --likelihood-mode fast \
  --chain-output outputs/eridanus_ii_nautilus_chain.csv \
  --checkpoint-output outputs/diagnostics/eridanus_ii_nautilus_sampler.h5
```

只作图：

```bash
python scripts/run_galaxy_nautilus.py \
  --mode plot \
  --galaxy eridanus2 \
  --chain-output outputs/eridanus_ii_nautilus_chain.csv \
  --figure-output outputs/figures/eridanus_ii_nautilus_density_profile.png \
  --burn 0
```

### 2.2 likelihood modes

当前支持：

```text
fast
mge
validation-halo
validation-strict
validation-convergence
```

含义：

- `fast`：默认快速模式，使用 `HaloForceGrid + RZMomentGrid`。
- `validation-halo`：使用 RZMomentGrid，但不使用 HaloForceGrid。
- `validation-convergence`：较高分辨率的 force/moment grid。
- `validation-strict`：最慢 reference projector，直接调用 `build_hernquist_projector(...)`。
- `mge`：MGE/JAM 近似路径，目前仍有数值偏差排查记录；但 halo 替换接口设计应优先考虑它需要的能力。

### 2.3 run_galaxy_nautilus.py 中最重要的函数

建议按下列顺序读：

1. 顶部常量：

```python
LIKELIHOOD_MODES
PARAMETER_NAMES
CHAIN_COLUMNS
```

2. 参数容器：

```python
SlowParams
```

3. 单点 sigma 计算：

```python
compute_sigma_unit(...)
```

这是 halo model 接入的关键函数。它目前会把 `SlowParams` 转为 `HayashiParameters`，然后根据 `likelihood_mode` 分发：

- `mge` -> `hayashi_jeans.mge`
- `validation-strict` -> `build_hernquist_projector(...)`
- `validation-halo` -> `RZMomentGridProjector`
- `fast` / `validation-convergence` -> `RZMomentGridWithForceGrid`

4. likelihood：

```python
log_likelihood_from_unit_sigma(...)
full_log_probability_vector(...)
```

其中 `rho0` 通过线性缩放进入：

```python
sigma_los2 = sigma_unit * 10.0**log10_rho0_msun_pc3
```

5. prior：

```python
log_prior_slow(...)
log_prior_fast(...)
log_prior_vector(...)
```

如果新 halo model 有不同参数，必须同步改这里。

6. vector/row schema：

```python
initial_vector(...)
vector_to_row(...)
row_to_vector(...)
```

新 halo model 的参数列、初值、输出行都要在这里统一。

7. Nautilus adapter：

搜索这些函数：

```text
run_nautilus_sample
make_nautilus_prior_transform
evaluate_log_probability
```

具体函数名可能随近期修改略有变化，但可用 `rg "nautilus|prior_transform|full_log_probability"` 定位。

8. plotting：

```python
plot_density_profile(...)
plot_corner(...)
```

density profile plot 目前假设 chain 里有 generalized Hernquist 参数：

```text
q_halo
log10_b_halo_pc
log10_rho0_msun_pc3
alpha
beta
gamma
```

加入新 halo model 后，plotting 也需要通过 halo factory 计算 \(\rho(r)\)，不能继续硬编码 `GeneralizedHernquistHalo`。

## 3. 加入可替换 halo model 的推荐设计

核心目标：不要为 `validation-strict`、`fast`、`mge` 分别做三套 halo 接入。应先抽象出统一的 halo/model 描述，然后让各 likelihood mode 自然复用它。

推荐的数据流是：

```text
sample vector
  -> model parameter object
  -> halo/model factory
  -> common model components
       - HaloModel.density / potential_gradients
       - optional MGE density descriptor/decomposer
       - plotting density evaluator
  -> compute_sigma_unit(..., likelihood_mode=...)
```

这样新增 halo model 时，主要工作应该集中在“如何描述和构造该 halo”，而不是每条 likelihood path 手动改一遍。

### 3.1 先定义 halo registry/factory

建议新增或扩展：

```text
src/hayashi_jeans/halos.py
```

保留：

```python
HaloModel
GeneralizedHernquistHalo
```

新增类似：

```python
def build_halo_model(model_name: str, params) -> HaloModel:
    ...
```

或者更结构化：

```python
HALO_MODEL_REGISTRY = {
    "generalized_hernquist": GeneralizedHernquistHalo,
}
```

为了避免 sampler 脚本到处知道每个 halo 的构造细节，推荐把“从参数对象构造 halo”的逻辑放在 `src/hayashi_jeans` 内，而不是散落在 `scripts/` 里。

### 3.2 参数 schema 不要硬塞进 HayashiParameters

当前 `HayashiParameters` 是 Hayashi generalized Hernquist 专用。若加入新模型，有两种路线：

路线 A：最小改动。

- 保留 `HayashiParameters`。
- 新增 `halo_model` 字符串。
- 对第一批新模型尽量复用已有字段，例如只换 density law 但仍使用 `q_halo`, `b_halo_pc`, `rho0`。
- 缺点是模型特异参数容易变得混乱。

路线 B：更干净。

- 新增通用参数容器，例如：

```python
@dataclass(frozen=True)
class JeansModelParameters:
    halo_model: str
    halo_params: dict[str, float]
    beta_z: float
    inclination_rad: float
    systemic_velocity_kms: float
```

- 保留 `HayashiParameters` 作为 generalized Hernquist 的 legacy wrapper。
- 修改 `compute_sigma_unit(...)` 和 plotting 使用 `JeansModelParameters` 或 halo factory。

如果目标是长期支持多个 halo model，推荐路线 B。但第一版可以先用路线 A 做小步验证。

### 3.3 统一接入，而不是按路径分叉

由于 `AxisymmetricJeansProjector` 已经接收：

```python
halo: HaloModel
```

strict/reference 路径已经天然接近可插拔；它可以作为 smoke test，但不应成为唯一或特殊的接入路径。接口设计应同时回答：

```text
1. strict/reference 如何拿到 HaloModel；
2. fast grid 如何拿到同一个 HaloModel 或同一个势梯度接口；
3. MGE 如何拿到该 halo 的 intrinsic density 表达或 MGE decomposition 所需函数；
4. plotting 如何用同一个 halo density evaluator 作 density profile。
```

最小目标应改为：

1. 新 halo class 实现 `density()` 和 `potential_gradients()`。
2. 新增统一 factory，例如 `build_halo_model(...)`。
3. `compute_sigma_unit(...)` 不直接硬编码 generalized Hernquist 参数，而是通过统一参数对象/factory 得到模型组件。
4. `validation-strict`、`fast`、`mge` 至少在接口层都能选择 halo model；若某个具体新模型暂时不能用于某个 mode，应在 factory 或 mode dispatch 中显式报错，而不是静默回退。
5. 单点 likelihood smoke test 覆盖 `validation-strict`、`fast`、`mge` 三类路径，至少确认接口一致。

### 3.4 fast grid 路径的注意点

当前 fast 模式使用实验模块：

```text
experiments/test_rz_moment_grid_interpolation.py
experiments/test_halo_force_grid_interpolation.py
```

核心类：

```python
RZMomentGridProjector
RZMomentGridWithForceGrid
```

这两个类目前接收的是 `HayashiParameters`，并可能直接假设 generalized Hernquist halo。新对话需要检查它们内部是否调用：

```python
GeneralizedHernquistHalo
build_hernquist_projector
params.halo_kwargs()
params.q_halo / params.alpha / params.beta / params.gamma
```

如果有硬编码，需要把它们改为接收已经构造好的 `HaloModel`，或接收通用 `model_params` + halo factory。

推荐目标接口：

```python
RZMomentGridProjector(
    galaxy=galaxy,
    tracer=...,
    halo=halo_model,
    beta_z=...,
    inclination_rad=...,
    ...
)
```

而不是在 grid 类内部构造特定 halo。

### 3.5 MGE 应优先纳入接口设计

当前有 MGE/JAM 相关代码：

```text
src/hayashi_jeans/mge.py
docs/mge_los_second_moment_derivation.md
docs/mge_physicality_guard_for_ai.md
```

当前解析 MGE 版本相对 `validation-strict` 仍有较大偏差，但这不意味着 halo 替换设计可以暂时忽略 MGE。相反，MGE 对 halo model 的要求更强，应该在接口设计阶段优先考虑：

```text
1. 新 halo model 是否能提供 intrinsic 3D density 函数 rho(m) 或 rho(R,z)；
2. MGE decomposition 应该直接分解哪个函数；
3. halo flattening、scale radius、normalization 在 MGE 振幅中如何定义；
4. signed MGE amplitudes 是否需要 physicality guard；
5. MGE density profile 与 HaloModel.density 是否能在同一组半径上互相验证。
```

建议新对话把 MGE 作为接口验收项之一：

```text
新增 halo model 后，至少应能调用 mge mode 的构造路径；若数值尚未可信，应输出明确 warning/diagnostic，而不是让 mge 路径继续默认假设 generalized Hernquist。
```

## 4. 如何阅读代码

建议新对话按以下顺序阅读。

### Step 1：读项目状态

先读：

```text
PROJECT_MEMORY.md
docs/run_galaxy_nautilus_commands.md
docs/nautilus_sampler_handoff.md
```

目的：

- 了解低风险星系选择；
- 了解 Hayashi Table 1 参数来源；
- 了解 Nautilus 命令；
- 了解现有 chain schema。

### Step 2：读数据层

读：

```text
src/hayashi_jeans/data.py
data/galaxies/README.md
data/processed/README.md
```

重点确认：

- `load_galaxy_data(...)` 如何选星；
- `x_pc`, `y_pc` 如何由 RA/Dec 得到；
- `GalaxyData` 给 likelihood 暴露了哪些数组。

### Step 3：读 halo 接口

读：

```text
src/hayashi_jeans/halos.py
```

重点看：

- `HaloModel` 抽象接口；
- `GeneralizedHernquistHalo.density(...)`；
- `GeneralizedHernquistHalo.potential_gradients(...)`；
- `force_integral_method="unit_interval"` 和 `force_quadrature_order`。

这是新 halo model 最应模仿的文件。

### Step 4：读 Jeans projector

读：

```text
src/hayashi_jeans/projection.py
```

重点看：

- `AxisymmetricJeansProjector.__init__(..., halo: HaloModel, ...)`
- `sigma_los2(...)`
- `jeans_moments(...)`
- `_vertical_pressure(...)`
- `_d_vertical_pressure_dr(...)`

理解数据流：

```text
HaloModel.potential_gradients
  -> vertical pressure
  -> Jeans moments
  -> LOS projection
  -> sigma_los2 for each star
```

### Step 5：读 likelihood

读：

```text
src/hayashi_jeans/likelihood.py
```

确认 likelihood 只依赖：

```text
sigma_los2
systemic_velocity_kms
observed velocities
velocity errors
```

暗晕替换不应改这里。

### Step 6：读核心 sampler

读：

```text
scripts/run_galaxy_nautilus.py
```

建议用搜索定位：

```bash
rg -n "PARAMETER_NAMES|CHAIN_COLUMNS|class SlowParams|def compute_sigma_unit|def full_log_probability_vector|def log_prior|def vector_to_row|def row_to_vector|nautilus|plot_density_profile" scripts/run_galaxy_nautilus.py
```

阅读重点：

- 参数向量如何定义；
- prior 如何限制参数；
- `compute_sigma_unit(...)` 如何根据 likelihood mode 调用不同 projector；
- `rho0` 如何线性缩放；
- Nautilus 如何把 unit cube prior transform 到物理参数；
- chain output 如何写出；
- plot 如何从 chain 重建 density profile。

### Step 7：读 fast grid 实验代码

读：

```text
experiments/test_rz_moment_grid_interpolation.py
experiments/test_halo_force_grid_interpolation.py
```

目标不是理解每行细节，而是找出 generalized Hernquist 的硬编码点。用：

```bash
rg -n "GeneralizedHernquistHalo|HayashiParameters|build_hernquist|halo_kwargs|alpha|gamma|potential_gradients" experiments/test_rz_moment_grid_interpolation.py experiments/test_halo_force_grid_interpolation.py
```

如果要让 `fast` 支持新 halo model，这两个文件通常要跟着改。

### Step 8：读 plotting

在：

```text
scripts/run_galaxy_nautilus.py
scripts/plot_eridanus2_density_profile_with_hayashi_overlay.py
scripts/plot_willman1_density_profile_with_hayashi_overlay.py
```

重点查找：

```bash
rg -n "GeneralizedHernquistHalo|density\\(|log10_b_halo|alpha|beta|gamma" scripts
```

当前 density profile plotting 基本假设 generalized Hernquist 参数列。新 halo model 接入后，plotting 应改为通过 halo factory 计算密度。

## 5. 建议的新对话任务拆分

建议按以下小步推进。

### Task A：只做接口设计，不改 sampler 行为

目标：

- 新增 halo factory；
- 保留 generalized Hernquist 默认行为；
- 跑现有 Eridanus II / Willman 1 smoke test，确认结果不变。

### Task B：让所有 likelihood mode 共用 halo_model 参数

目标：

- CLI 加：

```text
--halo-model generalized-hernquist
```

- 默认值保持现状；
- `compute_sigma_unit(...)` 通过统一参数对象/factory 构造或获取 halo/model 组件；
- `validation-strict`、`fast`、`validation-*`、`mge` 都不能继续暗中硬编码 generalized Hernquist；
- 跑 `--mode point` 或 `--mode smoke`，优先覆盖 `fast` 和 `mge`。

### Task C：扩展参数 schema

目标：

- 如果新 halo model 有新参数，新增 `PARAMETER_NAMES`、prior、prior transform、chain columns；
- 保证旧模型输出不变；
- 加 `schema-check`。

### Task D：让 fast grid 和 MGE 都通过统一接口支持新 halo

目标：

- 把 fast grid 内部的 halo 构造抽象出来；
- 确保 `RZMomentGridProjector` / `RZMomentGridWithForceGrid` 不依赖特定 density law；
- 确保 `src/hayashi_jeans/mge.py` 不再默认只接受 generalized Hernquist 参数，而是从统一 halo/MGE descriptor 获取 density；
- 和 `validation-strict` 做若干点误差测试；
- 对 `mge` 单独输出 density reconstruction / likelihood diagnostic，方便继续排查 MGE 数值误差。

### Task E：更新 plotting

目标：

- density profile plot 用 halo factory；
- 支持按 chain 中的 `halo_model` 计算 \(\rho(r)\)；
- 对 Hayashi generalized Hernquist 链保持原图一致。

## 6. 推荐给新对话的起始指令

可直接复制给新对话：

```text
请阅读 /Users/wangkaihao/Documents/New project 4/docs/halo_model_replacement_handoff.md，并以 scripts/run_galaxy_nautilus.py 为当前核心采样入口。目标是在不破坏现有 generalized Hernquist 默认行为的前提下，为管线加入可替换暗物质晕模型功能。

请先只做代码阅读和接口设计，重点检查：
1. src/hayashi_jeans/halos.py 的 HaloModel 接口；
2. src/hayashi_jeans/model.py 的 build_hernquist_projector；
3. scripts/run_galaxy_nautilus.py 中的 PARAMETER_NAMES、CHAIN_COLUMNS、SlowParams、compute_sigma_unit、log_prior_vector、full_log_probability_vector、vector_to_row、row_to_vector、plot_density_profile；
4. experiments/test_rz_moment_grid_interpolation.py 和 experiments/test_halo_force_grid_interpolation.py 是否硬编码 generalized Hernquist；
5. src/hayashi_jeans/mge.py 的 MGE density decomposition 入口如何改为使用统一 halo/MGE descriptor。MGE 应作为优先考虑项，不要把它排除在第一版接口设计之外。

请先给出最小改动方案和文件级修改计划，不要立即大改 sampler 或 likelihood。
```
