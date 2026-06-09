# SIDM 的 Shajib MGE 反演与复分支结构

本文记录当前 SIDM MGE 解析分解问题的完整判断。讨论对象是 Shajib
(2019) 的 Gaussian-kernel integral transform，以及当前
`src/hayashi_jeans/mge.py` 中采用的 Euler inverse-Laplace 离散实现。

主要参考：

- [Shajib (2019), Unified lensing and kinematic analysis for any elliptical mass profile](https://academic.oup.com/mnras/article/488/1/1387/5526256)
- [arXiv:1906.08263](https://arxiv.org/pdf/1906.08263)

## 1. 核心结论

SIDM 在正实半径上的密度轮廓可以被 Gaussian 和准确表示；当前问题不在
MGE 表示能力，也不在后续 MGE/JAM 矩积分，而在当前使用的复数 Euler
逆 Laplace 反演：

\[
\text{SIDM 复延拓的分支结构}
\longrightarrow
\text{有限 Euler 求和不能稳定表示逆变换}
\longrightarrow
\text{巨大交替权重放大误差}
\longrightarrow
\text{MGE 振幅失真}.
\]

更严格地说，需要区分三层：

1. Shajib 的连续 Gaussian 积分变换；
2. 论文采用的 Euler 离散逆 Laplace 算法；
3. 当前代码用 NumPy principal complex power 逐点计算 SIDM 复密度的实现。

实轴验证说明 SIDM 可以被有限 Gaussian 和准确逼近，但这本身不证明它满足
论文对精确连续变换的全部假设；即使第 1 层成立，也不代表第 2、3 层的当前
组合一定稳定。

## 2. Shajib 反演做了什么

目标是将径向轮廓 \(F(x)\) 写成

\[
F(x)\approx\sum_j A_j
\exp\left(-\frac{x^2}{2\sigma_j^2}\right).
\]

先定义连续 Gaussian 叠加：

\[
F(x)=
\int_0^\infty
\frac{f(\sigma)}{\sqrt{2\pi}\sigma}
\exp\left(-\frac{x^2}{2\sigma^2}\right)d\sigma.
\]

适当换元后，这是一个 Laplace 变换。其复积分反演为

\[
f(\sigma)=
\frac{1}{i\sigma^2}\sqrt{\frac{2}{\pi}}
\int_C zF(z)
\exp\left(\frac{z^2}{2\sigma^2}\right)dz.
\]

Shajib 使用 Euler inverse-Laplace algorithm 将它近似成

\[
f(\sigma)\approx
\sum_{n=0}^{2P}
\eta_n\operatorname{Re}[F(\sigma\chi_n)],
\]

\[
\chi_n=
\left[
\frac{2P\log 10}{3}+2\pi i n
\right]^{1/2}.
\]

因此算法需要在复半径

\[
z_n=\sigma\chi_n
\]

上计算 \(F(z_n)\)。物理密度只使用 \(r>0\)，但该数值反演主动离开正实轴。

得到 \(f(\sigma_j)\) 后，在对数等距的 \(\sigma_j\) 网格上做梯形离散：

\[
A_j=
\frac{w_j f(\sigma_j)\Delta\log\sigma}{\sqrt{2\pi}}.
\]

最终用这些 \(A_j,\sigma_j\) 重建有限项 MGE。

## 3. SIDM 的复分支点

当前 SIDM 密度的关键因子是

\[
F(z)\propto
\left(z^4+r_c^4\right)^{-\gamma/4}
\left[1+(z/r_s)^\alpha\right]^{-(3-\gamma)/\alpha}.
\]

它在正实轴上是平滑的，但非整数复数幂是多值函数，需要选择解析分支。

### 3.1 核项

核项

\[
F_{\rm core}(z)=
\left(z^4+r_c^4\right)^{-\gamma/4}
\]

的分支点由

\[
z^4+r_c^4=0
\]

给出：

\[
z_k=r_c
\exp\left[
i\left(\frac{\pi}{4}+\frac{k\pi}{2}\right)
\right].
\]

第一象限分支点位于

\[
z_k=r_c e^{i\pi/4}.
\]

在该点附近令 \(z=z_k+\delta z\)。一阶展开为

\[
z^4+r_c^4
=4z_k^3\delta z+O(\delta z^2),
\]

因此

\[
F_{\rm core}(z)
\sim C(\delta z)^{-\gamma/4}.
\]

这说明分支点附近的相对变化率为

\[
\frac{F'(z)}{F(z)}
\sim
-\frac{\gamma/4}{z-z_k}.
\]

当节点接近 \(z_k\) 时，同样大小的路径离散误差会引起更大的函数值和相位变化。

### 3.2 外层项

外层项为

\[
F_{\rm outer}(z)=
\left[1+(z/r_s)^\alpha\right]^{-b},
\qquad
b=\frac{3-\gamma}{\alpha}.
\]

除 \(z^\alpha\) 自身的分支结构外，括号的零点满足

\[
1+(z/r_s)^\alpha=0.
\]

第一象限可能出现

\[
z_k=r_s e^{i\pi/\alpha}.
\]

当 \(\alpha>2\) 时，

\[
0<\frac{\pi}{\alpha}<\frac{\pi}{2},
\]

所以该分支点进入第一象限。对于当前 SIDM `tau=1.08`，

\[
\alpha\simeq5.476,
\qquad
\arg z_k\simeq32.87^\circ.
\]

这位于当前 \(P=28\) Euler 节点约
\(0^\circ\)--\(41.52^\circ\) 的角度范围内。

## 4. 分支切线为何关键

### 4.1 分支点与分支切线不是同一件事

分支点是函数自身的多值结构。例如

\[
(z-z_k)^{-a}
\]

绕 \(z_k\) 一周后会多出因子

\[
e^{-2\pi ia},
\]

当 \(a\) 不是整数时不会回到原值。

为了把多值函数变成单值函数，必须切开复平面。NumPy 通常采用 principal
logarithm：

\[
-\pi<\arg w\leq\pi,
\]

即在中间变量 \(w\) 的负实轴上放置分支切线。

分支点由函数决定；分支切线的具体位置来自所选分支，但任何单值分支都必须
以某种方式阻止路径完整绕过分支点。

### 4.2 SIDM 核项的 principal branch cut

对

\[
w=z^4+r_c^4
\]

使用 principal power 时，\(w\) 平面的切线是 \(w\leq0\)。映射回第一象限的
\(z\) 平面，条件是

\[
z=re^{i\pi/4},\qquad r\geq r_c.
\]

因此切线从 \(r_c e^{i\pi/4}\) 开始，沿 \(45^\circ\) 射线向外延伸。

Euler 节点随 \(n\) 增大逐渐逼近 \(45^\circ\)，但有限 \(P\) 时仍位于其下方。
所以核项的主要问题不是节点必然穿过切线，而是反演路径不断逼近分支点和
收敛域边界，导致收敛变慢和求和病态。

### 4.3 SIDM 外层项的 principal branch cut

对

\[
w=1+(z/r_s)^\alpha
\]

使用 principal power 时，第一象限切线的一个分支为

\[
z=re^{i\pi/\alpha},\qquad r\geq r_s.
\]

当 \(\alpha>2\) 时，这条切线进入第一象限。`tau=1.08` 时它沿
\(32.87^\circ\) 向外延伸，因此当前 Euler 节点可以分布在切线两侧。

切线两侧令 \(w=|w|e^{\pm i\pi}\)，则

\[
w^{-b}=|w|^{-b}e^{\mp i\pi b}.
\]

两侧相差有限相位因子，而不是随着两侧点距离趋零而趋于同一个值。单独
\(w^{-b}\) 的两侧实部可能相同，但完整 SIDM 密度还乘有其他复因子，因此
\(\operatorname{Re}F(z)\) 一般也会受到相位跳变影响。

更根本的问题是：逐点调用 principal power 只根据当前点的 principal
argument 选值，不保证所有节点属于沿目标积分路径连续延拓得到的同一解析分支。

### 4.4 与 Shajib 反演路径的关系

论文附录的精确反演路径满足

\[
\operatorname{Re}(z^2)=d,
\]

并要求整条超曲线 \(C\) 位于 \(F(z)\) 的 region of convergence (ROC) 内。

核项 \(45^\circ\) 分支点满足

\[
\operatorname{Re}(z_k^2)=0,
\]

所以它位于常用右侧收敛扇区的边界。路径接近它时可能严重病态，但不等于
每条合法 \(d>0\) 路径都穿过切线。

晚期 SIDM 外层分支点则可能满足

\[
\operatorname{Re}(z_k^2)>0.
\]

这意味着它进入了通常的右半 \(z^2\) 平面。此时若仍使用当前未经平移或
变形的 Euler 节点，就不能自动保证对应路径位于正确 ROC 内；必须重新检查
路径是否跨越切线，并在发生 contour deformation 时判断是否需要加入
branch-cut contribution。

## 5. 误差如何产生和放大

需要避免一个误解：问题不一定是 NumPy 将每个
\(F(\sigma\chi_n)\) 都算错了。单点 principal-branch 值可能接近机器精度，
但它们的有限组合仍可能不是目标连续复积分的准确离散。

误差依次来自：

1. 分支点附近 \(F(z)\sim(z-z_k)^{-a}\)，函数值和相位快速变化；
2. 有限 Euler 节点可能错过这种窄的复积分结构，或让单个节点过度代表它；
3. 若节点跨越 principal branch cut，逐点值可能切换到另一解析分支；
4. 因此有限求和与目标 Bromwich/超曲线积分之间出现离散和分支选择误差；
5. Euler 权重巨大且正负交替，最终答案依赖高精度相消；
6. 小的中间误差被放大成巨大的 \(f(\sigma)\) 和 MGE 振幅误差。

令

\[
T_n=\eta_n\operatorname{Re}F(\sigma\chi_n).
\]

求和的病态程度可用

\[
\kappa_{\rm sum}
=
\frac{\sum_n|T_n|}
{\left|\sum_nT_n\right|}
\]

表示。当前 `tau=0.5` 的一个实际节点组中：

\[
\sum_n|T_n|\simeq1.34\times10^{12},
\qquad
\left|\sum_nT_n\right|\simeq0.215,
\]

所以

\[
\kappa_{\rm sum}\simeq6.2\times10^{12}.
\]

若希望最终结果保留约 \(1\%\) 相对精度，中间贡献的综合相对误差需低于约
\(1.6\times10^{-15}\)。这已经接近双精度极限，还没有包含有限节点对分支结构
的积分近似误差。

完整传播链为

\[
\begin{aligned}
&\text{复分支点/切线靠近或进入 Euler 采样区域}\\
&\Downarrow\\
&\text{有限节点不能稳定表示连续逆变换}\\
&\Downarrow\\
&\text{交替大权重破坏精细相消}\\
&\Downarrow\\
&f(\sigma)\text{ 和 }A_j\text{ 失真}\\
&\Downarrow\\
&\rho_{\rm MGE},\,M(r),\,\sigma_{\rm los}^2,\,
\log\mathcal L\text{ 产生误差}.
\end{aligned}
\]

## 6. 与标准 Hernquist profile 的对照

标准三维 Hernquist profile 可写成

\[
\rho_{\rm H}(z)
\propto
\frac{1}{(z/b)(1+z/b)^3}.
\]

它是有理函数，在 \(z=0\) 和 \(z=-b\) 有极点，但没有非整数幂产生的
有限半径分支点或分支切线。Shajib Euler 节点使用 principal square-root
后位于右半 \(z\) 平面，因此不会穿过 \(z=-b\)，而 \(z=0\) 位于
\(\operatorname{Re}(z^2)=0\) 边界。

使用当前相同设置：

- `n_gauss=60`
- `P=28`
- `r_min=1e-2 pc`
- `r_max=2e5 pc`
- 验证半径 `0.1--5e4 pc`
- 标度半径 `1600 pc`

得到密度重建误差：

| Profile | median relative error | p95 relative error | max relative error |
|---|---:|---:|---:|
| Standard Hernquist | \(1.56\times10^{-7}\) | \(4.34\times10^{-5}\) | \(5.70\times10^{-4}\) |
| SIDM, \(\tau=0\) | \(3.42\times10^{-7}\) | \(6.83\times10^{-4}\) | \(4.56\times10^{-3}\) |
| SIDM, \(\tau=0.1\) | \(3.41\times10^{-1}\) | \(7.46\times10^{-1}\) | \(9.39\times10^{-1}\) |
| SIDM, \(\tau=0.5\) | \(1.98\) | \(6.16\) | \(6.17\) |
| SIDM, \(\tau=1.08\) | \(4.02\times10^2\) | \(2.95\times10^6\) | \(3.35\times10^6\) |

对照解释：

- 标准 Hernquist 在 Euler 采样扇区内没有 SIDM 型分支切线，分解稳定；
- `tau=0` 为 NFW-like 极限，没有有限 \(r_c\)，也保持良好精度；
- `tau>0` 新增 \(45^\circ\) 核分支点后误差显著上升；
- `tau=1.08` 的外层分支点和切线进入节点角度范围，误差爆发；
- 实轴 signed-MGE 仍能准确拟合同一 SIDM 轮廓，说明失败的是当前复反演，不是
  Gaussian 表示或 JAM 公式。

这里的控制组是标准 Hernquist。Generalized Hernquist 若使用非整数指数，
尤其 \(\alpha>2\)，同样可能在第一象限产生分支结构，不能仅凭
“Hernquist family”名称假定解析 Euler 反演总是稳定。

## 7. Shajib 方法的限制条件

论文附录给出的数学条件为：

1. 逆权重 \(f(\sigma)\) 分段连续；
2. 当 \(\sigma\to0\) 时，
   \[
   f(\sigma)=O\left(e^{c/(2\sigma^2)}\right),\quad c\geq0;
   \]
3. 当 \(\sigma\to\infty\) 时，
   \[
   f(\sigma)=O(\sigma^\lambda),\quad\lambda<0;
   \]
4. 在这些条件下，变换存在于
   \[
   \operatorname{Re}(z^2)>c;
   \]
5. 唯一性定理额外假设 \(f(\sigma)\) 连续；
6. 反演定理假设 \(F(z)\) 确实是该 \(f(\sigma)\) 的变换，并且
   \(C:\operatorname{Re}(z^2)=d\) 完整位于 ROC 内。

特别需要注意：论文证明的是“若 \(F\) 是某个满足条件的 \(f\) 的变换，则可
按该公式反演”，并没有证明任意正实轴上光滑的 \(F(x)\) 都必然满足这些条件。
实轴光滑性不能替代复平面 ROC 和解析延拓条件。

因此，以下情况下不能直接套用当前 Shajib/Euler 实现：

- 只知道实轴轮廓，却没有可在全部复节点上连续、一致求值的复延拓；
- 目标 contour 穿过极点、分支点或分支切线，且没有做 contour deformation
  或 branch-cut contribution 修正；
- 使用的 principal branch 不是沿目标 contour 连续的同一解析分支；
- \(F(z)\) 在采样区域内不满足相应 ROC/解析性条件；
- 增加 \(P\)、Gaussian 数量或浮点精度后结果不收敛；
- 有限 \(\sigma\) 范围截断了重要宽度，或 Gaussian 数量不足；
- 二维轮廓不具椭圆对称性，或者轴比随半径变化，却仍尝试用一次一维变换和
  单一共同轴比表示。

论文推导本身并未因此“错误”。错误的是在不满足上述假设时，仍把方程 (4)
或其 Euler 近似方程 (5) 当作无条件成立，并沿用其理想误差估计
\(O(10^{-0.6P})\)。

此外，有限精度本身构成限制：论文指出双精度下 \(P\gtrsim27\) 后通常不会
继续改善。若求和已经高度病态，增大 \(P\) 还可能因权重增大而恶化结果。

## 8. 这是否意味着 Shajib 方法不适用于 SIDM

不能笼统地说 Shajib 的 Gaussian integral transform 不适用于 SIDM。
SIDM 的实轴轮廓仍可被有限 MGE 准确逼近，实轴拟合结果已经证明这一点；
但这不等同于证明论文的精确连续反演对该闭式复延拓无条件成立。

应当采用更窄、更准确的结论：

> 当前未经 contour/branch 修正的固定 Euler 公式，加上逐点 principal
> complex-power 的 SIDM 延拓，不适合作为 SIDM 全部 \(\tau\) 范围内的
> 无条件解析分解器。

具体区分如下：

- 早中期 SIDM：核分支结构主要位于 ROC 边界附近，理论反演未必不存在，但
  当前双精度 Euler 求和可能严重病态；
- 晚期 SIDM：外层分支点可进入右侧 \(z^2\) 平面，当前未平移的 Euler contour
  对部分 \(\sigma\) 不再自动满足论文反演条件，必须重新处理 ROC、contour 和
  解析分支；
- 若能选择位于所有奇点右侧的合法 contour、沿 contour 构造一致解析分支，
  并正确计入 branch-cut contribution，Shajib/Laplace 思路原则上仍可能使用；
- 在这些修改完成并验证前，应将实轴 signed-MGE 作为 SIDM 的可靠实现，并对
  密度、质量、速度弥散和 likelihood 做端到端误差检查。

## 9. 当前工程判断

现阶段推荐：

1. 标准 Hernquist/generalized Hernquist 继续保留现有解析入口，但每组参数仍做
   实轴重建验证；
2. SIDM 默认使用经过验证的 real-axis signed-MGE；
3. 解析 Euler 分解仅作为诊断或实验路径，不应在失败后静默进入 likelihood；
4. 后续若恢复 SIDM 解析反演，应优先研究 shifted/deformed Bromwich contour、
   连续 branch tracking 和 branch-cut contribution，而不是单纯增大 \(P\)；
5. 验收必须同时覆盖 density、enclosed mass、\(\sigma_{\rm los}^2\) 和
   \(\log\mathcal L\)，不能只检查复节点函数值或 MGE 振幅是否有限。

## 10. High-accuracy signed-MGE configuration

当前 SIDM sampler 默认采用：

```text
n_gauss_halo = 60
n_gauss_tracer = 60
real_fit_radii = 1600
real_lstsq_rcond = 1e-12
n_u = 96
```

generalized Hernquist 的 sampler 默认仍为 `45/45` 个 halo/tracer Gaussian，
因此该精度调整不会改变其原有默认分解规模。

在 Willman 1、\(r_{s0}=1600\ {\rm pc}\)、\(\rho_{s0}=0.1\
{\rm M_\odot\,pc^{-3}}\)、\(\beta_z=1-10^{-0.3}\)、\(i=75^\circ\) 的测试点：

| \(\tau\) | signed-MGE \(\log L\) | converged strict \(\log L\) | \(\Delta\log L\) | median time |
|---:|---:|---:|---:|---:|
| 0.50 | -139.515141841 | -139.515255741 | \(+1.14\times10^{-4}\) | 0.453 s |
| 1.08 | -171.066680968 | -171.066083785 | \(-5.97\times10^{-4}\) | 0.436 s |

全 \(\tau=0\)--1.08 验证网格上的最大 p95 density relative error 为
\(2.45\times10^{-4}\)。

这里的 strict reference 需要同时检查积分容差和有限积分边界。特别是
\(\tau=0.5\) 时，将 `zmax_factor=los_factor` 从 20 扩大到 60 会使
\(\log L\) 变化约 \(1.3\times10^{-3}\)。因此，用 factor 20 的旧 reference
衡量无限域 MGE/JAM 公式时，会把积分边界误差误认为 MGE 分解误差。
