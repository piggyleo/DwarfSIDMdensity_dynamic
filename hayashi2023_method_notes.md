# Hayashi et al. (2023) 方法笔记

这份笔记整理了论文 *Dark Matter Halo Properties of the Galactic Dwarf Satellites* 的动力学建模流程，重点是：

- 作者如何定义似然函数；
- 如何由观测数据和模型参数计算理论视向速度弥散；
- 代码实现时应拆分成哪些步骤与接口。

本文档主要依据：

- Hayashi et al. (2023), ApJ 953, 185
- Hayashi et al. (2020), ApJ 904, 45
- Hayashi & Chiba (2012, 2015) 的轴对称 Jeans 建模框架

## 1. 总体思路

作者对每个矮星系做的是**逐星、非分箱**拟合。

输入是成员星的观测数据：

- 天球平面位置 $(x_i, y_i)$
- 视向速度 $u_i$
- 速度测量误差 $\delta_{v,i}$

目标是给定一组模型参数

$$
\theta = (Q,\ b_{\rm halo},\ \rho_0,\ \beta_z,\ \alpha,\ \beta,\ \gamma,\ i)
$$

计算每颗星位置处的理论视向速度弥散

$$
\sigma_{\rm los}^2(x_i, y_i \mid \theta)
$$

再代回高斯似然，得到整组观测数据在该模型下的似然值。

## 2. 模型假设

作者的动力学框架基于以下假设：

- 星系处于稳态；
- 系统轴对称；
- 以暗物质势为主，恒星自引力忽略；
- 恒星是 tracer population；
- 速度椭球与柱坐标 $(R,\phi,z)$ 对齐；
- 交叉项 $\overline{v_R v_z} = 0$；
- 垂直各向异性参数 $\beta_z$ 为常数；
- streaming motion 可忽略，因此二阶矩主要由速度弥散贡献。

其中各向异性定义为

$$
\beta_z = 1 - \frac{\overline{v_z^2}}{\overline{v_R^2}}.
$$

因此有

$$
\overline{v_R^2} = \frac{\overline{v_z^2}}{1 - \beta_z}.
$$

## 3. 恒星 tracer 密度模型

作者对恒星分布采用轴对称 Plummer 模型：

$$
\nu(R,z) \propto \left(1 + \frac{m_*^2}{b_*^2}\right)^{-5/2},
$$

其中

$$
m_*^2 = R^2 + \frac{z^2}{q^2}.
$$

这里：

- $\nu(R,z)$ 是三维恒星数密度；
- $b_*$ 是沿主轴定义的半光半径；
- $q$ 是恒星分布的内禀轴比。

观测给出的通常是投影轴比 $q'$，而不是 $q$。二者通过倾角 $i$ 联系：

$$
q'^2 = \cos^2 i + q^2 \sin^2 i,
$$

从而

$$
q^2 = \frac{q'^2 - \cos^2 i}{\sin^2 i}.
$$

因此，给定 $q'$ 和自由参数 $i$，就能确定用于 Jeans 方程的内禀恒星密度形状。

## 4. 暗晕密度模型

2023 论文采用的是轴对称广义 Hernquist（也可看作 Zhao 型）密度轮廓：

$$
\rho_{\rm DM}(R,z) =
\rho_0
\left(\frac{m}{b_{\rm halo}}\right)^{-\gamma}
\left[1 + \left(\frac{m}{b_{\rm halo}}\right)^\alpha\right]^{-(\beta-\gamma)/\alpha},
$$

其中

$$
m^2 = R^2 + \frac{z^2}{Q^2}.
$$

参数含义：

- $\rho_0$：尺度密度
- $b_{\rm halo}$：尺度半径
- $\gamma$：内层斜率
- $\beta$：外层斜率
- $\alpha$：内外层过渡的锐利程度
- $Q$：暗晕内禀轴比

后续如果我们要替换成新的暗晕模型，最好不要把代码接口写死成 “Hernquist 专用”，而是抽象成：

- 给定 $(R,z)$ 和参数返回 $\rho_{\rm DM}$；
- 或更进一步，直接返回引力势梯度 $\partial\Phi/\partial R$ 与 $\partial\Phi/\partial z$。

因为 Jeans 方程真正直接使用的是势梯度，而不是密度本身。

## 5. 轴对称 Jeans 方程

在上述假设下，作者使用的 Jeans 方程可以写成：

$$
\nu \overline{v_z^2}(R,z)
=
\int_z^\infty \nu(R,z') \frac{\partial \Phi_{\rm DM}}{\partial z'}(R,z')\,dz',
$$

以及

$$
\overline{v_\phi^2}(R,z)
=
\frac{1}{1-\beta_z}
\left[
\overline{v_z^2}
+
\frac{R}{\nu}
\frac{\partial\bigl(\nu \overline{v_z^2}\bigr)}{\partial R}
\right]
+
R\frac{\partial \Phi_{\rm DM}}{\partial R}.
$$

再结合

$$
\overline{v_R^2} = \frac{\overline{v_z^2}}{1-\beta_z},
$$

就可以得到三个位于柱坐标系中的内禀二阶矩：

- $\overline{v_R^2}(R,z)$
- $\overline{v_\phi^2}(R,z)$
- $\overline{v_z^2}(R,z)$

这一步是整个计算链条的核心。

## 6. 从内禀二阶矩到视向速度弥散

### 6.1 视线几何变换

对天球平面上的一点 $(x,y)$，沿视线方向记积分变量为 $\ell$。则星系本征柱坐标中的 $(R,z)$ 可由 $(x,y,\ell,i)$ 决定。

常用的几何关系可写为：

$$
R^2 = x^2 + (y\cos i + \ell \sin i)^2,
$$

$$
z = y\sin i - \ell \cos i.
$$

这意味着：对某颗成员星，只要知道它在天球上的投影位置 $(x_i,y_i)$，我们就能对整条视线积分，得到该位置的理论 LOS 二阶矩。

### 6.2 先合成盘面内的二阶矩

在盘面内，与 LOS 投影相关的组合量为

$$
\overline{v_*^2}
=
\overline{v_\phi^2}\frac{x^2}{R^2}
+
\overline{v_R^2}\left(1-\frac{x^2}{R^2}\right).
$$

### 6.3 投影到视线方向

对应的视向二阶矩为

$$
\overline{v_\ell^2}
=
\overline{v_*^2}\sin^2 i
+
\overline{v_z^2}\cos^2 i.
$$

### 6.4 沿视线加权积分

最终理论视向速度弥散为：

$$
\sigma_{\rm los}^2(x,y)
=
\frac{1}{I(x,y)}
\int_{-\infty}^{+\infty}
\nu(R,z)\,\overline{v_\ell^2}(R,z)\,d\ell,
$$

其中 $I(x,y)$ 是表面数密度（或表面亮度）：

$$
I(x,y)
=
\int_{-\infty}^{+\infty} \nu(R,z)\,d\ell.
$$

因此，$\sigma_{\rm los}^2(x,y)$ 是一个**投影后的、沿视线加权平均的二阶矩**。

## 7. 似然函数的定义

对第 $i$ 颗成员星，观测量为：

- 视向速度 $u_i$
- 速度测量误差 $\delta_{v,i}$
- 天球位置 $(x_i, y_i)$

模型先在该位置给出理论值

$$
\sigma_{{\rm los},i}^2
=
\sigma_{\rm los}^2(x_i, y_i \mid \theta).
$$

再把观测误差加进去：

$$
s_i^2 = \sigma_{{\rm los},i}^2 + \delta_{v,i}^2.
$$

作者假设每颗成员星的 LOS 速度分布是以星系系统速度 $\langle u \rangle$ 为均值的高斯分布，因此整组数据的对数似然为

$$
\ln \mathcal{L}(\theta, \langle u \rangle)
=
-\frac{1}{2}
\sum_i
\left[
\frac{(u_i-\langle u\rangle)^2}{s_i^2}
+
\ln\left(2\pi s_i^2\right)
\right].
$$

这里：

- $\theta$ 是动力学模型参数；
- $\langle u \rangle$ 是 nuisance parameter；
- 成员资格不在该似然中重新估计，而是直接沿用原始光谱论文的成员样本。

## 8. 观测数据如何进入计算

这部分对于后续整理 `.csv` 和写代码很重要。

### 8.1 逐星运动学数据

每颗成员星至少需要：

- `x`
- `y`
- `v_los`
- `v_los_err`

它们分别用于：

- `x, y`：决定计算哪一个 $\sigma_{\rm los}(x,y)$
- `v_los`：进入残差项 $(u_i-\langle u\rangle)^2 / s_i^2$
- `v_los_err`：进入总方差 $s_i^2$

### 8.2 星系整体结构参数

每个星系还需要：

- 距离 $D_e$
- 半光半径 $b_*$
- 投影轴比 $q'$

作用分别是：

- $D_e$：把天球角距离换算成物理长度（pc）
- $b_*$：决定 tracer 密度尺度
- $q'$：结合倾角 $i$ 推出内禀轴比 $q$

### 8.3 其他量

如总光度、金属丰度等，多用于论文后续的物理讨论与相关图表，不直接进入上面的逐星速度似然。

## 9. 从参数到似然的计算流程

给定一个星系和一组参数，计算流程可以整理为：

1. 读取该星系的结构参数：$D_e$, $b_*$, $q'$
2. 读取成员星数据：$(x_i, y_i, u_i, \delta_{v,i})$
3. 由 $q'$ 和 $i$ 计算恒星分布内禀轴比 $q$
4. 构造 tracer 密度 $\nu(R,z)$
5. 构造暗晕密度 $\rho_{\rm DM}(R,z \mid \theta)$
6. 由暗晕密度求势梯度 $\partial\Phi/\partial R$ 和 $\partial\Phi/\partial z$
7. 解 Jeans 方程，得到 $\overline{v_z^2}$、$\overline{v_R^2}$、$\overline{v_\phi^2}$
8. 对每颗星的位置 $(x_i,y_i)$：
9. 通过视线几何把 $(x_i,y_i,\ell)$ 映射到 $(R,z)$
10. 沿视线积分，得到 $\sigma_{\rm los}^2(x_i,y_i)$
11. 构造 $s_i^2 = \sigma_{{\rm los},i}^2 + \delta_{v,i}^2$
12. 累加高斯对数似然
13. 与先验相乘，进入 MCMC

## 10. 代码实现建议

为了以后替换任意暗晕模型，建议把代码结构拆成下面几层：

### 10.1 `tracer.py`

负责：

- `q_from_qprime_i(qprime, i)`
- `nu_plummer(R, z, b_star, q)`
- `surface_density_I(x, y, ...)`

### 10.2 `halo_models.py`

负责统一暗晕接口，例如：

```python
rho_dm(R, z, params)
force_R(R, z, params)
force_z(R, z, params)
```

如果某个新模型没有解析力场，就需要额外的势求解器。

### 10.3 `jeans.py`

负责：

- 计算 $\overline{v_z^2}$
- 计算 $\overline{v_R^2}$
- 计算 $\overline{v_\phi^2}$

### 10.4 `projection.py`

负责：

- $(x,y,\ell,i) \to (R,z)$
- 计算 $\overline{v_\ell^2}$
- 沿 LOS 积分得到 $\sigma_{\rm los}^2(x,y)$

### 10.5 `likelihood.py`

负责：

- `log_likelihood(theta, u_mean, galaxy_data)`

### 10.6 `sampler.py`

负责：

- 先验
- 后验
- MCMC 调用

## 11. 后续实现时需要特别注意的点

- 倾角 $i$ 的取值不能让 $q^2 < 0$，这也是论文先验下限与 $q'$ 相联系的原因。
- 逐星位置必须统一到同一物理单位，推荐全部转成 `pc`。
- 速度误差不能省略，必须与理论弥散相加。
- 论文是 unbinned likelihood，不是先把星按半径分箱再拟合。
- 如果以后换新暗晕模型，最稳妥的是保留 `force_R/force_z` 这一层接口。
- 数值积分会是实现中的主要难点，尤其是 LOS 积分与 Jeans 方程积分的稳定性。

## 12. 对我们项目的直接结论

如果我们要复现作者流程并替换暗晕模型，那么代码最小闭环应是：

1. 先为一个星系读入逐星运动学数据；
2. 用论文的 Hernquist 模型实现 `sigma_los(x,y)`；
3. 验证能否得到与论文相容的后验或轮廓；
4. 再把 Hernquist 替换成新的暗晕模型，而不改 likelihood 和 projection 层。

只要这条链条成立，后面换暗晕模型就是一个相对局部的改动。
