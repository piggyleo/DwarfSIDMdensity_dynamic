# MGE 视向速度弥散推导记录

本文档记录当前尝试用于 `likelihood_mode="mge"` 的视向二阶矩计算思路。目标是把 Hayashi 2023 的 axisymmetric Jeans 计算替换为基于 3D Multi-Gaussian Expansion 的解析/半解析计算。

重要说明：本文档记录的是目前正在检查的推导链条。当前解析 MGE 实现相对 `validation-strict` 仍存在很大偏差，因此以下公式也用于人工核对，不应视为已经完成数值验证。

## 1. 坐标与参数

观测平面坐标记为 \((X,Y)\)，星系本征柱坐标记为 \((R,z,\phi)\)。倾角为 \(i\)，其中 \(i=90^\circ\) 为 edge-on。

投影扁率 \(q'\) 与本征 stellar 扁率 \(q_\star\) 的关系为

\[
q'^2=\cos^2 i+q_\star^2\sin^2 i ,
\]

因此

\[
q_\star^2=\frac{q'^2-\cos^2 i}{\sin^2 i}.
\]

暗晕本征扁率记为 \(q_h\)。Hayashi 模型中的速度各向异性参数为

\[
\beta_z = 1-\frac{\overline{v_z^2}}{\overline{v_R^2}},
\]

等价地

\[
\overline{v_R^2}=b\,\overline{v_z^2},\qquad
b=\frac{1}{1-\beta_z}.
\]

在采样代码中直接采样

\[
q_\beta=-\log_{10}(1-\beta_z).
\]

## 2. 直接对 3D density 做 MGE

这里不再先拟合 2D surface density 再转换到 3D density，而是直接将本征 3D tracer density 与本征 3D halo density 展开为 Gaussian 和：

\[
\nu(R,z)
=
\sum_{k=1}^{N_\star}
\nu_k
\exp\left[
-\frac{1}{2s_k^2}
\left(R^2+\frac{z^2}{q_\star^2}\right)
\right],
\]

\[
\rho(R,z)
=
\rho_0
\sum_{j=1}^{N_h}
\rho_j
\exp\left[
-\frac{1}{2a_j^2}
\left(R^2+\frac{z^2}{q_h^2}\right)
\right].
\]

其中 \(\rho_0\) 是 Hayashi halo density normalization，MGE likelihood 中通常先计算 unit-density 结果，即令 \(\rho_0=1\)，最后利用线性缩放：

\[
\sigma_{\rm los}^2(\rho_0)=\rho_0\,\sigma_{\rm los,unit}^2.
\]

stellar tracer 当前采用 Plummer-like 3D density：

\[
\nu(m_\star)\propto
\left(1+\frac{m_\star^2}{b_\star^2}\right)^{-5/2},
\qquad
m_\star^2=R^2+\frac{z^2}{q_\star^2}.
\]

halo 采用 generalized Hernquist-like 3D density：

\[
\rho(m_h)=\rho_0
x^{-\gamma}
\left(1+x^\alpha\right)^{-(\beta-\gamma)/\alpha},
\qquad
x=\frac{m_h}{b_h},
\]

\[
m_h^2=R^2+\frac{z^2}{q_h^2}.
\]

## 3. Shajib/AutoGalaxy 风格的解析 MGE 分解核

对任意 3D radial density \(f(r)\)，在一组对数均匀的 Gaussian 宽度 \(\sigma_l\) 上计算 MGE 振幅。设

\[
\Delta\log_{10}\sigma
=
\frac{\log_{10}\sigma_{\max}-\log_{10}\sigma_{\min}}{N_{\rm G}-1}.
\]

对给定整数 \(p\)，定义复数采样点

\[
\xi_n
=
\left(
\frac{2p\ln 10}{3}+2\pi i n
\right)^{1/2},
\qquad
n=0,1,\dots,2p.
\]

对应权重记为 \(\eta_n\)。当前代码使用 AutoGalaxy 中的构造方式：

\[
\eta_n
=
(-1)^n
2\sqrt{2\pi}\,10^{p/3}\,c_n,
\]

其中 \(c_n\) 是由二项式递推得到的一组系数。于是 3D density 的 Gaussian 振幅写为

\[
A_l
=
\Delta\log_{10}\sigma\,
\sigma_l\,
\operatorname{Re}
\left[
\sum_{n=0}^{2p}
\eta_n\,f(\sigma_l\xi_n)
\right].
\]

端点采用 trapezoid 权重：

\[
A_1\rightarrow \frac{1}{2}A_1,\qquad
A_{N_{\rm G}}\rightarrow \frac{1}{2}A_{N_{\rm G}}.
\]

注意：该解析分解一般会产生 signed amplitudes。后续所有 MGE 求和必须保留振幅符号，不能像 NNLS 正系数拟合那样丢弃负项。

## 4. 投影 tracer surface density

对单个本征 3D Gaussian tracer：

\[
\nu_k(R,z)
=
\nu_k
\exp\left[
-\frac{1}{2s_k^2}
\left(R^2+\frac{z^2}{q_\star^2}\right)
\right],
\]

投影到观测平面后的 Gaussian 扁率为

\[
q_k'^2=\cos^2 i+q_\star^2\sin^2 i.
\]

因为当前 tracer 各 Gaussian 共用同一个 \(q_\star\)，所以 \(q_k'=q'\)。投影 surface density 为

\[
I(X,Y)
=
\sum_k
\sqrt{2\pi}\,
\frac{q_\star s_k}{q'}
\nu_k
\exp\left[
-\frac{1}{2s_k^2}
\left(
X^2+\frac{Y^2}{q'^2}
\right)
\right].
\]

最终 likelihood 需要的是 luminosity-weighted LOS second moment：

\[
\overline{v_{\rm los}^2}(X,Y)
=
\frac{I\,\overline{v_{\rm los}^2}(X,Y)}{I(X,Y)}.
\]

这里的 \(I\,\overline{v_{\rm los}^2}\) 是下面 MGE/JAM 公式给出的 numerator。

## 5. Axisymmetric Jeans 方程

Hayashi 路径使用的本征 Jeans 关系可以写为

\[
\frac{\partial(\nu\overline{v_z^2})}{\partial z}
=
-\nu\frac{\partial\Phi}{\partial z},
\]

\[
\frac{\partial(\nu\overline{v_R^2})}{\partial R}
+
\nu
\left(
\frac{\overline{v_R^2}-\overline{v_\phi^2}}{R}
+\frac{\partial\Phi}{\partial R}
\right)
=0.
\]

令

\[
P_z(R,z)=\nu(R,z)\overline{v_z^2}(R,z),
\]

则

\[
P_z(R,z)
=
\int_z^\infty
\nu(R,z')
\frac{\partial\Phi}{\partial z'}(R,z')\,dz'.
\]

由 \(\overline{v_R^2}=b\overline{v_z^2}\)，可得

\[
\overline{v_\phi^2}
=
b\overline{v_z^2}
+
\frac{bR}{\nu}\frac{\partial P_z}{\partial R}
+
R\frac{\partial\Phi}{\partial R}.
\]

这也是当前 MGE physicality guard 检查局部 \(\overline{v_\phi^2}\) 是否显著为负的依据。

## 6. Gaussian halo 的势梯度

对每个 axisymmetric Gaussian halo 项，势梯度可以写成一个 \(u\in[0,1]\) 的一维积分。当前代码中对应的径向力形式为

\[
\frac{\partial\Phi}{\partial R}(R,z)
=
2\pi G q_h R
\sum_j \rho_j
\int_0^1
\frac{2u^2}{\sqrt{S_h(u)}}
\exp\left[
-\frac{u^2}{2\sigma_{h,j}^2}
\left(
R^2+\frac{z^2}{S_h(u)}
\right)
\right]du,
\]

其中

\[
S_h(u)=1+(q_h^2-1)u^2.
\]

当前代码中有一个需要重点核查的 convention：

\[
\sigma_{h,j}=q_h a_j
\]

被传入了 `mge_los_second_moment()` 的 halo Gaussian 宽度。这里 \(a_j\) 是直接分解 intrinsic density 时记录的 major-axis sigma。该转换是否与上式中的 Gaussian 宽度定义完全一致，是当前解析 MGE 偏差的重点可疑来源之一。

## 7. MGE/JAM 形式的 LOS numerator

对每一对 halo Gaussian \(j\) 与 tracer Gaussian \(k\)，当前实现采用如下 \(u\in[0,1]\) 核函数。记

\[
u_2=u^2,
\qquad
S_h(u)=1-(1-q_h^2)u^2,
\]

\[
a(u)
=
\frac{1}{2}
\left(
\frac{u^2q_h^2}{\sigma_{h,j}^2}
+
\frac{q_\star^2}{\sigma_{t,k}^2}
\right),
\]

\[
b_u(u)
=
\frac{1}{2}
\left[
\frac{1-q_\star^2}{\sigma_{t,k}^2}
+
\frac{q_h^2(1-q_h^2)u^4}
{\sigma_{h,j}^2S_h(u)}
\right],
\]

\[
c
=
1-q_h^2
-
\frac{q_h^2\sigma_{t,k}^2}{\sigma_{h,j}^2},
\]

\[
d(u)
=
1-bq_\star^2
-
\left[
(1-b)c
+
(1-q_h^2)b
\right]u^2.
\]

其中

\[
b=\frac{1}{1-\beta_z}.
\]

当前代码中的 Gaussian 宽度约定为

\[
\sigma_{t,k}=q_\star s_k,
\qquad
\sigma_{h,j}=q_h a_j.
\]

定义

\[
D(u)
=
\left(1-cu^2\right)
\sqrt{
\left[a(u)+b_u(u)\cos^2 i\right]S_h(u)
}.
\]

指数项为

\[
E(u;X,Y)
=
-a(u)
\left[
X^2
+
\frac{a(u)+b_u(u)}
{a(u)+b_u(u)\cos^2 i}
Y^2
\right].
\]

速度因子为

\[
V(u;X)
=
\sigma_{t,k}^2
\left(
\cos^2 i+b\sin^2 i
\right)
+
d(u)X^2\sin^2 i.
\]

于是 numerator 为

\[
I\overline{v_{\rm los}^2}(X,Y)
=
4\pi^{3/2}Gq_h
\sum_j\sum_k
\rho_j\nu_k
\int_0^1
\frac{
u^2 V(u;X)
}{
D(u)
}
\exp[E(u;X,Y)]\,du.
\]

最后

\[
\overline{v_{\rm los}^2}(X,Y)
=
\frac{
I\overline{v_{\rm los}^2}(X,Y)
}{
I(X,Y)
}.
\]

在 unit-density 模式下，上式中的 \(\rho_j\) 不含全局 \(\rho_0\)。实际 likelihood 使用

\[
\sigma_{\rm los}^2(X,Y)
=
10^{\log_{10}\rho_0}
\overline{v_{\rm los,unit}^2}(X,Y).
\]

## 8. Likelihood

对第 \(n\) 颗成员星，观测视向速度为 \(v_n\)，速度误差为 \(\delta v_n\)，系统速度为 \(u_{\rm sys}\)。模型总方差为

\[
s_n^2
=
\sigma_{{\rm los},n}^2+\delta v_n^2.
\]

Gaussian velocity likelihood 为

\[
\ln\mathcal{L}
=
-\frac{1}{2}
\sum_n
\left[
\frac{(v_n-u_{\rm sys})^2}{s_n^2}
+
\ln(2\pi s_n^2)
\right].
\]

因此 `mge` likelihood 的参数依赖大致分为：

- slow shape parameters：\(q_h,b_h,\alpha,\beta,\gamma,\beta_z,i\)，决定 MGE 与 \(\sigma_{\rm los,unit}^2\)；
- density normalization：\(\rho_0\)，线性缩放 \(\sigma_{\rm los,unit}^2\)；
- systemic velocity：\(u_{\rm sys}\)，只进入速度残差。

## 9. Physicality guard

当前保留的非物理解检查逻辑如下：

1. 用同一组 MGE 展开计算局部 tracer density \(\nu(R,z)\)、vertical pressure \(P_z(R,z)\)、\(\partial_R P_z\)、以及 \(\partial_R\Phi\)。
2. 计算

\[
\overline{v_z^2}=\frac{P_z}{\nu},
\qquad
\overline{v_R^2}=\frac{\overline{v_z^2}}{1-\beta_z},
\]

\[
\overline{v_\phi^2}
=
\overline{v_R^2}
+
\frac{R}{\nu}
\frac{\partial_R P_z}{1-\beta_z}
+
R\frac{\partial\Phi}{\partial R}.
\]

3. 如果 LOS 上抽样点的 \(\overline{v_\phi^2}\) 出现显著负值，或最终任意观测星的 \(\sigma_{\rm los,unit}^2\) 非有限/非正，则该参数点返回 \(-\infty\)。

这个 guard 不依赖 \(\rho_0\) 和 \(u_{\rm sys}\)，因为 \(\rho_0\) 只整体线性缩放所有二阶矩，系统速度不影响二阶矩正负性。

## 10. 当前最需要人工核查的点

上一轮固定五组 Eridanus II 样本复算显示，解析 MGE 版本虽然很快，但相对 `validation-strict` 偏差极大：

- 通过解析 MGE guard 的三组样本，\(|\Delta\log L|\) 约为 \(162\) 到 \(268\)；
- \(\sigma_{\rm los}^2\) 的 median relative error 可达到 \(10^2\) 到 \(10^3\) 量级；
- 另外两组样本被解析 MGE guard 判为非物理。

因此后续应优先检查：

1. 3D analytic MGE 分解的振幅归一化，尤其是 \(A_l=\Delta\log\sigma\,\sigma_l\,\mathrm{Re}\sum\eta_n f(\sigma_l\xi_n)\) 中 \(\sigma_l\)、\(\Delta\log\sigma\)、以及 \(2\pi\) 因子的 convention。
2. 直接 intrinsic 3D density 分解后，传入 JAM/MGE second-moment 公式的 Gaussian sigma 到底应是 major-axis sigma、minor-axis sigma，还是 projected convention 下的 sigma。
3. `mge_los_second_moment()` 中

\[
\sigma_{t,k}=q_\star s_k,\qquad
\sigma_{h,j}=q_h a_j
\]

是否与第 7 节公式的定义一致。
4. Projected tracer denominator

\[
I(X,Y)
=
\sum_k
\sqrt{2\pi}\frac{q_\star s_k}{q'}\nu_k
\exp\left[-\frac{X^2+Y^2/q'^2}{2s_k^2}\right]
\]

是否与 numerator 中 tracer Gaussian convention 完全匹配。
5. Signed MGE amplitudes 的振荡是否在当前 \([\sigma_{\min},\sigma_{\max}]\)、\(N_G\)、\(p\) 设置下造成严重 cancellation 或负 density 区域。

