#!/usr/bin/env python
"""Draw a Chinese overview of the current emcee likelihood acceleration path."""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "outputs" / ".matplotlib"))

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT_ROOT / "outputs/figures/emcee_likelihood_optimization_overview.png"


def main() -> None:
    configure_fonts()
    fig, ax = plt.subplots(figsize=(18.5, 11.2), dpi=180)
    ax.set_xlim(0, 18.5)
    ax.set_ylim(0, 11.2)
    ax.axis("off")
    fig.patch.set_facecolor("#f7f8f5")
    ax.set_facecolor("#f7f8f5")

    ax.text(0.55, 10.72, "当前 emcee 版本：likelihood 计算流程与网格加速思路", fontsize=21, weight="bold", color="#20252b")
    ax.text(
        0.58,
        10.36,
        "每个 walker 的 proposal 都重新计算完整 likelihood；当前加速来自 HaloForceGrid 与 RZMomentGrid，不包含三层参数空间缓存。",
        fontsize=10.8,
        color="#4b5563",
    )

    draw_left_tree(ax)
    draw_bottlenecks(ax)
    draw_grids(ax)
    draw_links(ax)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"wrote {OUTPUT}")


def configure_fonts() -> None:
    candidates = [
        "Arial Unicode MS",
        "PingFang SC",
        "Heiti SC",
        "STHeiti",
        "Songti SC",
        "Noto Sans CJK SC",
        "SimHei",
    ]
    available = {font.name for font in fm.fontManager.ttflist}
    for candidate in candidates:
        if candidate in available:
            plt.rcParams["font.sans-serif"] = [candidate, "DejaVu Sans"]
            break
    plt.rcParams["axes.unicode_minus"] = False


def rounded_box(ax, x, y, w, h, title, body, *, face, edge, title_size=12.2, body_size=8.8):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.035,rounding_size=0.07",
        linewidth=1.25,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(x + 0.16, y + h - 0.26, title, fontsize=title_size, weight="bold", color="#111827", va="top")
    ax.text(x + 0.16, y + h - 0.66, body, fontsize=body_size, color="#172033", va="top", linespacing=1.28)


def draw_left_tree(ax) -> None:
    ax.text(0.65, 9.78, "likelihood 计算树", fontsize=13.5, weight="bold", color="#233238")

    rounded_box(
        ax,
        0.55,
        8.35,
        5.25,
        1.05,
        "1. emcee proposal",
        r"$\theta=(Q,\log b_h,\log\rho_0,q_\beta,\alpha,\beta,\gamma,i,\langle u\rangle)$" "\n"
        r"$q_\beta=-\log_{10}(1-\beta_z)$",
        face="#e8f1f2",
        edge="#36636b",
        body_size=7.9,
    )
    rounded_box(
        ax,
        0.55,
        6.82,
        5.25,
        1.1,
        "2. halo force",
        r"在许多 $(R,z)$ 点计算 $F_R,F_z=-\nabla\Phi_{\rm halo}$" "\n"
        r"旧路径：对 $\tau\in[0,\infty)$ 做数值积分",
        face="#e8f1f2",
        edge="#36636b",
        body_size=8.0,
    )
    rounded_box(
        ax,
        0.55,
        4.95,
        5.25,
        1.52,
        "3. Jeans moments",
        r"$P_z(R,z)=\int_z^\infty \nu(R,z')F_z(R,z')\,dz'$" "\n"
        r"$\overline{v_z^2}=P_z/\nu$" "\n"
        r"$\overline{v_R^2}=\overline{v_z^2}/(1-\beta_z)$；"
        r"$\overline{v_\phi^2}$ 需要 $\partial P_z/\partial R$",
        face="#e8f1f2",
        edge="#36636b",
        body_size=7.75,
    )
    rounded_box(
        ax,
        0.55,
        3.28,
        5.25,
        1.25,
        "4. LOS projection",
        r"对每颗成员星 $(x_i,y_i)$：" "\n"
        r"$\sigma_{\rm los}^2(x_i,y_i|\theta)=\int \nu v_{\rm los}^2d\ell/\int\nu d\ell$",
        face="#e8f1f2",
        edge="#36636b",
        body_size=7.75,
    )
    rounded_box(
        ax,
        0.55,
        1.28,
        5.25,
        1.55,
        "5. Gaussian velocity likelihood",
        r"$s_i^2=\sigma_{\rm los}^2(x_i,y_i|\theta)+\delta v_i^2$" "\n"
        r"$\ln L=-\frac{1}{2}\sum_i[\ln(2\pi s_i^2)+(v_i-\langle u\rangle)^2/s_i^2]$" "\n"
        r"$\ln p(\theta|D)=\ln L+\ln p_{\rm Hayashi}(\theta)$",
        face="#e8f1f2",
        edge="#36636b",
        body_size=7.6,
    )

    for start_y, end_y in [(8.35, 7.92), (6.82, 6.47), (4.95, 4.53), (3.28, 2.83)]:
        arrow(ax, (3.18, start_y), (3.18, end_y), "#36636b", lw=1.55)


def draw_bottlenecks(ax) -> None:
    ax.text(7.05, 9.78, "速度瓶颈", fontsize=13.5, weight="bold", color="#4d3711")
    rounded_box(
        ax,
        6.75,
        6.95,
        3.75,
        1.36,
        "瓶颈 A",
        "halo force 是所有 moments 的底层输入。\n"
        "直接路径会在大量 $(R,z)$ 点反复做\n"
        r"$-\nabla\Phi_{\rm halo}$ 积分，且无穷区间可能慢收敛。",
        face="#fff4db",
        edge="#9a6a1f",
        body_size=7.8,
    )
    rounded_box(
        ax,
        6.75,
        4.35,
        3.75,
        2.0,
        "瓶颈 B",
        "LOS projection 的每个 quadrature 节点\n"
        r"都需要 $\nu$ 与 Jeans moments。" "\n"
        "若逐点直接求解，会反复执行垂向积分、\n"
        r"$\partial P_z/\partial R$ 求导和 moment 组合。" "\n"
        "这会把昂贵的 moment 计算嵌套进 LOS 积分。",
        face="#fff4db",
        edge="#9a6a1f",
        body_size=7.8,
    )


def draw_grids(ax) -> None:
    ax.text(12.4, 9.78, "两个 Grid", fontsize=13.5, weight="bold", color="#1f3a2b")
    rounded_box(
        ax,
        12.05,
        6.66,
        5.75,
        2.02,
        "HaloForceGrid",
        "网格坐标：$(R,|z|)$，通常比 moment grid 更粗\n\n"
        "计算内容：\n"
        r"$F_R(R_j,z_k),\,F_z(R_j,z_k)$" "\n\n"
        "作用：只在粗网格点上直接做 halo-force 积分；\n"
        "随后用 linear interpolation 给 RZMomentGrid 提供 force。",
        face="#e9edff",
        edge="#4d5faf",
        body_size=8.4,
    )
    rounded_box(
        ax,
        12.05,
        3.04,
        5.75,
        2.52,
        "RZMomentGrid",
        "网格坐标：$(R,|z|)$，覆盖所有成员星 LOS path\n\n"
        "计算内容：\n"
        r"$\nu,\ P_z,\ \overline{v_R^2},\ \overline{v_\phi^2},\ \overline{v_z^2}$" "\n"
        r"$\partial P_z/\partial R$ 用 cubic spline derivative" "\n\n"
        "作用：先在规则网格上一次性构造 moments；\n"
        "LOS projection 时只做 interpolation，避免把\n"
        "Jeans moments 的求解嵌套在每条 LOS 积分里。",
        face="#e6f5e8",
        edge="#337443",
        body_size=8.4,
    )
    rounded_box(
        ax,
        12.05,
        0.98,
        5.75,
        1.34,
        "未实践的优化思路：三层快慢参数空间",
        "Layer 1：halo density 参数更新时重建 HaloForceGrid。\n"
        "Layer 2：Jeans / projection 参数更新时重建 RZMomentGrid。\n"
        "Layer 3：仅更新 rho0 与系统速度，复用前两层网格。",
        face="#f0e7f6",
        edge="#76539a",
        body_size=7.9,
    )

    arrow(ax, (14.92, 6.66), (14.92, 5.62), "#4d5faf", lw=1.8)
    ax.text(15.08, 5.93, "插值后的 $F_R,F_z$\n输入 RZMomentGrid", fontsize=8.4, color="#3653a3", va="center")


def draw_links(ax) -> None:
    # Left steps to bottlenecks.
    arrow(ax, (5.83, 7.36), (6.7, 7.62), "#9a6a1f", lw=1.55)
    arrow(ax, (5.83, 5.72), (6.7, 5.45), "#9a6a1f", lw=1.55)
    arrow(ax, (5.83, 3.9), (6.7, 4.7), "#9a6a1f", lw=1.55)

    # Grids solve bottlenecks.
    arrow(ax, (10.55, 7.64), (12.0, 7.64), "#4d5faf", lw=1.75)
    ax.text(11.92, 7.94, "粗网格积分 + interpolation", fontsize=8.2, color="#4d5faf", ha="right")

    arrow(ax, (10.55, 5.05), (12.0, 4.45), "#337443", lw=1.75, connectionstyle="arc3,rad=-0.06")
    ax.text(
        11.28,
        5.22,
        "moments 表格化\nLOS 中插值",
        fontsize=7.8,
        color="#337443",
        ha="center",
        va="center",
    )


def arrow(ax, start, end, color, *, lw=1.7, connectionstyle="arc3,rad=0.0"):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=lw,
            color=color,
            connectionstyle=connectionstyle,
        )
    )


if __name__ == "__main__":
    main()
