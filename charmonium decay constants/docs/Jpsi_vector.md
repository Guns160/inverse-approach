# J/ψ 矢量流反问题

**中文** | [English](Jpsi_vector.en.md)

## 1. 流、关联函数与衰变常数

定义局域矢量流

$$
j_\mu(x)=\bar c(x)\gamma_\mu c(x),
\qquad
\langle0|j_\mu|V_i(p,\lambda)\rangle
=f_{V_i}m_{V_i}\epsilon_\mu^{(\lambda)}.
$$

它耦合于 $J^{PC}=1^{--}$ 态。矢量关联函数写为

$$
\Pi_{\mu\nu}(q)
=i\int d^4x\,e^{iq\cdot x}
\langle0|\operatorname{T}\{j_\mu(x)j_\nu^\dagger(0)\}|0\rangle
=(q_\mu q_\nu-q^2g_{\mu\nu})I_V(q^2).
$$

取 $z^2=0$，则 $z^\mu z^\nu\Pi_{\mu\nu}=(z\cdot q)^2I_V$。这是局域矢量流的投影。LCDA 矩流在 $n=0$ 时回到这一关联函数，且 $\langle\xi^0\rangle=1$ [1]。

谱密度为 $\rho_V(s)=\operatorname{Im}I_V(s+i0)/\pi$。窄态在谱中的贡献为 $R_i\delta(s-m_i^2)$，其中

$$
R_i=f_{V_i}^2.
$$

主程序提取 $J/\psi$ 与 $\psi(2S)$ 的两个留数，并据此同时求出两个矢量衰变常数。$\rho_V$ 无量纲，$R_i$ 的单位为 $\mathrm{GeV}^2$。

## 2. OPE

质量、耦合、凝聚数值和记号见[共用输入与反演方法](common_inputs_and_methods.md)。所有 Borel 点使用固定的 $\mu=2\ \mathrm{GeV}$ 输入，默认采用四圈 RG 跑动。

### 2.1 微扰谱

取 $m=m_c(\mu)$、$v=\sqrt{1-4m^2/s}$。当前 $n=0$ 微扰谱为

$$
\rho_V^{\mathrm{pert}}(s)
=\frac{v(3-v^2)}{8\pi^2}
\left[1+\frac{\alpha_s(\mu)}{\pi}\right],
\qquad s\geq4m^2.
$$

低于夸克阈值时取零。这里保留研究中沿用的 $1+\alpha_s/\pi$ 修正处方。

### 2.2 凝聚项

使用 $u=x(1-x)$、$p_k=x^k+(1-x)^k$ 和 $E=\exp[-m^2/(uM^2)]$，有

$$
C_4^V(M^2)=\frac{G_2}{6\pi}\int_0^1dx\,E
\left[\frac{1}{2M^4}-\frac{m^2p_3}{2M^6u^2}\right],
$$

$$
C_6^V(M^2)=\frac{G_3}{(4\pi)^2}\int_0^1dx\,E
\left[\frac{A_3}{2M^6u^3}
+\frac{B_3^V}{6M^8u^4}
+\frac{C_3}{24M^{10}u^5}\right],
$$

其中

$$
A_3=-\frac{45}{8}up_3-\frac{69}{72}u^2,
\qquad C_3=\frac25m^4up_5,
$$

$$
B_3^V=-\frac34m^2up_4
+m^2u^2\left[-\frac{23}{12}p_2-\frac{11}{3}u+2\right].
$$

以上是主程序采用的 $n=0$ 表达式。VV 采用局域矢量流；凝聚项整理所涉及的重赝标介子文献 [2] 不改变其物理通道。特别是，$B_3^V$ 中 $u$ 项系数为 $-11/3$，应与 AA 程序中的 $+1/3$ 区分。

## 3. 双极点 Borel 反问题

谱模型取

$$
\rho_V(s)=R_1\delta(s-m_{J/\psi}^2)
+R_2\delta(s-m_{\psi(2S)}^2)
+\rho_V^{\mathrm{cont}}(s)\mathbf1_{[s_h,\Lambda]}(s)
+\rho_V^{\mathrm{pert}}(s)\mathbf1_{(\Lambda,\infty)}(s).
$$

其中

$$
m_{J/\psi}=3.096900\ \mathrm{GeV},\qquad
m_{\psi(2S)}=3.686097\ \mathrm{GeV},
$$

$$
s_h=(2m_{D^0})^2=13.9105153024\ \mathrm{GeV}^2.
$$

做 Borel 变换，并减去 $s>\Lambda$ 的共同微扰尾：

$$
\begin{aligned}
\frac{1}{M^2}\biggl[
R_1e^{-m_{J/\psi}^2/M^2}
+R_2e^{-m_{\psi(2S)}^2/M^2}
+\int_{s_h}^{\Lambda}ds\,e^{-s/M^2}\rho_V^{\mathrm{cont}}(s)
\biggr]
={}&\frac{1}{M^2}\int_{4m^2}^{\Lambda}ds\,
e^{-s/M^2}\rho_V^{\mathrm{pert}}(s)\\
&+C_4^V(M^2)+C_6^V(M^2).
\end{aligned}
$$

未知量为 $R_1$、$R_2$ 和连续谱函数。$\Lambda>s_h$ 是未知谱与微扰尾的分界，单位为 $\mathrm{GeV}^2$。边界条件为

$$
\rho_V^{\mathrm{cont}}(s_h)=0,
\qquad \rho_V^{\mathrm{cont}}(\Lambda)=\rho_V^{\mathrm{pert}}(\Lambda).
$$

## 4. 强子谱

下表按谱位置排列。显式极点采用程序质量，其余共振采用 PDG 2026 中心值；二体阈值由[共用质量表](common_inputs_and_methods.md)计算。多体通道列出的是可组成 $1^{--}$ 的相应分波和电荷共轭组合。

| 态或通道 | 质量或阈值能量 / GeV | 谱位置 $s$ / $\mathrm{GeV}^2$ | 在当前模型中的处理 |
|---|---:|---:|---|
| $J/\psi(1S)$ | 3.096900 | 9.5908 | 显式极点 $R_1$ |
| $\psi(2S)$ | 3.686097 | 13.5873 | 显式极点 $R_2$ |
| $D^0\bar D^0$ | 3.72968 | 13.9105 | 连续谱起点 $s_h$；P 波 |
| $D^+D^-$ | 3.73932 | 13.9825 | 带电开粲阈值；P 波 |
| $\psi(3770)$ | 3.7737 | 14.2408 | 共振；未设独立极点 |
| $D^0\bar D^{*0}+\mathrm{c.c.}$ | 3.87170 | 14.9901 | 开粲阈值 |
| $D^+D^{*-}+\mathrm{c.c.}$ | 3.87993 | 15.0539 | 开粲阈值 |
| $D_s^+D_s^-$ | 3.93670 | 15.4976 | 开粲阈值 |
| $D^{*0}\bar D^{*0}$ | 4.01372 | 16.1099 | 中性开粲阈值 |
| $D^{*+}D^{*-}$ | 4.02054 | 16.1647 | 带电开粲阈值 |
| $\psi(4040)$ | 4.040 | 16.3216 | 共振；未设独立极点 |
| $D_s\bar D_s^*+\mathrm{c.c.}$ | 4.08055 | 16.6509 | 开粲阈值 |
| $\psi(4160)$ | 4.191 | 17.5645 | 共振；未设独立极点 |
| $\psi(4230)$ | 4.2217 | 17.8228 | 共振；未设独立极点 |
| $D_s^{*+}D_s^{*-}$ | 4.2244 | 17.8456 | 开粲阈值 |
| $\psi(4360)$ | 4.373 | 19.1231 | 共振；未设独立极点 |
| $\psi(4415)$ | 4.415 | 19.4922 | 共振；未设独立极点 |
| $\psi(4660)$ | 4.623 | 21.3721 | 共振；未设独立极点 |

$J/\psi$ 与 $\psi(2S)$ 位于所采用的开粲起点以下，因此作为两个独立留数保留。更高结构若落在 $[s_h,\Lambda]$ 内，其总贡献由待反演连续谱表示；$\Lambda$ 以上使用微扰对偶近似。程序没有为表中每个共振分别拟合留数或宽度。

这个开粲起点是一项谱建模约定。完整矢量谱还可含更低的隐藏粲多强子贡献，当前双极点模型没有单独描述它们。

## 5. 反演与正则化选参

连续谱采用帽函数展开，固定端点后通过 compact GSVD 求解

$$
\min\left\|W^{1/2}(K\boldsymbol x-\boldsymbol g)\right\|_2^2
+\alpha\boldsymbol c^{\mathsf T}H\boldsymbol c.
$$

仅连续谱受到平滑惩罚，两个留数保持自由。端点消元、GSVD 滤波以及 Q/MQ/HR 规则的完整表达式见[共用方法第 3–5 节](common_inputs_and_methods.md)。

每个 $(\Lambda,M_{\min}^2,M_{\max}^2)$ 窗口独立选择 $\alpha_*$。对数网格可向小 $\alpha$ 延伸至 $10^{-80}$，并记录小参数结构的覆盖情况。恢复留数与连续谱后，分别计算两极点和连续谱对 Borel 方程的贡献。

对于正留数，主程序输出

$$
f_{J/\psi}=\sqrt{R_1},\qquad
f_{\psi(2S)}=\sqrt{R_2}.
$$

留数单位为 $\mathrm{GeV}^2$，所以衰变常数单位为 GeV。非正留数不取绝对值，也不换算成实数衰变常数。

## 6. 平台

需要同时观察 $f_{J/\psi}$ 和 $f_{\psi(2S)}$ 对 $\alpha$、$\Lambda$、$M_{\min}^2$、$M_{\max}^2$ 的依赖。两个极点的指数核较接近，单个留数可能通过另一个留数和连续谱进行补偿，所以总回代良好仍可能伴随留数不稳定。

可报告的解应位于一族邻近稳定窗口中，两个留数均为正，且连续谱没有依赖极窄参数区间才出现的剧烈振荡。OPE 的微扰项与 $C_4^V$、$C_6^V$ 应分开观察，以判断窗口对高维项是否过于敏感。

## 参考文献

1. H.-B. Fu et al., *Longitudinal leading-twist distribution amplitude of the J/ψ meson within the background field theory*, Phys. Rev. D **97**, 074025 (2018). arXiv:1801.06832. DOI: 10.1103/PhysRevD.97.074025.
2. T. Zhong, X.-G. Wu and T. Huang, *Heavy Pseudoscalar Leading-Twist Distribution Amplitudes within QCD Theory in Background Fields*, Eur. Phys. J. C **75**, 45 (2015). arXiv:1408.2297. DOI: 10.1140/epjc/s10052-015-3271-6.

RG 与强子质量来源见[共用输入与反演方法](common_inputs_and_methods.md)。
