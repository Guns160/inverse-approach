# ηc 轴矢量流反问题

**中文** | [English](eta_c_axial_vector.en.md)

## 1. 流、关联函数与衰变常数

定义局域轴矢量流

$$
j_\mu^A(x)=\bar c(x)\gamma_\mu\gamma_5c(x).
$$

它的两点关联函数分解为

$$
\begin{aligned}
\Pi_{\mu\nu}^{AA}(q)
&=i\int d^4x\,e^{iq\cdot x}
\langle0|\operatorname{T}\{j_\mu^A(x)j_\nu^{A\dagger}(0)\}|0\rangle\\
&=(q_\mu q_\nu-q^2g_{\mu\nu})\Pi_1(q^2)
+q_\mu q_\nu\Pi_0(q^2).
\end{aligned}
$$

$\Pi_0$ 含 $0^{-+}$ 态，$\Pi_1$ 含 $1^{++}$ 态。当前主程序对应 LCDA 局域矩流的 $n=0$ 光状投影 [1,2]。取 $z^2=0$ 后，

$$
z^\mu z^\nu\Pi_{\mu\nu}^{AA}
=(z\cdot q)^2I_A(q^2),
\qquad I_A=\Pi_1+\Pi_0.
$$

因此这一不变量同时保留赝标量和轴矢量态。

相应矩阵元定义为

$$
\langle0|j_\mu^A|P_i(p)\rangle=if_{P_i}p_\mu,
\qquad
\langle0|j_\mu^A|A_i(p,\lambda)\rangle
=f_{A_i}m_{A_i}\epsilon_\mu^{(\lambda)}.
$$

对轴矢量态，偏振求和满足

$$
z^\mu z^\nu\sum_\lambda
\epsilon_\mu^{(\lambda)}\epsilon_\nu^{(\lambda)*}
=\frac{(z\cdot p)^2}{m_A^2}.
$$

除去 $(z\cdot q)^2$ 后，赝标量与轴矢量窄态的留数均为对应的 $f_i^2$。主程序保留 $\eta_c(1S)$ 与 $\chi_{c1}(1P)$ 两个极点：

$$
R_1=f_{\eta_c}^2,\qquad R_2=f_{\chi_{c1}}^2.
$$

谱密度 $\rho_A(s)=\operatorname{Im}I_A(s+i0)/\pi$ 无量纲，两个留数的单位为 $\mathrm{GeV}^2$。

## 2. OPE

质量、耦合、凝聚数值和记号见[共用输入与反演方法](common_inputs_and_methods.md)。默认固定 $\mu=2\ \mathrm{GeV}$，由四圈 RG 给出 $m=m_c(\mu)$ 与 $\alpha_s(\mu)$。

### 2.1 微扰谱

取 $v=\sqrt{1-4m^2/s}$，当前 AA 程序保留

$$
\rho_A^{\mathrm{pert}}(s)=\frac{v(3-v^2)}{8\pi^2},
\qquad s\geq4m^2.
$$

该式对应上述 $I_A=\Pi_1+\Pi_0$ 投影。低于夸克阈值时取零。这里没有乘上 VV 通道的 $1+\alpha_s/\pi$ 因子。

### 2.2 凝聚项

令 $u=x(1-x)$、$p_k=x^k+(1-x)^k$、$E=\exp[-m^2/(uM^2)]$，有

$$
C_4^A(M^2)=\frac{G_2}{6\pi}\int_0^1dx\,E
\left[\frac{1}{2M^4}-\frac{m^2p_3}{2M^6u^2}\right],
$$

$$
C_6^A(M^2)=\frac{G_3}{(4\pi)^2}\int_0^1dx\,E
\left[\frac{A_3}{2M^6u^3}
+\frac{B_3^A}{6M^8u^4}
+\frac{C_3}{24M^{10}u^5}\right],
$$

其中

$$
A_3=-\frac{45}{8}up_3-\frac{69}{72}u^2,
\qquad C_3=\frac25m^4up_5,
$$

$$
B_3^A=-\frac34m^2up_4
+m^2u^2\left[-\frac{23}{12}p_2+\frac13u+2\right].
$$

AA 主程序的交叉系数为 $+1/3$，VV 为 $-11/3$。上述公式分别记录两个主程序采用的 OPE，不能因其形式相近而互换。

## 3. 双极点 Borel 反问题

取

$$
\rho_A(s)=R_1\delta(s-m_{\eta_c}^2)
+R_2\delta(s-m_{\chi_{c1}}^2)
+\rho_A^{\mathrm{cont}}(s)\mathbf1_{[s_h,\Lambda]}(s)
+\rho_A^{\mathrm{pert}}(s)\mathbf1_{(\Lambda,\infty)}(s),
$$

其中

$$
m_{\eta_c}=2.9841\ \mathrm{GeV},\qquad
m_{\chi_{c1}}=3.51067\ \mathrm{GeV},\qquad
s_h=14.99\ \mathrm{GeV}^2.
$$

Borel 变换并减去高能微扰尾后，

$$
\begin{aligned}
\frac{1}{M^2}\biggl[
R_1e^{-m_{\eta_c}^2/M^2}
+R_2e^{-m_{\chi_{c1}}^2/M^2}
+\int_{s_h}^{\Lambda}ds\,e^{-s/M^2}\rho_A^{\mathrm{cont}}(s)
\biggr]
={}&\frac{1}{M^2}\int_{4m^2}^{\Lambda}ds\,
e^{-s/M^2}\rho_A^{\mathrm{pert}}(s)\\
&+C_4^A(M^2)+C_6^A(M^2).
\end{aligned}
$$

在给定 $(\Lambda,M_{\min}^2,M_{\max}^2)$ 的窗口中，同时求两个留数和连续谱。连续谱满足

$$
\rho_A^{\mathrm{cont}}(s_h)=0,
\qquad \rho_A^{\mathrm{cont}}(\Lambda)=\rho_A^{\mathrm{pert}}(\Lambda).
$$

## 4. 强子谱

本通道的光状投影同时含 $0^{-+}$ 和 $1^{++}$，因此谱表中需要同时列出 $\eta_c$ 系列和轴矢量结构。下表列出主要实测态及开粲阈值；显式极点采用程序质量，其余质量来源见[共用输入](common_inputs_and_methods.md)。

| 态或通道 | $J^{PC}$ | 质量或阈值能量 / GeV | 谱位置 $s$ / $\mathrm{GeV}^2$ | 在当前模型中的处理 |
|---|---|---:|---:|---|
| $\eta_c(1S)$ | $0^{-+}$ | 2.9841 | 8.9049 | 显式极点 $R_1$ |
| $\chi_{c1}(1P)$ | $1^{++}$ | 3.51067 | 12.3248 | 显式极点 $R_2$ |
| $\eta_c(2S)$ | $0^{-+}$ | 3.6378 | 13.2336 | 位于 $s_h$ 以下，当前未显式纳入 |
| $\chi_{c1}(3872)$ | $1^{++}$ | 3.87164 | 14.9896 | 紧邻 $s_h$，当前未设独立极点 |
| $D^0\bar D^{*0}+\mathrm{c.c.}$ | $0^{-+}$ 或 $1^{++}$ | 3.87170 | 14.9901 | 对应模型的名义开粲起点 |
| $D^+D^{*-}+\mathrm{c.c.}$ | $0^{-+}$ 或 $1^{++}$ | 3.87993 | 15.0539 | 带电开粲阈值 |
| $D^{*0}\bar D^{*0}$ | 可含 $0^{-+}$ | 4.01372 | 16.1099 | 开粲阈值 |
| $D^{*+}D^{*-}$ | 可含 $0^{-+}$ | 4.02054 | 16.1647 | 开粲阈值 |
| $D_s\bar D_s^*+\mathrm{c.c.}$ | $0^{-+}$ 或 $1^{++}$ | 4.08055 | 16.6509 | 开粲阈值 |
| $\chi_{c1}(4140)$ | $1^{++}$ | 4.1465 | 17.1935 | 共振；未设独立极点 |
| $\chi_{c1}(4274)$ | $1^{++}$ | 4.290 | 18.4041 | 共振；未设独立极点 |
| $\chi_{c1}(4685)$ | $1^{++}$ | 4.677 | 21.8743 | 共振；未设独立极点 |

$D\bar D^*$ 的相应电荷共轭组合可在 $^3P_0$ 分波形成 $0^{-+}$，或在 $^3S_1$ 分波形成 $1^{++}$。单列一个二体通道名称并不代表它只有一种总量子数。

$\chi_{c1}(3872)$ 的质量采用 PDG 的 $J/\psi X$ 模式平均值 $3871.64\pm0.06$ MeV，其位置与名义起点极接近。当前程序没有为它单独设置留数。更高 $\eta_c$ 径向激发的模型位置见[赝标量流谱表](eta_c_pseudoscalar.md)；这些态也可进入轴流的 $0^{-+}$ 部分。

AA 的双极点是一项低位态截断：$\eta_c(2S)$ 在 $s_h$ 以下，不能被从 $s_h$ 才开始的自由连续谱显式表示。与 PP 仅保留两个已知赝标量低位态的情形相比，AA 存在更多低位谱结构。解释所得 $f_{\eta_c}$ 与 $f_{\chi_{c1}}$ 时，需要考虑遗漏谱强度向已保留留数的转移。

## 5. 反演与正则化选参

采用[共用方法](common_inputs_and_methods.md)中的帽函数连续谱、Borel 梯形权重、$H^1$ 正则化和 compact GSVD。两个极点不直接受到平滑惩罚，连续谱的固定端点通过完整的交叉项平移消元。

在每个窗口上，使用 Q/MQ/HR 综合规则独立选择 $\alpha_*$；对数网格可向小 $\alpha$ 扩展至 $10^{-80}$。主程序输出两极点、连续谱及其分别对应的 Borel 回代贡献。

对正留数，

$$
f_{\eta_c}=\sqrt{R_1},\qquad
f_{\chi_{c1}}=\sqrt{R_2}.
$$

AA 直接使用轴流矩阵元定义，因而不需要在上述平方根之后再加入 PP 的赝标密度到轴流的换算项。

## 6. 平台

对 $\alpha$、$\Lambda$、$M_{\min}^2$ 和 $M_{\max}^2$ 的稳定性，应同时观察两个衰变常数。由于谱中还存在未显式保留的低位结构，即使得到平滑连续谱和较小回代残差，也需检查两个留数是否随窗口发生明显迁移。

稳定窗口首先说明当前双极点模型下的数值提取稳定；它本身不能证明低位态截断已经充分。判断 AA 结果时，应把这一模型条件与 OPE 分项的收敛性一起考虑。

## 参考文献

1. T. Zhong, X.-G. Wu and T. Huang, *Heavy Pseudoscalar Leading-Twist Distribution Amplitudes within QCD Theory in Background Fields*, Eur. Phys. J. C **75**, 45 (2015). arXiv:1408.2297. DOI: 10.1140/epjc/s10052-015-3271-6.
2. L. Zeng et al., *ηc leading-twist distribution amplitude and the Bc → ηc ℓ ν̄ℓ semileptonic decays using QCD sum rules*, Phys. Rev. D **111** (2025). arXiv:2503.17957. DOI: 10.1103/csrw-9hwp.

RG 与强子质量来源见[共用输入与反演方法](common_inputs_and_methods.md)。
