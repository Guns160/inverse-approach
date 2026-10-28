# 共用输入与反演方法

**中文** | [English](common_inputs_and_methods.en.md)

## 1. 物理输入与记号

三个通道分别见 [J/ψ 矢量流](Jpsi_vector.md)、[ηc 轴矢量流](eta_c_axial_vector.md) 和 [ηc 赝标量流](eta_c_pseudoscalar.md)。三个主程序共用下列 OPE 输入与数值方法。

| 输入 | 数值 | 定义 |
|---|---:|---|
| $\overline m_c(\overline m_c)$ | $1.275\ \mathrm{GeV}$ | $\overline{\mathrm{MS}}$ 粲夸克自质量 |
| $\alpha_s^{(5)}(M_Z)$ | $0.1184$ | 五味理论的耦合常数锚点 |
| $M_Z$ | $91.1876\ \mathrm{GeV}$ | 高能锚点 |
| $\overline m_b(\overline m_b)$ | $4.18\ \mathrm{GeV}$ | 底阈值匹配点 |
| $G_2$ | $0.038\ \mathrm{GeV}^4$ | $\langle\alpha_sG^a_{\mu\nu}G^{a\mu\nu}\rangle$ |
| $G_3$ | $0.013\ \mathrm{GeV}^6$ | $\langle g_s^3f^{abc}G^a_{\mu\nu}G^{b\nu}{}_{\rho}G^{c\rho\mu}\rangle$ |
| $\mu$ | $2\ \mathrm{GeV}$ | 默认重整化能标，在整个 Borel 窗口内固定 |

$G_2$、$G_3$ 在当前计算中固定，没有另行实施凝聚项的能标演化。$m_c(\mu)$ 与 $\alpha_s(\mu)$ 由共同的 RG 方程得到。

全文直接使用 Borel 参数 $M^2$。强子质量记为 $m_H$，避免与 Borel 参数混淆。$s$、$s_h$、$\Lambda$ 和 $M^2$ 的单位均为 $\mathrm{GeV}^2$；$\mu$ 与质量的单位为 GeV。谱表中的“能量”是 $m_H$ 或多体阈值的质量和，“谱位置”是其平方。

各通道的等质量 OPE 使用

$$
\bar x=1-x,\qquad u=x\bar x,\qquad
p_k=x^k+\bar x^k,\qquad
E(x,M^2)=\exp\left[-\frac{m_c^2(\mu)}{uM^2}\right].
$$

下面的 $m$ 均表示 $m_c(\mu)$。$C_4$ 与 $C_6$ 分别表示已经完成 Borel 变换的 $D=4$、$D=6$ 凝聚贡献，区别于真空凝聚输入 $G_2$、$G_3$。

## 2. 质量和耦合的跑动

采用 RunDec 的 $\overline{\mathrm{MS}}$ 约定 [1]。令 $a=\alpha_s/\pi$、$t=\ln\mu^2$，则 $N$ 圈方程为

$$
\frac{da}{dt}=-\sum_{k=0}^{N-1}\beta_k a^{k+2},
\qquad
\frac{d\ln m_c}{dt}=-\sum_{k=0}^{N-1}\gamma_k a^{k+1}.
$$

程序允许选择 $N=1,2,3,4$，默认 $N=4$。两组系数均截断到所选阶数。以 $t=\ln\mu^2$ 为自变量时，不再额外乘以 2。

默认输入依次经过以下步骤：

1. 在 $n_f=5$ 理论中，将 $\alpha_s(M_Z)$ 跑到 $\mu=\overline m_b$。
2. 在底夸克自质量处，将耦合由五味理论匹配到四味理论。
3. 在 $n_f=4$ 理论中，将耦合跑到 $\mu=\overline m_c$。
4. 施加 $m_c^{(4)}(\overline m_c)=\overline m_c$，将质量和耦合一起跑到目标能标。

在底夸克自质量处，匹配对数为零。用高能侧耦合 $a_h=\alpha_s^{(5)}/\pi$，三圈匹配写为

$$
\alpha_s^{(4)}=\alpha_s^{(5)}
\left[1+\frac{11}{72}a_h^2+d_3a_h^3\right],
$$

$$
d_3=\frac{564731}{124416}
-\frac{82043}{27648}\zeta(3)
-\frac{2633}{31104}n_l,
\qquad n_l=4.
$$

$N$ 圈跑动配合 $N-1$ 圈匹配。四圈跑动使用上式的全部项；三圈保留 $a_h^2$ 项；一、二圈在该匹配点均取连续耦合。粲夸克保持活跃，没有在粲质量锚点进行四味到三味的解耦。当前输入处方支持 $\overline m_c\leq\mu\leq\overline m_b$，不让粲质量跨越底阈值。

默认四圈输入为

$$
\alpha_s^{(4)}(2\ \mathrm{GeV})=0.3043984044,
\qquad
m_c^{(4)}(2\ \mathrm{GeV})=1.0932768472\ \mathrm{GeV}.
$$

这里上标 $(4)$ 表示四个活跃味。这组数值由有限阶 RG 微分方程直接积分得到。四圈描述的是质量与耦合的跑动阶数；三个通道各自的 OPE 截断阶数见对应 note。

## 3. Borel 核与连续谱离散化

采用

$$
\mathcal B_{M^2}\!\left[\frac{1}{s+Q^2}\right]
=\frac{1}{M^2}e^{-s/M^2},
\qquad Q^2=-q^2>0.
$$

各通道均保留两个固定质量的极点，在 $s_h\leq s\leq\Lambda$ 反演连续谱，并在 $s>\Lambda$ 使用微扰对偶尾。极点、连续谱和微扰积分使用同一个 $1/M^2$ 归一化。

在谱节点 $s_j$ 上用帽函数展开

$$
\rho_{\mathrm{cont}}(s)=\sum_jc_j\phi_j(s),
\qquad \phi_j(s_k)=\delta_{jk}.
$$

对 Borel 节点 $M_a^2$，矩阵的极点列与连续谱列分别为

$$
K_{ai}^{\mathrm{pole}}=
\frac{e^{-m_i^2/M_a^2}}{M_a^2},
\qquad
K_{aj}^{\mathrm{cont}}=
\int_{s_h}^{\Lambda}\frac{ds}{M_a^2}
e^{-s/M_a^2}\phi_j(s).
$$

连续谱列由分段线性基函数的解析积分构造。未知向量包含两个极点留数和谱节点系数，离散方程为 $K\boldsymbol x\simeq\boldsymbol g$。

默认谱间隔为 $\Delta s=0.1\ \mathrm{GeV}^2$，Borel 间隔为 $\Delta M^2=0.01\ \mathrm{GeV}^2$。两种网格均包含实际区间端点，最后一段允许短于默认间隔。

所有通道固定

$$
c_0=\rho_{\mathrm{cont}}(s_h)=0,
\qquad
c_N=\rho_{\mathrm{cont}}(\Lambda)=\rho_{\mathrm{pert}}(\Lambda).
$$

这些边界条件进入反问题的定义。谱在匹配点接上微扰尾是预先施加的条件；判断解的合理性还需要观察靠近边界的谱形及窗口稳定性。

## 4. Tikhonov 正则化与 compact GSVD

### 4.1 正则化泛函

Borel 核会平滑谱函数，因而某些谱形变化只对 OPE 回代产生很小的影响。为抑制这些方向上的振荡，求解

$$
\min_{\boldsymbol x}
\left\|W^{1/2}(K\boldsymbol x-\boldsymbol g)\right\|_2^2
+\alpha\boldsymbol c^{\mathsf T}H\boldsymbol c.
$$

$W$ 是 Borel 网格上的梯形积分权重，使第一项近似 $\int dM^2\,|K\boldsymbol x-g|^2$。这些权重来自积分离散化。

正则矩阵 $H$ 是连续谱帽函数的质量矩阵与刚度矩阵之和，约束谱的大小和一阶变化。以 $\widehat s=s/(1\ \mathrm{GeV}^2)$ 为固定数值坐标，正则项对应

$$
\int d\widehat s\left[
\rho_{\mathrm{cont}}^2+
\left(\frac{d\rho_{\mathrm{cont}}}{d\widehat s}\right)^2
\right].
$$

谱区间没有随 $\Lambda$ 重新缩放到 $[0,1]$。两个极点留数不受直接惩罚。参数以 $\alpha$ 而非 $\alpha^2$ 进入目标函数；其数值依赖谱、坐标和残差的归一化。

### 4.2 固定端点的消元

将谱节点分成内部 $I$ 与固定端点 $F$，则正则项含交叉项 $2\boldsymbol c_I^{\mathsf T}H_{IF}\boldsymbol c_F$。定义

$$
\boldsymbol q=H_{II}^{-1}H_{IF}\boldsymbol c_F,
\qquad \boldsymbol c_I=\boldsymbol w-\boldsymbol q,
$$

得到 $\boldsymbol w^{\mathsf T}H_{II}\boldsymbol w$ 加一个与待求系数无关的常数。把固定端点的谱贡献减到右端，同时加入内部平移产生的项：

$$
\boldsymbol b=
W^{1/2}\left(\boldsymbol g-K_F\boldsymbol c_F+K_I\boldsymbol q\right).
$$

令 $R_H^{\mathsf T}R_H=H_{II}$，则平移后的正则算子为

$$
L=\begin{pmatrix}0_{n_I\times2}&R_H\end{pmatrix}.
$$

这一消元保留了原正则泛函，包括端点与内部谱之间的耦合。

### 4.3 广义滤波与留数提取

对自由未知量对应的加权矩阵 $A=W^{1/2}K_{\mathrm{free}}$ 和 $L$ 做 compact GSVD。写成 $AX=UC$、$LX=VD$，记正则化方向上的对角因子为 $c_k,d_k$，数据投影为 $b_k=\boldsymbol u_k^{\mathsf T}\boldsymbol b$，则

$$
\boldsymbol y_\alpha
=\sum_{k\in\mathrm{reg}}
\frac{c_k b_k}{c_k^2+\alpha d_k^2}\boldsymbol x_k
+\sum_{k\in\mathrm{null}}
\frac{b_k}{c_k}\boldsymbol x_k.
$$

第二项对应 $L$ 的零空间。解出后再将内部谱平移回 $\boldsymbol c_I=\boldsymbol w-\boldsymbol q$，并恢复端点。Python 和 MATLAB 主程序均求解这一广义问题。

两个留数与连续谱都不施加非负约束。只有正留数按各流的定义换算为衰变常数；非正留数保留原值，对应的衰变常数记为未定义。

## 5. Q/MQ/HR 正则化参数选择

令 $\widetilde\sigma_k=c_k/d_k$，并记 $\delta_\perp^2$ 为平移后数据在算子列空间正交补上的平方范数。定义

$$
d_{\mathrm{MD}}^2(\alpha)
=\sum_k\frac{\alpha^3b_k^2}{(\widetilde\sigma_k^2+\alpha)^3}
+\delta_\perp^2,
$$

$$
\psi_Q^2(\alpha)
=\sum_k\frac{\alpha^2\widetilde\sigma_k^2b_k^2}
{(\widetilde\sigma_k^2+\alpha)^4},
\qquad
\psi_{HR}(\alpha)=\frac{d_{\mathrm{MD}}(\alpha)}{\sqrt\alpha}.
$$

$\psi_Q$ 衡量标准形式正则化解随 $\log\alpha$ 的变化，并非几何曲率本身。Q 曲线以 $d_{\mathrm{MD}}$ 和 $\psi_Q$ 为两个诊断量。

在递增网格 $\alpha_1<\cdots<\alpha_N$ 上，取 $p_0=1$，从大 $\alpha$ 向小 $\alpha$ 递推

$$
\psi_{MQ}(\alpha_N)=\psi_Q(\alpha_N),
$$

$$
\psi_{MQ}(\alpha_j)=
\max\left\{\psi_Q(\alpha_j),
\left[\frac{d_{\mathrm{MD}}(\alpha_j)}
{d_{\mathrm{MD}}(\alpha_{j+1})}\right]^{2/3}
\psi_{MQ}(\alpha_{j+1})\right\}.
$$

从 $\psi_Q$ 的局部极小值集合中，分别选出使 $\psi_Q$、$\psi_{MQ}$ 最小的 $\alpha_Q$、$\alpha_{MQ}$。$\alpha_{HR}$ 则使 $\psi_{HR}$ 在 $\alpha\geq\max(\alpha_1,\widetilde\sigma_{\min}^2)$ 的网格区间中最小。令 $r_Q=1/\alpha_Q$ 等，组合为

$$
r_H=\min\{r_Q,\max(r_{MQ},r_{HR})\},
\qquad \alpha_* = 1/r_H.
$$

若 $\psi_Q$ 没有局部极小值，先取 $\psi_{HR}$ 最低的局部极小值；若后者也不存在，再取上述受限区间中的 HR 最小值。主程序记录这类替代选择。

每个 $(\Lambda,M_{\min}^2,M_{\max}^2)$ 窗口独立选取 $\alpha_*$。VV、AA 的对数网格可向小 $\alpha$ 延伸至 $10^{-80}$；PP 默认使用 $[10^{-20},1]$ 上的 2001 个对数等距点。三个通道共享规则，各自保留上述搜索范围。

## 6. 平台判据

主程序给出指定窗口的解、正则化诊断和 OPE 回代。物理上需要进一步判断留数是否在一族邻近窗口中稳定。

首先，在固定 $(\Lambda,M_{\min}^2,M_{\max}^2)$ 时观察 $\alpha_*$ 附近两个留数的变化，判断是否存在有宽度的稳定区。然后改变 $\Lambda$、$M_{\min}^2$、$M_{\max}^2$，考察两个衰变常数是否在相邻窗口中同时稳定，并检查选参是否发生不连续跳变。

小 $M^2$ 增强低位态权重，也可能放大高维凝聚项；较大 $M^2$ 增强高能谱的影响，并降低相近极点的分辨能力。因此窗口应结合微扰项、$C_4$、$C_6$ 的相对大小及谱形确定。小回代残差说明积分方程得到较好满足，但不能单独保证两个留数已经被可靠分离。

## 7. 强子谱数值来源

各通道的显式极点采用主程序中的固定质量。其余实测态与开粲粒子质量取 PDG 2026 的 Summary Table 中 OUR AVERAGE/FIT 中心值 [2]；表内取适当小数位后计算谱位置，不把这些中心值自动作为新增程序输入。

开粲阈值由以下质量计算：

| 粒子 | 质量 / GeV |
|---|---:|
| $D^0$ | 1.86484 |
| $D^+$ | 1.86966 |
| $D^{*0}$ | 2.00686 |
| $D^{*+}$ | 2.01027 |
| $D_s^+$ | 1.96835 |
| $D_s^{*+}$ | 2.1122 |

二体阈值取 $s_{\mathrm{th}}=(m_a+m_b)^2$。表格只列与当前谱模型有关的主要粲偶素、类粲偶素结构及开粲阈值；完整谱还可包含未显式建模的多强子贡献。具有相同量子数只说明耦合允许，不确定局域 $\bar c\Gamma c$ 流对该态的耦合强度。

## 参考文献

1. K. G. Chetyrkin, J. H. Kühn and M. Steinhauser, *RunDec: A Mathematica package for running and decoupling of the strong coupling and quark masses*, Comput. Phys. Commun. **133**, 43–65 (2000). arXiv:hep-ph/0004189. DOI: 10.1016/S0010-4655(00)00155-7.
2. F. Takahashi et al. (Particle Data Group), Int. J. Mod. Phys. A **41**, 2630011 (2026). 本文质量查询使用官方 PDG 2026.0 数据库，数据发布标记为 2026-05-30 13:41:24 PDT。
