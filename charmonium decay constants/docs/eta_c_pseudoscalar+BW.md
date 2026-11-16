# ηc 赝标量谱中 BW 峰的引入

本文说明如何在原来的双极点谱模型中显式加入一个具有有限宽度的激发态。流的定义、OPE、共用输入及反问题方法沿用 [ηc 赝标量流笔记](eta_c_pseudoscalar.md) 和 [共用输入与反演方法](common_inputs_and_methods.md)，这里只讨论新增 BW 谱项及其物理构造。

记号沿用仓库：$M^2$ 为 Borel 参数，$s$ 为谱变量，$s_h$ 与 $\Lambda$ 为连续谱区间的端点；强子质量统一用小写 $m$ 并标明态名。对每次选定的 $\eta_c(nS)$，记 $m_{\mathrm{BW}}=m_{\eta_c(nS)}$，线形记为 $\mathrm{BW}(s)$，系数记为 $R_{\mathrm{BW}}$。剩余连续谱使用 $\rho_P^{\mathrm{cont}}(s)$，其帽函数展开系数仍为 $c_j$。

## 1. 从连续谱中分离一个有限宽度的激发态

### 1.1 新增的谱项

原模型将两个低位态以上的谱强度交给连续谱描述。显式引入一个 $\eta_c(nS)$ 共振后，在 $s_h\leq s\leq\Lambda$ 内将这部分谱写成

$$
\rho_P(s)=R_{\mathrm{BW}}\,\mathrm{BW}(s)
+\rho_P^{\mathrm{cont}}(s),
\qquad s_h\leq s\leq\Lambda.
$$

这里 $\rho_P^{\mathrm{cont}}(s)$ 表示分离该共振后的剩余连续谱。区间以下的两个低位窄极点和 $\Lambda$ 以上的微扰尾保持原来的形式。每次只显式加入一个选定径向态的 BW 峰。

采用运行宽度的 Breit–Wigner 形式。按文献 [1] 第 50.3.1 节式 (50.30) 的分母约定，其谱形可写成

$$
\boxed{
\begin{aligned}
\mathrm{BW}(s)
&=\frac{1}{\pi}\operatorname{Im}
\frac{1}{m_{\mathrm{BW}}^2-s-i m_{\mathrm{BW}}\Gamma(s)}\\
&=\frac{m_{\mathrm{BW}}\Gamma(s)}
{\pi\bigl[(s-m_{\mathrm{BW}}^2)^2+m_{\mathrm{BW}}^2\Gamma^2(s)\bigr]}.
\end{aligned}
}
$$

$m_{\mathrm{BW}}$ 为该态的 BW 质量参数，$\Gamma(s)$ 描述随可用衰变相空间变化的总宽度。在常宽 $\Gamma(s)=\Gamma_0$ 的窄宽极限下，

$$
\lim_{\Gamma_0\to0^+}
\frac{m_{\mathrm{BW}}\Gamma_0}
{\pi\bigl[(s-m_{\mathrm{BW}}^2)^2+m_{\mathrm{BW}}^2\Gamma_0^2\bigr]}
=\delta(s-m_{\mathrm{BW}}^2),
$$

从而恢复一个窄极点的谱表示。有限宽度则将谱强度分布在质量参数附近，并由各衰变道的阈值行为改变峰形。

$\mathrm{BW}(s)$ 的量纲为 $\mathrm{GeV}^{-2}$，其系数 $R_{\mathrm{BW}}$ 的量纲为 $\mathrm{GeV}^{4}$，因此 $R_{\mathrm{BW}}\mathrm{BW}(s)$ 与 $\rho_P^{\mathrm{cont}}(s)$ 均具有原 PP 谱的量纲 $\mathrm{GeV}^{2}$。

### 1.2 BW 对 Borel 关系的贡献

沿用仓库的 Borel 节点 $M_a^2$ 和核矩阵记号，新增的一列为

$$
K_a^{\mathrm{BW}}=
\int_{s_h}^{\Lambda}\frac{ds}{M_a^2}
 e^{-s/M_a^2}\mathrm{BW}(s),
\qquad
\Delta g_a=K_a^{\mathrm{BW}}R_{\mathrm{BW}}.
$$

它与原来的 $K_{ai}^{\mathrm{pole}}$、$K_{aj}^{\mathrm{cont}}$ 使用同一个 $1/M_a^2$ 归一化。这里直接对有限宽度谱积分，保留峰内及谱尾的贡献。若 $\Lambda<m_{\mathrm{BW}}^2$，该项对应区间内保留下来的低能谱尾。

有限区间中的 BW 面积记为

$$
A_{\mathrm{BW}}(s_h,\Lambda)
=\int_{s_h}^{\Lambda}\mathrm{BW}(s)\,ds.
$$

本模型保留原始 BW，不用 $A_{\mathrm{BW}}$ 对其重新归一化。因此，拟合系数为 $R_{\mathrm{BW}}$，区间内该共振的积分谱强度为 $R_{\mathrm{BW}}A_{\mathrm{BW}}$。固定这一区别可避免将积分截断吸收到共振系数的定义中。

### 1.3 总谱的端点匹配

分离共振后，原来的端点条件应由 BW 与剩余连续谱共同满足：

$$
\begin{aligned}
\rho_P^{\mathrm{cont}}(s_h)+R_{\mathrm{BW}}\mathrm{BW}(s_h)&=0,\\
\rho_P^{\mathrm{cont}}(\Lambda)+R_{\mathrm{BW}}\mathrm{BW}(\Lambda)
&=\rho_P^{\mathrm{pert}}(\Lambda).
\end{aligned}
$$

按仓库的 $c_j=\rho_P^{\mathrm{cont}}(s_j)$ 记号，两个端点系数相应变为

$$
\begin{aligned}
c_0&=-R_{\mathrm{BW}}\mathrm{BW}(s_h),\\
c_N&=\rho_P^{\mathrm{pert}}(\Lambda)
-R_{\mathrm{BW}}\mathrm{BW}(\Lambda).
\end{aligned}
$$

这使平滑总谱在低端满足原边界条件，在高端衔接微扰谱。对于正的 $R_{\mathrm{BW}}$，低端剩余谱承担抵消 BW 谱尾的作用，因此允许 $c_0<0$。

沿用对剩余谱系数 $\boldsymbol c$ 的平滑惩罚时，$R_{\mathrm{BW}}$ 会通过这两个端点进入罚项，与剩余谱产生耦合。两个原有窄极点仍对应无罚方向；新增 BW 的系数需要连同上述端点变化一并处理。这是本项目在原 PP 端点约束下加入 BW 的匹配方式。

## 2. 由分道宽度构造运行总宽度

将总宽度写为

$$
\Gamma(s)=\Gamma_{\mathrm{rest}}
+\sum_\nu\Gamma_\nu(s),
\qquad
\Gamma_{\mathrm{rest}}=\Gamma_{gg}+\Gamma_{\gamma\gamma}.
$$

$\nu$ 标记所保留的开粲衰变道，与 Borel 节点指标 $a$ 区分。开粲分宽取自 Wang、Yan 和 Ping [2] 的 Table 3 中 $\Gamma_b$ 一列，该列使用理论母态质量和实验末态质量。选定径向态的质量以及 $gg$、$\gamma\gamma$ 宽度沿用 Bokade 和 Bhaghyesh [3]；其中湮灭宽度取其含修正因子的结果。$gg$ 和 $\gamma\gamma$ 部分在本模型中近似为常数，开粲部分则随阈值、轨道角动量和末态相空间变化。

源表分宽记为 $\Gamma_\nu^{\mathrm{ref}}$，对应的母态参考质量记为 $m_{\mathrm{ref}}$。为了将该输入用于所选质量 $m_{\mathrm{BW}}$，本项目先用相空间函数 $\Phi_\nu$ 作质量迁移：

$$
\widetilde\Gamma_\nu
=\Gamma_\nu^{\mathrm{ref}}
\frac{\Phi_\nu(m_{\mathrm{BW}}^2)}{\Phi_\nu(m_{\mathrm{ref}}^2)}.
$$

随后保持选定的开粲总宽度 $\Gamma_{\mathrm{open}}$，对各道共同缩放：

$$
\Gamma_\nu^0\equiv\Gamma_\nu(m_{\mathrm{BW}}^2)
=\Gamma_{\mathrm{open}}
\frac{\widetilde\Gamma_\nu}{\sum_\kappa\widetilde\Gamma_\kappa}.
$$

运行分宽和总宽度因此为

$$
\boxed{
\begin{aligned}
\Gamma_\nu(s)&=\Gamma_\nu^0
\frac{\Phi_\nu(s)}{\Phi_\nu(m_{\mathrm{BW}}^2)},\\
\Gamma(s)&=\Gamma_{\mathrm{rest}}
+\sum_\nu\Gamma_\nu^0
\frac{\Phi_\nu(s)}{\Phi_\nu(m_{\mathrm{BW}}^2)}.
\end{aligned}
}
$$

在 $s=m_{\mathrm{BW}}^2$ 处，有

$$
\sum_\nu\Gamma_\nu^0=\Gamma_{\mathrm{open}},
\qquad
\Gamma_0\equiv\Gamma(m_{\mathrm{BW}}^2)
=\Gamma_{\mathrm{rest}}+\Gamma_{\mathrm{open}}.
$$

分宽随相空间比运行的形式参照文献 [1] 式 (50.35)；前面的质量迁移与保持开粲总宽度的共同缩放，是本项目连接两组模型输入的约定。改变势垒半径或分波假设时，应重新计算各道权重与共同缩放因子。

## 3. 双体相空间、分波与势垒因子

### 3.1 阈值行为

对于质量分别为 $m_A,m_B$ 的窄末态，质心动量定义为

$$
q_{AB}(s)=
\begin{cases}
\displaystyle
\frac{\sqrt{[s-(m_A+m_B)^2][s-(m_A-m_B)^2]}}{2\sqrt{s}},
& s\geq(m_A+m_B)^2,\\[6pt]
0,&s<(m_A+m_B)^2.
\end{cases}
$$

这与文献 [2] 式 (15) 的双体运动学一致。本文用 $\ell$ 表示衰变道轨道角动量，以区别于共用方法中的正则算子 $L$。相空间函数取为

$$
\Phi_{AB}^{(\ell)}(s)
=\frac{q_{AB}(s)^{2\ell+1}}{\sqrt{s}}
 B_\ell^2\!\left(q_{AB}(s)r_{\mathrm{BW}}\right).
$$

其中一个 $q$ 来自双体相空间，$q^{2\ell}$ 表示分波振幅在阈值附近的角动量抑制。S、P、D 波的阈值幂次分别为 $q$、$q^3$、$q^5$；这一构造对应文献 [1] 式 (50.11)、(50.32)–(50.35) 的相空间与分波因子。

$r_{\mathrm{BW}}$ 表示势垒半径，与共振系数 $R_{\mathrm{BW}}$ 区分。采用的 Blatt–Weisskopf 因子为 [1,4]

$$
B_0^2(z)=1,\qquad
B_1^2(z)=\frac{1}{1+z^2},\qquad
B_2^2(z)=\frac{1}{9+3z^2+z^4},
\qquad z=q r_{\mathrm{BW}}.
$$

上述未归一形式对应文献 [1] 式 (50.34)。对于单一窄双体道，文献 [4] 式 (12) 中在共振质量处归一的振幅因子 $X_\ell$ 满足

$$
X_\ell^2(qr_{\mathrm{BW}})
=\frac{B_\ell^2(qr_{\mathrm{BW}})}
{B_\ell^2(q_0r_{\mathrm{BW}})},
\qquad q_0=q_{AB}(m_{\mathrm{BW}}^2).
$$

将这一关系代入完整的 $\Phi(s)/\Phi(m_{\mathrm{BW}}^2)$，即可得到文献 [4] 式 (15) 中的 $q^{2\ell+1}$、$m_{\mathrm{BW}}/\sqrt{s}$ 和归一势垒因子。本文取 $r_{\mathrm{BW}}=4\ \mathrm{GeV}^{-1}$ 为现象学基准；文献 [4] 的振幅分析也使用这一半径。P、D 波势垒减缓远离阈值后的动量增长。

### 3.2 由量子数选择分波

母态具有 $J^{PC}=0^{-+}$。末态总自旋 $S$ 与轨道角动量 $\ell$ 耦合到 $J=0$ 时需要 $\ell=S$，宇称守恒还要求

$$
P_A P_B(-1)^\ell=-1.
$$

所用角动量耦合与双体衰变分波约定可参照文献 [2] 式 (13)–(15)。对 5S 保留的八个开粲道，分波取值如下。道的集合取自文献 [2] Table 3，表中简写包含相应反粒子及电荷共轭组合。

| 衰变道 | 所取 $\ell$ | 分波 |
|---|---:|---|
| $DD^*$ | 1 | P 波 |
| $D_sD_s^*$ | 1 | P 波 |
| $DD_0^*$ | 0 | S 波 |
| $D_sD_{s0}^*$ | 0 | S 波 |
| $DD_2^*$ | 2 | D 波 |
| $D^*D^*$ | 1 | P 波 |
| $D_s^*D_s^*$ | 1 | P 波 |
| $D^*D_1(2420)$ | 0；另考察 2 | 纯 S 波基准；纯 D 波对照 |

例如，$DD_0^*$ 的末态自旋均为零，取 $\ell=0$；$DD_2^*$ 的总自旋为二，取 $\ell=2$。$D^*D_1(2420)$ 可通过 $S=\ell=0$ 或 $S=\ell=2$ 形成 $0^-$。文献 [2] Table 3 给出该道的合计分宽，未列 S／D 波拆分，故本模型以纯 S 波作为基准，并单独考察纯 D 波假设。

### 3.3 电荷通道的平均

对同时包含中性和带电末态的道，本模型先在相空间层等权平均，例如

$$
\Phi_{DD^*}(s)
=\frac12\left[
\Phi_{D^0\bar D^{*0}}^{(1)}(s)
+\Phi_{D^+D^{*-}}^{(1)}(s)
\right].
$$

随后将这个整体 $\Phi_{DD^*}$ 用于分宽的质量迁移和运行。文献 [2] Table 3 已将相应电荷与共轭组合列入同一分宽，整体归一化使用该合计值。

## 4. 宽末态的谱卷积

### 4.1 将固定末态质量替换为质量谱

$DD_0^*$ 与 $D^*D_1(2420)$ 中的 $D_0^*$、$D_1$ 本身具有有限宽度。尤其在产生阈值附近，可用相空间对其质量分布敏感，因此对末态质量平方作谱卷积。用不稳定粒子的质量谱折叠双体衰变宽度的处理可参照 Roca、Oset 和 Singh [5] 第 3.3 节式 (24)。

沿用 $s$ 表示质量平方的记法，将宽末态 $B$ 的积分变量记为 $s_B$，以区别于 OPE 中的无量纲变量 $u=x(1-x)$。对质量为 $m_A$ 的旁观粒子和宽共振 $B$，采用

$$
\rho_B(s_B)=
\frac{m_B\Gamma_B}
{\pi\bigl[(s_B-m_B^2)^2+m_B^2\Gamma_B^2\bigr]}.
$$

其中 $\Gamma_B$ 为末态共振的常宽输入。本模型将双体相空间替换为

$$
\boxed{
\Phi_{AB}^{(\ell)}(s)=
\int_{s_{B,\mathrm{th}}}^{(\sqrt{s}-m_A)^2}
 ds_B\,\rho_B(s_B)
 \frac{q(s;m_A,\sqrt{s_B})^{2\ell+1}}{\sqrt{s}}
 B_\ell^2\!\left(q(s;m_A,\sqrt{s_B})r_{\mathrm{BW}}\right).
}
$$

$q(s;m_A,\sqrt{s_B})$ 由第 3.1 节的双体动量公式将 $m_B$ 替换为 $\sqrt{s_B}$ 得到。$s_{B,\mathrm{th}}$ 是 $B$ 的物理衰变阈值：$D_0^*\to D\pi$ 对应 $(m_D+m_\pi)^2$，$D_1(2420)\to D^*\pi$ 对应 $(m_{D^*}+m_\pi)^2$。当 $\sqrt{s}\leq m_A+\sqrt{s_{B,\mathrm{th}}}$ 时，该道相空间取零。

文献 [5] 在名义质量附近的有限区间对母态与末态同时折叠；这里固定母态的 $s$，只对宽末态积分，并以物理衰变阈值和运动学上限确定积分区间。卷积保留宽共振低质量部分的贡献，使产生阈值由实际多体末态决定。势垒因子也随 $s_B$ 改变，因此与动量幂次一起放在积分内部。

### 4.2 子衰变权重与电荷平均

在同位旋对称近似下，$D_0^*\to D\pi$ 和 $D_1\to D^*\pi$ 的 $\tfrac12\to\tfrac12\otimes1$ Clebsch–Gordan 系数平方给出带电与中性 $\pi$ 的 $2:1$ 权重。对每个中性或带电产生组，取

$$
\Phi_{\nu,k}(s)=
\frac23\,\Phi_{\nu,k}^{(\pi^{\pm})}(s)
+\frac13\,\Phi_{\nu,k}^{(\pi^0)}(s),
\qquad
\Phi_\nu(s)=\frac12\sum_{k=0,+}\Phi_{\nu,k}(s).
$$

每一项均使用对应的 $D^{(*)}\pi$ 物理阈值。末态质量与中性态宽度输入取自 PDG [6]。中性与带电伙伴使用各自质量；带电 $D_0^*$、$D_1$ 的宽度取对应中性态宽度，作为本模型的同位旋近似。

完成谱卷积和电荷平均后，再将所得整体 $\Phi_\nu$ 代入第 2 节的分宽公式。末态谱在有限可用区间内的积分面积保持原值：只有可达的末态质量部分参与该衰变道。

## 5. 激发态的选择与模型范围

各径向态的质量来源与原 PP 笔记一致，采用文献 [3] 的屏蔽势模型结果；所保留的开粲道以文献 [2] 的分宽表为基础。

**5S 为主要方案。** 使用第 3 节的八道运行宽度，其中 $DD_0^*$ 和 $D^*D_1(2420)$ 包含宽末态卷积；$D^*D_1(2420)$ 以纯 S 波为基准。

**4S 为备选方案。** 保留 $DD^*$、$D_sD_s^*$、$D^*D^*$、$D_s^*D_s^*$ 四道，均取 P 波。$DD_0^*$ 缺少相应 4S 耦合输入，在此四道模型中省略。

**3S 保留作对照。** 当前仅保留 $DD^*$ 道，峰高能侧开启的 $D^*D^*$ 道仍省略。3S 与前两个极点较近，经 Borel 核平滑后各谱成分的响应相近，留数分离较敏感，已有测试中的提取效果较差。

**6S 暂不构造运行宽度。** 其分道输入尚不完整：部分源表分宽与参考质量处的运动学阈值不相容，另有 $DD^*(2S)$、$D^*D_1'$ 的宽末态参数缺失。补齐这些输入后，才能按同一分道方法构造 6S 的 $\Gamma(s)$。

## 参考文献

1. D. M. Asner, C. Hanhart and M. Mikhasenko, *Resonances*, in R. L. Workman et al. (Particle Data Group), *Review of Particle Physics*, Prog. Theor. Exp. Phys. **2022**, 083C01 (2022), 2023 update. [官方综述](https://pdg.lbl.gov/2023/reviews/rpp2022-rev-resonances.pdf)。本文使用修订于 2023 年 8 月的第 50 章：第 50.3.1 节及式 (50.30)–(50.35) 给出 BW 分母、分道运行宽度与 Blatt–Weisskopf 因子。
2. H. Wang, Z. Yan and J. Ping, *Radially excited states of $\eta_c$*, Eur. Phys. J. C **75**, 196 (2015). [DOI: 10.1140/epjc/s10052-015-3418-5](https://doi.org/10.1140/epjc/s10052-015-3418-5). [arXiv:1412.7068](https://arxiv.org/abs/1412.7068)。Table 2 给出参考质量；[Table 3](https://link.springer.com/article/10.1140/epjc/s10052-015-3418-5/tables/3) 的 $\Gamma_b$ 列给出所用开粲分宽；式 (13)–(15) 给出分波衰变及双体动量约定。
3. C. A. Bokade and Bhaghyesh, *Conventional and XYZ charmonium states in a relativistic screened potential model*, Phys. Rev. D **111**, 014030 (2025). [DOI: 10.1103/PhysRevD.111.014030](https://doi.org/10.1103/PhysRevD.111.014030). [arXiv:2408.06759](https://arxiv.org/abs/2408.06759)。输入定位按 [arXiv v1](https://arxiv.org/abs/2408.06759v1)：Table II 的质量，以及 Tables XIV、XVI 中含修正因子的 $\Gamma_{\mathrm{cf}}$ 列，分别对应 $\gamma\gamma$、$gg$ 湮灭宽度。
4. R. Aaij et al. (LHCb Collaboration), *Amplitude analysis of $B^-\to D^+\pi^-\pi^-$ decays*, Phys. Rev. D **94**, 072001 (2016). [DOI: 10.1103/PhysRevD.94.072001](https://doi.org/10.1103/PhysRevD.94.072001). [arXiv:1608.01289](https://arxiv.org/abs/1608.01289)。式 (12)、(14)、(15) 分别给出归一的势垒因子、相对论 BW 振幅与质量依赖宽度，并使用 $r_{\mathrm{BW}}=4\ \mathrm{GeV}^{-1}$。
5. L. Roca, E. Oset and J. Singh, *Low lying axial-vector mesons as dynamically generated resonances*, Phys. Rev. D **72**, 014002 (2005). [DOI: 10.1103/PhysRevD.72.014002](https://doi.org/10.1103/PhysRevD.72.014002). [arXiv:hep-ph/0503273](https://arxiv.org/abs/hep-ph/0503273)。第 3.3 节式 (23)、(24) 展示双体衰变宽度及有限宽度粒子的质量谱卷积；本文采用的末态积分区间与势垒处方见第 4 节。
6. F. Takahashi et al. (Particle Data Group), *Review of Particle Physics*, Int. J. Mod. Phys. A **41**, 2630011 (2026). [PDG 2026 官方入口](https://pdg.lbl.gov/2026/)。用于末态强子质量和宽度，与仓库共用输入笔记采用同一版 PDG。
