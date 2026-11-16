# ηc 赝标量谱中的 Breit–Wigner 共振

**中文** | [English](eta_c_pseudoscalar+BW.en.md)

本文在原有两个低位窄极点和连续谱的基础上，显式加入一个具有有限宽度的 $\eta_c(nS)$ 激发态。流的定义、OPE 和反演方法沿用 [ηc 赝标量流笔记](eta_c_pseudoscalar.md) 与 [共用输入与反演方法](common_inputs_and_methods.md)。以下采用自然单位 $\hbar=c=1$。

记 $s$ 为不变质量平方，$M^2$ 为 Borel 参数，$[s_h,\Lambda]$ 为重建谱的区间。所选态的质量记为 $m_{\mathrm{BW}}=m_{\eta_c(nS)}$，共振谱形及系数分别记为 $\mathrm{BW}(s)$ 和 $R_{\mathrm{BW}}$。势垒半径 $r_{\mathrm{BW}}$ 是另一个独立参数。

## 1. Breit–Wigner 谱形

### 1.1 有限寿命与共振宽度

不稳定粒子的存活概率随时间指数衰减：

$$
P(t)=e^{-\Gamma_0t},\qquad \tau=\frac{1}{\Gamma_0}.
$$

对应振幅为 $\psi(t)\propto e^{-iE_0t}e^{-\Gamma_0t/2}$。作 Fourier 变换，得到

$$
\widetilde\psi(E)\propto\int_0^\infty dt\,e^{i(E-E_0)t-\Gamma_0t/2}
\propto\frac{1}{E-E_0+i\Gamma_0/2}.
$$

因此能量空间中的归一化谱形为

$$
\rho_E(E)=\frac{1}{\pi}\frac{\Gamma_0/2}
{(E-E_0)^2+\Gamma_0^2/4}.
$$

峰的半高全宽为 $\Gamma_0$：寿命越短，宽度越大，峰越宽。

### 1.2 相对论传播子与谱密度

在不变质量平方 $s$ 上，共振附近的常用传播子形式为 [1,4]

$$
D_{\mathrm{BW}}(s)=\frac{1}{m_{\mathrm{BW}}^2-s-i m_{\mathrm{BW}}\Gamma(s)}.
$$

在峰附近取 $E=\sqrt{s}\simeq m_{\mathrm{BW}}+\delta E$，则 $s-m_{\mathrm{BW}}^2\simeq2m_{\mathrm{BW}}\delta E$，可还原上面的 $1/(\delta E+i\Gamma_0/2)$ 形式。$\Gamma(s)$ 是随能量变化的总衰变宽度。

赝标量流两点关联函数的谱表示为

$$
\Pi_P(q^2)=\int ds'\,\frac{\rho_P(s')}{s'-q^2-i0}
+\text{subtractions},\qquad
\rho_P(s)=\frac1\pi\operatorname{Im}\Pi_P(s+i0).
$$

虚部与谱密度的关系可由
$1/(x-i0)=\operatorname{PV}(1/x)+i\pi\delta(x)$ 理解。将共振部分近似为 $D_{\mathrm{BW}}$，得到本文采用的 BW 谱形：

$$
\boxed{
\begin{aligned}
\mathrm{BW}(s)
&=\frac{1}{\pi}\operatorname{Im}
\frac{1}{m_{\mathrm{BW}}^2-s-i m_{\mathrm{BW}}\Gamma(s)}\\
&=\frac{m_{\mathrm{BW}}\Gamma(s)}
{\pi\big[(s-m_{\mathrm{BW}}^2)^2+m_{\mathrm{BW}}^2\Gamma^2(s)\big]}.
\end{aligned}}
$$

这里取传播子的虚部是因为计算对象是**关联函数的谱密度**；振幅模平方虽有相同分母，但不具有同样的分子。在常宽窄宽极限下，

$$
\lim_{\Gamma_0\to0^+}\mathrm{BW}(s)=\delta(s-m_{\mathrm{BW}}^2),
$$

即恢复窄极点。

## 2. 运行宽度与势垒因子

### 2.1 两体相空间和分波阈值

由 Fermi 黄金律，分道衰变宽度取决于跃迁振幅和末态相空间：

$$
\Gamma_\nu(s)=\frac{1}{2\sqrt{s}}\int d\Phi_\nu\,
\overline{|\mathcal M_\nu(s)|^2}.
$$

对于窄末态的双体道 $A B$，质心动量为

$$
\boxed{
q_{AB}(s)=
\begin{cases}
\displaystyle\frac{\sqrt{[s-(m_A+m_B)^2][s-(m_A-m_B)^2]}}{2\sqrt{s}},
&s\geq(m_A+m_B)^2,\\[6pt]
0,&s<(m_A+m_B)^2.
\end{cases}}
$$

两体相空间提供一个 $q$，而轨道角动量为 $\ell$ 的部分波振幅在阈值附近满足 $\mathcal M_\ell\sim q^\ell$，因此

$$
\boxed{\Gamma_\ell(s)\propto q^{2\ell+1}.}
$$

S、P、D 波分别为 $q$、$q^3$、$q^5$。轨道角动量越高，阈值附近的抑制越强。

### 2.2 Blatt–Weisskopf 势垒因子

为描述部分波在较高动量处的有限作用范围效应，引入 Blatt–Weisskopf 因子 [1,4]：

$$
\boxed{
\Phi_{AB}^{(\ell)}(s)
=\frac{q_{AB}(s)^{2\ell+1}}{\sqrt{s}}
B_\ell^2\!\left(q_{AB}(s)r_{\mathrm{BW}}\right).
}
$$

本文采用

$$
B_0^2(z)=1,\qquad
B_1^2(z)=\frac{1}{1+z^2},\qquad
B_2^2(z)=\frac{1}{9+3z^2+z^4},
\qquad z=qr_{\mathrm{BW}}.
$$

其中 $q^{2\ell+1}$ 给出阈值行为，$B_\ell^2(qr_{\mathrm{BW}})$ 调节远离阈值后的动量依赖。例如 P 波的因子是 $q^3/(1+q^2r_{\mathrm{BW}}^2)$；S 波 $B_0^2=1$。

令 $q_0=q_{AB}(m_{\mathrm{BW}}^2)$，窄两体道的运行宽度所需比值为

$$
\boxed{
\frac{\Phi_{AB}^{(\ell)}(s)}{\Phi_{AB}^{(\ell)}(m_{\mathrm{BW}}^2)}
=\left(\frac{q(s)}{q_0}\right)^{2\ell+1}
\frac{m_{\mathrm{BW}}}{\sqrt{s}}
\frac{B_\ell^2(q(s)r_{\mathrm{BW}})}{B_\ell^2(q_0r_{\mathrm{BW}})}.
}
$$

此式与文献 [4] 的归一势垒因子约定对应。

### 2.3 势垒半径的选取

本文以

$$
\boxed{r_{\mathrm{BW}}=4\,\mathrm{GeV}^{-1}\simeq0.79\,\mathrm{fm}}
$$

为基准。它是**现象学作用范围参数**，并非直接测得的 $\eta_c(5S)$ 几何半径。LHCb 的 $B^-\to D^+\pi^-\pi^-$ 振幅分析采用 $4\,\mathrm{GeV}^{-1}$，并在系统误差研究中考察 $3$–$5\,\mathrm{GeV}^{-1}$ [4]。因此取 4 有文献先例，适合作为本模型的默认值。

$r_{\mathrm{BW}}=2,3,5\,\mathrm{GeV}^{-1}$ 同样可以使用。它们分别对应约 $0.39,0.59,0.99\,\mathrm{fm}$。改变 $r_{\mathrm{BW}}$ 会改变 P、D 波的势垒因子，从而改变运行分宽和 BW 谱尾；对 S 波没有直接的势垒因子影响。建议以 $3,5\,\mathrm{GeV}^{-1}$ 作常规敏感性检验，以 $2\,\mathrm{GeV}^{-1}$ 作扩展检验。

更换半径后需重新计算分道相空间、质量迁移及共同归一化，而非只在最终谱形中替换一个参数。

### 2.4 $\eta_c(nS)$ 的分波选择

母态量子数为 $J^{PC}=0^{-+}$。两体末态总自旋 $S$ 与轨道角动量 $\ell$ 耦合到 $J=0$，要求 $\ell=S$；宇称守恒还要求

$$
P_A P_B(-1)^\ell=-1.
$$

对于 5S 使用的八类开粲道，取 [2]：

| 衰变道 | $\ell$ | 分波 |
|---|---:|---|
| $DD^*$ | 1 | P 波 |
| $D_sD_s^*$ | 1 | P 波 |
| $DD_0^*$ | 0 | S 波 |
| $D_sD_{s0}^*$ | 0 | S 波 |
| $DD_2^*$ | 2 | D 波 |
| $D^*D^*$ | 1 | P 波 |
| $D_s^*D_s^*$ | 1 | P 波 |
| $D^*D_1(2420)$ | 0；另取 2 | S 波基准；D 波对照 |

$D^*D_1(2420)$ 可经 $S=\ell=0$ 或 $S=\ell=2$ 形成 $0^-$；文献 [2] 的分宽表未区分这两个分波，故分别计算作比较。表中通道简写包含相应反粒子及电荷共轭组合。

## 3. 分道运行宽度

### 3.1 分宽输入与质量迁移

总宽度写作

$$
\boxed{
\Gamma(s)=\Gamma_{\mathrm{rest}}+\sum_\nu\Gamma_\nu(s),
\qquad \Gamma_{\mathrm{rest}}=\Gamma_{gg}+\Gamma_{\gamma\gamma}.
}
$$

开粲分宽 $\Gamma_\nu^{\mathrm{ref}}$ 使用 Wang、Yan 和 Ping [2] Table 3 的 $\Gamma_b$ 列；$m_{\mathrm{BW}}$、$\Gamma_{gg}$ 和 $\Gamma_{\gamma\gamma}$ 使用 Bokade 和 Bhaghyesh [3]，湮灭宽度取含修正因子的结果。湮灭宽度在此近似为常数。

两组输入使用的母态质量不同。记 [2] 的参考质量为 $m_{\mathrm{ref}}$，先按相空间比迁移分宽：

$$
\boxed{
\widetilde\Gamma_\nu
=\Gamma_\nu^{\mathrm{ref}}
\frac{\Phi_\nu(m_{\mathrm{BW}}^2)}{\Phi_\nu(m_{\mathrm{ref}}^2)}.
}
$$

再保持所选开粲总宽度 $\Gamma_{\mathrm{open}}$，对各道共同缩放：

$$
\boxed{
\Gamma_\nu^0\equiv\Gamma_\nu(m_{\mathrm{BW}}^2)
=\Gamma_{\mathrm{open}}
\frac{\widetilde\Gamma_\nu}{\sum_\kappa\widetilde\Gamma_\kappa}.
}
$$

因此 $\sum_\nu\Gamma_\nu^0=\Gamma_{\mathrm{open}}$，且 $\Gamma_0=\Gamma_{\mathrm{rest}}+\Gamma_{\mathrm{open}}$。

### 3.2 运行宽度

最终每个道按相空间比随 $s$ 运行：

$$
\boxed{
\begin{aligned}
\Gamma_\nu(s)&=\Gamma_\nu^0
\frac{\Phi_\nu(s)}{\Phi_\nu(m_{\mathrm{BW}}^2)},\\
\Gamma(s)&=\Gamma_{\mathrm{rest}}
+\sum_\nu\Gamma_\nu^0
\frac{\Phi_\nu(s)}{\Phi_\nu(m_{\mathrm{BW}}^2)}.
\end{aligned}}
$$

分宽的运行形式参考 [1]；质量迁移与固定 $\Gamma_{\mathrm{open}}$ 的共同缩放是本项目的模型约定。变更势垒半径或分波时，应重新计算 $\Phi_\nu$ 和 $\Gamma_\nu^0$。

### 3.3 电荷通道平均

对同时包含中性和带电末态的道，先对相空间函数等权平均，例如

$$
\Phi_{DD^*}(s)=\frac12\left[
\Phi_{D^0\bar D^{*0}}^{(1)}(s)
+\Phi_{D^+D^{*-}}^{(1)}(s)\right].
$$

随后使用整个 $\Phi_{DD^*}$ 作质量迁移和运行。[2] 的分宽已经包含相应电荷与共轭组合，归一化时使用其合计值。

## 4. 宽末态的谱卷积

### 4.1 卷积相空间

$DD_0^*$ 和 $D^*D_1(2420)$ 中的 $D_0^*$、$D_1$ 本身具有有限宽度。靠近阈值时，其较低质量的谱成分也可能参与衰变，因此不再把末态质量固定为 $m_B$，而是对 $s_B$ 积分 [5]。

采用末态常宽谱分布

$$
\rho_B(s_B)=\frac{m_B\Gamma_B}
{\pi[(s_B-m_B^2)^2+m_B^2\Gamma_B^2]}.
$$

相空间改写为

$$
\boxed{
\Phi_{AB}^{(\ell)}(s)=
\int_{s_{B,\mathrm{th}}}^{(\sqrt{s}-m_A)^2}
 ds_B\,\rho_B(s_B)
\frac{q(s;m_A,\sqrt{s_B})^{2\ell+1}}{\sqrt{s}}
B_\ell^2\!\left(q(s;m_A,\sqrt{s_B})r_{\mathrm{BW}}\right).
}
$$

这里 $q(s;m_A,\sqrt{s_B})$ 是将双体动量公式中的 $m_B$ 替换为 $\sqrt{s_B}$ 所得的动量。子衰变阈值分别为

$$
D_0^*\to D\pi:\quad s_{B,\mathrm{th}}=(m_D+m_\pi)^2,
\qquad
D_1(2420)\to D^*\pi:\quad s_{B,\mathrm{th}}=(m_{D^*}+m_\pi)^2.
$$

当 $\sqrt{s}\leq m_A+\sqrt{s_{B,\mathrm{th}}}$ 时相空间取零。势垒因子随 $s_B$ 变化，故保留在积分内部。本项目固定母态的 $s$，对末态的物理可达质量区间积分，且不重新归一化被截取的末态谱。

### 4.2 子衰变权重

同位旋近似下，$D_0^*\to D\pi$、$D_1\to D^*\pi$ 的带电与中性 $\pi$ 权重为 $2:1$，因此

$$
\Phi_{\nu,k}(s)=\frac23\Phi_{\nu,k}^{(\pi^\pm)}(s)
+\frac13\Phi_{\nu,k}^{(\pi^0)}(s),
\qquad
\Phi_\nu(s)=\frac12\sum_{k=0,+}\Phi_{\nu,k}(s).
$$

各项使用相应的物理阈值。末态质量和中性态宽度取自 PDG [6]；带电 $D_0^*$、$D_1$ 的宽度用对应中性态宽度近似。完成卷积与电荷平均后，再代入第 3 节的分宽公式。

## 5. BW 谱项与 Borel 反问题

### 5.1 共振谱项和 Borel 核

在 $s_h\leq s\leq\Lambda$ 内，显式分离一个激发态：

$$
\boxed{
\rho_P(s)=R_{\mathrm{BW}}\,\mathrm{BW}(s)
+\rho_P^{\mathrm{cont}}(s).
}
$$

两个低位窄极点及 $\Lambda$ 以上的微扰尾不变。$\mathrm{BW}(s)$、$R_{\mathrm{BW}}$、$\rho_P$ 的量纲分别为 $\mathrm{GeV}^{-2}$、$\mathrm{GeV}^{4}$、$\mathrm{GeV}^{2}$。

新增的 Borel 核矩阵列为

$$
\boxed{
K_a^{\mathrm{BW}}=\int_{s_h}^{\Lambda}\frac{ds}{M_a^2}
 e^{-s/M_a^2}\mathrm{BW}(s),
\qquad \Delta g_a=K_a^{\mathrm{BW}}R_{\mathrm{BW}}.
}
$$

它使用与原 $K_{ai}^{\mathrm{pole}}$、$K_{aj}^{\mathrm{cont}}$ 相同的 $1/M_a^2$ 归一化。有限宽度谱在不同 $s$ 处受到不同 Borel 权重；若 $\Lambda<m_{\mathrm{BW}}^2$，积分只包括区间内的低能谱尾。

区间内 BW 的面积为

$$
A_{\mathrm{BW}}(s_h,\Lambda)=\int_{s_h}^{\Lambda}\mathrm{BW}(s)\,ds.
$$

本文不按此面积重归一化。拟合系数仍为 $R_{\mathrm{BW}}$，区间内共振的积分谱强度为 $R_{\mathrm{BW}}A_{\mathrm{BW}}$。

### 5.2 端点匹配

端点条件由 BW 与剩余连续谱之和满足：

$$
\boxed{
\begin{aligned}
\rho_P^{\mathrm{cont}}(s_h)+R_{\mathrm{BW}}\mathrm{BW}(s_h)&=0,\\
\rho_P^{\mathrm{cont}}(\Lambda)+R_{\mathrm{BW}}\mathrm{BW}(\Lambda)
&=\rho_P^{\mathrm{pert}}(\Lambda).
\end{aligned}}
$$

若 $c_j=\rho_P^{\mathrm{cont}}(s_j)$，则

$$
\boxed{
\begin{aligned}
c_0&=-R_{\mathrm{BW}}\mathrm{BW}(s_h),\\
c_N&=\rho_P^{\mathrm{pert}}(\Lambda)-R_{\mathrm{BW}}\mathrm{BW}(\Lambda).
\end{aligned}}
$$

当 $R_{\mathrm{BW}}>0$ 时，低端剩余连续谱允许 $c_0<0$，用于抵消 BW 谱尾，使总谱满足边界条件。沿用原来的平滑惩罚时，$R_{\mathrm{BW}}$ 也通过端点条件与剩余连续谱系数耦合；两个原有窄极点仍为无罚方向。

## 6. 激发态与模型检验

**5S（主要方案）**：使用第 2.4 节八个开粲道；$DD_0^*$ 和 $D^*D_1(2420)$ 做宽末态卷积；后者取 S 波基准、D 波对照。

**4S（备选）**：保留 $DD^*$、$D_sD_s^*$、$D^*D^*$、$D_s^*D_s^*$ 四个 P 波道。当前缺少 $DD_0^*$ 的 4S 耦合输入，暂不纳入。

**3S（对照）**：目前只取 $DD^*$ 道，暂不包含峰高能侧的 $D^*D^*$ 道。其位置靠近已有低位极点，经 Borel 平滑后留数分离较敏感。

**6S**：相关分道输入尚不完整，部分参考分宽与参考质量下的阈值不相容，另有宽末态参数缺失，暂不构造运行宽度。

模型敏感性检验以 $r_{\mathrm{BW}}=4\,\mathrm{GeV}^{-1}$、5S 八道及 $D^*D_1(2420)$ 的 S 波为基准；分别改变势垒半径 $r_{\mathrm{BW}}=3,5$（扩展至 2）和该道的 S／D 波假设。每次从 $\Phi_\nu$ 开始完整重算，比较 $\Gamma(s)$、$\mathrm{BW}(s)$、$K_a^{\mathrm{BW}}$ 和反演得到的 $R_{\mathrm{BW}}$。

## 参考文献

1. D. M. Asner, C. Hanhart and M. Mikhasenko, *Resonances*, in R. L. Workman et al. (Particle Data Group), *Review of Particle Physics*, Prog. Theor. Exp. Phys. **2022**, 083C01 (2022), 2023 update. [官方综述](https://pdg.lbl.gov/2023/reviews/rpp2022-rev-resonances.pdf)。使用第 50.3.1 节及式 (50.30)–(50.35) 的 BW、运行宽度与势垒因子约定。
2. H. Wang, Z. Yan and J. Ping, *Radially excited states of $\eta_c$*, Eur. Phys. J. C **75**, 196 (2015). [DOI](https://doi.org/10.1140/epjc/s10052-015-3418-5)，[arXiv:1412.7068](https://arxiv.org/abs/1412.7068)。Table 2 的参考质量，Table 3 的 $\Gamma_b$ 分宽，式 (13)–(15) 的分波与运动学。
3. C. A. Bokade and Bhaghyesh, *Conventional and XYZ charmonium states in a relativistic screened potential model*, Phys. Rev. D **111**, 014030 (2025). [DOI](https://doi.org/10.1103/PhysRevD.111.014030)，[arXiv:2408.06759](https://arxiv.org/abs/2408.06759)。使用 arXiv v1 的 Table II 质量、Tables XIV、XVI 中含修正因子的 $\Gamma_{\mathrm{cf}}$ 输入。
4. R. Aaij et al. (LHCb Collaboration), *Amplitude analysis of $B^-\to D^+\pi^-\pi^-$ decays*, Phys. Rev. D **94**, 072001 (2016). [DOI](https://doi.org/10.1103/PhysRevD.94.072001)，[arXiv:1608.01289](https://arxiv.org/abs/1608.01289)。式 (12)、(14)、(15) 的势垒因子、BW 振幅、运行宽度；使用 $r_{\mathrm{BW}}=4\,\mathrm{GeV}^{-1}$，并考察 $3$–$5\,\mathrm{GeV}^{-1}$ 的变化。
5. L. Roca, E. Oset and J. Singh, *Low lying axial-vector mesons as dynamically generated resonances*, Phys. Rev. D **72**, 014002 (2005). [DOI](https://doi.org/10.1103/PhysRevD.72.014002)，[arXiv:hep-ph/0503273](https://arxiv.org/abs/hep-ph/0503273)。第 3.3 节式 (23)、(24) 的有限宽末态谱折叠。
6. F. Takahashi et al. (Particle Data Group), *Review of Particle Physics*, Int. J. Mod. Phys. A **41**, 2630011 (2026). [PDG 2026](https://pdg.lbl.gov/2026/)。使用末态质量与宽度。
