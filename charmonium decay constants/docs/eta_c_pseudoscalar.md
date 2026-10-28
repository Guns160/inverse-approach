# ηc 赝标量流反问题

**中文** | [English](eta_c_pseudoscalar.en.md)

## 1. 流、关联函数与衰变常数

定义局域赝标量流

$$
P(x;\mu)=\bar c(x)i\gamma_5c(x),
$$

其中 $\mu$ 表示重整化能标。该流具有 $J^{PC}=0^{-+}$，关联函数与谱密度为

$$
\Pi_P(q^2)=i\int d^4x\,e^{iq\cdot x}
\langle0|\operatorname{T}\{P(x;\mu)P^\dagger(0;\mu)\}|0\rangle,
\qquad
\rho_P(s)=\frac1\pi\operatorname{Im}\Pi_P(s+i0).
$$

对 $\eta_c$ 径向态定义

$$
h_i(\mu)=\langle0|P(\mu)|\eta_c(iS)\rangle,
\qquad R_i(\mu)=|h_i(\mu)|^2.
$$

反问题直接确定 $R_i$。$\rho_P$、$h_i$ 的单位为 $\mathrm{GeV}^2$，$R_i$ 的单位为 $\mathrm{GeV}^4$。

### 1.1 轴反常与衰变常数换算

衰变常数由轴流矩阵元定义：

$$
\langle0|\bar c\gamma_\mu\gamma_5c|\eta_c(iS)(p)\rangle
=if_i p_\mu.
$$

含反常项的 Ward 恒等式为

$$
\partial^\mu(\bar c\gamma_\mu\gamma_5c)
=2m_cP+\frac{\alpha_s}{4\pi}G^a_{\rho\sigma}\widetilde G^{a\rho\sigma}.
$$

取真空到介子态的矩阵元，令

$$
a_i=\left\langle0\left|
\frac{\alpha_s}{4\pi}G^a_{\rho\sigma}\widetilde G^{a\rho\sigma}
\right|\eta_c(iS)\right\rangle,
$$

得到

$$
f_i m_i^2=2m_c(\mu)h_i(\mu)+a_i.
$$

在选定 $h_i=+\sqrt{R_i}$ 且 $a_i$ 采用相同态相位的约定下，

$$
\boxed{f_i=\frac{2m_c(\mu)\sqrt{R_i(\mu)}+a_i}{m_i^2}}.
$$

主程序先计算并保留由留数确定的赝标密度贡献

$$
f_i^{(P)}\equiv\frac{2m_c(\mu)\sqrt{R_i(\mu)}}{m_i^2}.
$$

随后在同一个主程序内使用 $f_i=f_i^{(P)}+a_i/m_i^2$，返回加入已给定反常项的结果，同时保留 $f_i^{(P)}$ 便于核对。该换算用于选中解及整族正则化解。

对 $\eta_c(1S)$，项目采用文献 [2] 的模型输入

$$
a_1=-0.079\ \mathrm{GeV}^3,
\qquad
\frac{a_1}{m_{\eta_c}^2}\simeq-8.872\ \mathrm{MeV}.
$$

这里分母是强子质量平方，和 Borel 参数 $M^2$ 无关。对 $\eta_c(2S)$，默认没有对应的 $a_2$ 输入，保留未经该项修正的 $f_2^{(P)}$ 并标记其定义；这不等同于证明 $a_2=0$。两个反常矩阵元均可显式指定，也可省略某一态的输入以保留其赝标密度贡献。

## 2. OPE

OPE 保留 LO、$O(\alpha_s)$ 微扰修正以及 $D=4,6$ 凝聚项 [1]。共用质量、耦合和凝聚数值见[共用输入与反演方法](common_inputs_and_methods.md)。默认在固定的 $\mu=2\ \mathrm{GeV}$ 使用四圈跑动后的输入。

### 2.1 微扰谱

令 $m=m_c(\mu)$、$v=\sqrt{1-4m^2/s}$，则

$$
\rho_P^{\mathrm{pert}}(s)
=\rho_{\mathrm{LO}}(s)+\rho_{\mathrm{NLO}}(s),
$$

$$
\rho_{\mathrm{LO}}(s)=\frac{3sv}{8\pi^2},
\qquad
\rho_{\mathrm{NLO}}(s)=\frac{s\alpha_s(\mu)}{2\pi^3}F(v).
$$

定义

$$
L(v)=\ln\frac{1+v}{1-v},\qquad
\ell(v)=\ln\frac{1+v}{2}.
$$

有限函数为

$$
\begin{aligned}
F(v)={}&(1+v^2)\biggl[
\frac{\pi^2}{6}+L(v)\ell(v)
+2\operatorname{Li}_2\!\left(\frac{1-v}{1+v}\right)
-4\operatorname{Li}_2(v)+\operatorname{Li}_2(v^2)\\
&\hspace{25mm}
+2\operatorname{Li}_2\!\left(\frac{1+v}{2}\right)
-2\operatorname{Li}_2\!\left(\frac{1-v}{2}\right)\biggr]\\
&+\left(\frac{19}{16}-3v+\frac{v^2}{8}
+\frac{3v^4}{16}\right)L(v)
+\frac{29v}{8}-\frac{3v^2}{8}
+6v\ell(v)-4v\ln v.
\end{aligned}
$$

$\operatorname{Li}_2$ 为二重对数函数。在阈值处取连续极限

$$
\rho_{\mathrm{LO}}(4m^2)=0,
\qquad
\rho_{\mathrm{NLO}}(4m^2)=\frac{4m^2\alpha_s(\mu)}{4\pi}.
$$

### 2.2 凝聚项

对等质量粲夸克，令 $u=x(1-x)$、$E=\exp[-m^2/(uM^2)]$，可将凝聚项写成

$$
C_4^P(M^2)=\frac{G_2}{48\pi}\int_0^1dx\,E
\left[-\frac{2m^2}{M^4u^3}
+\frac{m^4(1-3u)}{M^6u^4}\right],
$$

$$
\begin{aligned}
C_6^P(M^2)=\frac{G_3}{32(4\pi)^2}\int_0^1dx\,E
\biggl\{&
\frac{199-\frac{416}{3}u-\frac{45}{u}}{M^4u^2}\\
&+\frac{m^2}{2M^6u^3}
\left(\frac{39}{u}-\frac{826}{9}-\frac{188}{3}u\right)\\
&+\frac{m^4}{6M^8u^4}
\left(\frac{46}{5u}-\frac{70}{3}-\frac{56}{3}u\right)\\
&-\frac{2m^6(1-5u+5u^2)}{15M^{10}u^6}
\biggr\}.
\end{aligned}
$$

这与主程序中保留 $x^k+(1-x)^k$ 的表达式等价，例如 $p_3=1-3u$、$p_4=1-4u+2u^2$、$p_5=1-5u+5u^2$。这里的简化使各次 Borel 幂次更容易辨认。

以上谱和凝聚项对应 $P=\bar c i\gamma_5c$。若改用质量加权流 $j_5=2m_cP$，谱、留数以及 OPE 每一项都要同时乘以 $4m_c^2$。

## 3. 双极点 Borel 反问题

谱模型为

$$
\rho_P(s)=R_1\delta(s-m_{\eta_c(1S)}^2)
+R_2\delta(s-m_{\eta_c(2S)}^2)
+\rho_P^{\mathrm{cont}}(s)\mathbf1_{[s_h,\Lambda]}(s)
+\rho_P^{\mathrm{pert}}(s)\mathbf1_{(\Lambda,\infty)}(s),
$$

采用

$$
m_1=2.9841\ \mathrm{GeV},\qquad
m_2=3.6378\ \mathrm{GeV},\qquad
s_h=14.99\ \mathrm{GeV}^2.
$$

Borel 变换并减去 $s>\Lambda$ 的微扰尾后，

$$
\begin{aligned}
\frac{1}{M^2}\biggl[
R_1e^{-m_1^2/M^2}+R_2e^{-m_2^2/M^2}
+\int_{s_h}^{\Lambda}ds\,e^{-s/M^2}\rho_P^{\mathrm{cont}}(s)
\biggr]
={}&\frac{1}{M^2}\int_{4m^2}^{\Lambda}ds\,
e^{-s/M^2}\rho_P^{\mathrm{pert}}(s)\\
&+C_4^P(M^2)+C_6^P(M^2).
\end{aligned}
$$

两个 $R_i=|h_i|^2$ 与连续谱为未知量。反常矩阵元 $a_i$ 用于由 $h_i$ 换算轴流衰变常数，不作为额外的项重复加入上述 PP OPE。

连续谱满足

$$
\rho_P^{\mathrm{cont}}(s_h)=0,
\qquad \rho_P^{\mathrm{cont}}(\Lambda)=\rho_P^{\mathrm{pert}}(\Lambda).
$$

## 4. 强子谱

PP 只含 $0^{-+}$ 谱。两个自旋零的 $D\bar D$ 若组成 $J=0$，必须取 $L=0$，此时宇称为正，因此不能形成本通道的 $0^{-+}$。相应 $D\bar D^*$ 电荷共轭组合可在 $^3P_0$ 分波形成 $0^{-+}$，故模型将连续谱起点取在中性 $D\bar D^*$ 阈值附近。

| 态或通道 | 质量或阈值能量 / GeV | 谱位置 $s$ / $\mathrm{GeV}^2$ | 性质与模型处理 |
|---|---:|---:|---|
| $\eta_c(1S)$ | 2.9841 | 8.9049 | 实测态；显式极点 $R_1$ |
| $\eta_c(2S)$ | 3.6378 | 13.2336 | 实测态；显式极点 $R_2$ |
| $D^0\bar D^{*0}+\mathrm{c.c.}$ | 3.87170 | 14.9901 | $^3P_0$ 阈值；对应名义起点 $s_h$ |
| $D^+D^{*-}+\mathrm{c.c.}$ | 3.87993 | 15.0539 | 带电 $^3P_0$ 阈值 |
| $\eta_c(3S)$ | $\sim4.0119$ | $\sim16.0953$ | 势模型预言；未设独立极点 |
| $D^{*0}\bar D^{*0}$ | 4.01372 | 16.1099 | $^3P_0$ 阈值 |
| $D^{*+}D^{*-}$ | 4.02054 | 16.1647 | 带电 $^3P_0$ 阈值 |
| $D_s\bar D_s^*+\mathrm{c.c.}$ | 4.08055 | 16.6509 | $^3P_0$ 阈值 |
| $D_s^{*+}D_s^{*-}$ | 4.2244 | 17.8456 | $^3P_0$ 阈值 |
| $\eta_c(4S)$ | $\sim4.2693$ | $\sim18.2269$ | 势模型预言；未设独立极点 |
| $\eta_c(5S)$ | $\sim4.4508$ | $\sim19.8096$ | 势模型预言；未设独立极点 |
| $\eta_c(6S)$ | $\sim4.5777$ | $\sim20.9553$ | 势模型预言；未设独立极点 |

$\eta_c(3S)$ 至 $\eta_c(6S)$ 的数值取 Bokade–Bhaghyesh 屏蔽势模型的 Table II“ Ours ”列 [3]，用于标示可能的高能谱结构，不能当作已确定的实验质量。主程序未将这些预言质量作为拟合输入。

开粲阈值由[共用质量表](common_inputs_and_methods.md)计算。程序采用固定的名义 $s_h=14.9900\ \mathrm{GeV}^2$，表中的质量和给出约 $14.9901\ \mathrm{GeV}^2$。位于 $[s_h,\Lambda]$ 的更高谱强度由连续函数表达，$\Lambda$ 以上使用微扰尾。

PP 的显式极点不包含 $\chi_{c1}$，因为其 $1^{++}$ 量子数不属于本流。完整 QCD 的赝标谱还可包含混合和多强子贡献；当前模型集中描述所列的两个低位粲偶素极点及开粲后的连续谱。

## 5. 反演与正则化选参

连续谱离散、端点消元和 compact GSVD 使用[共用方法](common_inputs_and_methods.md)。正则项只作用于连续谱，两个留数不直接受罚。广义滤波保留分母 $c_k^2+\alpha d_k^2$。

PP 在 $[10^{-20},1]$ 上的 2001 个对数等距 $\alpha$ 点构造整族解，并以 Q/MQ/HR 综合规则选择 $\alpha_*$。除选中解外，还保留各 $\alpha$ 对应的留数和连续谱，以观察参数依赖。

对正留数，主程序先取 $h_i=\sqrt{R_i}$ 并计算 $f_i^{(P)}$，再按第 1.1 节使用相应的反常输入得到 $f_i$。反常换算不改变已反演的谱和留数。对于固定的 $a_1$ 与强子质量，$f_1$ 相对 $f_1^{(P)}$ 只是一个与 $M^2$、$\Lambda$、$\alpha$ 无关的常数平移。后续绘图直接使用修正后的 $f_i$，不再重复加入这一项。

## 6. 平台

需要同时观察基态和径向激发态的留数。对一个给定窗口，先检查 $\alpha_*$ 附近的变化；再比较邻近 $\Lambda$、$M_{\min}^2$、$M_{\max}^2$ 下的结果。PP 保存的整族正则化解可直接显示候选附近是否存在稳定区。

OPE 中应分开考察 LO、NLO、$C_4^P$ 和 $C_6^P$，避免由不同项的较大抵消掩盖窗口对输入的敏感性。正留数、较小回代残差和稳定平台分别约束不同方面，需结合判断。

## 参考文献

1. T. Zhong, X.-G. Wu and T. Huang, *Heavy Pseudoscalar Leading-Twist Distribution Amplitudes within QCD Theory in Background Fields*, Eur. Phys. J. C **75**, 45 (2015). arXiv:1408.2297. DOI: 10.1140/epjc/s10052-015-3271-6.
2. W. Qin, Q. Zhao and X.-H. Zhong, *Revisiting the pseudoscalar meson and glueball mixing and key issues in the search for a pseudoscalar glueball state*, Phys. Rev. D **97**, 096002 (2018), Eq. (22). arXiv:1712.02550. DOI: 10.1103/PhysRevD.97.096002.
3. C. A. Bokade and Bhaghyesh, *Conventional and XYZ charmonium states in a relativistic screened potential model*, Phys. Rev. D **111**, 014030 (2025), Table II. arXiv:2408.06759. DOI: 10.1103/PhysRevD.111.014030.

RG 与实测强子质量来源见[共用输入与反演方法](common_inputs_and_methods.md)。
