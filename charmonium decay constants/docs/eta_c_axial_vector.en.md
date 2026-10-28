# ηc axial-vector-current inverse problem

[中文](eta_c_axial_vector.md) | **English**

## 1. Current, correlator, and decay constants

Define the local axial-vector current as

$$
j_\mu^A(x)=\bar c(x)\gamma_\mu\gamma_5c(x).
$$

Its two-point correlator is decomposed as

$$
\begin{aligned}
\Pi_{\mu\nu}^{AA}(q)
&=i\int d^4x\,e^{iq\cdot x}
\langle0|\operatorname{T}\{j_\mu^A(x)j_\nu^{A\dagger}(0)\}|0\rangle\\
&=(q_\mu q_\nu-q^2g_{\mu\nu})\Pi_1(q^2)
+q_\mu q_\nu\Pi_0(q^2).
\end{aligned}
$$

The invariant $\Pi_0$ contains $0^{-+}$ states, whereas $\Pi_1$ contains $1^{++}$ states. The main solver uses the $n=0$ lightlike projection of the local moment current for the light-cone distribution amplitude (LCDA) [1,2]. For $z^2=0$,

$$
z^\mu z^\nu\Pi_{\mu\nu}^{AA}
=(z\cdot q)^2I_A(q^2),
\qquad I_A=\Pi_1+\Pi_0.
$$

This invariant therefore retains both pseudoscalar and axial-vector states.

The corresponding vacuum-to-meson matrix elements are defined by

$$
\langle0|j_\mu^A|P_i(p)\rangle=if_{P_i}p_\mu,
\qquad
\langle0|j_\mu^A|A_i(p,\lambda)\rangle
=f_{A_i}m_{A_i}\epsilon_\mu^{(\lambda)}.
$$

For an axial-vector state, the polarization sum satisfies

$$
z^\mu z^\nu\sum_\lambda
\epsilon_\mu^{(\lambda)}\epsilon_\nu^{(\lambda)*}
=\frac{(z\cdot p)^2}{m_A^2}.
$$

After removing $(z\cdot q)^2$, each narrow pseudoscalar or axial-vector state has residue $f_i^2$. The main solver retains the two poles $\eta_c(1S)$ and $\chi_{c1}(1P)$:

$$
R_1=f_{\eta_c}^2,\qquad R_2=f_{\chi_{c1}}^2.
$$

The spectral density $\rho_A(s)=\operatorname{Im}I_A(s+i0)/\pi$ is dimensionless; both residues have units of $\mathrm{GeV}^2$.

## 2. Operator product expansion (OPE)

Masses, coupling, condensate inputs, and notation are specified in [Shared inputs and inversion methods](common_inputs_and_methods.en.md). The default scale is fixed at $\mu=2\ \mathrm{GeV}$, with $m=m_c(\mu)$ and $\alpha_s(\mu)$ obtained from four-loop renormalization-group (RG) running.

### 2.1 Perturbative spectral density

For $v=\sqrt{1-4m^2/s}$, the current axial–axial (AA) solver retains

$$
\rho_A^{\mathrm{pert}}(s)=\frac{v(3-v^2)}{8\pi^2},
\qquad s\geq4m^2.
$$

This expression corresponds to the projection $I_A=\Pi_1+\Pi_0$ defined above. The density vanishes below the quark threshold. The factor $1+\alpha_s/\pi$ used in the vector–vector (VV) channel is not included here.

### 2.2 Condensate contributions

With $u=x(1-x)$, $p_k=x^k+(1-x)^k$, and $E=\exp[-m^2/(uM^2)]$, the condensate contributions are

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

where

$$
A_3=-\frac{45}{8}up_3-\frac{69}{72}u^2,
\qquad C_3=\frac25m^4up_5,
$$

$$
B_3^A=-\frac34m^2up_4
+m^2u^2\left[-\frac{23}{12}p_2+\frac13u+2\right].
$$

The cross-term coefficient is $+1/3$ in the AA solver and $-11/3$ in VV. These expressions specify the OPE used by each solver; their similar forms do not make them interchangeable.

## 3. Two-pole Borel inverse problem

The spectral model is

$$
\rho_A(s)=R_1\delta(s-m_{\eta_c}^2)
+R_2\delta(s-m_{\chi_{c1}}^2)
+\rho_A^{\mathrm{cont}}(s)\mathbf1_{[s_h,\Lambda]}(s)
+\rho_A^{\mathrm{pert}}(s)\mathbf1_{(\Lambda,\infty)}(s),
$$

with

$$
m_{\eta_c}=2.9841\ \mathrm{GeV},\qquad
m_{\chi_{c1}}=3.51067\ \mathrm{GeV},\qquad
s_h=14.99\ \mathrm{GeV}^2.
$$

After the Borel transform and subtraction of the high-energy perturbative tail,

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

Both residues and the continuum are determined simultaneously within a specified $(\Lambda,M_{\min}^2,M_{\max}^2)$ window. The continuum satisfies

$$
\rho_A^{\mathrm{cont}}(s_h)=0,
\qquad \rho_A^{\mathrm{cont}}(\Lambda)=\rho_A^{\mathrm{pert}}(\Lambda).
$$

## 4. Hadronic spectrum

The lightlike projection contains both $0^{-+}$ and $1^{++}$ contributions, so the spectrum includes the $\eta_c$ family and axial-vector structures. The table lists the main observed states and open-charm thresholds. Explicit poles use the solver masses; other mass sources are given in [Shared inputs](common_inputs_and_methods.en.md).

| State or channel | $J^{PC}$ | Mass or threshold energy / GeV | Spectral position $s$ / $\mathrm{GeV}^2$ | Treatment in the current model |
|---|---|---:|---:|---|
| $\eta_c(1S)$ | $0^{-+}$ | 2.9841 | 8.9049 | Explicit pole $R_1$ |
| $\chi_{c1}(1P)$ | $1^{++}$ | 3.51067 | 12.3248 | Explicit pole $R_2$ |
| $\eta_c(2S)$ | $0^{-+}$ | 3.6378 | 13.2336 | Below $s_h$; not explicitly included |
| $\chi_{c1}(3872)$ | $1^{++}$ | 3.87164 | 14.9896 | Very close to $s_h$; no separate pole |
| $D^0\bar D^{*0}+\mathrm{c.c.}$ | $0^{-+}$ or $1^{++}$ | 3.87170 | 14.9901 | Corresponds to the nominal open-charm onset |
| $D^+D^{*-}+\mathrm{c.c.}$ | $0^{-+}$ or $1^{++}$ | 3.87993 | 15.0539 | Charged open-charm threshold |
| $D^{*0}\bar D^{*0}$ | Can contain $0^{-+}$ | 4.01372 | 16.1099 | Open-charm threshold |
| $D^{*+}D^{*-}$ | Can contain $0^{-+}$ | 4.02054 | 16.1647 | Open-charm threshold |
| $D_s\bar D_s^*+\mathrm{c.c.}$ | $0^{-+}$ or $1^{++}$ | 4.08055 | 16.6509 | Open-charm threshold |
| $\chi_{c1}(4140)$ | $1^{++}$ | 4.1465 | 17.1935 | Resonance; no separate pole |
| $\chi_{c1}(4274)$ | $1^{++}$ | 4.290 | 18.4041 | Resonance; no separate pole |
| $\chi_{c1}(4685)$ | $1^{++}$ | 4.677 | 21.8743 | Resonance; no separate pole |

The appropriate charge-conjugation combinations of $D\bar D^*$ can form $0^{-+}$ in the $^3P_0$ partial wave or $1^{++}$ in $^3S_1$. A two-body channel name alone does not specify a unique set of total quantum numbers.

For $\chi_{c1}(3872)$, the PDG average from $J/\psi X$ modes is $3871.64\pm0.06$ MeV, placing it very close to the nominal onset. The current solver assigns no separate residue to this state. Model positions of higher $\eta_c$ radial excitations are listed in the [pseudoscalar spectrum table](eta_c_pseudoscalar.en.md); these states can also enter the $0^{-+}$ axial contribution.

The AA two-pole model truncates the low-lying spectrum. The $\eta_c(2S)$ lies below $s_h$ and cannot be represented explicitly by a free continuum beginning at $s_h$. AA contains more low-lying structures than PP, whose explicit poles are the two known low-lying pseudoscalar states. Interpreting $f_{\eta_c}$ and $f_{\chi_{c1}}$ therefore requires considering how omitted spectral strength may shift the retained residues.

## 5. Inversion and regularization-parameter selection

The solver uses the hat-basis continuum, Borel trapezoidal weights, $H^1$ regularization, and compact generalized singular value decomposition (GSVD) described in the [shared methods](common_inputs_and_methods.en.md). Neither pole is directly penalized. Fixed continuum endpoints are eliminated through the shift that retains all cross terms.

The combined Q/MQ/HR rule selects $\alpha_*$ independently in each window. The logarithmic grid can extend toward smaller $\alpha$, down to $10^{-80}$. The main solver returns both poles, the continuum, and their separate contributions to the Borel equation.

For positive residues,

$$
f_{\eta_c}=\sqrt{R_1},\qquad
f_{\chi_{c1}}=\sqrt{R_2}.
$$

AA uses the axial-current matrix elements directly. No additional PP conversion from pseudoscalar density to axial-current decay constant is required after taking these square roots.

## 6. Stability plateaus

Stability against $\alpha$, $\Lambda$, $M_{\min}^2$, and $M_{\max}^2$ must be assessed for both decay constants. Additional low-lying structures are omitted from the explicit spectrum. Even with a smooth continuum and a small reconstruction residual, substantial migration of either residue across windows must therefore be checked.

A stable window establishes numerical stability within the current two-pole model. It does not by itself establish that the low-lying spectral truncation is adequate. AA results should be assessed under this model condition together with convergence of the individual OPE contributions.

## References

1. T. Zhong, X.-G. Wu and T. Huang, *Heavy Pseudoscalar Leading-Twist Distribution Amplitudes within QCD Theory in Background Fields*, Eur. Phys. J. C **75**, 45 (2015). arXiv:1408.2297. DOI: 10.1140/epjc/s10052-015-3271-6.
2. L. Zeng et al., *ηc leading-twist distribution amplitude and the Bc → ηc ℓ ν̄ℓ semileptonic decays using QCD sum rules*, Phys. Rev. D **111** (2025). arXiv:2503.17957. DOI: 10.1103/csrw-9hwp.

RG and hadron-mass sources are listed in [Shared inputs and inversion methods](common_inputs_and_methods.en.md).
