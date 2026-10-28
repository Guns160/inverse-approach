# J/ψ vector-current inverse problem

[中文](Jpsi_vector.md) | **English**

## 1. Current, correlator, and decay constants

Define the local vector current and its vacuum-to-meson matrix element by

$$
j_\mu(x)=\bar c(x)\gamma_\mu c(x),
\qquad
\langle0|j_\mu|V_i(p,\lambda)\rangle
=f_{V_i}m_{V_i}\epsilon_\mu^{(\lambda)}.
$$

This current couples to states with $J^{PC}=1^{--}$. The vector correlator is written as

$$
\Pi_{\mu\nu}(q)
=i\int d^4x\,e^{iq\cdot x}
\langle0|\operatorname{T}\{j_\mu(x)j_\nu^\dagger(0)\}|0\rangle
=(q_\mu q_\nu-q^2g_{\mu\nu})I_V(q^2).
$$

For $z^2=0$, the contraction gives $z^\mu z^\nu\Pi_{\mu\nu}=(z\cdot q)^2I_V$, which projects the local vector correlator. The local moment current for the light-cone distribution amplitude (LCDA) reduces to this correlator at $n=0$, with $\langle\xi^0\rangle=1$ [1].

The spectral density is $\rho_V(s)=\operatorname{Im}I_V(s+i0)/\pi$. A narrow state contributes $R_i\delta(s-m_i^2)$ to the spectrum, where

$$
R_i=f_{V_i}^2.
$$

The main solver extracts the two residues of $J/\psi$ and $\psi(2S)$, determining both vector decay constants simultaneously. The spectral density $\rho_V$ is dimensionless, whereas $R_i$ has units of $\mathrm{GeV}^2$.

## 2. Operator product expansion (OPE)

Masses, coupling, condensate inputs, and notation are specified in [Shared inputs and inversion methods](common_inputs_and_methods.en.md). All Borel points use inputs at fixed $\mu=2\ \mathrm{GeV}$, with four-loop renormalization-group (RG) running by default.

### 2.1 Perturbative spectral density

Let $m=m_c(\mu)$ and $v=\sqrt{1-4m^2/s}$. The implemented $n=0$ perturbative spectral density is

$$
\rho_V^{\mathrm{pert}}(s)
=\frac{v(3-v^2)}{8\pi^2}
\left[1+\frac{\alpha_s(\mu)}{\pi}\right],
\qquad s\geq4m^2.
$$

The density is set to zero below the quark threshold. This expression retains the $1+\alpha_s/\pi$ correction prescription adopted in the study.

### 2.2 Condensate contributions

Using $u=x(1-x)$, $p_k=x^k+(1-x)^k$, and $E=\exp[-m^2/(uM^2)]$, the condensate contributions are

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

where

$$
A_3=-\frac{45}{8}up_3-\frac{69}{72}u^2,
\qquad C_3=\frac25m^4up_5,
$$

$$
B_3^V=-\frac34m^2up_4
+m^2u^2\left[-\frac{23}{12}p_2-\frac{11}{3}u+2\right].
$$

These are the $n=0$ expressions used by the main solver. The vector–vector (VV) channel uses a local vector current. The heavy-pseudoscalar reference [2] used in compiling the condensate terms does not change this physical channel. In particular, the coefficient of $u$ in $B_3^V$ is $-11/3$, distinct from $+1/3$ in the axial–axial (AA) solver.

## 3. Two-pole Borel inverse problem

The spectral model is

$$
\rho_V(s)=R_1\delta(s-m_{J/\psi}^2)
+R_2\delta(s-m_{\psi(2S)}^2)
+\rho_V^{\mathrm{cont}}(s)\mathbf1_{[s_h,\Lambda]}(s)
+\rho_V^{\mathrm{pert}}(s)\mathbf1_{(\Lambda,\infty)}(s).
$$

with

$$
m_{J/\psi}=3.096900\ \mathrm{GeV},\qquad
m_{\psi(2S)}=3.686097\ \mathrm{GeV},
$$

$$
s_h=(2m_{D^0})^2=13.9105153024\ \mathrm{GeV}^2.
$$

After the Borel transform, subtracting the common perturbative tail at $s>\Lambda$ gives

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

The unknowns are $R_1$, $R_2$, and the continuum spectral density. The matching point $\Lambda>s_h$ separates the unknown spectrum from the perturbative tail and has units of $\mathrm{GeV}^2$. The boundary conditions are

$$
\rho_V^{\mathrm{cont}}(s_h)=0,
\qquad \rho_V^{\mathrm{cont}}(\Lambda)=\rho_V^{\mathrm{pert}}(\Lambda).
$$

## 4. Hadronic spectrum

The table is ordered by spectral position. Explicit poles use the solver masses; other resonances use PDG 2026 central values. Two-body thresholds are calculated from the [shared mass table](common_inputs_and_methods.en.md). Multihadron entries refer to the partial waves and charge-conjugation combinations that can form $1^{--}$.

| State or channel | Mass or threshold energy / GeV | Spectral position $s$ / $\mathrm{GeV}^2$ | Treatment in the current model |
|---|---:|---:|---|
| $J/\psi(1S)$ | 3.096900 | 9.5908 | Explicit pole $R_1$ |
| $\psi(2S)$ | 3.686097 | 13.5873 | Explicit pole $R_2$ |
| $D^0\bar D^0$ | 3.72968 | 13.9105 | Continuum onset $s_h$; P wave |
| $D^+D^-$ | 3.73932 | 13.9825 | Charged open-charm threshold; P wave |
| $\psi(3770)$ | 3.7737 | 14.2408 | Resonance; no separate pole |
| $D^0\bar D^{*0}+\mathrm{c.c.}$ | 3.87170 | 14.9901 | Open-charm threshold |
| $D^+D^{*-}+\mathrm{c.c.}$ | 3.87993 | 15.0539 | Open-charm threshold |
| $D_s^+D_s^-$ | 3.93670 | 15.4976 | Open-charm threshold |
| $D^{*0}\bar D^{*0}$ | 4.01372 | 16.1099 | Neutral open-charm threshold |
| $D^{*+}D^{*-}$ | 4.02054 | 16.1647 | Charged open-charm threshold |
| $\psi(4040)$ | 4.040 | 16.3216 | Resonance; no separate pole |
| $D_s\bar D_s^*+\mathrm{c.c.}$ | 4.08055 | 16.6509 | Open-charm threshold |
| $\psi(4160)$ | 4.191 | 17.5645 | Resonance; no separate pole |
| $\psi(4230)$ | 4.2217 | 17.8228 | Resonance; no separate pole |
| $D_s^{*+}D_s^{*-}$ | 4.2244 | 17.8456 | Open-charm threshold |
| $\psi(4360)$ | 4.373 | 19.1231 | Resonance; no separate pole |
| $\psi(4415)$ | 4.415 | 19.4922 | Resonance; no separate pole |
| $\psi(4660)$ | 4.623 | 21.3721 | Resonance; no separate pole |

Both $J/\psi$ and $\psi(2S)$ lie below the adopted open-charm onset and are retained as independent residues. Higher structures within $[s_h,\Lambda]$ contribute through the reconstructed continuum; above $\Lambda$, perturbative duality is assumed. The solver does not fit a separate residue or width for every listed resonance.

The open-charm onset is a spectral modeling assumption. The full vector spectrum can also contain lower hidden-charm multihadron contributions, which the current two-pole model does not describe separately.

## 5. Inversion and regularization-parameter selection

The continuum is expanded in hat basis functions. After fixing the endpoints, compact generalized singular value decomposition (GSVD) is used to solve

$$
\min\left\|W^{1/2}(K\boldsymbol x-\boldsymbol g)\right\|_2^2
+\alpha\boldsymbol c^{\mathsf T}H\boldsymbol c.
$$

Only the continuum receives a smoothness penalty; both residues remain unpenalized. Complete expressions for endpoint elimination, GSVD filtering, and the Q/MQ/HR rule appear in [Sections 3–5 of the shared methods](common_inputs_and_methods.en.md).

The parameter $\alpha_*$ is selected independently in each $(\Lambda,M_{\min}^2,M_{\max}^2)$ window. The logarithmic grid can extend toward smaller $\alpha$, down to $10^{-80}$, and coverage of small-parameter features is recorded. After reconstruction, the two pole contributions and the continuum contribution to the Borel equation are evaluated separately.

For positive residues, the main solver returns

$$
f_{J/\psi}=\sqrt{R_1},\qquad
f_{\psi(2S)}=\sqrt{R_2}.
$$

Residues have units of $\mathrm{GeV}^2$, so decay constants are expressed in GeV. Nonpositive residues are neither replaced by their absolute values nor converted into real decay constants.

## 6. Stability plateaus

The dependence of both $f_{J/\psi}$ and $f_{\psi(2S)}$ on $\alpha$, $\Lambda$, $M_{\min}^2$, and $M_{\max}^2$ must be examined. The exponential kernels of the two poles are similar. A change in one residue can therefore be compensated by the other residue and the continuum. Good agreement with the full Borel input can coexist with unstable residues.

A reportable solution should belong to a family of neighboring stable windows, with both residues positive. The continuum should not develop strong oscillations confined to an exceptionally narrow parameter range. The perturbative contribution, $C_4^V$, and $C_6^V$ should be examined separately to assess sensitivity to higher-dimensional terms.

## References

1. H.-B. Fu et al., *Longitudinal leading-twist distribution amplitude of the J/ψ meson within the background field theory*, Phys. Rev. D **97**, 074025 (2018). arXiv:1801.06832. DOI: 10.1103/PhysRevD.97.074025.
2. T. Zhong, X.-G. Wu and T. Huang, *Heavy Pseudoscalar Leading-Twist Distribution Amplitudes within QCD Theory in Background Fields*, Eur. Phys. J. C **75**, 45 (2015). arXiv:1408.2297. DOI: 10.1140/epjc/s10052-015-3271-6.

RG and hadron-mass sources are listed in [Shared inputs and inversion methods](common_inputs_and_methods.en.md).
