# ηc pseudoscalar-current inverse problem

[中文](eta_c_pseudoscalar.md) | **English**

## 1. Current, correlator, and decay constants

Define the local pseudoscalar current as

$$
P(x;\mu)=\bar c(x)i\gamma_5c(x),
$$

Here $\mu$ denotes the renormalization scale. The current has $J^{PC}=0^{-+}$, and its correlator and spectral density are

$$
\Pi_P(q^2)=i\int d^4x\,e^{iq\cdot x}
\langle0|\operatorname{T}\{P(x;\mu)P^\dagger(0;\mu)\}|0\rangle,
\qquad
\rho_P(s)=\frac1\pi\operatorname{Im}\Pi_P(s+i0).
$$

For radial $\eta_c$ states, define

$$
h_i(\mu)=\langle0|P(\mu)|\eta_c(iS)\rangle,
\qquad R_i(\mu)=|h_i(\mu)|^2.
$$

The inverse problem directly determines $R_i$. Both $\rho_P$ and $h_i$ have units of $\mathrm{GeV}^2$, whereas $R_i$ has units of $\mathrm{GeV}^4$.

### 1.1 Axial anomaly and decay-constant conversion

The decay constant is defined through the axial-current matrix element:

$$
\langle0|\bar c\gamma_\mu\gamma_5c|\eta_c(iS)(p)\rangle
=if_i p_\mu.
$$

The anomalous axial Ward identity is

$$
\partial^\mu(\bar c\gamma_\mu\gamma_5c)
=2m_cP+\frac{\alpha_s}{4\pi}G^a_{\rho\sigma}\widetilde G^{a\rho\sigma}.
$$

Taking the vacuum-to-meson matrix element and defining

$$
a_i=\left\langle0\left|
\frac{\alpha_s}{4\pi}G^a_{\rho\sigma}\widetilde G^{a\rho\sigma}
\right|\eta_c(iS)\right\rangle,
$$

gives

$$
f_i m_i^2=2m_c(\mu)h_i(\mu)+a_i.
$$

With the convention $h_i=+\sqrt{R_i}$ and a consistent state phase for $a_i$,

$$
\boxed{f_i=\frac{2m_c(\mu)\sqrt{R_i(\mu)}+a_i}{m_i^2}}.
$$

The main solver first computes and retains the pseudoscalar-density contribution determined by the residue:

$$
f_i^{(P)}\equiv\frac{2m_c(\mu)\sqrt{R_i(\mu)}}{m_i^2}.
$$

Within the same solver, it then applies $f_i=f_i^{(P)}+a_i/m_i^2$ using the supplied anomaly matrix elements. The density contribution $f_i^{(P)}$ is retained for comparison. This conversion is applied both to the selected solution and to the full family of regularized solutions.

For $\eta_c(1S)$, the project adopts the model input from Ref. [2]:

$$
a_1=-0.079\ \mathrm{GeV}^3,
\qquad
\frac{a_1}{m_{\eta_c}^2}\simeq-8.872\ \mathrm{MeV}.
$$

The denominator is the hadron mass squared, independent of the Borel parameter $M^2$. For $\eta_c(2S)$, no $a_2$ is supplied by default. Its uncorrected $f_2^{(P)}$ is retained with an explicit label identifying its definition; this does not establish that $a_2=0$. Both anomaly matrix elements can be specified explicitly. Either input can also be omitted to retain only that state's pseudoscalar-density contribution.

## 2. Operator product expansion (OPE)

The OPE retains the leading-order (LO) term, its next-to-leading-order (NLO) $O(\alpha_s)$ correction, and the $D=4,6$ condensate contributions [1]. Shared mass, coupling, and condensate values are given in [Shared inputs and inversion methods](common_inputs_and_methods.en.md). By default, four-loop running inputs are used at fixed $\mu=2\ \mathrm{GeV}$.

### 2.1 Perturbative spectral density

Let $m=m_c(\mu)$ and $v=\sqrt{1-4m^2/s}$. Then

$$
\rho_P^{\mathrm{pert}}(s)
=\rho_{\mathrm{LO}}(s)+\rho_{\mathrm{NLO}}(s),
$$

$$
\rho_{\mathrm{LO}}(s)=\frac{3sv}{8\pi^2},
\qquad
\rho_{\mathrm{NLO}}(s)=\frac{s\alpha_s(\mu)}{2\pi^3}F(v).
$$

Define

$$
L(v)=\ln\frac{1+v}{1-v},\qquad
\ell(v)=\ln\frac{1+v}{2}.
$$

The finite function is

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

Here $\operatorname{Li}_2$ denotes the dilogarithm. The continuous limits at threshold are

$$
\rho_{\mathrm{LO}}(4m^2)=0,
\qquad
\rho_{\mathrm{NLO}}(4m^2)=\frac{4m^2\alpha_s(\mu)}{4\pi}.
$$

### 2.2 Condensate contributions

For equal charm-quark masses, let $u=x(1-x)$ and $E=\exp[-m^2/(uM^2)]$. The condensate contributions can be written as

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

These expressions are equivalent to the solver's form retaining $x^k+(1-x)^k$, using $p_3=1-3u$, $p_4=1-4u+2u^2$, and $p_5=1-5u+5u^2$. The simplified form makes the powers of the Borel parameter easier to identify.

The spectral density and condensate contributions above correspond to $P=\bar c i\gamma_5c$. If the mass-weighted current $j_5=2m_cP$ is used, the spectral density, residues, and every OPE term must all be multiplied by $4m_c^2$.

## 3. Two-pole Borel inverse problem

The spectral model is

$$
\rho_P(s)=R_1\delta(s-m_{\eta_c(1S)}^2)
+R_2\delta(s-m_{\eta_c(2S)}^2)
+\rho_P^{\mathrm{cont}}(s)\mathbf1_{[s_h,\Lambda]}(s)
+\rho_P^{\mathrm{pert}}(s)\mathbf1_{(\Lambda,\infty)}(s),
$$

with

$$
m_1=2.9841\ \mathrm{GeV},\qquad
m_2=3.6378\ \mathrm{GeV},\qquad
s_h=14.99\ \mathrm{GeV}^2.
$$

After the Borel transform and subtraction of the perturbative tail at $s>\Lambda$,

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

The unknowns are the two residues $R_i=|h_i|^2$ and the continuum. The anomaly matrix elements $a_i$ convert $h_i$ into axial-current decay constants. They are not added again as extra terms in the PP OPE above.

The continuum satisfies

$$
\rho_P^{\mathrm{cont}}(s_h)=0,
\qquad \rho_P^{\mathrm{cont}}(\Lambda)=\rho_P^{\mathrm{pert}}(\Lambda).
$$

## 4. Hadronic spectrum

PP contains only the $0^{-+}$ spectrum. Two spin-zero mesons $D\bar D$ can form $J=0$ only with $L=0$, which has positive parity and cannot yield $0^{-+}$. The appropriate charge-conjugation combination of $D\bar D^*$ can form $0^{-+}$ in the $^3P_0$ partial wave. The model therefore places the continuum onset near the neutral $D\bar D^*$ threshold.

| State or channel | Mass or threshold energy / GeV | Spectral position $s$ / $\mathrm{GeV}^2$ | Status and model treatment |
|---|---:|---:|---|
| $\eta_c(1S)$ | 2.9841 | 8.9049 | Observed state; explicit pole $R_1$ |
| $\eta_c(2S)$ | 3.6378 | 13.2336 | Observed state; explicit pole $R_2$ |
| $D^0\bar D^{*0}+\mathrm{c.c.}$ | 3.87170 | 14.9901 | $^3P_0$ threshold; corresponds to the nominal onset $s_h$ |
| $D^+D^{*-}+\mathrm{c.c.}$ | 3.87993 | 15.0539 | Charged $^3P_0$ threshold |
| $\eta_c(3S)$ | $\sim4.0119$ | $\sim16.0953$ | Potential-model prediction; no separate pole |
| $D^{*0}\bar D^{*0}$ | 4.01372 | 16.1099 | $^3P_0$ threshold |
| $D^{*+}D^{*-}$ | 4.02054 | 16.1647 | Charged $^3P_0$ threshold |
| $D_s\bar D_s^*+\mathrm{c.c.}$ | 4.08055 | 16.6509 | $^3P_0$ threshold |
| $D_s^{*+}D_s^{*-}$ | 4.2244 | 17.8456 | $^3P_0$ threshold |
| $\eta_c(4S)$ | $\sim4.2693$ | $\sim18.2269$ | Potential-model prediction; no separate pole |
| $\eta_c(5S)$ | $\sim4.4508$ | $\sim19.8096$ | Potential-model prediction; no separate pole |
| $\eta_c(6S)$ | $\sim4.5777$ | $\sim20.9553$ | Potential-model prediction; no separate pole |

The $\eta_c(3S)$ through $\eta_c(6S)$ values come from the “Ours” column of Table II in the Bokade–Bhaghyesh screened-potential model [3]. They indicate possible higher-energy structures and are not established experimental masses. The main solver does not use these predicted masses as fit inputs.

Open-charm thresholds are calculated from the [shared mass table](common_inputs_and_methods.en.md). The solver uses the fixed nominal value $s_h=14.9900\ \mathrm{GeV}^2$; the listed mass sum gives approximately $14.9901\ \mathrm{GeV}^2$. Higher spectral strength within $[s_h,\Lambda]$ is represented by the continuum function, with the perturbative tail used above $\Lambda$.

The explicit PP poles exclude $\chi_{c1}$ because its $1^{++}$ quantum numbers do not match this current. The full QCD pseudoscalar spectrum can also contain mixing and multihadron contributions. The present model focuses on the two listed low-lying charmonium poles and the continuum above the open-charm onset.

## 5. Inversion and regularization-parameter selection

Continuum discretization, endpoint elimination, and compact generalized singular value decomposition (GSVD) follow the [shared methods](common_inputs_and_methods.en.md). The regularization penalty acts only on the continuum; neither residue is directly penalized. The generalized filter retains the denominator $c_k^2+\alpha d_k^2$.

PP constructs the full solution family on 2001 logarithmically spaced $\alpha$ points over $[10^{-20},1]$, then selects $\alpha_*$ using the combined Q/MQ/HR rule. In addition to the selected solution, residues and continuum spectra are retained at every $\alpha$ to examine parameter dependence.

For positive residues, the main solver takes $h_i=\sqrt{R_i}$, computes $f_i^{(P)}$, and uses the anomaly inputs in Section 1.1 to obtain $f_i$. This conversion leaves the reconstructed spectrum and residues unchanged. For fixed $a_1$ and hadron mass, $f_1$ differs from $f_1^{(P)}$ by a constant independent of $M^2$, $\Lambda$, and $\alpha$. Subsequent plots use the corrected $f_i$ directly, without adding the same term again.

## 6. Stability plateaus

The ground-state and radially excited-state residues should be examined together. For a given window, first inspect their variation near $\alpha_*$. Then compare results at neighboring $\Lambda$, $M_{\min}^2$, and $M_{\max}^2$. The stored PP solution family directly shows whether a stable region exists near the candidate.

The LO, NLO, $C_4^P$, and $C_6^P$ contributions should be examined separately. Large cancellations between terms can otherwise conceal the window's sensitivity to the inputs. Positive residues, small reconstruction residuals, and stable plateaus constrain different aspects of the extraction and should be assessed together.

## References

1. T. Zhong, X.-G. Wu and T. Huang, *Heavy Pseudoscalar Leading-Twist Distribution Amplitudes within QCD Theory in Background Fields*, Eur. Phys. J. C **75**, 45 (2015). arXiv:1408.2297. DOI: 10.1140/epjc/s10052-015-3271-6.
2. W. Qin, Q. Zhao and X.-H. Zhong, *Revisiting the pseudoscalar meson and glueball mixing and key issues in the search for a pseudoscalar glueball state*, Phys. Rev. D **97**, 096002 (2018), Eq. (22). arXiv:1712.02550. DOI: 10.1103/PhysRevD.97.096002.
3. C. A. Bokade and Bhaghyesh, *Conventional and XYZ charmonium states in a relativistic screened potential model*, Phys. Rev. D **111**, 014030 (2025), Table II. arXiv:2408.06759. DOI: 10.1103/PhysRevD.111.014030.

RG and measured hadron-mass sources are listed in [Shared inputs and inversion methods](common_inputs_and_methods.en.md).
