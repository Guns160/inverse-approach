# Breit–Wigner resonances in the ηc pseudoscalar spectrum

[中文](eta_c_pseudoscalar+BW.md) | **English**

This note adds one finite-width $\eta_c(nS)$ excitation to the original model of two low-lying narrow poles and a continuum. The current definition, OPE and inversion method follow the [ηc pseudoscalar-current note](eta_c_pseudoscalar.en.md) and [shared inputs and inversion methods](common_inputs_and_methods.en.md). Natural units, $\hbar=c=1$, are used throughout.

Let $s$ denote the invariant mass squared, $M^2$ the Borel parameter, and $[s_h,\Lambda]$ the interval over which the spectrum is reconstructed. The selected state's mass is $m_{\mathrm{BW}}=m_{\eta_c(nS)}$. Its resonance line shape and coefficient are denoted by $\mathrm{BW}(s)$ and $R_{\mathrm{BW}}$, respectively. The barrier radius $r_{\mathrm{BW}}$ is a separate parameter.

## 1. Breit–Wigner line shape

### 1.1 Finite lifetime and resonance width

The survival probability of an unstable particle decreases exponentially with time:

$$
P(t)=e^{-\Gamma_0t},\qquad \tau=\frac{1}{\Gamma_0}.
$$

The corresponding amplitude is $\psi(t)\propto e^{-iE_0t}e^{-\Gamma_0t/2}$. Its Fourier transform gives

$$
\widetilde\psi(E)\propto\int_0^\infty dt\,e^{i(E-E_0)t-\Gamma_0t/2}
\propto\frac{1}{E-E_0+i\Gamma_0/2}.
$$

The normalized line shape in energy space is therefore

$$
\rho_E(E)=\frac{1}{\pi}\frac{\Gamma_0/2}
{(E-E_0)^2+\Gamma_0^2/4}.
$$

The full width at half maximum is $\Gamma_0$. A shorter lifetime corresponds to a larger width and a broader peak.

### 1.2 Relativistic propagator and spectral density

As a function of the invariant mass squared $s$, a commonly used propagator near a resonance is [1,4]

$$
D_{\mathrm{BW}}(s)=\frac{1}{m_{\mathrm{BW}}^2-s-i m_{\mathrm{BW}}\Gamma(s)}.
$$

Near the peak, let $E=\sqrt{s}\simeq m_{\mathrm{BW}}+\delta E$. Then $s-m_{\mathrm{BW}}^2\simeq2m_{\mathrm{BW}}\delta E$, which recovers the form $1/(\delta E+i\Gamma_0/2)$ above. The total decay width $\Gamma(s)$ depends on energy.

The spectral representation of the pseudoscalar-current two-point correlator is

$$
\Pi_P(q^2)=\int ds'\,\frac{\rho_P(s')}{s'-q^2-i0}
+\text{subtractions},\qquad
\rho_P(s)=\frac1\pi\operatorname{Im}\Pi_P(s+i0).
$$

The relation between the imaginary part and the spectral density follows from
$1/(x-i0)=\operatorname{PV}(1/x)+i\pi\delta(x)$. Approximating the resonance contribution by $D_{\mathrm{BW}}$ gives the BW line shape used here:

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

The imaginary part of the propagator is used because the quantity of interest is the **spectral density of the correlator**. The squared amplitude has the same denominator but a different numerator. In the narrow-width limit with a constant width,

$$
\lim_{\Gamma_0\to0^+}\mathrm{BW}(s)=\delta(s-m_{\mathrm{BW}}^2),
$$

which recovers a narrow pole.

## 2. Running width and barrier factors

### 2.1 Two-body phase space and partial-wave thresholds

According to Fermi's golden rule, a partial decay width depends on the transition amplitude and the final-state phase space:

$$
\Gamma_\nu(s)=\frac{1}{2\sqrt{s}}\int d\Phi_\nu\,
\overline{|\mathcal M_\nu(s)|^2}.
$$

For a two-body channel $A B$ with narrow final-state particles, the center-of-mass momentum is

$$
\boxed{
q_{AB}(s)=
\begin{cases}
\displaystyle\frac{\sqrt{[s-(m_A+m_B)^2][s-(m_A-m_B)^2]}}{2\sqrt{s}},
&s\geq(m_A+m_B)^2,\\[6pt]
0,&s<(m_A+m_B)^2.
\end{cases}}
$$

Two-body phase space supplies one factor of $q$. Near threshold, a partial-wave amplitude with orbital angular momentum $\ell$ satisfies $\mathcal M_\ell\sim q^\ell$. Thus,

$$
\boxed{\Gamma_\ell(s)\propto q^{2\ell+1}.}
$$

The S-, P- and D-wave threshold factors are $q$, $q^3$ and $q^5$, respectively. Higher orbital angular momentum gives stronger suppression near threshold.

### 2.2 Blatt–Weisskopf barrier factors

Blatt–Weisskopf factors describe finite interaction-range effects in partial waves at higher momenta [1,4]:

$$
\boxed{
\Phi_{AB}^{(\ell)}(s)
=\frac{q_{AB}(s)^{2\ell+1}}{\sqrt{s}}
B_\ell^2\!\left(q_{AB}(s)r_{\mathrm{BW}}\right).
}
$$

The factors used here are

$$
B_0^2(z)=1,\qquad
B_1^2(z)=\frac{1}{1+z^2},\qquad
B_2^2(z)=\frac{1}{9+3z^2+z^4},
\qquad z=qr_{\mathrm{BW}}.
$$

The factor $q^{2\ell+1}$ determines the threshold behavior, while $B_\ell^2(qr_{\mathrm{BW}})$ modifies the momentum dependence away from threshold. For example, the P-wave factor is $q^3/(1+q^2r_{\mathrm{BW}}^2)$; for the S wave, $B_0^2=1$.

With $q_0=q_{AB}(m_{\mathrm{BW}}^2)$, the ratio needed for the running width of a narrow two-body channel is

$$
\boxed{
\frac{\Phi_{AB}^{(\ell)}(s)}{\Phi_{AB}^{(\ell)}(m_{\mathrm{BW}}^2)}
=\left(\frac{q(s)}{q_0}\right)^{2\ell+1}
\frac{m_{\mathrm{BW}}}{\sqrt{s}}
\frac{B_\ell^2(q(s)r_{\mathrm{BW}})}{B_\ell^2(q_0r_{\mathrm{BW}})}.
}
$$

This expression corresponds to the normalized barrier-factor convention in Ref. [4].

### 2.3 Choice of barrier radius

The baseline value is

$$
\boxed{r_{\mathrm{BW}}=4\,\mathrm{GeV}^{-1}\simeq0.79\,\mathrm{fm}}
$$

This is a **phenomenological interaction-range parameter**, rather than a directly measured geometric radius of $\eta_c(5S)$. The LHCb amplitude analysis of $B^-\to D^+\pi^-\pi^-$ uses $4\,\mathrm{GeV}^{-1}$ and varies it over $3$–$5\,\mathrm{GeV}^{-1}$ in its systematic studies [4]. The choice of 4 therefore has a literature precedent and serves as a suitable default for this model.

The values $r_{\mathrm{BW}}=2,3,5\,\mathrm{GeV}^{-1}$ can also be used. They correspond to approximately $0.39,0.59,0.99\,\mathrm{fm}$, respectively. Changing $r_{\mathrm{BW}}$ modifies the P- and D-wave barrier factors, and therefore the running partial widths and BW tails. It has no direct barrier-factor effect on the S wave. The values $3,5\,\mathrm{GeV}^{-1}$ are recommended for routine sensitivity checks, with $2\,\mathrm{GeV}^{-1}$ as an extended check.

After changing the radius, the channel phase spaces, mass shifts and common normalization must be recalculated. Replacing a parameter only in the final line shape is insufficient.

### 2.4 Partial-wave assignments for $\eta_c(nS)$

The parent state has quantum numbers $J^{PC}=0^{-+}$. Coupling the total final-state spin $S$ and orbital angular momentum $\ell$ to $J=0$ requires $\ell=S$. Parity conservation also requires

$$
P_A P_B(-1)^\ell=-1.
$$

For the eight classes of open-charm channels used for 5S, the assignments are [2]:

| Decay channel | $\ell$ | Partial wave |
|---|---:|---|
| $DD^*$ | 1 | P wave |
| $D_sD_s^*$ | 1 | P wave |
| $DD_0^*$ | 0 | S wave |
| $D_sD_{s0}^*$ | 0 | S wave |
| $DD_2^*$ | 2 | D wave |
| $D^*D^*$ | 1 | P wave |
| $D_s^*D_s^*$ | 1 | P wave |
| $D^*D_1(2420)$ | 0; alternatively 2 | S-wave baseline; D-wave comparison |

The $D^*D_1(2420)$ channel can form $0^-$ through either $S=\ell=0$ or $S=\ell=2$. The partial-width table in Ref. [2] does not separate these waves, so they are calculated separately for comparison. The channel abbreviations include the corresponding antiparticles and charge-conjugate combinations.

## 3. Channel-dependent running widths

### 3.1 Partial-width inputs and mass shifts

The total width is written as

$$
\boxed{
\Gamma(s)=\Gamma_{\mathrm{rest}}+\sum_\nu\Gamma_\nu(s),
\qquad \Gamma_{\mathrm{rest}}=\Gamma_{gg}+\Gamma_{\gamma\gamma}.
}
$$

The open-charm partial widths $\Gamma_\nu^{\mathrm{ref}}$ are taken from the $\Gamma_b$ column of Table 3 in Wang, Yan and Ping [2]. The values of $m_{\mathrm{BW}}$, $\Gamma_{gg}$ and $\Gamma_{\gamma\gamma}$ are taken from Bokade and Bhaghyesh [3], using annihilation widths that include the correction factors. The annihilation widths are approximated as constants here.

The two input sets use different parent masses. Denoting the reference mass in Ref. [2] by $m_{\mathrm{ref}}$, first shift each partial width using a phase-space ratio:

$$
\boxed{
\widetilde\Gamma_\nu
=\Gamma_\nu^{\mathrm{ref}}
\frac{\Phi_\nu(m_{\mathrm{BW}}^2)}{\Phi_\nu(m_{\mathrm{ref}}^2)}.
}
$$

Then rescale all channels together while preserving the selected total open-charm width $\Gamma_{\mathrm{open}}$:

$$
\boxed{
\Gamma_\nu^0\equiv\Gamma_\nu(m_{\mathrm{BW}}^2)
=\Gamma_{\mathrm{open}}
\frac{\widetilde\Gamma_\nu}{\sum_\kappa\widetilde\Gamma_\kappa}.
}
$$

Thus, $\sum_\nu\Gamma_\nu^0=\Gamma_{\mathrm{open}}$ and $\Gamma_0=\Gamma_{\mathrm{rest}}+\Gamma_{\mathrm{open}}$.

### 3.2 Running widths

Each channel then runs with $s$ according to its phase-space ratio:

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

The running-width form follows Ref. [1]. The mass shift and common rescaling at fixed $\Gamma_{\mathrm{open}}$ are model conventions adopted in this project. Changing the barrier radius or partial-wave assignment requires recalculating $\Phi_\nu$ and $\Gamma_\nu^0$.

### 3.3 Charge-channel averaging

For channels containing both neutral and charged final states, first average their phase-space functions with equal weights. For example,

$$
\Phi_{DD^*}(s)=\frac12\left[
\Phi_{D^0\bar D^{*0}}^{(1)}(s)
+\Phi_{D^+D^{*-}}^{(1)}(s)\right].
$$

The complete $\Phi_{DD^*}$ is then used for the mass shift and running width. The partial widths in Ref. [2] already include the corresponding charge and conjugate combinations, so normalization uses their combined value.

## 4. Spectral convolution for broad final-state particles

### 4.1 Convolved phase space

The $D_0^*$ and $D_1$ particles in $DD_0^*$ and $D^*D_1(2420)$ have finite widths themselves. Near threshold, their lower-mass spectral components may also contribute to the decay. The final-state mass is therefore integrated over $s_B$ rather than fixed at $m_B$ [5].

Use the constant-width final-state spectral distribution

$$
\rho_B(s_B)=\frac{m_B\Gamma_B}
{\pi[(s_B-m_B^2)^2+m_B^2\Gamma_B^2]}.
$$

The phase space becomes

$$
\boxed{
\Phi_{AB}^{(\ell)}(s)=
\int_{s_{B,\mathrm{th}}}^{(\sqrt{s}-m_A)^2}
 ds_B\,\rho_B(s_B)
\frac{q(s;m_A,\sqrt{s_B})^{2\ell+1}}{\sqrt{s}}
B_\ell^2\!\left(q(s;m_A,\sqrt{s_B})r_{\mathrm{BW}}\right).
}
$$

Here $q(s;m_A,\sqrt{s_B})$ is obtained by replacing $m_B$ with $\sqrt{s_B}$ in the two-body momentum formula. The daughter-decay thresholds are

$$
D_0^*\to D\pi:\quad s_{B,\mathrm{th}}=(m_D+m_\pi)^2,
\qquad
D_1(2420)\to D^*\pi:\quad s_{B,\mathrm{th}}=(m_{D^*}+m_\pi)^2.
$$

The phase space is zero when $\sqrt{s}\leq m_A+\sqrt{s_{B,\mathrm{th}}}$. The barrier factor depends on $s_B$ and therefore remains inside the integral. This project fixes the parent's $s$ and integrates over the physically accessible final-state mass interval, without renormalizing the truncated final-state spectrum.

### 4.2 Daughter-decay weights

In the isospin approximation, the charged-to-neutral $\pi$ weights in $D_0^*\to D\pi$ and $D_1\to D^*\pi$ are $2:1$. Thus,

$$
\Phi_{\nu,k}(s)=\frac23\Phi_{\nu,k}^{(\pi^\pm)}(s)
+\frac13\Phi_{\nu,k}^{(\pi^0)}(s),
\qquad
\Phi_\nu(s)=\frac12\sum_{k=0,+}\Phi_{\nu,k}(s).
$$

Each term uses its corresponding physical threshold. Final-state masses and neutral-state widths are taken from the PDG [6]. Charged $D_0^*$ and $D_1$ widths are approximated by their neutral counterparts. After convolution and charge averaging, the result enters the partial-width formulas in Section 3.

## 5. The BW spectral term and the Borel inverse problem

### 5.1 Resonance spectral term and Borel kernel

Explicitly separate one excited state over $s_h\leq s\leq\Lambda$:

$$
\boxed{
\rho_P(s)=R_{\mathrm{BW}}\,\mathrm{BW}(s)
+\rho_P^{\mathrm{cont}}(s).
}
$$

The two low-lying narrow poles and the perturbative tail above $\Lambda$ remain unchanged. The dimensions of $\mathrm{BW}(s)$, $R_{\mathrm{BW}}$ and $\rho_P$ are $\mathrm{GeV}^{-2}$, $\mathrm{GeV}^{4}$ and $\mathrm{GeV}^{2}$, respectively.

The additional Borel-kernel matrix column is

$$
\boxed{
K_a^{\mathrm{BW}}=\int_{s_h}^{\Lambda}\frac{ds}{M_a^2}
 e^{-s/M_a^2}\mathrm{BW}(s),
\qquad \Delta g_a=K_a^{\mathrm{BW}}R_{\mathrm{BW}}.
}
$$

It uses the same $1/M_a^2$ normalization as the original $K_{ai}^{\mathrm{pole}}$ and $K_{aj}^{\mathrm{cont}}$. The finite-width spectrum receives different Borel weights at different values of $s$. If $\Lambda<m_{\mathrm{BW}}^2$, the integral includes only the low-energy tail within the interval.

The BW area over this interval is

$$
A_{\mathrm{BW}}(s_h,\Lambda)=\int_{s_h}^{\Lambda}\mathrm{BW}(s)\,ds.
$$

This area is not used to renormalize the line shape. The fitted coefficient remains $R_{\mathrm{BW}}$, and the integrated resonance strength over the interval is $R_{\mathrm{BW}}A_{\mathrm{BW}}$.

### 5.2 Endpoint matching

The sum of the BW term and remaining continuum satisfies the endpoint conditions:

$$
\boxed{
\begin{aligned}
\rho_P^{\mathrm{cont}}(s_h)+R_{\mathrm{BW}}\mathrm{BW}(s_h)&=0,\\
\rho_P^{\mathrm{cont}}(\Lambda)+R_{\mathrm{BW}}\mathrm{BW}(\Lambda)
&=\rho_P^{\mathrm{pert}}(\Lambda).
\end{aligned}}
$$

With $c_j=\rho_P^{\mathrm{cont}}(s_j)$,

$$
\boxed{
\begin{aligned}
c_0&=-R_{\mathrm{BW}}\mathrm{BW}(s_h),\\
c_N&=\rho_P^{\mathrm{pert}}(\Lambda)-R_{\mathrm{BW}}\mathrm{BW}(\Lambda).
\end{aligned}}
$$

For $R_{\mathrm{BW}}>0$, the lower-end continuum coefficient may satisfy $c_0<0$ to cancel the BW tail and enforce the total-spectrum boundary condition. Under the original smoothness penalty, the endpoint conditions also couple $R_{\mathrm{BW}}$ to the remaining continuum coefficients. The two original narrow poles remain unpenalized directions.

## 6. Excited states and model checks

**5S (primary model):** Use the eight open-charm channels in Section 2.4. Apply broad-final-state convolutions to $DD_0^*$ and $D^*D_1(2420)$. For the latter, use an S-wave baseline and a D-wave comparison.

**4S (alternative):** Retain the four P-wave channels $DD^*$, $D_sD_s^*$, $D^*D^*$ and $D_s^*D_s^*$. The $DD_0^*$ channel is omitted because its 4S coupling input is currently unavailable.

**3S (comparison):** Currently retain only $DD^*$ and omit the $D^*D^*$ channel on the high-energy side of the peak. Its proximity to the existing low-lying poles makes residue separation sensitive after Borel smoothing.

**6S:** The channel inputs remain incomplete. Some reference partial widths are incompatible with thresholds at the reference mass, and some broad-final-state parameters are missing. No running width is constructed at present.

The baseline for model sensitivity checks uses $r_{\mathrm{BW}}=4\,\mathrm{GeV}^{-1}$, the eight 5S channels and an S wave for $D^*D_1(2420)$. Vary the radius to $r_{\mathrm{BW}}=3,5$ (extended to 2) and test the S-/D-wave alternatives separately. For each variation, recalculate from $\Phi_\nu$ onward. Compare $\Gamma(s)$, $\mathrm{BW}(s)$, $K_a^{\mathrm{BW}}$ and the inferred $R_{\mathrm{BW}}$.

## References

1. D. M. Asner, C. Hanhart and M. Mikhasenko, *Resonances*, in R. L. Workman et al. (Particle Data Group), *Review of Particle Physics*, Prog. Theor. Exp. Phys. **2022**, 083C01 (2022), 2023 update. [Official review](https://pdg.lbl.gov/2023/reviews/rpp2022-rev-resonances.pdf). Section 50.3.1 and Eqs. (50.30)–(50.35) provide the BW, running-width and barrier-factor conventions used here.
2. H. Wang, Z. Yan and J. Ping, *Radially excited states of $\eta_c$*, Eur. Phys. J. C **75**, 196 (2015). [DOI](https://doi.org/10.1140/epjc/s10052-015-3418-5), [arXiv:1412.7068](https://arxiv.org/abs/1412.7068). Table 2 provides the reference masses; Table 3 provides the $\Gamma_b$ partial widths; Eqs. (13)–(15) specify the partial-wave and kinematic conventions.
3. C. A. Bokade and Bhaghyesh, *Conventional and XYZ charmonium states in a relativistic screened potential model*, Phys. Rev. D **111**, 014030 (2025). [DOI](https://doi.org/10.1103/PhysRevD.111.014030), [arXiv:2408.06759](https://arxiv.org/abs/2408.06759). The inputs are taken from arXiv v1: Table II for masses and the corrected $\Gamma_{\mathrm{cf}}$ entries in Tables XIV and XVI.
4. R. Aaij et al. (LHCb Collaboration), *Amplitude analysis of $B^-\to D^+\pi^-\pi^-$ decays*, Phys. Rev. D **94**, 072001 (2016). [DOI](https://doi.org/10.1103/PhysRevD.94.072001), [arXiv:1608.01289](https://arxiv.org/abs/1608.01289). Equations (12), (14) and (15) give the barrier factors, BW amplitude and running width. The analysis uses $r_{\mathrm{BW}}=4\,\mathrm{GeV}^{-1}$ and studies variations over $3$–$5\,\mathrm{GeV}^{-1}$.
5. L. Roca, E. Oset and J. Singh, *Low lying axial-vector mesons as dynamically generated resonances*, Phys. Rev. D **72**, 014002 (2005). [DOI](https://doi.org/10.1103/PhysRevD.72.014002), [arXiv:hep-ph/0503273](https://arxiv.org/abs/hep-ph/0503273). Section 3.3, Eqs. (23) and (24), describes spectral folding for finite-width final states.
6. F. Takahashi et al. (Particle Data Group), *Review of Particle Physics*, Int. J. Mod. Phys. A **41**, 2630011 (2026). [PDG 2026](https://pdg.lbl.gov/2026/). Used for final-state masses and widths.
