# Shared inputs and inversion methods

[中文](common_inputs_and_methods.md) | **English**

## 1. Physical inputs and notation

The three channels are described in the [J/ψ vector-current](Jpsi_vector.en.md), [ηc axial-vector-current](eta_c_axial_vector.en.md), and [ηc pseudoscalar-current](eta_c_pseudoscalar.en.md) notes. Their main solvers share the following operator product expansion (OPE) inputs and numerical methods.

| Input | Value | Definition |
|---|---:|---|
| $\overline m_c(\overline m_c)$ | $1.275\ \mathrm{GeV}$ | $\overline{\mathrm{MS}}$ charm mass evaluated at its own scale |
| $\alpha_s^{(5)}(M_Z)$ | $0.1184$ | Coupling boundary value in the five-flavor theory |
| $M_Z$ | $91.1876\ \mathrm{GeV}$ | High-energy reference scale |
| $\overline m_b(\overline m_b)$ | $4.18\ \mathrm{GeV}$ | Bottom-threshold matching point |
| $G_2$ | $0.038\ \mathrm{GeV}^4$ | $\langle\alpha_sG^a_{\mu\nu}G^{a\mu\nu}\rangle$ |
| $G_3$ | $0.013\ \mathrm{GeV}^6$ | $\langle g_s^3f^{abc}G^a_{\mu\nu}G^{b\nu}{}_{\rho}G^{c\rho\mu}\rangle$ |
| $\mu$ | $2\ \mathrm{GeV}$ | Default renormalization scale, fixed throughout the Borel window |

The current calculation holds $G_2$ and $G_3$ fixed, without implementing separate scale evolution for the condensates. Both $m_c(\mu)$ and $\alpha_s(\mu)$ follow the same renormalization-group (RG) prescription.

The Borel parameter is written directly as $M^2$ throughout. Hadron masses are denoted by $m_H$ to distinguish them from the Borel parameter. The variables $s$, $s_h$, $\Lambda$, and $M^2$ have units of $\mathrm{GeV}^2$; $\mu$ and masses have units of GeV. In spectrum tables, “energy” denotes $m_H$ or the mass sum at a multihadron threshold; “spectral position” denotes its square.

The equal-mass OPE expressions use

$$
\bar x=1-x,\qquad u=x\bar x,\qquad
p_k=x^k+\bar x^k,\qquad
E(x,M^2)=\exp\left[-\frac{m_c^2(\mu)}{uM^2}\right].
$$

Below, $m$ denotes $m_c(\mu)$. The symbols $C_4$ and $C_6$ denote the Borel-transformed $D=4$ and $D=6$ contributions, respectively. They are distinct from the vacuum condensate inputs $G_2$ and $G_3$.

## 2. Running mass and coupling

We use the RunDec $\overline{\mathrm{MS}}$ convention [1]. With $a=\alpha_s/\pi$ and $t=\ln\mu^2$, the equations at $N$ loops are

$$
\frac{da}{dt}=-\sum_{k=0}^{N-1}\beta_k a^{k+2},
\qquad
\frac{d\ln m_c}{dt}=-\sum_{k=0}^{N-1}\gamma_k a^{k+1}.
$$

The selectable orders are $N=1,2,3,4$, with $N=4$ as the default. Both coefficient series are truncated at the selected order. Since the evolution variable is $t=\ln\mu^2$, no extra factor of 2 is included.

The default inputs are evolved in the following sequence:

1. Evolve $\alpha_s(M_Z)$ to $\mu=\overline m_b$ in the $n_f=5$ theory.
2. Match the five-flavor coupling to the four-flavor theory at the bottom mass evaluated at its own scale.
3. Evolve the coupling to $\mu=\overline m_c$ in the $n_f=4$ theory.
4. Impose $m_c^{(4)}(\overline m_c)=\overline m_c$, then evolve the mass and coupling together to the target scale.

The matching logarithm vanishes at the bottom mass evaluated at its own scale. Using the higher-flavor coupling $a_h=\alpha_s^{(5)}/\pi$, three-loop matching reads

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

Running at $N$ loops is paired with matching at $N-1$ loops. Four-loop running uses every term above; three-loop running retains the $a_h^2$ term. For one- and two-loop running, the coupling is continuous at this matching point. Charm remains active, with no four-to-three-flavor decoupling at the charm boundary. The current prescription supports $\overline m_c\leq\mu\leq\overline m_b$ and does not evolve the charm mass across the bottom threshold.

The default four-loop inputs are

$$
\alpha_s^{(4)}(2\ \mathrm{GeV})=0.3043984044,
\qquad
m_c^{(4)}(2\ \mathrm{GeV})=1.0932768472\ \mathrm{GeV}.
$$

The superscript $(4)$ denotes four active flavors. These values are obtained by directly integrating the finite-order RG differential equations. Four loops refers to the running of the mass and coupling; each channel's OPE truncation is specified in its own note.

## 3. Borel kernel and continuum discretization

We use the Borel-transform convention

$$
\mathcal B_{M^2}\!\left[\frac{1}{s+Q^2}\right]
=\frac{1}{M^2}e^{-s/M^2},
\qquad Q^2=-q^2>0.
$$

Each channel retains two poles at fixed masses, reconstructs the continuum over $s_h\leq s\leq\Lambda$, and uses a perturbative duality tail for $s>\Lambda$. The poles, continuum, and perturbative integral share the same $1/M^2$ normalization.

On spectral nodes $s_j$, expand the continuum in hat basis functions:

$$
\rho_{\mathrm{cont}}(s)=\sum_jc_j\phi_j(s),
\qquad \phi_j(s_k)=\delta_{jk}.
$$

For Borel nodes $M_a^2$, the pole and continuum columns of the matrix are

$$
K_{ai}^{\mathrm{pole}}=
\frac{e^{-m_i^2/M_a^2}}{M_a^2},
\qquad
K_{aj}^{\mathrm{cont}}=
\int_{s_h}^{\Lambda}\frac{ds}{M_a^2}
e^{-s/M_a^2}\phi_j(s).
$$

Continuum columns are constructed by analytically integrating the piecewise-linear basis functions. The unknown vector contains two pole residues and the continuum nodal coefficients, giving the discrete equation $K\boldsymbol x\simeq\boldsymbol g$.

The default spectral spacing is $\Delta s=0.1\ \mathrm{GeV}^2$, and the Borel spacing is $\Delta M^2=0.01\ \mathrm{GeV}^2$. Both grids include the exact interval endpoints; the final cell may be shorter than the default spacing.

All channels impose

$$
c_0=\rho_{\mathrm{cont}}(s_h)=0,
\qquad
c_N=\rho_{\mathrm{cont}}(\Lambda)=\rho_{\mathrm{pert}}(\Lambda).
$$

These boundary conditions form part of the inverse-problem definition. Matching the spectrum to the perturbative tail at the endpoint is imposed in advance. Assessing a solution also requires examining the spectrum near that boundary and its stability across windows.

## 4. Tikhonov regularization and compact GSVD

### 4.1 Regularization functional

The Borel kernel smooths the spectrum, so some spectral variations barely change the reconstructed OPE. To suppress oscillations along these directions, we solve

$$
\min_{\boldsymbol x}
\left\|W^{1/2}(K\boldsymbol x-\boldsymbol g)\right\|_2^2
+\alpha\boldsymbol c^{\mathsf T}H\boldsymbol c.
$$

The matrix $W$ contains trapezoidal weights on the Borel grid. The first term therefore approximates $\int dM^2\,|K\boldsymbol x-g|^2$; these weights arise from numerical integration.

The regularization matrix $H$ is the sum of the hat-basis mass and stiffness matrices, penalizing continuum magnitude and first-derivative variation. With the fixed numerical coordinate $\widehat s=s/(1\ \mathrm{GeV}^2)$, the penalty corresponds to

$$
\int d\widehat s\left[
\rho_{\mathrm{cont}}^2+
\left(\frac{d\rho_{\mathrm{cont}}}{d\widehat s}\right)^2
\right].
$$

The spectral interval is not rescaled to $[0,1]$ as $\Lambda$ changes. Neither pole residue is directly penalized. The parameter enters the objective as $\alpha$, rather than $\alpha^2$; its numerical value depends on the normalization of the spectrum, coordinate, and residual.

### 4.2 Elimination of fixed endpoints

Partition the spectral nodes into interior nodes $I$ and fixed endpoints $F$. The penalty contains the cross term $2\boldsymbol c_I^{\mathsf T}H_{IF}\boldsymbol c_F$. Defining

$$
\boldsymbol q=H_{II}^{-1}H_{IF}\boldsymbol c_F,
\qquad \boldsymbol c_I=\boldsymbol w-\boldsymbol q,
$$

reduces the penalty to $\boldsymbol w^{\mathsf T}H_{II}\boldsymbol w$ plus a constant independent of the unknown coefficients. Subtracting the fixed-endpoint contribution from the right-hand side and including the interior shift gives

$$
\boldsymbol b=
W^{1/2}\left(\boldsymbol g-K_F\boldsymbol c_F+K_I\boldsymbol q\right).
$$

With $R_H^{\mathsf T}R_H=H_{II}$, the shifted regularization operator is

$$
L=\begin{pmatrix}0_{n_I\times2}&R_H\end{pmatrix}.
$$

This elimination preserves the original functional, including the coupling between endpoint values and the interior spectrum.

### 4.3 Generalized filtering and residue extraction

Apply compact generalized singular value decomposition (GSVD) to the weighted matrix $A=W^{1/2}K_{\mathrm{free}}$ and the regularization operator $L$. Write $AX=UC$ and $LX=VD$, with diagonal factors $c_k,d_k$ on the regularized modes and data projections $b_k=\boldsymbol u_k^{\mathsf T}\boldsymbol b$. The solution is

$$
\boldsymbol y_\alpha
=\sum_{k\in\mathrm{reg}}
\frac{c_k b_k}{c_k^2+\alpha d_k^2}\boldsymbol x_k
+\sum_{k\in\mathrm{null}}
\frac{b_k}{c_k}\boldsymbol x_k.
$$

The second sum corresponds to the null space of $L$. After solving, recover the interior spectrum through $\boldsymbol c_I=\boldsymbol w-\boldsymbol q$ and restore the endpoints. Both Python and MATLAB main solvers solve this generalized problem.

Neither the two residues nor the continuum is constrained to be nonnegative. Only positive residues are converted into decay constants using the current-specific definitions. Nonpositive residues are retained, with their corresponding decay constants marked as undefined.

## 5. Q/MQ/HR regularization-parameter selection

Let $\widetilde\sigma_k=c_k/d_k$, and let $\delta_\perp^2$ denote the squared norm of the shifted data component orthogonal to the operator's column space. Define

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

The quantity $\psi_Q$ measures the variation of the standard-form regularized solution with $\log\alpha$; it is not the geometric curvature itself. The Q-curve uses $d_{\mathrm{MD}}$ and $\psi_Q$ as its two diagnostic coordinates.

On an increasing grid $\alpha_1<\cdots<\alpha_N$, set $p_0=1$ and recurse from larger $\alpha$ to smaller $\alpha$ using

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

Within the set of local minima of $\psi_Q$, choose $\alpha_Q$ and $\alpha_{MQ}$ to minimize $\psi_Q$ and $\psi_{MQ}$, respectively. Choose $\alpha_{HR}$ to minimize $\psi_{HR}$ over grid points satisfying $\alpha\geq\max(\alpha_1,\widetilde\sigma_{\min}^2)$. With $r_Q=1/\alpha_Q$ and analogous definitions for the other candidates, combine them as

$$
r_H=\min\{r_Q,\max(r_{MQ},r_{HR})\},
\qquad \alpha_* = 1/r_H.
$$

If $\psi_Q$ has no local minimum, use the lowest local minimum of $\psi_{HR}$. If none exists, use the HR minimum in the restricted interval above. The main solver records these fallback selections.

The parameter $\alpha_*$ is selected independently in every $(\Lambda,M_{\min}^2,M_{\max}^2)$ window. The VV and AA logarithmic grids can extend toward smaller $\alpha$, down to $10^{-80}$. PP uses 2001 logarithmically spaced points over $[10^{-20},1]$ by default. All three channels share the selection rule while retaining these search ranges.

## 6. Stability-plateau criterion

The main solver returns the solution, regularization diagnostics, and reconstructed OPE for a specified window. Physical interpretation also requires checking whether the residues remain stable over a family of neighboring windows.

First fix $(\Lambda,M_{\min}^2,M_{\max}^2)$ and examine both residues near $\alpha_*$ to identify a stable interval of nonzero width. Then vary $\Lambda$, $M_{\min}^2$, and $M_{\max}^2$ to assess whether both decay constants remain stable in neighboring windows. Check also for discontinuous jumps in the selected regularization parameter.

Small $M^2$ enhances low-lying states but can also amplify higher-dimensional condensate contributions. Larger $M^2$ increases sensitivity to the high-energy spectrum and reduces the ability to resolve nearby poles. Window selection should therefore consider the relative sizes of the perturbative term, $C_4$, and $C_6$, together with the spectral shape. A small reconstruction residual shows agreement with the integral equation but cannot alone establish reliable separation of the two residues.

## 7. Sources of hadronic-spectrum values

Explicit poles use the fixed masses in each main solver. Other observed-state and open-charm particle masses use PDG 2026 Summary Table central values labeled OUR AVERAGE/FIT [2]. Spectral positions are computed after rounding masses to the displayed precision. These tabulated values are not automatically added as solver inputs.

Open-charm thresholds are calculated from the following masses:

| Particle | Mass / GeV |
|---|---:|
| $D^0$ | 1.86484 |
| $D^+$ | 1.86966 |
| $D^{*0}$ | 2.00686 |
| $D^{*+}$ | 2.01027 |
| $D_s^+$ | 1.96835 |
| $D_s^{*+}$ | 2.1122 |

Two-body thresholds are calculated as $s_{\mathrm{th}}=(m_a+m_b)^2$. The tables list the main charmonium states, charmonium-like structures, and open-charm thresholds relevant to the current model. The full spectrum may contain additional multihadron contributions that are not modeled explicitly. Matching quantum numbers permits a coupling but does not determine the strength with which the local $\bar c\Gamma c$ current couples to a state.

## References

1. K. G. Chetyrkin, J. H. Kühn and M. Steinhauser, *RunDec: A Mathematica package for running and decoupling of the strong coupling and quark masses*, Comput. Phys. Commun. **133**, 43–65 (2000). arXiv:hep-ph/0004189. DOI: 10.1016/S0010-4655(00)00155-7.
2. F. Takahashi et al. (Particle Data Group), Int. J. Mod. Phys. A **41**, 2630011 (2026). Mass values were queried from the official PDG 2026.0 database, with release timestamp 2026-05-30 13:41:24 PDT.
