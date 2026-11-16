# Differences from the reference project

[中文](reference-comparison.md) | **English**

The reference project concerns the axial-vector-current channel of $\pi$ and $a_1(1260)$, whereas this project concerns charmonium. Their input lists overlap only partially.

## 1. Item-by-item comparison

| Reference project ($\pi$) | This project (charmonium) | Difference |
|---|---|---|
| $m_u\langle\bar uu\rangle+m_d\langle\bar dd\rangle$ | None | Light-quark condensates, absent from the charmonium OPE |
| $\langle\alpha_sGG\rangle$ | $G_2$ | Same name, different values (0.07 vs 0.038) |
| $\sum_{q}m_q\langle g_s\bar q\sigma TGq\rangle$ (dim-5) | None | Mixed condensates, not included in the charmonium OPE |
| $\langle g_s\bar qq\rangle^2$, $\langle g_s^2\bar qq\rangle^2$ | None | Four-quark condensates; charmonium retains the three-gluon term instead |
| $\kappa\equiv\langle\bar ss\rangle/\langle\bar qq\rangle$ | None | No strange-quark input |
| $M_{a_1}$, $\Gamma_{a_1}$ | `MASSES` (pole masses) | The reference project treats resonance parameters as prior inputs with experimental errors; in this project, pole masses are fixed model quantities |
| $\mu$ scan over $[1,3]$ GeV | Fixed $\mu$ | Different treatment |
| OPE truncation (dim-6 four-quark) | OPE truncation ($G_3$) | Same form, different omitted terms |
| None | $\bar m_c(\bar m_c)$ | **Unique to this project**, the starting point of heavy-quark running |
| None | $a_1$ (axial anomaly) | Pseudoscalar channel only |
| $G_3=8.2\,G_2$ (derived quantity) | Independent constant $G_3$ | In the reference project, the three-gluon condensate is derived from $G_2$ and is not independent; here it is independent |

## 2. Three conclusions

**First, this project adds the heavy-quark mass input.**

$\bar m_c$ has no counterpart in the reference project. The error from the light-quark mass combination in the reference project is only 0.002 MeV and is negligible; here $\bar m_c$ changes both $m_c(\mu)$ and the perturbative spectrum through the RG and is the dominant input source.

This project uses $\alpha_s$ already evolved to the calculation scale and assigns it no separate error, consistent with the treatment in the reference project.

**Second, this project omits all light-quark condensates.**

Of the ten input-error rows in the reference project, five concern light-quark or four-quark condensates; none appears in this project. The charmonium OPE retains only gluonic operators.

**Third, the independence of $G_3$ needs to be confirmed.**

The reference project uses vacuum saturation to write the three-gluon condensate as $G_3=8.2\,G_2$, so the error of $G_2$ automatically covers $G_3$. This project writes both as independent constants.

If the same vacuum-saturation relation is adopted, $G_3$ should not be treated as a separate error source, as this would double count the contribution. This point must be settled before deciding whether $G_3$ is an independent input or a quantity derived from $G_2$.
