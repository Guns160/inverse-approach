# Renormalization scale $\mu$

[中文](mu-scale.md) | **English**

$\mu$ is neither a sampling error nor a truncation error. It is a scheme choice.

> **This document records the method only.** The numerical values given here are current conventions, not final choices.

## 1. Convention in this project

**$\mu$ is fixed at 2 GeV and is not included in the error.**

Reason: $\mu$ is a scheme choice, not a measured value, and has no distribution to sample. Earlier project budgets also adopted this convention.

$\alpha_s(\mu)$ and $m_c(\mu)$ are calculated by RG evolution from their starting points to $\mu$ and automatically change with $\mu$. Fixing $\mu$ means that both are also fixed.

## 2. What changes when $\mu$ varies

Only two parameters run with $\mu$:

| Quantity | How it changes | Dependence |
|---|---|---|
| $\alpha_s(\mu)$ | Evolved by the RG from $\alpha_s^{(5)}(M_Z)$ to $\mu$ | Beta function |
| $m_c(\mu)$ | Evolved by the RG from $\overline m_c(\overline m_c)$ to $\mu$ | Mass anomalous dimension |

The scale dependence of $G_2$ and $G_3$ themselves is small and can be neglected.

Although $m_c(\mu)$ is a single parameter, it enters the OPE in multiple places: the perturbative term, the coefficients of the $G_2$ and $G_3$ terms, the perturbative spectral threshold, and the residue conversion in the pseudoscalar channel. Thus the actual effect of $\mu$ extends through the entire chain rather than changing just one factor.

This has two consequences:

**First, a scan yields the combined scheme dependence.** The span cannot be separated into an “$\alpha_s$ component” and an “$m_c$ component”, because both are driven by the same $\mu$ and are inherently correlated.

**Second, it cannot be added to the $\overline m_c$ error row.** A $\mu$ scan already includes the variation of $m_c(\mu)$; adding the sampling error of $\overline m_c$ would count the same effect twice. This is also one reason why this project fixes $\mu$ and samples only $\overline m_c$.

## 3. Difference from the reference project

The reference project scans $\mu$ from 1 to 3 GeV, with a central value of 2 GeV, obtaining 1.69 MeV ($f_\pi$) and 0.91 MeV ($f_{a_1}$).

The budget in this project cannot be compared item by item with that row. This was already stated in the project's earlier J/ψ budget.

## 4. Diagnostic scan

If an additional scale diagnostic is desired, use a scan, but do not include it in the main error.

```text
PROCEDURE mu_diagnostic(center_inputs, mu_values):
    ys ← []
    FOR each mu IN mu_values:
        y ← SOLVE(center_inputs, mu = mu)
        ys.append(y)
    RETURN (max(ys) − min(ys)) / 2
```

During the scan, $\alpha_s(\mu)$ is automatically recalculated by the RG, while $m_c(\mu)$ changes simultaneously; the two are not independent.

The diagnostic value only indicates the magnitude of scale dependence. It does not enter the quadrature sum or appear in the main error table.

## 5. Reporting

Do not include a $\mu$ row in the main error table.

If a diagnostic scan was performed, describe it in a separate paragraph, giving the scan interval, number of points, and resulting half-width, and label it “diagnostic value, not included in the total error”.
