# OPE truncation error

[中文](ope-truncation.md) | **English**

The OPE is retained only up to a finite operator dimension. The effect of the omitted higher-dimensional operators on the result is the truncation error.

> **This document records the method only.** The retained orders and windows given here are examples or current settings; the actual choices are determined by the specific calculation.

## 1. What is being estimated

The OPE in this repository consists of the perturbative term, $G_2$, and $G_3$; the highest retained term is $G_3$.

Higher-dimensional operators are not explicitly included. Their contributions to the output cannot be read directly from the existing calculation and can only be estimated from the sizes of the other terms.

## 2. Estimation method

Treat the highest-dimensional retained term as an upper bound on all higher-dimensional terms, and take the difference between results calculated “with” and “without” it.

Rationale: higher operator dimensions lead to stronger suppression by powers of $1/Q^2$. The highest retained term is the largest among the remaining terms, so using it as an upper bound is a conservative estimate.

```text
PROCEDURE ope_truncation(center_inputs):
    y_full  ← SOLVE(center_inputs, G3 = G3_center)
    y_trunc ← SOLVE(center_inputs, G3 = 0)
    RETURN |y_full − y_trunc|
```

The two calculations are identical except for $G_3$: the same inputs, window, and boundary conditions.

## 3. Checks

1. **Successive decrease.** The $G_3$ contribution should be smaller than the $G_2$ contribution. If $G_3$ is larger than $G_2$, the window lies outside the OPE convergence region and the truncation estimate is unreliable.
2. **Magnitude.** Truncation error should not dominate the total error. If it is the largest component, a higher OPE order is needed or the window choice is problematic.
3. **Sign.** Use the absolute value in the quadrature sum, but the signed response direction can be given in the report to help identify systematic tendencies.
4. **Consistent boundaries.** Both calculations must use the same rule for the boundary condition at $s=\Lambda$. If one uses the full OPE and the other uses the truncated OPE, the difference will include boundary contributions.

## 4. Difference from sampling

| Item | Gaussian sampling | OPE truncation |
|---|---|---|
| Is there a distribution? | Yes | No |
| Number of endpoints | Determined by the sample count | Only two |
| Included in the ensemble? | Yes | No, a separate number |
| Combination method | Calculate the standard deviation or range first, then combine in quadrature | Combine directly in quadrature |

There is no intermediate value such as “half of $G_3$”. $G_3$ is either included or omitted.

This component has no random numbers, and neither a standard deviation nor quantiles can be reported. Enter only one number in the table.

## 5. Reporting format

| Retained OPE order | Removed term | $\Delta y_1$ | $\Delta y_2$ |
|---|---|---|---|
| Perturbative $+G_2+G_3$ | $G_3$ | — | — |

Below the table, state the window, $\mu$, and the ratio of the $G_3$ and $G_2$ contributions.
