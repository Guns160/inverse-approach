# Plateau-parameter error

[中文](platform-parameters.md) | **English**

This file explains how to estimate plateau-parameter uncertainties in the three-channel charmonium inverse problem. The method follows the reference paper, with parameter names adapted to this repository's Borel formulation.

> **This document records the method only.** The specific values and intervals given here are examples, not the actual choices for this project. Plateau boundaries and scan ranges must be determined by the specific calculation and then entered in the report.

## 1. Parameter list

This problem has three plateau parameters, all in GeV².

| Parameter | Code name | Meaning |
|---|---|---|
| Lower bound on the Borel parameter | `M2_min` | Lower end of the Borel window |
| Upper bound on the Borel parameter | `M2_max` | Upper end of the Borel window |
| Continuum boundary | `Lambda` | Boundary between the nonperturbative and OPE regions, also the upper end of the $s$ grid |

Correspondence with the reference paper:

| Reference paper | This work | Explanation |
|---|---|---|
| $q^2_{\rm up}$ | `M2_max` | Upper end of the input interval |
| $q^2_{\rm low}$ | `M2_min` | Lower end of the input interval |
| $\Lambda$ | `Lambda` | Corresponding parameter |
| $\sigma_0$ (subtraction point) | None | The Borel transform eliminates the subtraction, so this problem has no $\sigma_0$ |

Of the four plateau parameters in the reference paper, only three remain in this problem.

## 2. Meaning of a plateau

A plateau is not a numerical trick but the intersection of physical conditions. Each parameter is constrained from both above and below; the plateau is the interval in which both constraints are satisfied.

### 2.1 `M2_min`

As $M^2$ decreases, the Borel weight $e^{-s/M^2}$ becomes more sensitive to low $s$, increasingly suppressing contributions from higher excited states and the continuum. However, if $M^2$ is too small, higher-order terms in the operator product expansion become relatively larger and the OPE itself becomes unreliable.

The lower bound is therefore determined by **OPE convergence**.

### 2.2 `M2_max`

As $M^2$ increases, OPE convergence improves, but the Borel weight suppresses large $s$ less strongly. Higher excited states and the continuum re-enter, and the ground state no longer dominates.

The upper bound is therefore determined by **ground-state dominance**.

Together, the two ends define a window $[M^2_{\min}, M^2_{\max}]$. Every $M^2$ within the window satisfies both conditions, and the result should not change appreciably as the window endpoints move within their plateaus.

### 2.3 `Lambda`

$\Lambda$ is the boundary between the nonperturbative and OPE regions. The spectrum for $s<\Lambda$ is obtained by solving the inverse problem, while the OPE is used directly for $s>\Lambda$.

- $\Lambda$ too small: the solution interval is narrow, the information is insufficient, and the solution is sensitive to boundary conditions.
- $\Lambda$ too large: the OPE is already invalid near $\Lambda$, so the boundary condition itself is biased.

The plateau of $\Lambda$ is therefore determined by both ends simultaneously.

## 3. Error definition

For each plateau parameter $p$, scan within its plateau $\mathcal{P}$, record the maximum and minimum of the output $y$, and take

$$
\Delta_{\rm plat}(p) \;=\; \frac{\max_{p\in\mathcal{P}} y(p) \;-\; \min_{p\in\mathcal{P}} y(p)}{2}
$$

Plateau parameters are treated as **uniform**, not Gaussian. The reason is that they have no measured value to serve as a center, only an admissible interval. There are two implementations, both uniform: scan point by point over a uniform grid within the plateau, or sample randomly and uniformly within it.

The distribution type for a plateau does not change between single-parameter and joint modes.

$y$ may be a ground-state decay constant, an excited-state decay constant, or a spectral functional. Calculate a separate row for each output type.

## 4. Three conventions

**First, use the output range, not the width of the parameter interval.**

The $\Lambda$ plateau might be $[20,30]$ GeV² (illustrative values), an interval of width 10 GeV², whereas $f_{J/\psi}$ might vary by only 2 MeV within it. Use half of the latter range, not $10/2$.

**Second, vary only one parameter at a time, holding the others at their baseline values.**

When scanning `Lambda`, keep `M2_min` and `M2_max` fixed at their baseline values. When scanning `M2_min`, keep `M2_max` and `Lambda` fixed at their baseline values.

The cost of this approach is that it neglects correlations between parameters. Both endpoints of the two-dimensional window can actually move together; separate scans cover only two edges of the rectangular region, not its corners. For a more rigorous treatment, see Section 7.

**Third, plateau boundaries must be justified and must not be adjusted afterward.**

The boundaries are determined by the physical conditions in Section 2. Do not narrow the interval to reduce the error or widen it to appear robust. The criteria used to choose the boundaries must be stated in the report.

## 5. Steps

Overall procedure:

```text
PROCEDURE platform_error(p, P_range, base, n_scan):
    y ← empty list
    FOR v IN linspace(P_range, n_scan):
        args ← base with p ← v          # 其余两个参数保持基准值
        (f1, f2, alpha, flags) ← SOLVE(args)
        IF alpha on boundary OR f1 not finite:
            record (v, reason); continue
        y.append((v, f1, f2))
    IF count(y) < 2:
        RETURN undefined
    Δ(p) ← (max(y) − min(y)) / 2
    RETURN Δ(p), count(y)
```

Call this once per parameter to obtain three $\Delta$ values. The following five steps expand this procedure.

### 5.1 Step 1: determine the baseline point and plateaus

First find a baseline window in the solver for each of the three channels. Take the baseline point at the plateau center, not at an endpoint.

```python
base = dict(M2_min=M2_MIN_CENTER, M2_max=M2_MAX_CENTER, Lambda=LAMBDA_CENTER)
```

Set the scan interval for each of the three parameters. Each interval must lie within the plateau and have a physical justification.

### 5.2 Step 2: single-parameter scans

Call `solve()` independently at each scan point. $\alpha$ is automatically reselected internally, so the scan code does not need to handle it.

```python
import numpy as np
from Jpsi_vector import solve

def scan_parameter(name, values, base):
    """扫描一个平台参数，其余两个固定在 base 中。"""
    rows = []
    for v in values:
        args = dict(base)
        args[name] = v
        res = solve(
            args["M2_min"], args["M2_max"], args["Lambda"],
            channel="2pole", mu=2.0, rg_loops=4,
        )
        ch = res["channels"]["2pole"]
        rows.append(dict(
            value=v,
            f1=float(ch["f"][0]),
            f2=float(ch["f"][1]),
            alpha=float(ch["alpha_star"]),
            on_boundary=bool(ch["alpha_diagnostics"]["on_boundary"]),
            rel_res=float(ch["relative_weighted_residual"]),
            neg_rho=int(ch["negative_rho_count"]),
        ))
    return rows
```

Example call:

```python
rows = scan_parameter("Lambda", np.arange(LAMBDA_LO, LAMBDA_HI + 1e-9, LAMBDA_STEP), base)
```

`eta_c_axial_vector` and `eta_c_pseudoscalar` have the same interface, differing only in the import path. In the pseudoscalar channel, distinguish `f` from `f_raw`: `f` includes the axial-anomaly correction, whereas `f_raw` does not. Use only one definition within a given error table.

### 5.3 Step 3: calculate the range

```python
def plateau_error(rows, key):
    values = np.array([r[key] for r in rows], float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        return np.nan, values.size
    return 0.5 * (values.max() - values.min()), values.size
```

`np.isfinite` filters out `NaN` values caused by negative pole residues. **The sample count after filtering must be recorded**; see Section 9.

### 5.4 Step 4: checks

Check the following three items at every scan point. A point that fails any one of them must not enter the range statistics.

```python
ok = [
    r for r in rows
    if not r["on_boundary"]
    and np.isfinite(r["f1"])
    and r["rel_res"] < RESIDUAL_TOL
]
```

Set `RESIDUAL_TOL` to a multiple of the baseline residual, for example 3 times that value.

### 5.5 Step 5: recording

For each scan, record at least the following quantities in the results table:

| Field | Meaning |
|---|---|
| Parameter name and baseline values | Which parameter is scanned and the values of the other two |
| Scan interval and number of points | Plateau boundaries and their justification |
| $\alpha_{\rm opt}$ at each scan point | Used to determine whether it lies within the stable window |
| `on_boundary` | Whether $\alpha$ lies at an endpoint of the candidate grid |
| `relative_weighted_residual` | Quality of the data fit |
| Number of valid samples | Number of points remaining after filtering |
| $\Delta_{\rm plat}$ | Final error |

## 6. Checks

### 6.1 Whether $\alpha$ lies at an endpoint

If `alpha_diagnostics.on_boundary` returned by `solve()` is true, $\alpha$ at that point was selected at one end of the candidate sequence. That point is unusable.

This is particularly important when comparing points: different $\Lambda$ values correspond to different solution intervals, and the stable window for $\alpha$ also moves. If some points leave that window, the scan has exceeded the method's domain of applicability.

### 6.2 Whether the residual is stable

`relative_weighted_residual` should remain of the same order of magnitude throughout the scan interval. A sudden increase by an order of magnitude at a point indicates a problem with the match between its OPE input and solution.

### 6.3 Whether the solution is unphysical

`negative_rho_count` counts grid points at which the reconstructed spectrum is negative. A small number of negative values is acceptable and should remain of the same order of magnitude as the parameter varies. A sudden increase in this count at a particular parameter value indicates that the point lies outside the plateau.

### 6.4 Whether the response is monotonic

The output should be nearly flat within the plateau. If it changes monotonically with the parameter and has an appreciable slope, the scan interval is not entirely within the plateau or the plateau was identified incorrectly.

Plotting $y$ against $p$ makes this problem easier to identify than inspecting numbers alone. Fig. 3 of the reference paper serves this purpose.

### 6.5 Whether pole residues are negative

When a residue is negative, `f` is set to `NaN`, not to zero or its absolute value. This affects the range statistics, so the number of excluded samples must be recorded.

## 7. Treatment of the two-dimensional window

`M2_min` and `M2_max` jointly define a two-dimensional window, and both endpoints can move simultaneously. The single-parameter scans in Section 5 cover only two edges of the rectangle, not its corners.

Three approaches, ordered by rigor:

**Approach 1: single-parameter scans (the approach of the reference paper).**

Scan each edge separately and combine the two $\Delta$ values in quadrature. This has the lowest computational cost, but neglects the corners and implicitly assumes that the two endpoints are independent.

**Approach 2: a two-dimensional grid.**

Place a grid within the rectangle and perform a calculation at every point.

```python
for m_lo in np.arange(lo_plateau[0], lo_plateau[1] + 1e-9, 1.0):
    for m_hi in np.arange(hi_plateau[0], hi_plateau[1] + 1e-9, 1.0):
        if m_hi <= m_lo:
            continue
        ...
```

Take half of the two-dimensional range as the final error:

$$
\Delta_{\rm window} = \frac{\max y - \min y}{2}
$$

This quantity already contains the joint effect of the two endpoints. **Do not add it to the single-parameter results**, which would double count the effect.

**Approach 3: random sampling within the plateau.**

Draw points uniformly within the rectangle and take the range. With too few points, the range underestimates the true span, so a sufficiently large sample is needed.

The three approaches should give similar conclusions. If approach 2 gives a substantially larger result than approach 1, corner effects are significant; the report must state which approach was used.

## 8. Reporting format

Use a separate table for each channel.

| Parameter | Plateau | Basis | $\Delta f_1$ [MeV] | $\Delta f_2$ [MeV] | Valid points |
|---|---|---|---|---|---|
| `M2_min` | $[\ldots,\ldots]$ GeV² | OPE convergence | — | — | — |
| `M2_max` | $[\ldots,\ldots]$ GeV² | Ground-state dominance | — | — | — |
| `Lambda` | $[\ldots,\ldots]$ GeV² | OPE matching | — | — | — |

Fill the blanks with the calculated values. Below the table, state:

- The baseline point $(M^2_{\min}, M^2_{\max}, \Lambda)$.
- That $\alpha$ is reselected at every point, together with the selection criterion.
- The excluded points and reasons for exclusion.
- Which window-treatment approach from Section 7 was used.

Provide separate tables for the three channels (`Jpsi_vector`, `eta_c_axial_vector`, and `eta_c_pseudoscalar`); do not combine them into one row. Their pole models and residue conversions differ, so their errors are not interchangeable.

## 9. Points to note

**1. Parameter-interval width is not the error.** See the first convention in Section 4. This is the most common mistake.

**2. Plateau boundaries must not be adjusted afterward.** Determine boundaries using physical criteria first, then calculate the error. Reversing this order makes the error meaningless.

**3. Filtering `NaN` values introduces selection bias.** Once points with negative residues are excluded, the remaining points form a conditional sample. Report both the valid and original point counts.

**4. $\alpha$ is outside the scope of this file.** This file covers only plateau parameters; the error associated with $\alpha$ is discussed separately.

**5. The two $\Delta$ values must not be counted twice.** If a two-dimensional grid is used to obtain the joint span directly, do not add the single-parameter results.

**6. Scan step size affects the range.** The range is $\max-\min$; coarser steps are more likely to miss the true extrema. Use a sufficiently fine step size and state it in the report.
