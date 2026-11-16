# Input-sampling error

[中文](input-sampling.md) | **English**

Input error originates outside the method. It does not depend on the inverse-problem solution procedure itself and remains present even if a different inverse solver is used.

This file covers only **Gaussian sampling**. The other three error types have separate files:

| Error type | File |
|---|---|
| Plateau parameters | [platform-parameters.md](platform-parameters.en.md) |
| Gaussian input sampling | This file |
| OPE truncation | [ope-truncation.md](ope-truncation.en.md) |
| Renormalization scale $\mu$ | [mu-scale.md](mu-scale.en.md) |

> **This document records the method only.** The intervals in the text and numerical values in the code are examples or current settings, not final choices. The input central values and $1\sigma$ errors have been confirmed and can be used directly.

## 1. Sampling list

Five inputs enter Gaussian sampling in this project.

| Input | Code location | Central value | $1\sigma$ | Source |
|---|---|---|---|---|
| $G_2$ | `physics.G2` | $0.038\ \mathrm{GeV}^4$ | 0.011 | 1801.06832 |
| $G_3$ | `physics.G3` | $0.013\ \mathrm{GeV}^6$ | 0.007 | 1801.06832 |
| $\bar m_c(\bar m_c)$ | `RGInputs.mbar_c` | 1.275 GeV | 0.020 GeV | Project input table |
| $M_{1S}$ | `MASSES` | See table below | See table below | PDG |
| $M_{2S}$ | `MASSES` | See table below | See table below | PDG |

Pole masses differ by channel:

| Channel | $M_{1S}$ | $\sigma$ | $M_{2S}$ | $\sigma$ |
|---|---|---|---|---|
| `Jpsi_vector` | $J/\psi$ 3.096900 | $6\times10^{-6}$ GeV | $\psi(2S)$ 3.686097 | $6\times10^{-5}$ GeV |
| `eta_c_axial_vector` | $\eta_c$ 2.9841 | $3\times10^{-4}$ GeV | $\chi_{c1}$ 3.51067 | $5\times10^{-5}$ GeV |
| `eta_c_pseudoscalar` | $\eta_c$ 2.9841 | $3\times10^{-4}$ GeV | $\eta_c(2S)$ 3.6378 | $6\times10^{-4}$ GeV |

Variations in $\bar m_c$ propagate through the RG and change both $m_c(\mu)$ and the perturbative spectrum, not merely the explicit mass factors in the OPE. The previous budget omitted the latter contribution, underestimating the input error. When recalculating, ensure that every location using $m_c$ receives the updated value.

Inputs not included in the sampling covered by this file:

| Input | Treatment | File |
|---|---|---|
| $\alpha_s(2\ \mathrm{GeV})=0.3029$ | Fixed, with no error assigned | This file, see below |
| $M_Z$, $\bar m_b$ | RG internal endpoints only | This file |
| $\mu$ | Fixed or scanned diagnostically | [mu-scale.md](mu-scale.en.md) |
| $a_1=-0.079$ GeV³ | No error assigned for now | This file |
| OPE truncation | Difference | [ope-truncation.md](ope-truncation.en.md) |

$\alpha_s$ has already been evolved by the RG to the calculation scale, and the project convention assigns it no independent error. The source of $a_1$ must be established before it can be sampled.

## 2. Distributions

### 2.1 Single parameter

Given a central value $\theta_0$ and error $\sigma$, draw

$$
\theta = \theta_0 + \sigma z, \qquad z \sim \mathcal{N}(0,1)
$$

The error is interpreted as $1\sigma$. This is the standard meaning of “the measured value follows a normal distribution”.

### 2.2 Asymmetric errors

When the upper and lower errors differ, use a split normal: first draw $z$, then select the error according to its sign.

$$
\theta = \theta_0 + z \times
\begin{cases}
\sigma_{-}, & z < 0\\
\sigma_{+}, & z \ge 0
\end{cases}
$$

Each side contains half of the probability mass. The central value is therefore the median, not the mean. Both central values must be clearly stated in the report.

All five current inputs in this project have symmetric upper and lower errors, so this provision is not needed at present. It is retained for future asymmetric inputs.

### 2.3 Correlations

When multiple inputs are correlated, use a correlation matrix $C$. First perform the eigendecomposition $C=V\sqrt{D}$, then draw

$$
\boldsymbol\theta = \boldsymbol\theta_0 + V\sqrt{D}\,\boldsymbol z,
\qquad \boldsymbol z\sim\mathcal{N}(\mathbf 0,I)
$$

$C$ must be symmetric and positive semidefinite, with unit diagonal entries. Correlations must have a justified source rather than being assigned intuitively.

**This project currently treats the inputs as independent.** The report must explicitly state “independence assumed”.

## 3. Two sampling modes

The same inputs can be organized in two ways. Both use Gaussian sampling; they differ in how many inputs are varied at once.

**Distribution type is independent of sampling mode.** Single-parameter and joint sampling are only organizational choices: all five inputs remain Gaussian, regardless of the mode. Plateau parameters always remain uniform; see [platform-parameters.md](platform-parameters.en.md).

### 3.1 Mode 1: single-parameter sampling (source by source)

Vary only one input at a time, holding the others at their central values. Each input produces its own output distribution, and the components are then combined in quadrature.

```text
PROCEDURE per_source(inputs, N):
    FOR each input i:
        values ← []
        FOR k IN 1..N:
            θ ← inputs with only θ_i resampled
            y ← SOLVE(θ)
            IF finite(y): values.append(y)
        Δ_i ← std(values)
    Δ_total ← sqrt(Σ_i Δ_i²)
    RETURN {Δ_i}, Δ_total
```

Features:

- Produces a component-wise table showing which source dominates.
- Implicitly assumes independent inputs and approximately linear responses.
- Computational cost increases linearly with the number of inputs.

### 3.2 Mode 2: joint sampling

Draw all inputs simultaneously. Each sample provides a complete input set for one inversion.

```text
PROCEDURE joint(inputs, N):
    values ← []
    FOR k IN 1..N:
        θ ← resample every input at once
        y ← SOLVE(θ)
        IF finite(y): values.append(y)
    Δ_total ← std(values)
    RETURN Δ_total
```

Features:

- Provides only a total uncertainty, without a component-wise breakdown.
- Automatically includes correlations between inputs and nonlinear responses.
- Computational cost is independent of the number of inputs and depends only on the sample count.

### 3.3 Relation between the two modes

| Situation | Single-parameter result | Joint result |
|---|---|---|
| Independent, linear | $\sqrt{\sum\Delta_i^2}$ | Agrees with the result on the left |
| Correlated | Quadrature sum can be too large or too small | Directly gives the joint distribution |
| Clearly nonlinear | Only a first-order approximation | Closer to the true distribution |

The modes therefore serve different purposes rather than replacing one another.

### 3.4 Choosing a mode

- Use mode 1 to produce an error budget listing individual sources.
- Use mode 2 to report a total uncertainty or when inputs are correlated.
- Use both: mode 1 for the component-wise table and mode 2 to check the total. A large discrepancy means that the linearity or independence assumption fails, which must be explained in the report.

In this project, the $\bar m_c$ responses of 1S and 2S have opposite directions and differ in magnitude by more than a factor of two. This is a clearly nonlinear case, so using both modes is recommended.

## 4. Algorithm

```text
PROCEDURE input_error(theta0, sigma_minus, sigma_plus, n_samples, seed):
    rng.seed(seed)
    values ← empty list
    FOR i IN 1..n_samples:
        z ← N(0,1)
        theta ← theta0 + z * (sigma_minus if z < 0 else sigma_plus)
        IF theta violates physical constraint:
            record and resample
        y ← SOLVE(theta)
        IF y is finite:
            values.append(y)
    RETURN mean(values), std(values), quantiles(values, [0.025,0.1587,0.5,0.8413,0.975])
```

Key points:

- `SOLVE` automatically reselects $\alpha$ internally on every call. The sampling code does not intervene in $\alpha$ selection.
- Constraints are checked only when inputs are generated. Samples that fail after entering the inversion are not replaced.
- Use the sample standard deviation (denominator $n-1$), not the standard error.
- Record both the number of valid samples and the number of rejected samples.

## 5. Statistics

### 5.1 Output statistics

Calculate separately for each output quantity:

```text
mean     ← arithmetic mean of finite values
sd       ← sample standard deviation (n-1)
q_lo     ← 15.87th percentile
q_med    ← median
q_hi     ← 84.13rd percentile
err_minus ← q_med - q_lo
err_plus  ← q_hi - q_med
```

Report the sample standard deviation, not the standard error. The aim is to measure “how much the result varies when the inputs change”, not “the precision of the mean”. Increasing the sample count does not reduce this quantity.

### 5.2 Checks

Check four items after running each parameter:

1. **Number of valid samples.** If it is substantially smaller than `n_samples`, the reported statistics represent only a subset.
2. **Distribution shape.** Check whether the histogram is unimodal and whether it is strongly skewed. Quantiles are more reliable than the standard deviation when skewness is pronounced.
3. **Unphysical outputs.** Count the `NaN` values caused by negative pole residues. This count itself varies with the parameter and must be recorded.
4. **Whether $\alpha$ leaves the plateau.** If $\alpha$ is selected at an endpoint for some samples, their reliability is reduced.

## 6. Reporting format

Use a separate table for each channel and each output quantity.

| Input | Central value | $1\sigma$ | $\Delta y_1$ | $\Delta y_2$ | Valid samples |
|---|---|---|---|---|---|
| $G_2$ | $0.038\ \mathrm{GeV}^4$ | 0.011 | — | — | — |
| $G_3$ | $0.013\ \mathrm{GeV}^6$ | 0.007 | — | — | — |
| $\bar m_c(\bar m_c)$ | 1.275 GeV | 0.020 GeV | — | — | — |
| $M_{1S}$ | See Section 1 | See Section 1 | — | — | — |
| $M_{2S}$ | See Section 1 | See Section 1 | — | — | — |

Below the table, state the sampling mode, sample count, random seed, distribution types, correlation assumptions, constraint-rejection rate, invalid-sample count, output units, and whether the individual components were combined in quadrature.

## 7. Points to note

**1. The number of valid samples must be reported.**

Samples with negative pole residues are set to `NaN` and excluded. The statistics after exclusion describe a conditional sample.

**2. Standard deviation is not standard error.**

The quantity of interest is the spread after input propagation. Do not divide by $\sqrt{n}$.

**3. Report two values for asymmetric errors.**

Report $\Delta_-$ and $\Delta_+$ separately relative to the median.

**4. Do not introduce correlations without justification.**

Assuming independence is acceptable, but must be stated. If a correlation matrix is used, the source of its correlation coefficients must be explained.

**5. Combining in quadrature requires independence and approximate linearity.**

Before combining single-parameter results in quadrature, check that the response is approximately linear. If a parameter causes the output to cross a sign change or threshold, the quadrature sum is not valid.

**6. The coupling between $S_H$ and $\bar m_c$ remains to be confirmed.**

Currently, `S_H` is a hard-coded constant in all three channels and does not change with `mbar_c`. If the physical definition of $S_H$ is $4m_c^2$, it should also move when $\bar m_c$ changes. This must be confirmed before deciding whether to treat $S_H$ as an independent input.
