# Error analysis for inverse problems

[中文](README.md) | **English**

This folder contains methodological documents on error analysis. No general-purpose program is provided: the operator, grid, and outputs differ across inverse problems, so each project implements its own code; these documents provide steps, criteria, and code examples.

## 1. Two categories of error

Errors are divided into two categories by **source**, not into statistical and systematic errors.

| Category | Source | Includes |
|---|---|---|
| Input-induced | External to the method | Condensates, heavy-quark mass, pole masses, OPE truncation ($\mu$ is fixed, with no error assigned) |
| Method-intrinsic | The method itself | Plateau parameters `M2_min`, `M2_max`, and `Lambda`; regularization approximation $\alpha$ |

In experimental and lattice analyses, the statistical/systematic classification depends on whether an error is random or bias-dominated. The present classification depends on whether the error originates outside or inside the method. These are different perspectives, without a one-to-one correspondence, and the two categories should not be mapped to statistical/systematic errors.

## 2. Basic inequality for error estimation

The difference between the regularized solution and the true solution satisfies

$$\|f_\alpha^\delta - f\| \;\le\; \delta\,\|R_\alpha\| \;+\; \|R_\alpha K f - f\|$$

- The first term is determined by the input error $\delta$ and corresponds to input-induced error.
- The second term is determined by the finite regularization parameter $\alpha$ and corresponds to the regularization-approximation component of method-intrinsic error.

This inequality provides the mathematical basis for the classification and a way to estimate each term. The bound need not be tight in practice, but each term has a corresponding computable quantity.

## 3. Prerequisite criteria

Two criteria must be met before performing error analysis:

1. **Plateau criterion.** Each of the three method parameters (`M2_min`, `M2_max`, and `Lambda`) has a stable plateau.
2. **Endpoint criterion.** The solution approaches the OPE spectrum at the large-$s$ endpoint.

If either criterion fails, the solution cannot be treated as a physical result, and error analysis loses its meaning. Even if numerical values can be obtained, the errors then do not cover the true value.

## 4. Distribution assignments

The distribution assigned to each error source is determined by the nature of that source.

| Error source | Distribution | Reason |
|---|---|---|
| $G_2$, $G_3$, $\bar m_c$, $M_{1S}$, $M_{2S}$ | Gaussian ($1\sigma$) | The literature or PDG provides central values and errors |
| `M2_min`, `M2_max`, `Lambda` | Uniform (within the plateau) | No unique physical value, only an admissible interval |
| $\mu$ | Not sampled (fixed) | A scheme choice, not a measurement |
| OPE truncation | Not sampled (difference) | No distribution to sample |
| $\alpha$ | Not sampled (forward–inverse evaluation of the mean solution) | Determined by the regularization operator itself, not a measurement distribution |

**These assignments are independent of the single-parameter/joint sampling modes.** The two modes merely organize the same distributions into different sampling schemes; they do not change the distribution type of any source.

For details, see [input-sampling.md](input-sampling.en.md), [platform-parameters.md](platform-parameters.en.md), [ope-truncation.md](ope-truncation.en.md), and [mu-scale.md](mu-scale.en.md).

## 5. Document index

One file per error type.

| Document | Error type | Category | Status |
|---|---|---|---|
| [platform-parameters.md](platform-parameters.en.md) | Plateau parameters `M2_min`, `M2_max`, and `Lambda` | Method-intrinsic | Written |
| [alpha-regularization.md](alpha-regularization.en.md) | Regularization approximation $\alpha$ (forward–inverse evaluation of the mean solution) | Method-intrinsic | Written |
| [input-sampling.md](input-sampling.en.md) | Gaussian input sampling (single-parameter/joint modes) | Input-induced | Written |
| [ope-truncation.md](ope-truncation.en.md) | OPE truncation | Input-induced | Written |
| [mu-scale.md](mu-scale.en.md) | Renormalization scale $\mu$ (fixed, with no error assigned) | Input-induced | Written |
| [reference-comparison.md](reference-comparison.en.md) | Differences in input lists relative to the reference project | — | Written |
| reporting.md | Budget-table format, combination rules, and required reporting fields | — | To be written |
