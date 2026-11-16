# Regularization-approximation error $\alpha$

[中文](alpha-regularization.md) | **English**

A finite $\alpha$ makes the regularization operator differ from the exact inverse operator. This section explains how to estimate that bias.

This is the other component of method-intrinsic error besides the plateau parameters. For plateau parameters, see [platform-parameters.md](platform-parameters.en.md).

> **This document records the method only.** The numerical values given here are current conventions, not final choices.

## 1. Quantity to calculate

$$
\Delta_\alpha \;=\; \|R_\alpha K f - f\|
$$

$R_\alpha$ is the regularized inverse operator, $K$ is the forward operator, and $f$ is the true solution.

As $\alpha\to0$, $R_\alpha K\to I$, and this term vanishes. The larger $\alpha$ is, the more small components regularization suppresses, and the larger this term becomes. It therefore measures the “bias introduced by finite $\alpha$”.

The true solution $f$ is unknown and is replaced by the numerical solution $\hat f$:

$$
\Delta_\alpha \;\approx\; \|R_\alpha K \hat f - \hat f\|
$$

## 2. Steps

Five steps.

**Step 1: sample to obtain an ensemble of solutions.**

Sample according to [input-sampling.md](input-sampling.en.md). Run a complete inversion for each sample to obtain one solution.

**Step 2: take the mean solution $\hat f$ of this ensemble.**

Take the mean in solution-vector space, not in output-quantity space. For example, if a decay constant is $\sqrt{\Pi}$ multiplied by 1000, first average $\Pi$ and then take the square root; do not take square roots first and then average.

**Step 3: apply the forward problem to $\hat f$ to obtain synthetic data.**

$$
\text{合成数据} = K\hat f
$$

Physical meaning: assume a dataset exists whose exact solution is precisely $\hat f$, with no measurement error.

**Step 4: invert the synthetic data once using the same inverse solver and the same $\alpha$.**

```
新解 = R_alpha · 合成数据
```

This step checks whether the regularized inverse solver reproduces the solution unchanged when fed data that this solution can already explain exactly.

**Step 5: take the difference between the two solutions.**

$$
\Delta_\alpha = |\text{新解} - \hat f|
$$

The part that cannot be recovered is the bias introduced by $\alpha$.

## 3. Algorithm

```text
PROCEDURE alpha_error(samples, alpha_opt, K, INVERT):
    f_hat    ← MEAN(samples)              # 在解向量空间取均值
    合成数据 ← K · f_hat                  # 走一遍正问题
    新解     ← INVERT(合成数据, alpha_opt)  # 同一个 alpha、同一套边界条件
    RETURN |新解 − f_hat|
```

Project the difference in the solution vector onto the reported physical quantities:

```text
FOR each output q IN {y_1, y_2, spectrum}:
    Δ_alpha(q) ← PROJECT(|新解 − f_hat|, q)
```

`PROJECT` must be consistent with the definition in the main program. If the output is $\sqrt{\Pi}\times1000$, take the difference using that same function; do not compare $\Pi$ in one place and $\sqrt{\Pi}$ in another.

## 4. Implementation details

**1. Use the mean of the solution vectors, not the solution at the central input parameters.**

These differ when the response is nonlinear. The former represents the “center of this ensemble”; the latter is a separate calculation.

**2. Use the $\alpha$ value employed in the actual run.**

Do not use the plateau center or perform a scan. $\Delta_\alpha$ varies with $\alpha$, so the value of $\alpha$ used must be recorded for reproducibility.

**3. The inversion procedure must be identical to that of the main program.**

Use the same forward operator, penalty matrix, boundary conditions, and $\alpha$ selection rule. Any difference between the settings of the two inversions introduces other contributions into the difference.

**4. Record the size and source of the sampled ensemble.**

The report must state which ensemble produced $\hat f$, how many samples it contains, and which distributions were used. A different ensemble gives a slightly different mean.

## 5. Checks

**1. Number of valid samples.** The mean should include only valid samples. Samples set to `NaN` do not enter the average.

**2. Magnitude of the difference.** Normally, this component should be much smaller than the input error. If it is instead large, either the $\alpha$ selection is problematic or the inversion procedure differs from that of the main program.

**3. Independent cross-check.** Also calculate the difference between the solution at the current $\alpha$ and the solution at a smaller $\alpha$, and compare it with the forward–inverse estimate. Their magnitudes should not differ greatly. A difference of an order of magnitude indicates that the forward–inverse estimate is unreliable.

## 6. Risk of bias

This estimate can be **systematically too small**, which must be stated in the report.

The reason is that $\hat f$ is itself already a regularized solution whose small components have been suppressed. Applying $R_\alpha K$ to it again removes less than applying the operator to an unregularized true solution. Thus

$$
\|R_\alpha K\hat f-\hat f\| \;\lesssim\; \|R_\alpha K f-f\|
$$

This difference is a lower-bound-type estimate, not a rigorous upper bound.

## 7. Reporting format

| Component | $\alpha$ value | $\Delta y_1$ | $\Delta y_2$ | Samples used for the mean |
|---|---|---|---|---|
| Regularization approximation | — | — | — | — |

Below the table, state which sampled ensemble produced $\hat f$ (mode, sample count, and seed), the $\alpha$ used, how consistency of the forward and inversion procedures with the main program was checked, and that this is a “lower-bound-type estimate”.
