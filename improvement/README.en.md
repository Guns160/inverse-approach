# improvement: examples of improvements to the equation-solving step (compact GSVD filtering vs per-alpha normal equations)

[中文](README.md) | **English**

This folder contains only minimal runnable examples of the **equation-solving step**, illustrating this repository's improvements over the original solver. The two examples (Python + MATLAB) construct **the same** synthetic inverse problem, with outputs that can be compared item by item.

## Original and improved methods

| | Original method | Improved method (shared by the three programs in this repository) |
|---|---|---|
| Source | Original code (the original equation-solving code is quoted in the header comments of `matlab/gsvd_example.m`) | Shared solver core in [`../charmonium decay constants/`](../charmonium%20decay%20constants/) |
| Assembly | Normal equations `A2=A'WA`, `D2=A'Wb`, and penalty `H2` | Retain the tall matrix `A=sqrt(w)*K_free`, penalty rows `B=[0, chol(H_ii)]`, and eliminate spectral endpoints |
| Solution | Reassemble and fully factorize `(A2+alpha*H2)` for each alpha (earlier versions used `inv(A2+alpha*H2)*D2`) | **One** economy QR + compact GSVD, followed by scalar filtering for any alpha |
| Full alpha grid | Solve point by point, with nAlpha full factorizations (the L-curve / quasi-optimality criterion requires solutions over the entire grid) | One factorization + vectorized filtering, with `X'\[coef; un]` solving for the entire grid in a single back-substitution |
| Diagnostics | Recompute residuals / penalties by back-substitution at each point | Analytical expressions for the residual norm `norm(C*y-beta)^2+null` and penalty norm `norm(S*y)` |

Filtering formula (the same as in `solveInverse.m` and `reconstruction.py`):

$$
y_i(\alpha)=\frac{\sigma_i\,\beta_i}{\sigma_i^2+\alpha\,\mu_i^2}\ \ (i\le p),
\qquad
y_i=\frac{\beta_i}{\sigma_i}\ \ (i>p,\ \sigma_i=1),
\qquad
x_\alpha=X\,y,
$$

Here `(sigma_i, mu_i)` are the generalized singular values of the matrix pair `(A, B)`, and `beta=U'*(Q'*b)`. The normal equations `(A'A+alpha*B'B)x=A'b` solved by the old method correspond to **the same** Tikhonov functional; only the solver differs in the two examples, making the comparison fair.

Corresponding production implementations:

- Python: `../charmonium decay constants/inverse_common/gsvd.py` (LAPACKE `dggsvd3`), `factors.py` (assembly), and `reconstruction.py` (filtering over the full alpha grid)
- MATLAB: `../charmonium decay constants/matlab/+charmonium/compactGsvd.m`, `factorsOneWindow.m`, and `solveInverse.m`

## Why it is faster

- Selecting alpha with the original method requires solving at every point of the grid: 1000 points in the original code mean 1000 full O(n^3) factorizations (each also solving for the three right-hand sides `[D2, e1, en]` simultaneously).
- The improved method performs the GSVD only once; each subsequent alpha solution uses O(p) scalar filtering, and the entire grid is evaluated at once. At the example size (n=161, 1001 alpha values), the measured filtering time is only 1.8–3.7 ms, compared with 103–140 ms for point-by-point solves with the old method.

## Why it is more accurate

- Normal equations square the condition number: in the example, `cond(A)≈1e+18`, so theoretically `cond(A'A)≈cond(A)^2≈1e+36`; the smallest singular directions of `A'A` have already been overwhelmed by floating-point rounding. At small alpha, measured relative errors of the old method are 30–40%, while the new method remains at the 1e-9 level (see the table below).
- No explicit matrix inverse is used; the GSVD-filtered solution agrees with the direct augmented least-squares solution `[A; sqrt(alpha)*B]\[b; 0]` to the order of machine precision (this is precisely the `CheckSolution` check in the production program `solveInverse.m`).
- Analytical expressions are available for the residual / penalty norms; these are more stable than direct recomputation at very small alpha (at the smallest alpha in the example, direct recomputation has cancellation errors of order 4e-7–3e-6, whereas the analytical expressions do not have this problem).

## Measurements (2026-10-15)

Environment: Python 3.9.6 + NumPy 2.0.2 (GSVD calls the LP64 `LAPACKE_dggsvd3` in SciPy's bundled OpenBLAS); MATLAB R2026a (built-in `gsvd`). Problem size: 301 Borel data points, 161 unknowns (2 unpenalized poles + 159 interior hat functions), 159 penalty rows, and 1001 alpha grid points (logspace 1e-16…1e-2). Errors are relative errors with respect to the direct augmented least-squares solution:

| alpha | Old method (Python / MATLAB) | New method (Python / MATLAB) |
|---|---|---|
| 1e-16 | 3.9e-01 / 3.0e-01 | 8.3e-09 / 2.3e-09 |
| 3.2e-13 | 1.9e-04 / 1.0e-04 | 2.0e-10 / 5.8e-11 |
| 1e-09 | 6.3e-08 / 2.5e-08 | 1.3e-12 / 3.1e-12 |
| 3.2e-06 | 5.8e-10 / 2.0e-10 | 3.1e-13 / 4.7e-13 |
| 1e-02 | 2.1e-12 / 3.9e-12 | 2.3e-14 / 9.3e-15 |

Timing (the entire grid of 1001 alpha values):

| | Old method (point-by-point solves) | New method (one GSVD + full-grid filtering) | Speedup |
|---|---|---|---|
| Python (LAPACKE dggsvd3) | 103.3 ms | 56.1 ms + 1.8 ms | 1.8-fold |
| MATLAB (built-in gsvd) | 140.2 ms | 5.0 ms + 3.7 ms | 16.2-fold |

Note: factorization is a one-time cost in the new method, so amortization becomes more pronounced for longer alpha grids and larger problems (old-method cost ∝ nAlpha*n^3, new-method filtering ∝ n^2*nAlpha + one GSVD). The example is demonstration-sized; the production program defaults to 2001 alpha grid points, and the historical production spectral grid (ds=0.01) is larger.

Other self-checks:

- Analytical penalty norm vs direct recomputation: Python 6.0e-15, MATLAB 7.9e-16.
- Independent MATLAB cross-check using `cgsvd`+`fil_fac` from Hansen's Regularization Tools: relative difference in generalized singular values (sigma/mu) 0.0, maximum relative difference in solutions 3.8e-10.
- Both languages construct the same problem (grid values differ only in the last few rounding digits), and the outputs agree to rounding / order-of-magnitude precision.

## Running the examples

Python (requires NumPy; GSVD requires LP64 `LAPACKE_dggsvd3`, with automatic detection of system / NumPy / SciPy OpenBLAS, or an explicit path through the `JPSI_LAPACKE_LIBRARY` environment variable. SciPy does not expose a gsvd interface, so the example calls LAPACKE directly through ctypes, as does the production implementation):

```sh
python3 python/gsvd_example.py
```

MATLAB (only base MATLAB is required; the final cgsvd cross-check requires Regularization Tools and is automatically skipped if it is not installed):

```matlab
matlab -batch "run('matlab/gsvd_example.m')"   % 或在 MATLAB 中 cd 到 matlab/ 后运行 gsvd_example
```

## Relation to Regularization Tools

The ordering convention of the compact GSVD follows `cgsvd` in Hansen's *Regularization Tools* (as also noted in the documentation of the production implementation `gsvd.py`). Convention difference: that package writes the `'Tikh'` filter factor as `sigma^2/(sigma^2+lambda^2*mu^2)`; alpha in this repository corresponds to `lambda^2`.

## Scope

The examples construct synthetic data using a smooth illustrative kernel (`exp(-s*t/60)`) and demonstrate only the equation-solving step; for the production physics pipeline (analytical hatKernel, nonzero-boundary q shift, Raus / quasi-optimal alpha selection, three-channel physics and RG conventions), see [`../charmonium decay constants/`](../charmonium%20decay%20constants/).
