# Inverse approach

Research code for extracting charmonium pole residues and decay constants from Borel sum rules.

The programs are in [`charmonium decay constants/`](charmonium%20decay%20constants/):

| Program | Current | Explicit pole models |
|---|---|---|
| `Jpsi_vector` | Vector | J/psi + psi(2S) |
| `eta_c_axial_vector` | Axial vector | eta_c + chi_c1 |
| `eta_c_pseudoscalar` | Pseudoscalar density | eta_c + eta_c(2S) |

The second axial-vector pole is chi_c1, not eta_c(2S). The programs preserve the separate OPE and residue conventions of these three channels.
All three programs use the model identifier `2pole`.

Numerical outputs depend on the running-coupling prescription, spectral grid, Borel window, pole model and regularization rule. A successful solve by itself does not establish a physical stability window.

## Version history

- `ed20389`: original Python source upload.
- `06fb9a2`: particle/current directory layout and corrected entry-point imports.
- Later commits separate the common-core refactor, running-scale update and MATLAB implementation.

Historical results and parameter scans are not bundled with the source packages.

## Solver improvement example

[`improvement/`](improvement/) contains a minimal self-contained Python + MATLAB example of the solver improvement: one economy-QR + compact GSVD replaces the original per-alpha normal-equation solves, and the whole alpha grid is produced by vectorized generalized filters. The MATLAB example quotes the original per-alpha solve loop it replaces.

## Uncertainty analysis notes

[`uncertainty-analysis/`](uncertainty-analysis/) records the uncertainty-analysis method for this class of inverse problems. Errors are separated by where they come from — outside the method versus inside it — and each kind has its own note:

| Note | Error kind | Class |
|---|---|---|
| [`input-sampling.md`](uncertainty-analysis/input-sampling.md) | Gaussian sampling of physical inputs | Input-induced |
| [`ope-truncation.md`](uncertainty-analysis/ope-truncation.md) | Truncation of the OPE | Input-induced |
| [`mu-scale.md`](uncertainty-analysis/mu-scale.md) | Renormalization scale (fixed, not an error) | Input-induced |
| [`platform-parameters.md`](uncertainty-analysis/platform-parameters.md) | Platform parameters | Method-intrinsic |
| [`alpha-regularization.md`](uncertainty-analysis/alpha-regularization.md) | Regularization approximation | Method-intrinsic |
| [`reference-comparison.md`](uncertainty-analysis/reference-comparison.md) | Input-list differences from the reference study | — |

Each note gives the procedure, the checks and the reporting format. The implementation stays in each project, since the operator, grid and output quantities differ. [`README.md`](uncertainty-analysis/README.md) is the index.

The notes also record the distribution assigned to each source: Gaussian for inputs with a measured central value and uncertainty, uniform over the plateau for platform parameters, and no sampling at all for the renormalization scale, the OPE truncation and the regularization approximation. That assignment does not depend on whether the sampling is per-source or joint.
