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

## Reusable analysis Skill

[`uncertainty-analysis/SKILL.md`](uncertainty-analysis/SKILL.md) describes a reusable separated-ensemble uncertainty workflow for AI agents. Statistical inputs and systematic platform/window parameters are sampled in separate Monte Carlo runs and reported separately; the tool does not combine them into a total uncertainty. It does not prescribe the model's kinematic variables, condensates or number of observables. The Skill is stored in this repository and is not installed into a local agent directory.

### Uncertainty-analysis TODO

- Add an explicit regularization-parameter (`alpha`) ensemble or diagnostic; the current reference workflow keeps `alpha` fixed during both Monte Carlo ensembles.
- Add running-scale (`mu`) variation; the current charmonium workflow fixes `mu = 2.0` and does not sample it.
- Track further uncertainty-method changes in future pull requests rather than silently changing completed ensembles or reports.
