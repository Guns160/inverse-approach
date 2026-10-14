# Validation record: 2026-10-14

These checks validate implementations and selected numerical contracts. They do not establish new physical extraction windows or replace the existing error-budget analysis.

## Environment

- Python 3.9.6; NumPy 2.0.2; SciPy 1.13.1.
- Python GSVD: LP64 `LAPACKE_dggsvd3` from the existing SciPy OpenBLAS library.
- MATLAB R2026a Update 4, version 26.1.0.3312084, built-in `gsvd`.
- macOS host. The same Python code and library discovery have not been validated on other operating systems in this run.

## Separate changes

1. Original source/layout: repository commit `06fb9a2`.
2. Shared Python infrastructure: repository commit `bfe62ea`. Old input constants were retained. At `(M2_min,M2_max,Lambda)=(5,20,20)`, `ds=0.1`, the three double-pole models matched the original arrays exactly for M2, s, OPE, continuum, residues, f (including NaNs), and back substitution. Ten tests passed at this stage.
3. Selectable running and generalized P-P filter: repository commit `849f178`. This intentionally changes numerical results. It is not claimed to reproduce the original frozen-input physical outputs.
4. MATLAB implementation and cross-language tests are recorded in the subsequent commit.

## RG references

The four-loop mass and coupling use direct integration of the truncated beta and gamma functions. An independent 55-digit inverse-beta quadrature and gamma/beta integral gave:

| Quantity | Independent reference |
|---|---:|
| alpha_s^(5)(4.18 GeV) | 0.226141845816242854 |
| alpha_s^(4)(4.18 GeV) | 0.226374296730276341 |
| alpha_s^(4)(1.275 GeV) | 0.391619277160514745 |
| alpha_s^(4)(2 GeV) | 0.304398404402813591 |
| mc^(4)(2 GeV), GeV | 1.093276847184784040 |

Python and MATLAB ODE results agree with these references to better than 4e-13 in the displayed dimensionless/GeV units. Additional digits are numerical checks, not physical input precision. Tests also recover the prior three-loop values, verify a one-loop analytic solution, round trips, matching direction, moving mass/scale anchors, invalid inputs, and independent JSON CLI execution.

## Python checks

Forty-one Python tests pass with warnings treated as errors. The optional MATLAB test is separately gated. The tests cover physical back substitution, exact requested grid endpoints, nonzero boundary shifts, unpenalized poles, standard/generalized singular-value conversion, and QOC finite differences. A well-conditioned artificial problem with nontrivial GSVD mu values is compared against an independently assembled augmented least-squares solution; this catches the original P-P missing-mu-squared error.

The three `run.py` entry points were also run from a different working directory, with four-loop RG, JSON output and pickle-free NPZ output. All three used the same derived mass/coupling.

Additional default-grid runs, with the double-pole model and four-loop inputs, completed with finite spectra and residues:

| Package | M2_min | M2_max | Lambda | ds | Hat nodes |
|---|---:|---:|---:|---:|---:|
| Jpsi_vector | 18.5 | 80 | 68 | 0.1 | 542 |
| eta_c_axial_vector | 18.5 | 80 | 68 | 0.1 | 532 |
| eta_c_pseudoscalar | 15.5 | 80 | 67.5 | 0.1 | 527 |

These are execution checks, not a full scan or a re-estimation of the uncertainties.

## MATLAB and Python parity

At `(5,20,20)` and `ds=0.1`, all three double-pole models were compared at both three and four loops. The tests compare RG metadata, Borel grids, each OPE component, alpha selection, pole residues, spectra, f (including its NaN mask), and residuals. The alpha comparison requires `abs(log10(alpha_matlab/alpha_python)) <= 1e-8`.

The standalone MATLAB smoke test checks a nontrivial penalty, a nonzero continuum boundary and the identity `A=U*C*X'`, `B=V*S*X'`, plus direct first-order conditions. The final three double-pole cases use the default four-loop inputs. Direct-solution, scaled-gradient and factor-reconstruction checks each pass their 1e-10 tolerance. The toy system also has exactly two unpenalized pole modes. Extra direct-solution checks are enabled explicitly with `CheckSolution=true`.

To repeat the MATLAB comparison from the program directory:

```sh
RUN_MATLAB_TESTS=1 MATLAB_PARITY_LOOPS=4 python -B -W error -m unittest discover -s tests -p test_matlab_parity.py -v
RUN_MATLAB_TESTS=1 MATLAB_PARITY_LOOPS=3 python -B -W error -m unittest discover -s tests -p test_matlab_parity.py -v
```

Set `MATLAB_EXECUTABLE` if MATLAB is not on PATH and is not at the tested macOS location. The tests create their own temporary results and use the current Python implementation as the cross-language reference.
