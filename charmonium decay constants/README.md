# Charmonium decay constants

Three channel packages use common numerical infrastructure and a shared MSbar QCD running prescription. Inputs and outputs use GeV and powers of GeV.

| Package | Pole model names | Second pole | Residue conversion |
|---|---|---|---|
| `Jpsi_vector` | `2pole` | psi(2S) | `f = sqrt(R)` |
| `eta_c_axial_vector` | `2pole` | chi_c1 | `f = sqrt(R)` |
| `eta_c_pseudoscalar` | `2pole` | eta_c(2S) | `f = (2*mc(mu)*sqrt(R)+a_i)/m_i^2` when a_i is supplied |

Nonpositive pole residues produce `NaN` decay constants. Neither pole residues nor the reconstructed continuum are constrained to be positive.

## Python

Python 3.9+, NumPy and SciPy are required. Install the Python dependencies from `requirements.txt`. The GSVD backend additionally requires a library exporting the LP64 C interface `LAPACKE_dggsvd3`. It detects supported system, NumPy and SciPy OpenBLAS libraries; `JPSI_LAPACKE_LIBRARY` can specify the library explicitly. `ETA_C_LAPACKE_LIBRARY` remains a compatibility alias. There is no normal-equation or Gram-matrix fallback.

Run from this directory:

```sh
python -m Jpsi_vector 5 20 20 --ds 0.1 --mu 2 --rg-loops 4
python -m eta_c_axial_vector 5 20 20 --ds 0.1 --mu 2 --rg-loops 4
python -m eta_c_pseudoscalar 5 20 20 --ds 0.1 --mu 2 --rg-loops 4
```

The positional arguments are `M2_min`, `M2_max`, and `Lambda`, all in GeV². Every program computes its double-pole model, identified by `2pole`. These small windows demonstrate execution; they are not recommended physical extraction windows.

```python
from Jpsi_vector import solve
from inverse_common.rg import RGInputs

result = solve(5, 20, 20, ds=0.1,
               mu=2.0, rg_loops=4, rg_inputs=RGInputs())
channel = result["channels"]["2pole"]
print(result["metadata"]["rg"])
print(channel["f"], channel["alpha_star"])
```

The API returns data without writing files. CLI runs create unique directories containing `result.npz` and `summary.json`; `--output` selects the output root. NPZ values do not require pickle. JSON replaces nonfinite values with `null`.

## Pseudoscalar axial-anomaly conversion

The PP solver now includes the axial Ward conversion in the main result:
`f_i = f_raw_i + a_i/m_i^2`, with `f_raw_i = 2*mc(mu)*sqrt(R_i)/m_i^2`.
Here `m_i` is the hadron mass, not the Borel parameter. The default ground-state
matrix element is `a_1=-0.079 GeV^3` from Qin et al., Phys. Rev. D 97, 096002
(2018), Eq. (22). Its shift is approximately -8.872 MeV. No 2S matrix element
is assumed: by default the second entry remains density-only, as recorded in
`metadata.anomaly.applied` and `metadata.anomaly.f_definition`.

```python
from eta_c_pseudoscalar import solve
out = solve(5, 20, 20)  # a_1=-0.079; a_2 omitted
raw = solve(5, 20, 20, anomaly_1s=None)  # omit both anomaly inputs
custom = solve(5, 20, 20, anomaly_1s=-0.079, anomaly_2s=0.02)
```

The custom 2S value above only illustrates the API; it is not a physical input
recommended by this project. Supply matrix elements in GeV^3 with the same
state phase convention as `h_i=+sqrt(R_i)`. `None` means missing, while zero
explicitly supplies a zero matrix element.

```sh
python -m eta_c_pseudoscalar 5 20 20 --anomaly-1s none
python -m eta_c_pseudoscalar 5 20 20 --anomaly-1s -0.079 --anomaly-2s 0.02
```

PP output schema v3 uses `f`, `f_all`, and `f_eta_c` for values with supplied
anomaly terms included. `f_raw` and `f_raw_all` retain the density contribution;
`anomaly_shift_GeV` gives the additive shift per state. Summary JSON contains
both `f_GeV` and `f_raw_GeV`. Nonpositive residues remain undefined. The anomaly
conversion leaves the OPE, residues, spectrum, and alpha selection unchanged.
Plots reading `f` must not add the same correction again; old postprocessors
that perform the correction themselves should read `f_raw` instead.

## Running scale

Running is a separate program, `inverse_common/rg.py`. It can be used without running an inverse problem:

```sh
python -m inverse_common.rg --mu 2 --loops 4
python -m inverse_common.rg --mu 2 --loops 3
```

The output is JSON; omitting `--loops` selects four loops. The three solvers call this same module through `solve(..., rg_loops=4)` or `--rg-loops 4`.

Default inputs are `alpha_s^(5)(MZ)=0.1184`, `MZ=91.1876`, `mbbar(mbbar)=4.18`, `mcbar(mcbar)=1.275`, and `mu=2`. All three solvers use four-loop beta and mass-anomalous-dimension functions and three-loop bottom-threshold matching. The supported physical target domain is `mcbar <= mu <= mbbar`, with charm active (`nf=4`).

The same derived mass and coupling enter the OPE, perturbative boundary and pseudoscalar residue conversion. Changing `RGInputs.mbar_c` changes both the anchor mass and its scale. Calls do not mutate global physics parameters.

`rg_loops=1,2,3,4` selects consistent running orders with matching order `rg_loops-1`. The historical eta_c leading-log mass prescription is not equivalent to `rg_loops=1`; original versions remain in Git history.

See [RG conventions](docs/RG_CONVENTIONS.md) for formulas, matching direction and numerical references. Four-loop RG does not raise the OPE truncation order: VV retains its existing NLO factor, A-A retains LO plus condensates, and P-P retains its existing LO/NLO plus condensates.

## Numerical conventions

- The Borel step is 0.01 GeV²; requested endpoints are included exactly.
- The default spectral hat step is 0.1 GeV². The P-P historical production grid uses 0.01 GeV²; this is a different discretization with a substantially larger GSVD.
- The continuum starts at `S_H=13.9105153024` (VV) or `14.99` (A-A and P-P). These are channel model inputs.
- The objective is `||sqrt(W)*(K*z-g)||² + alpha*rho_hat'*H*rho_hat`, with `H` the mass-plus-stiffness matrix. Alpha is not squared. Pole residues are unpenalized.
- The two continuum endpoints are fixed at zero and the perturbative spectral density at Lambda. Endpoint cross terms are retained.
- VV and A-A use their adaptive alpha grid; P-P uses the historical 2001-point `10^[-20,0]` grid. Each uses the Raus Q/MQ/HR rule.
- The P-P full-alpha reconstruction uses the generalized filter `sigma/(sigma²+alpha*mu_gsvd²)`. Its selector uses `sigma/mu_gsvd`. This corrects the original standalone P-P filter convention.

Shared infrastructure lives in `inverse_common/`. The three package directories retain their own physical kernels and pole definitions. Keep the shared directory alongside the packages when copying or installing the code.

## MATLAB

The three MATLAB versions use the shared `matlab/+charmonium` package and MATLAB's built-in `gsvd`. No separate GSVD toolbox is required. From this directory:

```matlab
addpath("Jpsi_vector/matlab", "eta_c_axial_vector/matlab", ...
        "eta_c_pseudoscalar/matlab", "matlab");
rg = charmonium.runDec(2, "Loops", 4);
j = solve_jpsi_vector(5, 20, 20, "RGLoops", 4);
a = solve_eta_c_axial_vector(5, 20, 20, "RGLoops", 4);
p = solve_eta_c_pseudoscalar(5, 20, 20, "RGLoops", 4);
```

Omitting `Loops` or `RGLoops` selects four loops. Set it to 1, 2, or 3 to change the order. All three MATLAB programs return their double-pole result in `two_pole`, since field names cannot begin with a digit. For example, `p.two_pole.f` holds the two pseudoscalar decay constants. Common scale information is in `p.metadata.rg`.

See [MATLAB instructions](matlab/README.md) for options and validation. `CheckSolution=true` enables the extra augmented least-squares comparison; it is off for normal calculations.

## Validation

```sh
python -B -W error -m unittest discover -s tests -v
```

The tests cover RG coefficients and independent quadrature references, parameter propagation, nonzero boundary terms, stationarity, residue conventions, QOC derivatives and serialization. MATLAB parity is an explicit optional test and is documented with the MATLAB implementation.

See [validation record](docs/VALIDATION.md) for the tested environment, test windows and limits of this validation.
