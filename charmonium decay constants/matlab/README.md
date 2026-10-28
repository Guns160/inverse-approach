# MATLAB implementation

The shared numerical core is the `+charmonium` package in this directory. It implements the physical-s hat basis, H1 penalty, channel-specific OPE formulas, economy-QR compact GSVD, reconstruction, and Raus 2024 selector. MATLAB's built-in `gsvd` is required; no toolbox or external GSVD package is needed.

Add one public entry directory to the path before calling its solver. Each entry also adds this shared `matlab` directory:

```matlab
addpath("Jpsi_vector/matlab");
out = solve_jpsi_vector(5,20,20,"RGLoops",3);
```

The independent RG program can be called without running an inverse problem:

```matlab
addpath("matlab");
rg = charmonium.runDec(2,"Loops",3);   % alpha_s and mc at mu=2 GeV
```

The three public solvers default to four loops and select the order through `RGLoops`:

```matlab
addpath("Jpsi_vector/matlab");
addpath("eta_c_axial_vector/matlab");
addpath("eta_c_pseudoscalar/matlab");
j = solve_jpsi_vector(5,20,20,"RGLoops",3);
a = solve_eta_c_axial_vector(5,20,20,"RGLoops",3);
p = solve_eta_c_pseudoscalar(5,20,20,"RGLoops",3);
```

All three solvers run only `2pole` and return the result in `out.two_pole`. `Channel="2pole"` is the optional explicit spelling of this default. Metadata, Borel points, weights, and OPE components remain in `out.metadata`, `out.M2`, `out.weights`, and `out.OPE`. The two pole masses stay channel-specific:

| Solver | Pole masses (GeV) |
| --- | --- |
| `solve_jpsi_vector` | 3.096900, 3.686097 |
| `solve_eta_c_axial_vector` | 2.9841, 3.51067 |
| `solve_eta_c_pseudoscalar` | 2.9841, 3.6378 |

Common options are `Channel`, `DS`, `RGLoops`, `Mu`, `RGInputs`, `Alpha`, and `CheckSolution`. `RGLoops` accepts 1 through 4 and defaults to 4; the main solvers call the independent `charmonium.runDec` program with that order. Other defaults are `DS=0.1`, `Mu=2`, and `CheckSolution=false`. Setting `CheckSolution=true` adds an independent augmented least-squares regression and direct gradient check. `RGInputs` is the complete coupled struct with fields `mbar_c`, `alpha_mz`, `mz`, and `mb`; the charm mass and its anchor scale both use `mbar_c`.

## Pseudoscalar axial anomaly

Only the PP entry accepts `Anomaly1S` and `Anomaly2S` (GeV^3). Defaults are
`Anomaly1S=-0.079` and `Anomaly2S=[]`. An empty input leaves that state
density-only; zero explicitly supplies a zero matrix element.

```matlab
p = solve_eta_c_pseudoscalar(5,20,20); % ground-state anomaly included
raw = solve_eta_c_pseudoscalar(5,20,20,"Anomaly1S",[]);
```

`two_pole.f`, `f_all` and `f_eta_c` include supplied `a_i/m_i^2` terms.
`f_raw` and `f_raw_all` preserve `2*mc*sqrt(R_i)/m_i^2`, and
`anomaly_shift_GeV` records the shifts. `metadata.anomaly` records the inputs
and which states remain density-only. Do not add the correction again when
plotting `f`. The OPE and spectral inversion do not depend on anomaly inputs.

Run `test_pp_anomaly` after adding this directory to verify the Ward identity,
optional inputs, and unchanged spectral inversion. The Python/MATLAB parity
check also runs this test.

Run the self-contained numerical smoke test with:

```matlab
addpath("matlab");
result = test_smoke;
```

It checks a generalized filter with two unpenalized poles and a nonzero endpoint against a direct augmented least-squares solution. It then checks all three physical entries with default four-loop inputs, independent stationarity, GSVD factor identities, and the unique `2pole` result contract. The Python/MATLAB cross-language test in `tests/test_matlab_parity.py` is gated by `RUN_MATLAB_TESTS=1` because it launches MATLAB. Set `MATLAB_PARITY_LOOPS=3` or `MATLAB_PARITY_LOOPS=4` to choose the validation order.
