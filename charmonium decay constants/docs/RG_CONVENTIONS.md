# Shared four-loop QCD RG conventions

## Source and scope

Source: K. G. Chetyrkin, J. H. Kuhn and M. Steinhauser, **RunDec: a
Mathematica package for running and decoupling of the strong coupling and
quark masses**, hep-ph/0004189v1 (19 April 2000). The provided PDF was read
with its extracted text; Eq. (7) and Eqs. (20)-(21) were also checked on
rendered PDF pages 6 and 11 (printed pages 5 and 10).

The code implements Eqs. (1), (2), (6), (7), (19), and (20)/(21) at the
MSbar heavy-quark self scale. All masses and scales are in GeV.
This upgrades **input running**; it does not change the perturbative order
of the OPE, condensate evolution, spectral models, or inverse solver.

## API for Python solver integration

```python
from inverse_common.rg import RGInputs, derive, evolve_fixed_nf

r = derive()  # mu=2 GeV, loops=4
mc, alpha_s = r["mc"], r["alpha_s"]
r3 = derive(loops=3)  # compare with frozen three-loop physics
r_varied = derive(inputs=RGInputs(mbar_c=1.305))
# Low-level mathematical fixed-nf evolution, no automatic thresholds:
alpha_end, mass_end = evolve_fixed_nf(alpha_s, mc, 2.0, 3.0, nf=4, loops=4)
```

- `RGInputs(mbar_c=1.275, alpha_mz=0.1184, mz=91.1876, mb=4.18)` is
  a frozen, hashable dataclass. Inputs are validated and converted to floats.
- `derive(mu=2.0, loops=4, inputs=None) -> dict`: fields include `mu`,
  `mc`, `alpha_s`, `nf=4`, `loops`, `matching_loops`, `inputs`,
  `input_source`, `coefficient_source`, `mass_method="ode"`, `mass_anchor`, `threshold`, and
  integration tolerances. Each call returns fresh mutable metadata.
- `beta_coefficients(nf)` and `gamma_coefficients(nf)` return immutable
  tuples with indices 0 through 3 in the normalization below.
- `decouple_alpha_down(alpha_s, nf=5, matching_loops=3)` returns the
  lower-flavour coupling **only at** `mu=mhbar(mhbar)`.
- `evolve_fixed_nf(alpha_s, mass, start, end, nf=4, loops=4)` returns
  `(alpha_s(end), mass(end))`. It integrates precisely the requested beta
  and gamma orders, in either direction, with fixed explicit nf.

`derive` supports `mbar_c <= mu <= mb`. At `mu=mb`, it returns the nf=4
**low** side of the threshold. Targets outside the interval raise an error.
Light-mass matching across bottom is not implemented. Charm always stays
active: no nf=4 to nf=3 decoupling at the charm anchor. The low-level
fixed-nf API is a mathematical operation and performs no threshold inference.

The input chain is:

1. Evolve `alpha_s^(5)(MZ)=0.1184` from `MZ=91.1876` to `mb=4.18`, nf=5.
2. Match the coupling nf=5 to nf=4 at `mu=mbbar(mbbar)=4.18`.
3. Evolve nf=4 coupling from mb to the **moving** anchor `mu=mbar_c`.
4. Impose `mc^(4)(mbar_c)=mbar_c`, then evolve coupling and mass to target mu.

Steps 1 and 3 use a dummy mass of 1 in the fixed-nf routine; that mass is
discarded. No physical charm mass crosses the bottom threshold. Changing
`mbar_c` changes both its boundary mass and boundary scale. Changing `mb`
changes the MSbar bottom self mass/threshold together; it does not implement
arbitrary matching-scale variation with fixed bottom mass.

## Independent running program (Python and MATLAB)

One independent RG program supports orders 1, 2, 3, and 4. Callers select
an order without duplicating RG coefficients in each inverse solver.
The default is four loops. From the project directory:

```sh
python -m inverse_common.rg --mu 2 --loops 4
python -m inverse_common.rg --mu 2 --loops 3
python inverse_common/rg.py --mu 2 --loops 4
python -m inverse_common.rg --mu 2 --loops 4 --mcbar 1.305 --alpha-mz 0.119 --mb 4.2 --mz 91.2
```

Both entry points print one standard JSON object to stdout, including
physical results, input sources, loop/matching orders, and
`mass_method="ode"`. They write no result files. Anchor defaults match
`RGInputs`; `--mu` defaults to 2 GeV and `--loops` defaults to 4.
Use `--help` to list options. Invalid arguments/domains report a clear
stderr error and exit status 2. Integration failures exit status 1.
Errors emit no partial JSON result or traceback.

MATLAB already exposes the independent running entry point:

```matlab
r4 = charmonium.runDec(2, "Loops", 4);
r3 = charmonium.runDec(2, "Loops", 3);
r_default = charmonium.runDec(2); % four loops by default
```

The inverse-program caller selects the running order; parameter variations
may additionally supply optional anchors.

System-Python validation passes all 19 RG/CLI unittest methods. CLI checks
run real subprocesses for both entry points, loop orders 1-4, defaults,
custom anchors, JSON output, an empty temporary output directory, help,
invalid inputs (status 2), and integration failure (status 1).

## Precise coefficients for MATLAB

Let `a=alpha_s/pi`, `t=ln(mu^2)`, and `n=nf`.

```text
da/dt       = -beta0*a^2 - beta1*a^3 - beta2*a^4 - beta3*a^5
d ln(mc)/dt = -gamma0*a - gamma1*a^2 - gamma2*a^3 - gamma3*a^4

beta0 = (11 - 2*n/3)/4
beta1 = (102 - 38*n/3)/16
beta2 = (2857/2 - 5033*n/18 + 325*n^2/54)/64
beta3 = (149753/6 + 3564*zeta3
         + (-1078361/162 - 6508*zeta3/27)*n
         + (50065/162 + 6472*zeta3/81)*n^2
         + 1093*n^3/729)/256

gamma0 = 1
gamma1 = (202/3 - 20*n/9)/16
gamma2 = (1249 + (-2216/27 - 160*zeta3/3)*n - 140*n^2/81)/64
gamma3 = (4603055/162 + 135680*zeta3/27 - 8800*zeta5
          + (-91723/27 - 34192*zeta3/9 + 880*zeta4 + 18400*zeta5/9)*n
          + (5242/243 + 800*zeta3/9 - 160*zeta4/3)*n^2
          + (-332/243 + 64*zeta3/27)*n^3)/256

zeta3 = 1.2020569031595942854
zeta4 = pi^4/90
zeta5 = 1.0369277551433699263
```

No extra factor of 2 is applied to these coefficients: the independent
variable is already `ln(mu^2)`. If integrating in `ln(mu)`, multiply both
right-hand sides by 2. Integrate `ln(mc)` to keep positive masses.
At N loops truncate both coefficient vectors to indices `0..N-1`.

| nf | beta0 | beta1 | beta2 | beta3 |
|---|---:|---:|---:|---:|
| 4 | 2.08333333333333333 | 3.20833333333333333 | 6.34924768518518519 | 31.38744695230514049 |
| 5 | 1.91666666666666667 | 2.41666666666666667 | 2.82667824074074074 | 18.85217315933944027 |

| nf | gamma0 | gamma1 | gamma2 | gamma3 |
|---|---:|---:|---:|---:|
| 4 | 1 | 3.65277777777777778 | 9.94704026107295732 | 27.30293363421854576 |
| 5 | 1 | 3.51388888888888889 | 7.41986321214366579 | 11.03432124206805957 |

## Exact matching convention

RunDec Eq. (19) is `alpha_low = zeta_g^2 * alpha_high`. At the
MSbar self mass the logarithm `L=ln(mu^2/mhbar(mu)^2)=0`. Eqs. (20) and
(21) agree there:

```text
a_high = alpha_high/pi
nl = nf_high - 1
d1 = 0
d2 = 11/72
d3 = 564731/124416 - (82043/27648)*zeta3 - (2633/31104)*nl
alpha_low = alpha_high * (1 + d2*a_high^2 + d3*a_high^3)
```

At bottom `nl=4`, `d3=0.633450719548517290326579933996`. Use the
**high-side** coupling in the matching polynomial. This is not the pole-mass
matching formula Eq. (22). N-loop running uses N-1-loop matching:

| loops | matching_loops | terms in zeta_g squared at L=0 |
|---|---|---|
| 1 | 0 | 1 |
| 2 | 1 | 1 (d1=0) |
| 3 | 2 | 1 + d2*a_high^2 |
| 4 | 3 | 1 + d2*a_high^2 + d3*a_high^3 |

## Numerical results and validation

Default inputs, `mu=2 GeV`; these are RG outputs, not fitted decay constants:

| loops | matching_loops | alpha_s^(4)(2) | mc^(4)(2) GeV |
|---|---|---:|---:|
| 1 | 0 | 0.26976783425397743 | 1.171916344348932 |
| 2 | 1 | 0.2992873428686697 | 1.115375821681162 |
| 3 | 2 | 0.30294509881312925 | 1.0990536546335794 |
| 4 | 3 | 0.3043984044027081 | 1.0932768471848942 |

Production integration is SciPy `solve_ivp(method="DOP853")`,
`rtol=2e-12`, `atol=2e-14`. It integrates the finite-order beta/gamma
ODEs directly, rather than using the truncated mass-ratio expansion of
RunDec Eq. (10). An implementation based on that expansion can differ by
formally higher-order terms and is not the same numerical convention.
The result explicitly records `mass_method="ode"`. The independent physics
audit supplied on 2026-10-14 confirms the ODE reference with 55-digit
quadrature. For the same four-loop inputs, the truncated `c(a)` ratio gives
`mc(2)=1.092612428731451 GeV`, lower by `0.664418453333 MeV`. That value
is not the ODE regression target.

An independent 50-decimal mpmath calculation, using exact rational
coefficients and `mp.zeta`, solves instead in coupling space. With
`B(a)=sum(beta_i*a^(i+2))` and `G(a)=sum(gamma_i*a^(i+1))`, solve

```text
integral[a_start -> a_end] da/B(a) + 2*ln(mu_end/mu_start) = 0
m_end = m_start * exp(integral[a_start -> a_end] G(a)/B(a) da)
```

Reference results are:

```text
alpha_s^(5)(4.18)   = 0.226141845816242853648784456086645162252737161
alpha_s^(4)(4.18)   = 0.226374296730276340625973984089232024321228921
alpha_s^(4)(1.275)  = 0.391619277160514744523610869382628899904314963
alpha_s^(4)(2)      = 0.304398404402813591292861525129163528655322391
mc^(4)(2)          = 1.09327684718478403978880450988362917080470269 GeV
```

The DOP853/reference difference is below 2e-13 for each listed quantity.
These references come from independent quadrature, not an installed RunDec
binary. Installing the temporary `rundec==0.7` validation package failed
because its distribution lacked `CRunDec.3.1.h`; no RunDec-library cross-check
is claimed. Independent exact coefficients and numerical method provide
the four-loop reference.

Tests cover the frozen three-loop values (absolute tolerance 2e-11),
independent four-loop values, nf=4/5 coefficients, analytic one-loop
running, 1-4-loop round trips and composition, threshold correction/order,
the nf=4 side of mb, moving charm anchors 1.245/1.275/1.305, actual numerical
changes when loops changes, immutable input/cache isolation, invalid input
orders, unsupported target domains, and a nonperturbative integration guard.

```text
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 -m unittest discover -s tests -p test_rg.py -v
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 -m unittest discover -s tests -v
```

Validation environment: system `/usr/bin/python3` (Python 3.9.6),
NumPy 2.0.2, SciPy 1.13.1, and standard-library `unittest`; independent
reference arithmetic used mpmath 1.4.1. NumPy/SciPy are runtime dependencies.
The test suite requires no pytest or mpmath installation; its numeric
references are frozen and its live quadrature uses the runtime SciPy.
No dependency is installed into the project by this change.

## Caching and failure behavior

The internal bounded LRU caches store immutable tuples: 128 boundary
configurations `(RGInputs, loops)` and 512 fixed-nf evolutions. Cache keys
include every numerical input and loop order. There is no mutable cached
result dict. Nonfinite/nonpositive scales, masses and coupling are rejected;
bool and noninteger flavour/loop counts are rejected. Couplings `alpha_s>=pi`
are rejected as a permissive perturbative guard; this is not an accuracy
guarantee near that bound. Integration stops and raises when it reaches
`a=1`. Solver failure and invalid final states raise explicit errors.

## MATLAB RG implementation and measured validation

The separately authorized MATLAB RG repair uses these package files:
`matlab/+charmonium/runDec.m`, `qcdCoefficients.m`, `integrateRG.m`, and
`matchAlphaDownMS.m`. The incorrect `matchAlphaDownOS.m` was removed.

```matlab
r = charmonium.runDec(2, 'Loops', 4, 'AlphaMZ', .1184, ...
    'MZ', 91.1876, 'Mb', 4.18, 'McBar', 1.275);
[alpha_end, mass_end] = charmonium.integrateRG(.3, 1.1, 2, 3, 4, 4);
c = charmonium.qcdCoefficients(4, 4);
alpha_low = charmonium.matchAlphaDownMS(alpha_high, 5, 3);
```

`runDec` returns the same physical field names and conventions as Python,
including `mass_method='ode'`, `inputs`, `mass_anchor`, and `threshold`.
The target domain is `[McBar,Mb]`. There is no independent `MassScale`
option. `integrateRG` accepts/returns alpha_s, converts the internal state
to `[alpha_s/pi; ln(mass)]`, and integrates in `2*ln(mu)` with adaptive
`ode113`, `RelTol=2e-12`, `AbsTol=2e-14`.

Measured with MATLAB R2026a on 2026-10-14:

| loops | alpha_s^(4)(2) | mc^(4)(2) GeV |
|---|---:|---:|
| 3 | 0.30294509881300485 | 1.099053654633755 |
| 4 | 0.30439840440251487 | 1.0932768471851821 |

Four-loop bottom threshold: high-side `0.22614184581613855`, low-side
`0.22637429673017168`. Each four-loop reference discrepancy is below
4e-13. Temporary MATLAB validation passed the three/four-loop numeric
references at absolute tolerance 2e-11, nf=4/5 gamma3, 1-4-loop round
trips/composition, moving charm anchors, low-side mb, zero-log matching
orders, invalid scales/orders, and rejection of `MassScale`.

System-Python `unittest` discovered and passed all 13 RG test methods
(the previous parameterized cases are retained as subtests). An earlier
full-suite discovery found two unrelated failures: ndarray JSON
serialization in `test_contract` and a syntax error importing
`test_matlab_parity`. Those files belong to another implementation agent;
this stage does not claim that full-suite run passed. Git state and
production datasets were not operated on.
