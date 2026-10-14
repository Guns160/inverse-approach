# Inverse approach

Research code for extracting charmonium pole residues and decay constants from Borel sum rules.

The programs are in [`charmonium decay constants/`](charmonium%20decay%20constants/):

| Program | Current | Explicit pole models |
|---|---|---|
| `Jpsi_vector` | Vector | J/psi; J/psi + psi(2S) |
| `eta_c_axial_vector` | Axial vector | eta_c; eta_c + chi_c1 |
| `eta_c_pseudoscalar` | Pseudoscalar density | eta_c; eta_c + eta_c(2S) |

The second axial-vector pole is chi_c1, not eta_c(2S). The programs preserve the separate OPE and residue conventions of these three channels.

Numerical outputs depend on the running-coupling prescription, spectral grid, Borel window, pole model and regularization rule. A successful solve by itself does not establish a physical stability window.

## Version history

- `ed20389`: original Python source upload.
- `06fb9a2`: particle/current directory layout and corrected entry-point imports.
- Later commits separate the common-core refactor, running-scale update and MATLAB implementation.

Historical results and parameter scans are not bundled with the source packages.
