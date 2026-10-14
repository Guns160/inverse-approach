"""Independent eta_c pseudoscalar-current (P-P) spectral sum-rule solver, mu = 2 GeV."""
from .solver import solve, lapacke_info

__all__ = ["solve", "lapacke_info"]
