"""Reconstruction helpers for compact-GSVD factors."""
from __future__ import annotations

import numpy as np

from .numerics import _dot, trapw


def reconstruct_alpha_family(factors, boundary, *, n_poles, s_grid,
                             block_size=128):
    """Reconstruct the alpha family in blocks.

    The filter uses the generalized GSVD form with the original X columns:
    sigma/(sigma**2 + alpha*mu**2).  Derivatives with respect to log10(alpha)
    include alpha*ln(10) and the required mu**2 factor.
    """
    alpha = np.asarray(factors["alpha"], float)
    sigma = np.asarray(factors["sigma"], float)
    mu = np.asarray(factors["mu"], float)
    beta = np.asarray(factors["beta_c"], float)
    x_reg = factors["X"][:, : len(sigma)]
    x_unreg = factors["X"][:, len(sigma) :]
    unreg = factors["beta_n"] / factors["sigma_unreg"]
    pole_offset = _dot(x_unreg, unreg)
    n_unknown = factors["X"].shape[0]
    n_s = n_unknown - n_poles + 2
    residues = np.empty((n_poles, alpha.size), float)
    rho = np.empty((n_s, alpha.size), float)
    rho[0, :] = boundary[0]
    rho[-1, :] = boundary[1]
    qoc = np.empty(alpha.size, float)
    qoc_weights = trapw(s_grid)[1:-1]
    for start in range(0, alpha.size, block_size):
        stop = min(start + block_size, alpha.size)
        aa = alpha[start:stop]
        den = sigma[:, None] ** 2 + aa[None, :] * mu[:, None] ** 2
        coef = (sigma * beta)[:, None] / den
        shifted = _dot(x_reg, coef)
        shifted += pole_offset[:, None]
        residues[:, start:stop] = shifted[:n_poles]
        rho[1:-1, start:stop] = shifted[n_poles:] - factors["q"][:, None]
        dcoef = -coef * mu[:, None] ** 2 / den
        dshifted = _dot(x_reg, dcoef)[n_poles:]
        dshifted *= (np.log(10.0) * aa)[None, :]
        qoc[start:stop] = np.sum(qoc_weights[:, None] * dshifted**2, axis=0)
    if not all(np.all(np.isfinite(v)) for v in (residues, rho, qoc)):
        raise FloatingPointError("full-alpha reconstruction returned NaN or Inf")
    return residues, rho, qoc
