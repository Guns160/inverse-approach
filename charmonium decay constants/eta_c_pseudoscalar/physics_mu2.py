"""eta_c pseudoscalar-current OPE with explicit RG-derived inputs.

This file retains the accepted P-P formulas and integration grids, but no
frozen charm mass or coupling is bound at a default parameter or module global.
Every call must supply the values derived at the requested scale.
"""
from __future__ import annotations

import numpy as np

from inverse_common.numerics import _dot, trapw

G2 = 0.038
G3 = 0.013
S_H = 14.99
ALPHAS = np.logspace(-20, 0, 2001)


def pp_li2_series(z):
    z = np.asarray(z, dtype=float)
    term = z.copy()
    total = term.copy()
    for k in range(2, 10001):
        term *= z
        add = term / (k * k)
        updated = total + add
        if np.all(np.abs(add) <= 4 * np.spacing(np.maximum(1.0, np.abs(updated)))):
            return updated
        total = updated
    raise RuntimeError("Li2 series did not converge")


def pp_li2(x):
    values = np.asarray(x, dtype=float)
    flat = values.reshape(-1)
    out = np.empty_like(flat)
    az = flat == 0
    ao = flat == 1
    low = (flat > 0) & (flat <= 0.5)
    high = (flat > 0.5) & (flat < 1)
    out[az] = 0
    out[ao] = np.pi**2 / 6
    out[low] = pp_li2_series(flat[low])
    z = 1 - flat[high]
    out[high] = (
        np.pi**2 / 6
        - np.log(flat[high]) * np.log1p(-flat[high])
        - pp_li2_series(z)
    )
    return out.reshape(values.shape)


def pp_finite_function(v):
    lp = np.log1p(v) - np.log1p(-v)
    lh = np.log1p(v) - np.log(2.0)
    bracket = (
        np.pi**2 / 6 + lp * lh
        + 2 * pp_li2((1 - v) / (1 + v))
        - 4 * pp_li2(v)
        + pp_li2(v * v)
        + 2 * pp_li2((1 + v) / 2)
        - 2 * pp_li2((1 - v) / 2)
    )
    return (
        (1 + v * v) * bracket
        + (19 / 16 - 3 * v + v * v / 8 + 3 * v**4 / 16) * lp
        + 29 * v / 8 - 3 * v * v / 8 + 6 * v * lh - 4 * v * np.log(v)
    )


def pp_rho_pert_components(s, *, mc: float, alpha_s: float):
    """Return (LO, NLO) continuum components at the supplied charm inputs."""
    values = np.asarray(s, dtype=float)
    lo = np.zeros_like(values)
    nlo = np.zeros_like(values)
    threshold = 4 * mc * mc
    above = values > threshold
    at = values == threshold
    if np.any(above):
        sa = values[above]
        v = np.sqrt(1 - threshold / sa)
        lo[above] = 3 * sa * v / (8 * np.pi**2)
        nlo[above] = sa * alpha_s * pp_finite_function(v) / (2 * np.pi**3)
    nlo[at] = threshold * alpha_s / (4 * np.pi)
    return lo, nlo


def pp_rho_pert(s, *, mc: float, alpha_s: float):
    lo, nlo = pp_rho_pert_components(s, mc=mc, alpha_s=alpha_s)
    return lo + nlo


def pp_condensates(m2, *, mc: float, x_nodes=4001):
    """Return the dimension-four and dimension-six condensate contributions."""
    x = np.linspace(1e-8, 1 - 1e-8, x_nodes)
    xb = 1 - x
    q = x * xb
    wx = trapw(x)
    cubic = mc * mc * (x**3 + xb**3)
    a = mc * x + mc * xb
    x3 = x**3 + xb**3
    b4 = 64 - (416 / 3) * q - 45 * x3 / q
    b6 = (
        (45 * a * a / q - (2 / 3) * mc * mc) * x3
        - 6 * mc * mc * (x**4 + xb**4) / q
        - 12 * mc * mc * x3
        - (478 / 9) * mc * mc * q
        + (287 / 9) * a * a
        - (320 / 9) * q * mc * mc
    )
    b8 = (
        (6 * a * a / q + 28 * mc * mc) * mc * mc * (x**4 + xb**4)
        - 4 * mc**4 * (x**4 + xb**4)
        + (16 / 5) * mc**4 * (x**5 + xb**5) / q
        + (30 * a * a + 88 * mc * mc * q) * mc * mc * (x * x + xb * xb)
        - (112 / 3) * a**4
        + (64 / 3) * mc**4 * q
        + 128 * mc**4 * q * q
    )
    g2 = np.empty_like(m2)
    g3 = np.empty_like(m2)
    for start in range(0, len(m2), 128):
        mm = m2[start:start + 128, None]
        expo = np.exp(-mc * mc / (q * mm))
        term4_g2 = (-6 * mc * mc - 2 * cubic / q) / (mm**2 * x**2 * xb**2)
        term6_g2 = a * a * cubic / (mm**3 * x**4 * xb**4)
        i2 = expo * (term4_g2 + term6_g2)
        i3 = expo * (
            b4 / (mm**2 * x**2 * xb**2)
            + b6 / (2 * mm**3 * x**3 * xb**3)
            + b8 / (6 * mm**4 * x**4 * xb**4)
            - (2 / 15) * a * a * mc**4 * (x**5 + xb**5)
              / (mm**5 * x**6 * xb**6)
        )
        g2[start:start + len(mm)] = G2 * _dot(i2, wx) / (48 * np.pi)
        g3[start:start + len(mm)] = G3 * _dot(i3, wx) / (32 * (4 * np.pi)**2)
    return g2, g3


def pp_ope(upper, m2, *, mc: float, alpha_s: float,
           nv=16001, nx=4001, components=False):
    """Return P-P LO, NLO, D4, and D6 with explicit charm inputs."""
    m2 = np.asarray(m2, dtype=float)
    threshold = 4 * mc * mc
    velocity = np.linspace(0, np.sqrt(1 - threshold / upper), nv)
    s = threshold / (1 - velocity * velocity)
    ws = trapw(s)
    rho_lo, rho_nlo = pp_rho_pert_components(s, mc=mc, alpha_s=alpha_s)
    lo = np.empty_like(m2)
    nlo = np.empty_like(m2)
    for start in range(0, len(m2), 256):
        mm = m2[start:start + 256]
        weight = np.exp(-s[None, :] / mm[:, None]) / mm[:, None]
        lo[start:start + len(mm)] = _dot(weight * rho_lo, ws)
        nlo[start:start + len(mm)] = _dot(weight * rho_nlo, ws)
    g2, g3 = pp_condensates(m2, mc=mc, x_nodes=nx)
    if components:
        return lo, nlo, g2, g3
    return lo + nlo + g2 + g3


def sgrid_pp(upper):
    """Return the accepted production P-P grid with step 0.01 GeV^2."""
    n = int(np.floor((upper - S_H) / 0.01 + 1e-10))
    s = S_H + 0.01 * np.arange(n + 1)
    if upper - s[-1] > 1e-10:
        s = np.append(s, upper)
    else:
        s[-1] = upper
    if np.min(np.diff(s)) < 0.009999:
        raise AssertionError(
            f"P-P s grid contains a spurious small cell at Lambda={upper}"
        )
    return s


__all__ = [
    "pp_li2", "pp_li2_series", "pp_finite_function",
    "pp_rho_pert_components", "pp_rho_pert", "pp_condensates", "pp_ope",
    "sgrid_pp",
]
