"""eta_c axial-vector n=0 OPE with explicitly supplied RG-derived inputs."""
from __future__ import annotations

import numpy as np

from inverse_common.numerics import (
    _dot, borel_grid, h_mass_stiffness, hats, sgrid_hats, trapw,
)

G2 = 0.038
G3 = 0.013
S_H = 14.99
DS = 0.1
M2STEP = 0.01
MASSES = {
    "1delta": (2.9841,),
    "2delta": (2.9841, 3.51067),
}


def aa_rho_pert(s, *, mc: float, alpha_s: float):
    """Return the LO A-A continuum; alpha_s is accepted and intentionally unused."""
    _ = float(alpha_s)
    s = np.asarray(s, float)
    out = np.zeros_like(s)
    q = 4 * mc * mc
    m = s > q
    v = np.sqrt(1 - q / s[m])
    out[m] = 3 * v * (1 - v * v / 3) / (8 * np.pi**2)
    return out


def aa_ope(L, m2, *, mc: float, alpha_s: float,
           nv=16001, nx=4001, components=False):
    """Return LO, D4, and D6 with explicit charm inputs."""
    m2 = np.asarray(m2, float)
    q = 4 * mc * mc
    vel = np.linspace(0, np.sqrt(max(1.0 - q / L, 0.0)), nv)
    denom = np.maximum(1.0 - vel * vel, q / L)
    sg = q / denom
    sg[-1] = L
    ws = trapw(sg)
    rw = aa_rho_pert(sg, mc=mc, alpha_s=alpha_s)
    pert = np.empty_like(m2)
    for start in range(0, len(m2), 256):
        mm = m2[start:start + 256]
        kernel = np.exp(-sg[None, :] / mm[:, None]) / mm[:, None] * rw
        pert[start:start + len(mm)] = _dot(kernel, ws)

    x = np.linspace(1e-8, 1 - 1e-8, nx)
    xb = 1 - x
    wx = trapw(x)
    a3 = -(45 / 8) * x * xb * (x**3 + xb**3) - (69 / 72) * x**2 * xb**2
    b3 = (
        -(3 / 4) * mc**2 * x * xb * (x**4 + xb**4)
        + mc**2 * x**2 * xb**2 * (
            -(23 / 12) * (x**2 + xb**2) + x * xb / 3 + 2
        )
    )
    c3 = (2 / 5) * mc**4 * x * xb * (x**5 + xb**5)
    g2 = np.empty_like(m2)
    g3 = np.empty_like(m2)
    for start in range(0, len(m2), 128):
        mm = m2[start:start + 128, None]
        E = np.exp(-mc**2 / (x * xb * mm))
        i2 = E * (
            1 / (2 * mm**2)
            - mc**2 * (x**3 + xb**3) / (2 * mm**3 * x**2 * xb**2)
        )
        i3 = E * (
            a3 / (2 * mm**3 * x**3 * xb**3)
            + b3 / (6 * mm**4 * x**4 * xb**4)
            + c3 / (24 * mm**5 * x**5 * xb**5)
        )
        g2[start:start + len(mm)] = G2 * _dot(i2, wx) / (6 * np.pi)
        g3[start:start + len(mm)] = G3 * _dot(i3, wx) / (4 * np.pi) ** 2
    if components:
        return pert, g2, g3
    return pert + g2 + g3


def _channel_grid(L, nodes=None, *, step=DS):
    return sgrid_hats(L, nodes, step=step, s_h=S_H)


def sgrid_hats(L, nodes=None, *, step=DS):
    return _channel_grid(L, nodes, step=step)


def sgrid_fine(L, *, step=DS):
    return _channel_grid(L, step=step)


__all__ = [
    "aa_ope", "aa_rho_pert", "borel_grid", "h_mass_stiffness", "hats",
    "sgrid_fine", "sgrid_hats", "trapw",
]
