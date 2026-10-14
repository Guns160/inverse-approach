"""Frozen eta_c n=0 OPE and hat kernels for the compact-GSVD run.

A-A n=0 LO + D4 + D6, no alpha_s term. Fixed mu=2, leading-log mass.
Physical coefficients from src/vendor/eta_aa_physics.py (2026-09-11 snapshot).
The inversion uses 0.1 GeV^2 hat spacing on [s_h, Lambda]; arbitrary Borel windows.
"""
from __future__ import annotations

import numpy as np

MC = 1.1334710909
MU = 2.0
G2 = 0.038
G3 = 0.013
S_H = 14.99
DS = 0.1
M2STEP = 0.01
MASSES = {
    "1delta": (2.9841,),
    "2delta": (2.9841, 3.51067),
}


def _dot(a, b):
    """1D/2D products without Accelerate matmul's spurious FP warnings.

    dot does not report real overflow either, so reject nonfinite results.
    """
    out = np.dot(a, b)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("matrix product returned NaN or Inf")
    return out


def trapw(x: np.ndarray) -> np.ndarray:
    d = np.diff(x)
    w = np.empty_like(x)
    w[0] = d[0] / 2
    w[-1] = d[-1] / 2
    w[1:-1] = (d[:-1] + d[1:]) / 2
    return w


def aa_rho_pert(s):
    s = np.asarray(s, float)
    out = np.zeros_like(s)
    q = 4 * MC * MC
    m = s > q
    v = np.sqrt(1 - q / s[m])
    out[m] = 3 * v * (1 - v * v / 3) / (8 * np.pi**2)
    return out


def aa_ope(L, m2, nv=16001, nx=4001, components=False):
    m2 = np.asarray(m2, float)
    q = 4 * MC * MC
    vel = np.linspace(0, np.sqrt(max(1.0 - q / L, 0.0)), nv)
    denom = np.maximum(1.0 - vel * vel, q / L)
    sg = q / denom
    sg[-1] = L
    ws = trapw(sg)
    rw = aa_rho_pert(sg)
    pert = np.empty_like(m2)
    for a in range(0, len(m2), 256):
        mm = m2[a : a + 256]
        pert[a : a + len(mm)] = _dot(np.exp(-sg[None, :] / mm[:, None]) / mm[:, None] * rw, ws)
    x = np.linspace(1e-8, 1 - 1e-8, nx)
    xb = 1 - x
    wx = trapw(x)
    a3 = -(45 / 8) * x * xb * (x**3 + xb**3) - (69 / 72) * x**2 * xb**2
    b3 = (
        -(3 / 4) * MC**2 * x * xb * (x**4 + xb**4)
        + MC**2 * x**2 * xb**2 * (-(23 / 12) * (x**2 + xb**2) + x * xb / 3 + 2)
    )
    c3 = (2 / 5) * MC**4 * x * xb * (x**5 + xb**5)
    g2 = np.empty_like(m2)
    g3 = np.empty_like(m2)
    for a in range(0, len(m2), 128):
        mm = m2[a : a + 128, None]
        E = np.exp(-MC**2 / (x * xb * mm))
        i2 = E * (1 / (2 * mm**2) - MC**2 * (x**3 + xb**3) / (2 * mm**3 * x**2 * xb**2))
        i3 = E * (
            a3 / (2 * mm**3 * x**3 * xb**3)
            + b3 / (6 * mm**4 * x**4 * xb**4)
            + c3 / (24 * mm**5 * x**5 * xb**5)
        )
        g2[a : a + len(mm)] = G2 * _dot(i2, wx) / (6 * np.pi)
        g3[a : a + len(mm)] = G3 * _dot(i3, wx) / (4 * np.pi) ** 2
    if components:
        return pert, g2, g3
    return pert + g2 + g3


def sgrid_fine(L: float, *, step: float = DS) -> np.ndarray:
    """Output grid, identical to the requested inversion grid (no densification)."""
    return sgrid_hats(L, step=step)


def sgrid_hats(L: float, nodes: int | None = None, *, step: float = DS) -> np.ndarray:
    """Anchor at S_H, step by ``step``, retain exact Lambda and a final short interval.

    Explicit ``nodes`` selects the historical linspace for regression only.
    The default spectral step is 0.1 GeV^2; callers may pass another positive
    step when reproducing a separately frozen discretisation.
    """
    if not np.isfinite(L) or L <= S_H:
        raise ValueError("Lambda must be finite and exceed S_H")
    if isinstance(step, (bool, np.bool_)) or np.ndim(step) != 0:
        raise ValueError("step must be a finite positive scalar")
    try:
        step = float(step)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("step must be a finite positive scalar") from exc
    if not np.isfinite(step) or step <= 0:
        raise ValueError("step must be a finite positive scalar")
    if nodes is not None:
        if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)) or nodes < 3:
            raise ValueError("need at least 3 hat nodes (integer nodes)")
        return np.linspace(S_H, L, nodes)
    n = int(np.floor((L - S_H) / step))
    s = S_H + step * np.arange(n + 1)
    s = s[s < L]
    tol = 8 * np.finfo(float).eps * max(1., abs(S_H), abs(L))
    if len(s) > 1 and L - s[-1] <= tol:
        s[-1] = L
    else:
        s = np.r_[s, L]
    if len(s) < 3:
        raise ValueError("need at least 3 hat nodes on the physical-s grid; increase Lambda")
    if not np.all(np.diff(s) > 0):
        raise ValueError("hat grid must have strictly increasing float64 values")
    return s


def hats(m2, s):
    K = np.zeros((len(m2), len(s)))
    for j, h in enumerate(np.diff(s)):
        q = h / m2
        st = np.exp(-s[j] / m2)
        small = np.abs(q) < 1e-4
        asc = np.empty_like(q)
        des = np.empty_like(q)
        z = q[small]
        asc[small] = 0.5 - z / 3 + z * z / 8 - z**3 / 30 + z**4 / 144
        des[small] = 0.5 - z / 6 + z * z / 24 - z**3 / 120 + z**4 / 720
        z = q[~small]
        e = np.exp(-z)
        asc[~small] = (1 - (1 + z) * e) / z**2
        des[~small] = (z - 1 + e) / z**2
        fac = h / m2
        K[:, j] += st * fac * des
        K[:, j + 1] += st * fac * asc
    return K


def h_mass_stiffness(s: np.ndarray) -> np.ndarray:
    """H = M + M1 on the physical-s hat grid."""
    h = np.diff(s)
    n = len(s)
    H = np.zeros((n, n))
    for i, hi in enumerate(h):
        H[i, i] += hi / 3 + 1 / hi
        H[i + 1, i + 1] += hi / 3 + 1 / hi
        H[i, i + 1] += hi / 6 - 1 / hi
        H[i + 1, i] += hi / 6 - 1 / hi
    return H

