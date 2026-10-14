"""Numerical and geometric helpers shared by all three channels."""
from __future__ import annotations

import numpy as np


def _dot(a, b):
    """Use BLAS dot and reject a non-finite result."""
    out = np.dot(a, b)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("matrix product returned NaN or Inf")
    return out


def trapw(x):
    x = np.asarray(x, float)
    if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)) or not np.all(np.diff(x) > 0):
        raise ValueError("trapezoid input must be a finite strictly increasing vector")
    d = np.diff(x)
    w = np.empty_like(x)
    w[0] = d[0] / 2
    w[-1] = d[-1] / 2
    w[1:-1] = (d[:-1] + d[1:]) / 2
    return w


def sgrid_hats(L, nodes=None, *, step=0.1, s_h):
    """Construct the physical-s hat grid with an exact Lambda endpoint."""
    if not np.isfinite(L) or L <= s_h:
        raise ValueError("Lambda must be finite and exceed the spectral threshold")
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
        return np.linspace(s_h, L, nodes)
    n = int(np.floor((L - s_h) / step))
    s = s_h + step * np.arange(n + 1)
    s = s[s < L]
    tol = 8 * np.finfo(float).eps * max(1.0, abs(s_h), abs(L))
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
    """Analytic hat-basis kernel on the physical-s grid."""
    m2 = np.asarray(m2, float)
    s = np.asarray(s, float)
    if m2.ndim != 1 or s.ndim != 1:
        raise ValueError("hats requires one-dimensional M2 and s grids")
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


def h_mass_stiffness(s):
    """Return H = M + M1 on a physical-s hat grid."""
    s = np.asarray(s, float)
    if s.ndim != 1 or s.size < 3 or not np.all(np.diff(s) > 0):
        raise ValueError("H1 assembly requires a strictly increasing hat grid")
    h = np.diff(s)
    n = len(s)
    H = np.zeros((n, n))
    for i, hi in enumerate(h):
        H[i, i] += hi / 3 + 1 / hi
        H[i + 1, i + 1] += hi / 3 + 1 / hi
        H[i, i + 1] += hi / 6 - 1 / hi
        H[i + 1, i] += hi / 6 - 1 / hi
    return H


def borel_grid(lo, hi, *, step=0.01):
    """Retain the historical anchored Borel grid and exact user endpoints."""
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        raise ValueError("Borel bounds must be finite and increasing")
    if isinstance(step, (bool, np.bool_)) or np.ndim(step) != 0:
        raise ValueError("Borel step must be a finite positive scalar")
    try:
        step = float(step)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("Borel step must be a finite positive scalar") from exc
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Borel step must be a finite positive scalar")
    index = round((lo - 1.5) / step)
    anchor = 1.5 + step * index
    tol = 8 * np.finfo(float).eps * max(1.0, abs(lo), abs(hi))
    n = int(np.floor((hi - lo) / step))
    if abs(anchor - lo) <= tol:
        grid = 1.5 + step * (index + np.arange(n + 1))
    else:
        grid = lo + step * np.arange(n + 1)
    grid[0] = lo
    grid = grid[grid < hi]
    if len(grid) > 1 and hi - grid[-1] <= tol:
        grid[-1] = hi
    else:
        grid = np.r_[grid, hi]
    if not np.all(np.diff(grid) > 0):
        raise ValueError("Borel grid cannot be represented with increasing float64 values")
    return grid
