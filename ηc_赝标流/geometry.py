"""Shared geometry kernels: Borel trapezoid weights, the frozen physical-s hat
basis, and H = M + M1.

Copied unchanged from the A-A package's ``physics.py`` minus its A-A OPE.  This
package carries no OPE of its own: the P-P mu = 2 GeV physics is the accepted
``physics_mu2.py``, byte-identical to the 赝标流 production source.
"""
from __future__ import annotations

import numpy as np

S_H = 14.99
DS = 0.1
M2STEP = 0.01


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


def sgrid_hats(L: float, nodes: int | None = None, *, step: float = DS) -> np.ndarray:
    """Anchor at S_H, step by `step`, retain exact Lambda and a final short interval.

    Near-coincident endpoints are merged at float64 roundoff, not at step scale.
    The default step DS is this package's default spectral step; pass a step to
    reproduce another discretisation (0.01 for the accepted production chain).
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
    """Analytic hat-basis kernel on the physical-s grid."""
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
