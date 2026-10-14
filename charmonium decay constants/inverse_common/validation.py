"""Shared input validation for a one-window spectral inversion."""
from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np

from .numerics import borel_grid, sgrid_hats


def positive_scalar(value, name):
    if isinstance(value, (bool, str, bytes)) or np.ndim(value) != 0:
        raise ValueError(f"{name} must be a finite positive scalar")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite positive scalar") from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive scalar")
    return value


def validate_window(M2_min, M2_max, Lambda, *, masses, channel, s_h, ds=0.1, dm2=0.01):
    """Validate solve arguments and construct the Borel and spectral grids."""
    lo = positive_scalar(M2_min, "M2_min")
    hi = positive_scalar(M2_max, "M2_max")
    lam = positive_scalar(Lambda, "Lambda")
    ds = positive_scalar(ds, "ds")
    dm2 = positive_scalar(dm2, "dM2")
    if hi <= lo:
        raise ValueError("M2_max must exceed M2_min")
    if channel is not None and (not isinstance(channel, str) or channel not in masses):
        choices = ", ".join(repr(name) for name in masses)
        raise ValueError(f"channel must be None or one of {choices}")
    channels = tuple(masses) if channel is None else (channel,)
    if lam <= s_h or any(lam <= m * m for c in channels for m in masses[c]):
        raise ValueError("Lambda must exceed S_H and all included pole masses squared")
    s_hat = sgrid_hats(lam, step=ds, s_h=s_h)
    t = borel_grid(lo, hi, step=dm2)
    for c in channels:
        free = len(s_hat) - 2 + len(masses[c])
        if len(t) < free:
            raise ValueError(f"{c} needs at least {free} Borel points (free unknowns); got {len(t)}")
    return lo, hi, lam, ds, channels, t, s_hat
