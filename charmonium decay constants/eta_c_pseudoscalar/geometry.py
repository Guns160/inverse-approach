"""Compatibility wrappers around the shared physical-s geometry kernels."""
from inverse_common.numerics import (  # noqa: F401
    _dot, borel_grid, h_mass_stiffness, hats, trapw,
)
from . import physics_mu2 as _physics

S_H = _physics.S_H
DS = 0.1
M2STEP = 0.01


def sgrid_hats(L, nodes=None, *, step=DS):
    from inverse_common.numerics import sgrid_hats as _sgrid_hats
    return _sgrid_hats(L, nodes, step=step, s_h=S_H)
