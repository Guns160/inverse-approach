"""Standalone eta_c pseudoscalar-current (P-P) inversion at mu = 2 GeV.

Same factorisation, endpoint elimination and GSVD path as the A-A package
(``charmonium/ηc/solver.py``): only the physics input and the alpha rule
differ.  The frozen P-P OPE (``physics_mu2.py``) and the accepted Raus p.53
selector (``raus2024_pp.py``) are byte-for-byte copies of the published
赝标流 program; they are not rederived here.
"""
from __future__ import annotations

import numpy as np

from . import geometry as G
from . import physics_mu2 as W
from . import raus2024_pp as R
from .factors import factors_one_window
from inverse_common.gsvd import lapacke_info as _common_lapacke_info
from inverse_common.numerics import borel_grid as _common_borel_grid
from inverse_common.reconstruction import reconstruct_alpha_family
from inverse_common.validation import positive_scalar, validate_window

# Accepted P-P pole content (select_mu2.py / thin_full_alpha.py).
MASSES = {"1pole": (2.9841,), "2pole": (2.9841, 3.6378)}
SCHEMA = "eta_c-pp-standalone-v1"
OPE_LABEL = "PP n=0 LO+NLO(alpha_s)+G2+G3"
RULE = "raus2024_pp p.53 heuristic"
# Spectral hat step: DS_DEFAULT matches the A-A package; DS_PRODUCTION is the
# step of the accepted 赝标流 mu=2 chain (physics_mu2.sgrid_pp).
DS_DEFAULT = 0.1
DS_PRODUCTION = 0.01


def lapacke_info():
    """Probe installed libraries, without installing anything."""
    return _common_lapacke_info()


def _positive(value, name):
    return positive_scalar(value, name)


def borel_grid(lo, hi):
    """Include exact requested endpoints; retain historical anchored arithmetic.

    On 1.5 + .01*k use that expression for interior points. Else start at lo.
    A last partial step is allowed. A near-endpoint sample is replaced, not doubled.
    """
    return _common_borel_grid(lo, hi, step=G.M2STEP)


def _pp_f(residues, masses):
    """Accepted P-P converter f_i = 2 mc sqrt(R_i) / m_i^2 (GeV).

    Non-positive residues give NaN, never an absolute value or a clip.
    """
    residues = np.asarray(residues, float)
    positive = residues > 0
    safe = np.where(positive, residues, 1.0)
    return np.where(positive, 2.0 * W.MC * np.sqrt(safe) / masses[:, None] ** 2, np.nan)


def solve(M2_min, M2_max, Lambda, *, channel=None, ds=DS_DEFAULT):
    """Return {metadata, M2, weights, OPE, channels}; channel=None means both.

    All dimensionful input values are in GeV^2.  No files are written.
    channel may be '1pole' or '2pole'.  Each channel has independent alpha,
    selected on the accepted logspace(-20, 0, 2001) grid.
    ds is the spectral hat step; pass DS_PRODUCTION = 0.01 to reproduce the
    accepted production chain, whose s grid is physics_mu2.sgrid_pp.
    """
    lo, hi, lam, ds, channels, t, s_hat = validate_window(
        M2_min, M2_max, Lambda, masses=MASSES, channel=channel,
        s_h=W.S_H, ds=ds, dm2=G.M2STEP)
    n_hats = len(s_hat)
    info = lapacke_info()
    if not info["available"]:
        raise RuntimeError("LAPACKE_dggsvd3 unavailable. Point JPSI_LAPACKE_LIBRARY to an existing LP64 LAPACKE library; ETA_C_LAPACKE_LIBRARY is accepted as a legacy alias; no numerical fallback is used.")
    s = s_hat
    g = W.pp_ope(lam, t)
    hat_kernel = G.hats(t, s_hat)
    H = G.h_mass_stiffness(s_hat)
    boundary = np.array([0., float(W.pp_rho_pert(np.array([lam]))[0])])
    weights = G.trapw(t)
    if not all(np.all(np.isfinite(a)) for a in (g, hat_kernel, H)):
        raise RuntimeError("Nonfinite physical operator/OPE for this window")
    alpha = np.asarray(W.ALPHAS, float)
    results = {}
    for c in channels:
        masses = np.array(MASSES[c]); p = len(masses)
        poles = np.column_stack([np.exp(-m * m / t) / t for m in masses])
        K = np.column_stack((poles, hat_kernel))
        fac = factors_one_window(K, g, t, H, p, boundary, "lapacke")
        sig = fac["sigma"]; n = len(sig)
        fac["alpha"] = alpha
        residues, cont, qoc = reconstruct_alpha_family(
            fac, boundary, n_poles=p, s_grid=s)
        # Accepted selector: sigma/beta are flat arrays with one contiguous
        # block per window (here a single window), null_sq one entry per window.
        sel = R.select_row(sig, fac["beta_c"], np.array([0, n], dtype=np.int64),
                           np.array([fac["null"]], float), alpha, qoc[None, :])
        if sel["status"][0] != "defined":
            raise RuntimeError(f"{c}: alpha selection {sel['status'][0]}")
        k = int(sel["index"][0])
        a = float(alpha[k])
        backsub = G._dot(K, np.r_[residues[:, k], cont[:, k]])
        residual = backsub - g
        wr = float(np.sqrt(np.sum(weights * residual * residual)))
        f_all = _pp_f(residues, masses)
        d_md, psi = R.quantities(sig, fac["beta_c"], fac["null"], alpha)
        if not all(np.all(np.isfinite(v)) for v in (residues, cont, residual)):
            raise RuntimeError(f"{c}: nonfinite reconstructed solution")
        pen = float(np.sqrt(max(0., float(G._dot(cont[:, k], G._dot(H, cont[:, k]))))))
        results[c] = dict(
            masses=masses, s=s, s_hat=s_hat, boundary=boundary,
            alpha=alpha, residues_all=residues, f_all=f_all, rho_all=cont, qoc=qoc,
            residues=residues[:, k], f=f_all[:, k], f_eta_c=float(f_all[0, k]),
            rho=cont[:, k], hat=cont[:, k], alpha_star=a, selection_index=k,
            alpha_diagnostics=dict(
                index=k, index_q_local=int(sel["index_q_local"][0]),
                index_mq=int(sel["index_mq"][0]), index_hr=int(sel["index_hr"][0]),
                source=str(sel["source"][0]), fallback=str(sel["fallback"][0]),
                status=str(sel["status"][0]), n_local_minima=int(sel["n_local_minima"][0]),
                rule_b_index=int(sel["rule_b_index"][0]),
                qoc_argmin_index=int(sel["qoc_argmin_index"][0]),
                on_boundary=bool(k in (0, alpha.size - 1)),
                grid=alpha, d_MD=d_md, psi_Q=psi,
                psi_HR=d_md / np.sqrt(alpha), sigma_min_sq=float(np.min(sig) ** 2)),
            backsub=backsub, residual=residual,
            weighted_residual=wr,
            relative_weighted_residual=wr / max(fac["data_norm"], np.finfo(float).tiny),
            unweighted_residual=float(np.linalg.norm(residual)),
            penalty=pen,
            objective=wr**2+a*pen**2,
            min_rho=float(np.min(cont[:, k])),
            negative_rho_count=int(np.count_nonzero(cont[:, k] < 0)),
            factors=fac)
    return dict(metadata=dict(schema=SCHEMA, M2_min=lo, M2_max=hi,
                Lambda=lam, n_hats=n_hats, ds=float(ds), dM2=G.M2STEP,
                mu=2.0, mc=W.MC, alpha_s=W.ALPHA_S,
                ope_order=OPE_LABEL, G2=W.G2, G3=W.G3, S_H=W.S_H,
                rule=RULE, lapacke=info),
                M2=t, weights=weights, OPE=dict(total=g), channels=results)
