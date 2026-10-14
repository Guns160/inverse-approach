"""Standalone J/psi n=0 spectral inversion; no cached physics or workspace imports."""
from __future__ import annotations

import numpy as np
from . import physics as P
from .alpha_grid import selection_grid
from .factors import factors_one_window
from . import raus2024 as R
from inverse_common.gsvd import lapacke_info as _common_lapacke_info
from inverse_common.numerics import borel_grid as _common_borel_grid
from inverse_common.rg import derive
from inverse_common.validation import positive_scalar, validate_window


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
    return _common_borel_grid(lo, hi, step=P.M2STEP)


def solve(M2_min, M2_max, Lambda, *, channel=None, ds=P.DS,
          mu=2.0, rg_loops=4, rg_inputs=None):
    """Return {metadata, M2, weights, OPE, channels}; channel=None means both.

    All dimensionful input values are in GeV^2. No files are written.
    channel may be '1delta' or '2delta'. Each channel has independent alpha.
    mu and rg_inputs are passed to the shared RG derivation on every call;
    no frozen charm mass or coupling is rebound as a default.
    """
    lo, hi, lam, ds, channels, t, s_hat = validate_window(
        M2_min, M2_max, Lambda, masses=P.MASSES, channel=channel,
        s_h=P.S_H, ds=ds, dm2=P.M2STEP)
    n_hats = len(s_hat)
    rg_state = derive(mu=mu, loops=rg_loops, inputs=rg_inputs)
    mc = float(rg_state["mc"])
    alpha_s = float(rg_state["alpha_s"])
    info = lapacke_info()
    if not info["available"]:
        raise RuntimeError("LAPACKE_dggsvd3 unavailable. Point JPSI_LAPACKE_LIBRARY to an existing LP64 LAPACKE library; no numerical fallback is used.")
    s = s_hat
    pert, d4, d6 = P.aa_ope(
        lam, t, mc=mc, alpha_s=alpha_s, components=True)
    g = pert + d4 + d6
    hat_kernel = P.hats(t, s_hat)
    H = P.h_mass_stiffness(s_hat)
    boundary = np.array(
        [0., float(P.aa_rho_pert(np.array([lam]), mc=mc, alpha_s=alpha_s)[0])])
    weights = P.trapw(t)
    if not all(np.all(np.isfinite(a)) for a in (g, hat_kernel, H)):
        raise RuntimeError("Nonfinite physical operator/OPE for this window")
    results = {}
    for c in channels:
        masses = np.array(P.MASSES[c]); p = len(masses)
        poles = np.column_stack([np.exp(-m*m/t)/t for m in masses])
        K = np.column_stack((poles, hat_kernel))
        fac = factors_one_window(K, g, t, H, p, boundary, "lapacke")
        alpha, covered = selection_grid(fac["sigma_tilde"], fac["beta_c"], fac["null"])
        sel = R.select(fac["sigma_tilde"], fac["beta_c"], fac["null"], alpha)
        if sel["status"] != "defined":
            raise RuntimeError(f"{c}: alpha selection {sel['status']}")
        a = sel["alpha"]; sig = fac["sigma"]; gsvd_mu = fac["mu"]; n = len(sig)
        coef = sig*fac["beta_c"]/(sig*sig+a*gsvd_mu*gsvd_mu)
        un = fac["beta_n"]/fac["sigma_unreg"]
        shifted = P._dot(fac["X"][:, :n], coef) + P._dot(fac["X"][:, n:], un)
        residues = shifted[:p]
        cont = np.r_[boundary[0], shifted[p:]-fac["q"], boundary[1]]
        rho = np.interp(s, s_hat, cont)
        backsub = P._dot(K, np.r_[residues, cont])
        residual = backsub-g
        wr = float(np.sqrt(np.sum(weights*residual*residual)))
        pen = float(np.sqrt(max(0., float(P._dot(cont, P._dot(H, cont))))))
        f = np.full(p, np.nan)
        positive = residues > 0
        f[positive] = np.sqrt(residues[positive])
        d_md, psi = R.quantities(fac["sigma_tilde"], fac["beta_c"], fac["null"], alpha)
        d_factor = float(np.sqrt(np.sum((a/(fac["sigma_tilde"]**2+a)*fac["beta_c"])**2)+fac["null"]))
        if not all(np.all(np.isfinite(v)) for v in (residues, cont, residual, rho)):
            raise RuntimeError(f"{c}: nonfinite reconstructed solution")
        results[c] = dict(
            masses=masses, s=s, rho=rho, s_hat=s_hat, hat=cont,
            alpha=alpha, selection_index=int(sel["index"]),
            residues=residues, f=f, f_jpsi=float(f[0]), boundary=boundary,
            alpha_star=a, alpha_diagnostics=dict(
                **sel, grid=alpha, d_MD=d_md, psi_Q=psi,
                psi_MQ=R.psi_mq_tikhonov(psi, d_md), psi_HR=d_md/np.sqrt(alpha),
                on_boundary=sel["index"] in (0,len(alpha)-1),
                small_alpha_peak_covered=covered, sigma_min_sq=float(np.min(fac["sigma_tilde"])**2)),
            backsub=backsub, backsub_poles=P._dot(poles, residues),
            backsub_continuum=P._dot(hat_kernel, cont), residual=residual,
            weighted_residual=wr, relative_weighted_residual=wr/max(fac["data_norm"], np.finfo(float).tiny),
            unweighted_residual=float(np.linalg.norm(residual)), penalty=pen,
            objective=wr**2+a*pen**2, residual_factor=d_factor,
            residual_factor_discrepancy=abs(d_factor-wr), min_rho=float(np.min(rho)),
            negative_rho_count=int(np.count_nonzero(rho<0)), factors=fac)
    return dict(metadata=dict(schema="jpsi-standalone-v2", M2_min=lo, M2_max=hi,
                Lambda=lam, n_hats=n_hats, ds=ds, dM2=P.M2STEP,
                mu=float(rg_state["mu"]), mc=mc, alpha_s=alpha_s,
                rg=rg_state, G2=P.G2, G3=P.G3, S_H=P.S_H,
                ope_order="VV n=0 LO*(1+alpha_s/pi)+G2+G3",
                rule="raus2024_etac_equiv_p0_1", lapacke=info),
                M2=t, weights=weights, OPE=dict(pert=pert, D4=d4, D6=d6, total=g), channels=results)
