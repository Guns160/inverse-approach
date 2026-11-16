"""P-P inversion with two delta poles and one running-width BW state.

The original P-P OPE, RG, grids and compact-GSVD backend are reused.
This core returns the requested alpha family without selecting a platform.
"""
from __future__ import annotations

import numpy as np

from eta_c_pseudoscalar import physics_mu2 as W
from eta_c_pseudoscalar.solver import (
    DEFAULT_ANOMALY_1S, DS_DEFAULT, DS_PRODUCTION, MASSES, OPE_LABEL,
    _anomaly_value, _pp_f, lapacke_info,
)
from inverse_common.numerics import hats, h_mass_stiffness, trapw
from inverse_common.rg import derive
from inverse_common.validation import validate_window

from .factors import assemble, diagnose_family, factors_one_window, reconstruct_alpha_family
from .physics import RunningWidth

SCHEMA = "eta_c-pp-bw-standalone-v1"


def _alpha_grid(alpha):
    if alpha is None:
        return np.asarray(W.ALPHAS, float).copy()
    value = np.asarray(alpha)
    if value.dtype.kind not in "iuf" or value.ndim > 1 or value.size == 0:
        raise ValueError("alpha must be a positive scalar or a strictly increasing real vector")
    value = np.atleast_1d(value).astype(float)
    if not np.all(np.isfinite(value)) or np.any(value <= 0) or np.any(np.diff(value) <= 0):
        raise ValueError("alpha must be positive, finite and strictly increasing")
    return value


def solve(M2_min, M2_max, Lambda, *, channel=None, ds=DS_DEFAULT,
          alpha=None, mu=2.0, rg_loops=4, rg_inputs=None,
          anomaly_1s=DEFAULT_ANOMALY_1S, anomaly_2s=None,
          radius=4.0, dstar_d1_l=0, quadrature_nodes=1024):
    """Return the original P-P result structure extended by one BW state.

    channel chooses 3S, 4S, 5S (default), or the unavailable 6S input slot.
    Only one BW state is added per solve. 5S is primary, 4S is alternative,
    and the nearby 3S state is retained for diagnostic comparisons.
    alpha=None returns the original PP logspace(-20,0,2001) family, without
    automatic alpha selection. A scalar alpha also exposes single-solution
    fields for the shared CLI. No files are written by this function.

    s, Lambda, ds and M2 use GeV^2; mu uses GeV; radius uses GeV^-1.
    c is the remaining continuum in GeV^2. rho=c/s is dimensionless.
    The BW term has its own raw, unnormalized line shape in GeV^-2 and a
    fitted residue in GeV^4. Endpoints constrain c+R_BW*BW, not c alone.
    """
    channel = "5S" if channel is None else channel
    model = RunningWidth(channel, radius=radius, dstar_d1_l=dstar_d1_l,
                         quadrature_nodes=quadrature_nodes)
    alpha = _alpha_grid(alpha)
    anomaly_values = (_anomaly_value(anomaly_1s, "anomaly_1s"),
                      _anomaly_value(anomaly_2s, "anomaly_2s"), None)
    lo, hi, lam, ds, _, t, s = validate_window(
        M2_min, M2_max, Lambda, masses={channel: MASSES["2pole"]},
        channel=channel, s_h=W.S_H, ds=ds, dm2=0.01)
    if len(t) < len(s) + 1:
        raise ValueError(f"BW model needs at least {len(s) + 1} Borel points; got {len(t)}")
    rg = derive(mu=mu, loops=rg_loops, inputs=rg_inputs)
    mc, alpha_s = float(rg["mc"]), float(rg["alpha_s"])
    info = lapacke_info()
    if not info["available"]:
        raise RuntimeError("LAPACKE_dggsvd3 unavailable. Point JPSI_LAPACKE_LIBRARY to an existing LP64 LAPACKE library; no numerical fallback is used.")
    lo_ope, nlo_ope, d4, d6 = W.pp_ope(lam, t, mc=mc, alpha_s=alpha_s, components=True)
    pert = lo_ope + nlo_ope
    g = pert + d4 + d6
    G, H = hats(t, s), h_mass_stiffness(s)
    weights = trapw(t)
    boundary = np.array([0.0, float(W.pp_rho_pert(np.array([lam]), mc=mc, alpha_s=alpha_s)[0])])
    masses = np.array((*MASSES["2pole"], model.mass))
    poles = np.column_stack([np.exp(-m*m/t)/t for m in masses[:2]])
    k_B, integration = model.borel_column(t, W.S_H, lam)
    bw_shape = model.bw(s)
    if not all(np.all(np.isfinite(v)) for v in (g, G, H, boundary, k_B, bw_shape)):
        raise RuntimeError("Nonfinite physical operator/OPE for this window")
    system = assemble(poles, G, k_B, g, weights, H, bw_shape[[0, -1]], boundary)
    fac = factors_one_window(system, g, weights)
    residues, cont = reconstruct_alpha_family(fac, system, alpha)
    diagnostics = diagnose_family(poles, G, k_B, g, weights, H, system,
                                  residues, cont, alpha)
    f_raw = _pp_f(residues, masses, mc=mc)
    shifts = np.array([0.0 if a is None else a for a in anomaly_values]) / masses**2
    f = f_raw + shifts[:, None]
    bw_hat = bw_shape[:, None] * residues[2]
    total_hat = cont + bw_hat
    endpoint_errors = total_hat[[0, -1]] - boundary[:, None]
    positive = np.all(residues > 0, axis=0)
    endpoint_pass = np.max(np.abs(endpoint_errors), axis=0) <= 1e-10*(1 + abs(boundary[-1]))
    stationarity_pass = np.max(diagnostics["stationarity_relative_all"], axis=0) < 1e-8
    result = dict(
        masses=masses, s=s, s_hat=s, alpha=alpha, residues_all=residues,
        f_all=f, f_raw_all=f_raw, anomaly_shift_GeV=shifts,
        c_all=cont, hat_all=cont, rho_all=cont/s[:, None],
        bw_shape=bw_shape, gamma=model.gamma(s), bw_hat_all=bw_hat,
        total_hat_all=total_hat, rho_total_all=total_hat/s[:, None],
        k_B=k_B, boundary_total=boundary, positive_all=positive,
        endpoint_error_all=endpoint_errors,
        validation_pass_all=endpoint_pass & stationarity_pass,
        negative_c_count_all=np.count_nonzero(cont < 0, axis=0),
        factors=fac, **diagnostics)
    if len(alpha) == 1:
        result.update(
            alpha_star=float(alpha[0]), selection_index=0, residues=residues[:, 0],
            f=f[:, 0], f_raw=f_raw[:, 0], f_eta_c=float(f[0, 0]),
            c=cont[:, 0], hat=cont[:, 0], rho=cont[:, 0]/s,
            total_hat=total_hat[:, 0], rho_total=total_hat[:, 0]/s,
            backsub=diagnostics["backsub_all"][:, 0], residual=diagnostics["residual_all"][:, 0],
            weighted_residual=float(diagnostics["weighted_residual_all"][0]),
            relative_weighted_residual=float(diagnostics["relative_weighted_residual_all"][0]),
            penalty=float(diagnostics["penalty_all"][0]), objective=float(diagnostics["objective_all"][0]),
            alpha_diagnostics=dict(source="user_supplied", fallback="none", on_boundary=False))
    anomaly = dict(eta_c_1s_GeV3=anomaly_values[0], eta_c_2s_GeV3=anomaly_values[1],
                   bw_GeV3=None, applied=[a is not None for a in anomaly_values],
                   f_definition=["axial_with_supplied_anomaly" if a is not None
                                 else "pseudoscalar_density_only" for a in anomaly_values],
                   convention="h_i=+sqrt(R_i); supplied a_i uses the same state phase")
    metadata = dict(schema=SCHEMA, M2_min=lo, M2_max=hi, Lambda=lam, n_hats=len(s),
                    ds=ds, dM2=0.01, mu=float(rg["mu"]), mc=mc, alpha_s=alpha_s, rg=rg,
                    ope_order=OPE_LABEL, G2=W.G2, G3=W.G3, S_H=W.S_H,
                    channel=channel, alpha_selection="none", lapacke=info,
                    anomaly=anomaly, bw=model.metadata(), integration=integration,
                    residue_unit="GeV^4", continuum_hat_unit="GeV^2", rho_unit="dimensionless",
                    penalty="c^T H c; remaining continuum only",
                    endpoints="c(S_H)+R_BW*BW(S_H)=0; c(Lambda)+R_BW*BW(Lambda)=rho_pert(Lambda)")
    return dict(metadata=metadata, M2=t, weights=weights,
                OPE=dict(LO=lo_ope, NLO=nlo_ope, pert=pert, D4=d4, D6=d6, total=g),
                channels={channel: result})
