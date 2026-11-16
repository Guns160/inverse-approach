"""Affine total-spectrum endpoints with a penalty on the remaining continuum."""
from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, qr, solve_triangular

from inverse_common.gsvd import compact_gsvd
from inverse_common.numerics import _dot


def assemble(poles, hats, bw_column, g, weights, H, endpoint_bw, boundary):
    """Eliminate c endpoints using c=a+E v, with v=(R_BW,c_interior)."""
    n_data, n_s = hats.shape
    if n_data < n_s + 1:
        raise ValueError(f"BW model needs at least {n_s + 1} Borel points")
    a = np.zeros(n_s)
    a[[0, -1]] = boundary
    d = np.zeros(n_s)
    d[[0, -1]] = -endpoint_bw
    E = np.zeros((n_s, n_s - 1))
    E[:, 0] = d
    E[1:-1, 1:] = np.eye(n_s - 2)
    L_H = cholesky(H, lower=False, check_finite=True)
    L_E = _dot(L_H, E)
    Q, T = qr(L_E, mode="economic", check_finite=True)
    if np.any(np.diag(T) == 0) or not np.all(np.isfinite(T)):
        raise RuntimeError("rank(L_H E) must equal n_hats-1")
    L_a = _dot(L_H, a)
    projection = _dot(Q.T, L_a)
    q = solve_triangular(T, projection, lower=False)
    orthogonal = L_a - _dot(Q, projection)
    F = np.column_stack((bw_column + _dot(hats, d), hats[:, 1:-1]))
    sqrt_w = np.sqrt(weights)
    A = sqrt_w[:, None] * np.column_stack((poles, F))
    # R_BW changes both continuum endpoints and is therefore penalized through E.
    B = np.column_stack((np.zeros((n_s - 1, 2)), T))
    b0 = sqrt_w * (g - _dot(hats, a) + _dot(F, q))
    return dict(A=A, B=B, b0=b0, a=a, d=d, E=E, q=q, T=T,
                kappa=float(_dot(orthogonal, orthogonal)),
                qr_relative_error=float(np.linalg.norm(L_E - _dot(Q, T)) /
                                        np.linalg.norm(L_E)))


def factors_one_window(system, g, weights):
    """Use the repository's LP64 compact GSVD, with exactly two zero columns."""
    A, B, b0 = system["A"], system["B"], system["b0"]
    decomp = compact_gsvd(A, B, method="lapacke", qr_reduce=True)
    n_reg = B.shape[0]
    if decomp.p != n_reg or decomp.k != 2:
        raise RuntimeError(f"GSVD ranks p={decomp.p} k={decomp.k}, expected p={n_reg} k=2")
    if np.any(decomp.mu <= np.finfo(float).tiny):
        raise RuntimeError("GSVD mu has non-positive compact modes")
    beta = _dot(decomp.U.T, b0)
    remainder = b0 - _dot(decomp.U, beta)
    target_a = decomp.U * np.r_[decomp.sigma, decomp.sigma_unregularized]
    target_b = np.column_stack((decomp.V * decomp.mu, np.zeros((n_reg, 2))))
    error_a = float(np.linalg.norm(_dot(A, decomp.X) - target_a) / np.linalg.norm(target_a))
    error_b = float(np.linalg.norm(_dot(B, decomp.X) - target_b) / np.linalg.norm(target_b))
    if max(error_a, error_b) > 1e-8:
        raise RuntimeError("BW GSVD reconstruction error exceeds 1e-8")
    return dict(sigma=decomp.sigma, mu=decomp.mu,
                sigma_unreg=decomp.sigma_unregularized,
                sigma_tilde=decomp.sigma / decomp.mu,
                beta_c=beta[:n_reg], beta_n=beta[n_reg:],
                null=float(_dot(remainder, remainder)), X=decomp.X,
                q=system["q"], T=system["T"], kappa=system["kappa"],
                data_norm=float(np.sqrt(np.sum(weights * g*g))), n_data=np.int64(len(g)),
                n_regularized=n_reg, n_unregularized=2,
                AX_relative_error=error_a, BX_relative_error=error_b,
                qr_relative_error=system["qr_relative_error"],
                source=decomp.source, library=decomp.library)


def reconstruct_alpha_family(factors, system, alpha, *, block_size=128):
    """Apply sigma/(sigma^2+alpha*mu^2); keep signed residues and continuum."""
    sigma, mu = factors["sigma"], factors["mu"]
    n_reg = len(sigma)
    X = factors["X"]
    offset = _dot(X[:, n_reg:], factors["beta_n"] / factors["sigma_unreg"])
    residues = np.empty((3, len(alpha)))
    cont = np.empty((len(system["a"]), len(alpha)))
    for start in range(0, len(alpha), block_size):
        stop = min(start + block_size, len(alpha))
        den = sigma[:, None]**2 + alpha[None, start:stop] * mu[:, None]**2
        coef = (sigma * factors["beta_c"])[:, None] / den
        shifted = _dot(X[:, :n_reg], coef) + offset[:, None]
        v = shifted[2:] - system["q"][:, None]
        residues[:2, start:stop] = shifted[:2]
        residues[2, start:stop] = v[0]
        cont[:, start:stop] = system["a"][:, None] + _dot(system["E"], v)
    if not all(np.all(np.isfinite(v)) for v in (residues, cont)):
        raise FloatingPointError("BW reconstruction returned NaN or Inf")
    return residues, cont


def diagnose_family(poles, hats, bw_column, g, weights, H, system, residues, cont, alpha):
    """Evaluate the original model and its three independent gradient blocks."""
    prediction = _dot(poles, residues[:2]) + bw_column[:, None]*residues[2] + _dot(hats, cont)
    residual = prediction - g[:, None]
    weighted = weights[:, None] * residual
    Hc = _dot(H, cont)
    rss = np.sum(weights[:, None] * residual**2, axis=0)
    penalty_squared = np.sum(cont * Hc, axis=0)
    vc = _dot(hats.T, weighted) + Hc*alpha
    delta_gradient = _dot(poles.T, weighted)
    bw_gradient = _dot(bw_column, weighted) + _dot(system["d"], vc)
    vc_scale = (np.linalg.norm(_dot(hats.T, weights[:, None]*prediction), axis=0) +
                np.linalg.norm(_dot(hats.T, weights*g)) + alpha*np.linalg.norm(Hc, axis=0))
    delta_scale = (np.linalg.norm(_dot(poles.T, weights[:, None]*prediction), axis=0) +
                   np.linalg.norm(_dot(poles.T, weights*g)))
    bw_scale = (np.abs(_dot(bw_column, weights[:, None]*prediction)) +
                abs(_dot(bw_column, weights*g)) + np.linalg.norm(system["d"])*vc_scale)
    tiny = np.finfo(float).tiny
    stationarity = np.vstack((
        np.linalg.norm(delta_gradient, axis=0) / np.maximum(delta_scale, tiny),
        np.linalg.norm(vc[1:-1], axis=0) / np.maximum(vc_scale, tiny),
        np.abs(bw_gradient) / np.maximum(bw_scale, tiny)))
    if not all(np.all(np.isfinite(v)) for v in (residual, penalty_squared, stationarity)):
        raise FloatingPointError("BW diagnostics returned NaN or Inf")
    return dict(backsub_all=prediction, residual_all=residual,
                weighted_residual_all=np.sqrt(rss),
                relative_weighted_residual_all=np.sqrt(rss) / max(np.sqrt(np.sum(weights*g*g)), tiny),
                penalty_all=np.sqrt(np.maximum(penalty_squared, 0.0)),
                objective_all=rss + alpha*penalty_squared,
                stationarity_relative_all=stationarity)
