"""Raus MMA2024 abstract p.53 heuristic parameter choice, applied to the
P-P one-pole/two-pole standard-form SVD factors.

Selection logic is a verbatim port of the eta_c-standard implementation used
by the J/psi platform viewer (`raus2024/src/raus2024_etac.py`), so the two
selections are directly comparable: p0 = 1, psi_Q without the leading r
factor, psi_HR = d_MD / sqrt(alpha), r = 1/alpha parametrisation,
r_H = min(r_Q, max(r_MQ, r_HR)) -> alpha_hat = 1/r_H.
"""
from __future__ import annotations

import numpy as np

TINY = np.finfo(np.float64).tiny
P0_TIKHONOV = 1.0


def quantities(sigma, beta, null, alpha):
    """d_MD and psi_Q on the alpha grid.

    Same values as the reference implementation, via the algebraically equal
    two-sum form (fewer full-size temporaries):
      d^2  = alpha^3 * sum(beta^2/(sigma^2+alpha)^3) + null
      psi^2= alpha^2 * sum(sigma^2 beta^2/(sigma^2+alpha)^4)
    """
    s = np.asarray(sigma, float)
    b = np.asarray(beta, float)
    a = np.asarray(alpha, float)
    u = s * s
    b2 = b * b
    b2u = b2 * u
    s3 = np.empty(a.size)
    s4 = np.empty(a.size)
    blk = 64
    for k0 in range(0, a.size, blk):
        k1 = min(k0 + blk, a.size)
        r = 1.0 / (u[None, :] + a[k0:k1, None])
        r2 = r * r
        s3[k0:k1] = (b2[None, :] * (r2 * r)).sum(axis=1)
        s4[k0:k1] = (b2u[None, :] * (r2 * r2)).sum(axis=1)
    d = np.sqrt(a ** 3 * s3 + null)
    psi = np.sqrt(a ** 2 * s4)
    return d, psi


def local_minima(values):
    values = np.asarray(values, float)
    return np.array(
        [i for i in range(1, len(values) - 1) if values[i] <= values[i - 1] and values[i] < values[i + 1]],
        dtype=np.int64,
    )


def psi_mq_tikhonov(psi_q, d_md, p0=P0_TIKHONOV):
    psi_q = np.asarray(psi_q, float)
    d_md = np.asarray(d_md, float)
    n = psi_q.size
    mq = np.empty(n)
    mq[-1] = psi_q[-1]
    expn = 2.0 / (2.0 + 1.0 / p0)
    for j in range(n - 2, -1, -1):
        floor = (d_md[j] / max(d_md[j + 1], TINY)) ** expn * mq[j + 1]
        mq[j] = max(psi_q[j], floor)
    return mq


def raus2024_rule(alpha, psi_q, d_md, sigma_min_sq, p0=P0_TIKHONOV):
    alpha = np.asarray(alpha, float)
    psi_q = np.asarray(psi_q, float)
    d_md = np.asarray(d_md, float)
    mq = psi_mq_tikhonov(psi_q, d_md, p0)
    psi_hr = d_md / np.sqrt(np.maximum(alpha, TINY))
    in_hq = alpha >= max(float(alpha[0]), float(sigma_min_sq))
    i_hr = int(np.argmin(np.where(in_hq, psi_hr, np.inf)))
    mins = local_minima(psi_q)
    fallback = ""
    if mins.size == 0:
        hr_mins = local_minima(psi_hr)
        if hr_mins.size:
            i_h = int(hr_mins[np.argmin(psi_hr[hr_mins])])
            fallback = "psiHR_local_min"
        else:
            i_h = i_hr
            fallback = "psiHR_global"
        i_q_loc = i_h
        i_mq = i_h
        source = fallback
    else:
        i_q_loc = int(mins[np.argmin(psi_q[mins])])
        i_mq = int(mins[np.argmin(mq[mins])])
        r_h = min(1.0 / alpha[i_q_loc], max(1.0 / alpha[i_mq], 1.0 / alpha[i_hr]))
        i_h = int(np.argmin(np.abs(alpha - 1.0 / r_h)))
        source = "q" if i_h == i_q_loc else ("mq" if i_h == i_mq else ("hr" if i_h == i_hr else "nearest"))
    return {
        "index": i_h,
        "index_q_local": int(i_q_loc),
        "index_mq": int(i_mq),
        "index_hr": int(i_hr),
        "source": source,
        "fallback": fallback,
        "n_psiQ_local_minima": int(mins.size),
    }


def rule_b_index(qoc):
    """Current published pp rule: smallest QOC among interior local minima."""
    q = np.asarray(qoc, float)
    local = (q[1:-1] <= q[:-2]) & (q[1:-1] < q[2:])
    idx = np.flatnonzero(local) + 1
    if idx.size == 0:
        return -1
    return int(idx[np.argmin(q[idx])])


def select_row(sigma_parts, beta_parts, offsets, null_sq, alpha, qoc):
    """Apply the rule to every M2min window of one (channel, Lambda) row."""
    n = len(null_sq)
    out = {
        "index": np.full(n, -1, dtype=np.int64),
        "index_q_local": np.full(n, -1, dtype=np.int64),
        "index_mq": np.full(n, -1, dtype=np.int64),
        "index_hr": np.full(n, -1, dtype=np.int64),
        "n_local_minima": np.zeros(n, dtype=np.int64),
        "status": np.empty(n, dtype="<U32"),
        "source": np.empty(n, dtype="<U32"),
        "fallback": np.empty(n, dtype="<U32"),
        "rule_b_index": np.full(n, -1, dtype=np.int64),
        "qoc_argmin_index": np.full(n, -1, dtype=np.int64),
    }
    alpha = np.asarray(alpha, float)
    for j in range(n):
        sig = np.asarray(sigma_parts[offsets[j]:offsets[j + 1]], float)
        bet = np.asarray(beta_parts[offsets[j]:offsets[j + 1]], float)
        d, q = quantities(sig, bet, float(null_sq[j]), alpha)
        r = raus2024_rule(alpha, q, d, float(np.min(sig) ** 2))
        i = int(r["index"])
        ok = np.isfinite(alpha[i]) and np.all(np.isfinite(d + q))
        out["index"][j] = i
        out["index_q_local"][j] = r["index_q_local"]
        out["index_mq"][j] = r["index_mq"]
        out["index_hr"][j] = r["index_hr"]
        out["n_local_minima"][j] = r["n_psiQ_local_minima"]
        out["status"][j] = "defined" if ok else "undefined_nonfinite"
        out["source"][j] = r["source"]
        out["fallback"][j] = r["fallback"]
        out["rule_b_index"][j] = rule_b_index(qoc[j])
        out["qoc_argmin_index"][j] = int(np.nanargmin(qoc[j]))
    return out
