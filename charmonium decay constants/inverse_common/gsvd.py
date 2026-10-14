"""Actual compact LAPACKE GSVD, copied without changing decomposition arithmetic.
Only the required library discovery and decomposition routines are retained.
"""
from __future__ import annotations
from dataclasses import dataclass
import ctypes
import ctypes.util
from functools import lru_cache
import glob
import os
from pathlib import Path
import sys
from typing import Any
import numpy as np
from .numerics import _dot
_LAPACK_COL_MAJOR = 102
_LAPACK_INT = ctypes.c_int
_DOUBLE_PTR = ctypes.POINTER(ctypes.c_double)
_INT_PTR = ctypes.POINTER(_LAPACK_INT)

@dataclass(frozen=True)
class CompactGSVD:
    """Compact GSVD in the ordering used by Regularization Tools ``cgsvd``.

    The first ``p`` columns are regularized modes.  The final ``k`` columns are
    the unregularized amplitude modes.  Thus ``A = U C W`` and ``B = V S W``
    with ``X = inv(W)``.
    """

    U: np.ndarray
    V: np.ndarray
    X: np.ndarray
    W: np.ndarray
    sigma: np.ndarray
    mu: np.ndarray
    sigma_unregularized: np.ndarray
    source: str
    library: str

    @property
    def p(self) -> int:
        return int(self.sigma.size)

    @property
    def k(self) -> int:
        return int(self.sigma_unregularized.size)


def _candidate_libraries() -> list[str]:
    """Return likely shared libraries that may export LAPACKE_dggsvd3."""

    candidates: list[str] = []
    explicit = os.environ.get("JPSI_LAPACKE_LIBRARY") or os.environ.get("ETA_C_LAPACKE_LIBRARY")
    if explicit:
        candidates.append(explicit)
    for name in ("lapacke", "lapack", "openblas", "blas"):
        found = ctypes.util.find_library(name)
        if found:
            candidates.append(found)
    for prefix in (Path(sys.prefix), Path("/usr"), Path("/usr/local"), Path("/opt"), Path("/opt/homebrew/opt/lapack"), Path("/usr/local/opt/lapack")):
        for directory in (prefix / "lib", prefix / "lib64"):
            for pattern in ("liblapacke.so*", "liblapack.so*", "libopenblas.so*", "libblas.so*", "liblapacke*.dylib", "liblapack*.dylib", "libopenblas*.dylib"):
                candidates.extend(sorted(glob.glob(str(directory / pattern))))
    numpy_root = Path(np.__file__).resolve().parent
    for directory in (numpy_root.parent / "numpy.libs", numpy_root / ".dylibs", numpy_root.parent / ".dylibs"):
        if directory.is_dir():
            for pattern in ("*lapack*", "*openblas*", "*blas*"):
                candidates.extend(sorted(glob.glob(str(directory / pattern))))
    for entry in sys.path:
        root = Path(entry)
        for directory in (root / "scipy" / ".dylibs", root / "scipy.libs"):
            if directory.is_dir():
                for pattern in ("*openblas*.dylib", "*openblas*.so*", "*lapack*.so*"):
                    candidates.extend(sorted(glob.glob(str(directory / pattern))))
    seen: set[str] = set()
    return [item for item in candidates if item and not (item in seen or seen.add(item))]


@lru_cache(maxsize=1)
def _load_dggsvd3() -> tuple[Any, str] | None:
    """Locate LAPACKE ``dggsvd3`` without relying on SciPy wrappers."""

    for candidate in _candidate_libraries():
        try:
            library = ctypes.CDLL(candidate)
            routine = library.LAPACKE_dggsvd3
        except (OSError, AttributeError):
            continue
        config = getattr(library, "openblas_get_config", None)
        if config is not None:
            config.restype = ctypes.c_char_p
            if b"USE64BITINT" in (config() or b""):
                if candidate == (os.environ.get("JPSI_LAPACKE_LIBRARY") or os.environ.get("ETA_C_LAPACKE_LIBRARY")):
                    raise RuntimeError("The specified OpenBLAS uses ILP64; this GSVD interface requires LP64.")
                continue
        routine.restype = _LAPACK_INT
        routine.argtypes = [
            _LAPACK_INT,
            ctypes.c_char,
            ctypes.c_char,
            ctypes.c_char,
            _LAPACK_INT,
            _LAPACK_INT,
            _LAPACK_INT,
            _INT_PTR,
            _INT_PTR,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _DOUBLE_PTR,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _INT_PTR,
        ]
        return routine, str(candidate)
    return None


def lapacke_gsvd_available() -> bool:
    """Return whether a LAPACKE implementation exporting ``dggsvd3`` is visible."""

    return _load_dggsvd3() is not None


def lapacke_info():
    """Probe installed libraries, without installing anything."""
    loaded = _load_dggsvd3()
    return {"available": loaded is not None, "library": loaded[1] if loaded else None}


def _compact_from_lapacke_direct(A: np.ndarray, B: np.ndarray) -> CompactGSVD:
    """Direct LAPACKE GSVD; retained for validation of the QR reduction."""
    loaded = _load_dggsvd3()
    if loaded is None:
        raise RuntimeError("LAPACKE_dggsvd3 was not found")
    routine, library_path = loaded
    m, n = A.shape
    p, n_b = B.shape
    if n_b != n or m < n or p >= n:
        raise ValueError("compact M3 GSVD requires m >= n > p")

    a_work = np.array(A, dtype=np.float64, order="F", copy=True)
    b_work = np.array(B, dtype=np.float64, order="F", copy=True)
    alpha = np.empty(n, dtype=np.float64)
    beta = np.empty(n, dtype=np.float64)
    u = np.empty((m, m), dtype=np.float64, order="F")
    v = np.empty((p, p), dtype=np.float64, order="F")
    q = np.empty((n, n), dtype=np.float64, order="F")
    iwork = np.empty(n, dtype=np.intc)
    k = _LAPACK_INT()
    ell = _LAPACK_INT()
    info = routine(
        _LAPACK_COL_MAJOR,
        b"U",
        b"V",
        b"Q",
        m,
        n,
        p,
        ctypes.byref(k),
        ctypes.byref(ell),
        a_work.ctypes.data_as(_DOUBLE_PTR),
        m,
        b_work.ctypes.data_as(_DOUBLE_PTR),
        p,
        alpha.ctypes.data_as(_DOUBLE_PTR),
        beta.ctypes.data_as(_DOUBLE_PTR),
        u.ctypes.data_as(_DOUBLE_PTR),
        m,
        v.ctypes.data_as(_DOUBLE_PTR),
        p,
        q.ctypes.data_as(_DOUBLE_PTR),
        n,
        iwork.ctypes.data_as(_INT_PTR),
    )
    if info != 0:
        raise RuntimeError(f"LAPACKE_dggsvd3 failed with INFO={info}")
    if k.value + ell.value != n or ell.value != p or k.value != n - p:
        raise RuntimeError(
            "unexpected dggsvd3 ranks: "
            f"k={k.value}, l={ell.value}, expected k={n - p}, l={p}"
        )
    if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(beta)):
        raise RuntimeError("LAPACKE_dggsvd3 returned non-finite generalized singular values")

    raw_r = np.triu(a_work[:n, :n])
    if np.any(np.abs(np.diag(raw_r)) <= np.finfo(float).tiny):
        raise RuntimeError("LAPACKE_dggsvd3 returned a singular triangular factor")
    raw_w = _dot(raw_r, q.T)
    compact_raw = np.arange(k.value, n)
    null_raw = np.arange(k.value)
    ordering = np.r_[compact_raw, null_raw]
    x_raw = _dot(q, np.linalg.solve(raw_r, np.eye(n)))
    return CompactGSVD(
        U=np.asarray(u[:, :n])[:, ordering],
        V=np.asarray(v),
        X=x_raw[:, ordering],
        W=raw_w[ordering, :],
        sigma=np.asarray(alpha[compact_raw]),
        mu=np.asarray(beta[compact_raw]),
        sigma_unregularized=np.asarray(alpha[null_raw]),
        source="LAPACKE_dggsvd3",
        library=library_path,
    )


def _compact_from_lapacke(A: np.ndarray, B: np.ndarray, *, qr_reduce: bool) -> CompactGSVD:
    """Actual LAPACKE GSVD, optionally after an exact economy-QR reduction.

    For the physical M3 system ``A=K0`` is tall.  If ``A=Q_A R_A`` with
    economy ``Q_A``, minimizing ``||A x-b||`` differs from minimizing
    ``||R_A x-Q_A.T b||`` only by the alpha-independent orthogonal residual.
    The GSVD of ``(R_A, B)`` is therefore exactly the GSVD needed for the
    Tikhonov solution.  Mapping ``U_R`` back as ``Q_A @ U_R`` preserves all
    active projections, reconstruction identities, and analytic derivatives,
    while avoiding LAPACK's dense m-by-m U allocation.
    """

    if not qr_reduce:
        return _compact_from_lapacke_direct(A, B)
    q_active, r_active = np.linalg.qr(A, mode="reduced")
    reduced = _compact_from_lapacke_direct(r_active, B)
    return CompactGSVD(
        U=_dot(q_active, reduced.U),
        V=reduced.V,
        X=reduced.X,
        W=reduced.W,
        sigma=reduced.sigma,
        mu=reduced.mu,
        sigma_unregularized=reduced.sigma_unregularized,
        source="LAPACKE_dggsvd3_economy_qr",
        library=reduced.library,
    )


def compact_gsvd(A, B, *, method="lapacke", qr_reduce=True):
    if method != "lapacke":
        raise ValueError("Physical solves require LAPACKE; no Gram-matrix fallback")
    return _compact_from_lapacke(np.asarray(A, dtype=np.float64), np.asarray(B, dtype=np.float64), qr_reduce=qr_reduce)
