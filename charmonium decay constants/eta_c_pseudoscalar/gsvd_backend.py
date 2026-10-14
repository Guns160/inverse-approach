"""Compatibility wrapper around the shared real-LAPACKE GSVD backend."""
from inverse_common.gsvd import (  # noqa: F401
    CompactGSVD, _load_dggsvd3, compact_gsvd, lapacke_gsvd_available, lapacke_info,
)
