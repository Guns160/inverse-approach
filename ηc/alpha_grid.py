"""Historical grid construction only; no historical selection rule is shipped."""
import math
import numpy as np

def make_alpha_grid(log10_min: float, log10_max: float, step: float) -> np.ndarray:
    if step <= 0.0:
        raise ValueError("log10_alpha_step must be positive")
    if log10_max < log10_min:
        raise ValueError("log10_alpha_max must be >= log10_alpha_min")
    n = int(round((log10_max - log10_min) / step)) + 1
    return 10.0 ** np.linspace(log10_min, log10_max, n)


def quantities(
    sigma: np.ndarray,
    beta: np.ndarray,
    null: float,
    alpha: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (d_D, d_MD, ψ_Q) on the given α grid."""
    sigma = np.asarray(sigma, float).reshape(-1, 1)
    beta = np.asarray(beta, float).reshape(-1, 1)
    if sigma.shape != beta.shape:
        raise ValueError(f"sigma and beta length mismatch: {sigma.size} vs {beta.size}")
    alpha = np.asarray(alpha, float).reshape(-1)
    if alpha.size == 0:
        raise ValueError("alpha grid is empty")
    if np.any(alpha <= 0.0):
        raise ValueError("alpha must be strictly positive")
    null = float(null)
    if null < 0.0:
        raise ValueError("orthogonal residual square must be >= 0")
    den = sigma**2 + alpha
    d_d = np.sqrt(np.sum((alpha * beta / den) ** 2, axis=0) + null)
    d_md = np.sqrt(np.sum((alpha**1.5 * beta / den**1.5) ** 2, axis=0) + null)
    psi_q = np.sqrt(np.sum((alpha * sigma * beta / den**2) ** 2, axis=0))
    return d_d, d_md, psi_q


def local_maxima(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    return np.array(
        [
            i
            for i in range(1, len(values) - 1)
            if values[i] >= values[i - 1] and values[i] > values[i + 1]
        ],
        dtype=np.int64,
    )


def alpha_grid_with_small_peak(sigma, beta, null, log10_max=0.0, step=0.01, floor=-80.0):
    sig_min = max(float(np.min(sigma)), 1e-300)
    log10_min = min(-50.0, 2.0 * math.log10(sig_min) - 8.0)
    log10_min = max(log10_min, floor)
    last = None
    for _ in range(10):
        alpha = make_alpha_grid(log10_min, log10_max, step)
        _, _, psi = quantities(sigma, beta, null, alpha)
        maxima = local_maxima(psi)
        has_left_peak = (
            maxima.size > 0
            and int(maxima[0]) >= 3
            and float(psi[0]) < float(psi[int(maxima[0])])
        )
        last = alpha
        if has_left_peak:
            return alpha
        if log10_min <= floor:
            break
        log10_min = max(floor, log10_min - 5.0)
    return last


def selection_grid(sigma, beta, null):
    old = alpha_grid_with_small_peak(sigma, beta, null)
    _, _, psi = quantities(sigma, beta, null, old)
    peaks = local_maxima(psi)
    covered = bool(peaks.size and int(peaks[0]) >= 3 and psi[0] < psi[int(peaks[0])])
    # Exactly watch_and_select.py: log10 of the stored first alpha, then linspace.
    actual = 10.0 ** np.linspace(float(np.log10(old[0])), 0.0, len(old))
    return actual, covered
