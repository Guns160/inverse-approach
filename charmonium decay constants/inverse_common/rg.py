"""MSbar QCD running for the shared charmonium inputs.

RunDec, hep-ph/0004189v1, Eqs. (1), (2), (6), (7), (19)-(21).
Conventions: a = alpha_s/pi, t = ln(mu**2); masses/scales in GeV.
The production API keeps charm active (nf=4) and supports mcbar <= mu <= mb.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import argparse
import json
import math
from numbers import Integral, Real

import numpy as np
from scipy.integrate import solve_ivp

__all__ = [
    "RGInputs", "derive", "evolve_fixed_nf", "beta_coefficients",
    "gamma_coefficients", "decouple_alpha_down",
]

ZETA3 = 1.2020569031595942854
ZETA4 = math.pi**4 / 90.0
ZETA5 = 1.0369277551433699263
COEFFICIENT_SOURCE = "Chetyrkin, Kuhn, Steinhauser, hep-ph/0004189v1, Eqs. (2), (7), (20)-(21)"
INPUT_SOURCE = "Project inputs: mcbar(mcbar)=1.275 GeV; alpha_s^(5)(MZ)=0.1184; MZ=91.1876 GeV; mbbar(mbbar)=4.18 GeV"
RTOL = 2e-12
ATOL = 2e-14


def _positive(name: str, value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite positive real number")
    out = float(value)
    if not math.isfinite(out) or out <= 0:
        raise ValueError(f"{name} must be a finite positive real number")
    return out


def _integer(name: str, value: int, lo: int, hi: int) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or not lo <= value <= hi:
        raise ValueError(f"{name} must be an integer in [{lo}, {hi}]")
    return int(value)


def _alpha(value: float) -> float:
    value = _positive("alpha_s", value)
    if value >= math.pi:
        raise ValueError("alpha_s must be below pi (perturbative guard a < 1)")
    return value


@dataclass(frozen=True)
class RGInputs:
    """Immutable boundary data; the charm anchor scale equals mbar_c.

    Varying mbar_c moves BOTH the mass and its nf=4 anchor scale.
    mb denotes the nf=5 MSbar bottom mass at its own scale, not a pole mass.
    """

    mbar_c: float = 1.275
    alpha_mz: float = 0.1184
    mz: float = 91.1876
    mb: float = 4.18

    def __post_init__(self) -> None:
        for name in ("mbar_c", "mz", "mb"):
            object.__setattr__(self, name, _positive(name, getattr(self, name)))
        object.__setattr__(self, "alpha_mz", _alpha(self.alpha_mz))
        if not self.mbar_c < self.mb < self.mz:
            raise ValueError("require mbar_c < mb < mz for nf=5 -> nf=4 running")


def beta_coefficients(nf: int) -> tuple[float, ...]:
    """Return beta_0..beta_3 in da/dt = -sum beta_i*a**(i+2)."""
    n = _integer("nf", nf, 1, 6)
    return (
        (11 - 2*n/3) / 4,
        (102 - 38*n/3) / 16,
        (2857/2 - 5033*n/18 + 325*n*n/54) / 64,
        (149753/6 + 3564*ZETA3
         + (-1078361/162 - 6508*ZETA3/27)*n
         + (50065/162 + 6472*ZETA3/81)*n*n
         + 1093*n**3/729) / 256,
    )


def gamma_coefficients(nf: int) -> tuple[float, ...]:
    """Return gamma_0..gamma_3 in d ln(m)/dt = -sum gamma_i*a**(i+1)."""
    n = _integer("nf", nf, 1, 6)
    return (
        1.0,
        (202/3 - 20*n/9) / 16,
        (1249 + (-2216/27 - 160*ZETA3/3)*n - 140*n*n/81) / 64,
        (4603055/162 + 135680*ZETA3/27 - 8800*ZETA5
         + (-91723/27 - 34192*ZETA3/9 + 880*ZETA4 + 18400*ZETA5/9)*n
         + (5242/243 + 800*ZETA3/9 - 160*ZETA4/3)*n*n
         + (-332/243 + 64*ZETA3/27)*n**3) / 256,
    )


def decouple_alpha_down(alpha_s: float, *, nf: int = 5, matching_loops: int = 3) -> float:
    """alpha^(nf) -> alpha^(nf-1) at mu = mhbar(mhbar) ONLY.

    The heavy mass is MSbar, so L=ln(mu**2/mhbar(mu)**2)=0.
    Eq. (20)/(21): alpha_low=alpha_high*(1+d2*a_high**2+d3*a_high**3).
    Zero/one-loop matching is identity at this scale. No mass matching is
    implemented: derive never evolves the charm mass across bottom.
    """
    alpha_s = _alpha(alpha_s)
    nf = _integer("nf", nf, 2, 6)
    order = _integer("matching_loops", matching_loops, 0, 3)
    a = alpha_s / math.pi
    factor = 1.0
    if order >= 2:
        factor += (11/72) * a*a
    if order >= 3:
        factor += (564731/124416 - 82043*ZETA3/27648 - 2633*(nf-1)/31104) * a**3
    return _alpha(alpha_s * factor)


@lru_cache(maxsize=512)
def _evolve_cached(alpha_s: float, mass: float, start: float, end: float,
                   nf: int, loops: int) -> tuple[float, float]:
    if start == end:
        return alpha_s, mass
    beta = beta_coefficients(nf)[:loops]
    gamma = gamma_coefficients(nf)[:loops]

    def rhs(_t, y):
        a = y[0]
        b = beta[-1]
        g = gamma[-1]
        for i in range(loops-2, -1, -1):
            b = beta[i] + a*b
            g = gamma[i] + a*g
        return [-a*a*b, -a*g]

    def strong_coupling(_t, y):
        return 1.0 - y[0]

    strong_coupling.terminal = True
    strong_coupling.direction = 0
    try:
        sol = solve_ivp(
            rhs, (2*math.log(start), 2*math.log(end)),
            [alpha_s/math.pi, math.log(mass)], method="DOP853",
            rtol=RTOL, atol=ATOL, events=strong_coupling,
        )
    except (OverflowError, FloatingPointError, ValueError) as exc:
        raise RuntimeError(f"QCD RG integration failed: {exc}") from exc
    if not sol.success or sol.status == 1:
        detail = "reached perturbative guard alpha_s=pi" if sol.status == 1 else sol.message
        raise RuntimeError(f"QCD RG integration failed: {detail}")
    a, log_mass = sol.y[:, -1]
    if not np.isfinite(a) or not np.isfinite(log_mass) or not 0 < a < 1:
        raise RuntimeError("QCD RG integration returned a nonfinite or nonperturbative state")
    try:
        final_mass = math.exp(float(log_mass))
    except OverflowError as exc:
        raise RuntimeError("QCD RG mass overflow") from exc
    if not math.isfinite(final_mass) or final_mass <= 0:
        raise RuntimeError("QCD RG integration returned a nonpositive/nonfinite mass")
    return float(a*math.pi), final_mass


def evolve_fixed_nf(alpha_s: float, mass: float, start: float, end: float,
                    nf: int = 4, loops: int = 4) -> tuple[float, float]:
    """Return (alpha_s(end), mass(end)) by integrating beta AND gamma.

    This low-level mathematical operation holds nf fixed; it performs no
    threshold inference or matching. Physical charmonium inputs use derive.
    At loops=N, beta and gamma include precisely indices 0..N-1.
    """
    return _evolve_cached(
        _alpha(alpha_s), _positive("mass", mass),
        _positive("start", start), _positive("end", end),
        _integer("nf", nf, 1, 6), _integer("loops", loops, 1, 4),
    )


@lru_cache(maxsize=128)
def _boundary(inputs: RGInputs, loops: int) -> tuple[float, float, float]:
    a5_mb, _ = evolve_fixed_nf(inputs.alpha_mz, 1.0, inputs.mz, inputs.mb, 5, loops)
    a4_mb = decouple_alpha_down(a5_mb, nf=5, matching_loops=loops-1)
    a4_anchor, _ = evolve_fixed_nf(a4_mb, 1.0, inputs.mb, inputs.mbar_c, 4, loops)
    return a5_mb, a4_mb, a4_anchor


def derive(mu: float = 2.0, loops: int = 4, inputs: RGInputs | None = None) -> dict:
    """Derive physical nf=4 charm inputs, with N-loop RG/(N-1)-loop matching.

    Only mbar_c <= mu <= mb is supported; mu=mb returns the LOW side (nf=4).
    The nf=4 mass anchor is mc(mbar_c)=mbar_c with no charm decoupling.
    Each call returns fresh metadata, so callers cannot corrupt cached states.
    """
    mu = _positive("mu", mu)
    loops = _integer("loops", loops, 1, 4)
    if inputs is None:
        inputs = RGInputs()
    if not isinstance(inputs, RGInputs):
        raise TypeError("inputs must be an RGInputs instance or None")
    if not inputs.mbar_c <= mu <= inputs.mb:
        raise ValueError(
            f"mu must lie in [mbar_c, mb]=[{inputs.mbar_c}, {inputs.mb}] GeV; "
            "charm stays active (nf=4); mass matching outside this domain is not implemented"
        )
    a5_mb, a4_mb, a4_anchor = _boundary(inputs, loops)
    alpha_s, mc = evolve_fixed_nf(a4_anchor, inputs.mbar_c, inputs.mbar_c, mu, 4, loops)
    return {
        "mu": mu, "mc": mc, "alpha_s": alpha_s, "nf": 4,
        "loops": loops, "matching_loops": loops-1, "inputs": asdict(inputs),
        "input_source": INPUT_SOURCE if inputs == RGInputs() else "Caller-supplied RGInputs",
        "coefficient_source": COEFFICIENT_SOURCE,
        "scheme": "MSbar", "mass_method": "ode",
        "a_convention": "alpha_s/pi", "t_convention": "ln(mu^2)",
        "mass_anchor": {"mu": inputs.mbar_c, "mc": inputs.mbar_c,
                        "alpha_s": a4_anchor, "nf": 4},
        "threshold": {"mu": inputs.mb, "heavy_mass_scheme": "MSbar at own scale",
                      "nf_high": 5, "nf_low": 4, "alpha_s_high": a5_mb,
                      "alpha_s_low": a4_mb, "matching_loops": loops-1},
        "integration": {"method": "DOP853", "rtol": RTOL, "atol": ATOL},
    }


def main(argv: list[str] | None = None) -> int:
    """Independent running CLI; emit one strict JSON result, write no files."""
    parser = argparse.ArgumentParser(
        description="MSbar charm QCD RG: 1-4 loop running with matching order loops-1. "
                    "Print JSON to stdout; target must lie between mcbar and mb.",
        allow_abbrev=False,
    )
    defaults = RGInputs()
    parser.add_argument("--mu", type=float, default=2.0, help="target scale in GeV (default: 2)")
    parser.add_argument("--loops", type=int, choices=range(1, 5), default=4,
                        help="running loop order (default: 4)")
    parser.add_argument("--mcbar", type=float, default=defaults.mbar_c,
                        help="nf=4 mcbar(mcbar), moving mass/scale anchor in GeV")
    parser.add_argument("--alpha-mz", type=float, default=defaults.alpha_mz,
                        help="nf=5 alpha_s at MZ")
    parser.add_argument("--mb", type=float, default=defaults.mb,
                        help="nf=5 MSbar bottom self mass/threshold in GeV")
    parser.add_argument("--mz", type=float, default=defaults.mz, help="MZ reference scale in GeV")
    args = parser.parse_args(argv)
    try:
        inputs = RGInputs(mbar_c=args.mcbar, alpha_mz=args.alpha_mz, mz=args.mz, mb=args.mb)
        result = derive(mu=args.mu, loops=args.loops, inputs=inputs)
        output = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)
    except (ValueError, TypeError) as exc:
        parser.error(str(exc))
    except RuntimeError as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
