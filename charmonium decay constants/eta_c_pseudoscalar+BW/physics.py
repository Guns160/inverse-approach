"""Running-width Breit-Wigner inputs for one excited eta_c state.

The tabulated partial widths are inclusive. Charge averages do not introduce
another charge-conjugate multiplicity. All masses and widths use GeV.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.integrate import quad_vec
from scipy.special import roots_legendre

from inverse_common.validation import positive_scalar

STATES = {
    "3S": dict(mass=4.0119, reference_mass=4.007,
               open_width=0.06458, rest_width=0.00916374,
               partial_widths=(0.06458,), channels=("DDstar",),
               gg_width=0.00916, gamma_gamma_width=0.00000374,
               scope="DDstar-only truncation; high-side DstarDstar channel omitted",
               status="diagnostic: too close to the two delta poles; poor extraction"),
    "4S": dict(mass=4.2693, reference_mass=4.276,
               open_width=0.02727, rest_width=0.00629257,
               partial_widths=(0.00012, 0.00077, 0.02610, 0.00028),
               channels=("DDstar", "DsDsstar", "DstarDstar", "DsstarDsstar"),
               gg_width=0.00629, gamma_gamma_width=0.00000257,
               scope="four-channel truncation; DD0star omitted without a 4S coupling input",
               status="alternative"),
    "5S": dict(mass=4.4508, reference_mass=4.470,
               open_width=0.07515, rest_width=0.00437178,
               partial_widths=(0.00292, 0.00021, 0.03846, 0.00107,
                               0.00092, 0.00366, 0.00013, 0.02778),
               channels=("DDstar", "DsDsstar", "DD0star", "DsDs0star",
                         "DD2star", "DstarDstar", "DsstarDsstar", "DstarD1"),
               gg_width=0.00437, gamma_gamma_width=0.00000178,
               scope="eight-channel running-width model",
               status="primary"),
    "6S": dict(mass=4.5777, reference_mass=4.612, open_width=None,
               rest_width=None, partial_widths=None, channels=None,
               known_open_subtotal=0.02273, gg_width=0.00295,
               gamma_gamma_width=0.00000120, known_subtotal=0.02568120,
               status="unavailable: channel partial widths are missing",
               blocker="DsstarDs1(2535) threshold exceeds the source mass; "
                       "DDstar(2S) and DstarD1prime lack daughter-width inputs"),
}
D = (1.86484, 1.86961)
DSTAR = (2.00696, 2.01026)
DS, DSSTAR, DS0STAR = 1.9683, 2.1121, 2.3177
D0STAR = (2.318, 2.403)
D1 = (2.4214, 2.4232)
D2STAR = (2.4626, 2.4643)
D0_WIDTH, D1_WIDTH = 0.2292725269223561, 0.03133891996221452
PI0, PIC = 0.1349768, 0.13957039
SOURCES = (
    "https://doi.org/10.1140/epjc/s10052-015-3418-5",
    "https://pdg.lbl.gov/2026/",
    "https://pdg.lbl.gov/2023/reviews/rpp2022-rev-resonances.pdf",
    "https://arxiv.org/abs/1608.01289",
)


@lru_cache(maxsize=8)
def _legendre(n):
    nodes, weights = roots_legendre(n)
    theta = (nodes + 1.0) * np.pi / 4.0
    return np.sin(theta)**2, np.sin(2.0 * theta) * np.pi / 4.0, weights


def _s_values(s):
    if np.iscomplexobj(s):
        raise ValueError("s must be real, finite and positive")
    s = np.asarray(s, dtype=float)
    if not np.all(np.isfinite(s)) or np.any(s <= 0):
        raise ValueError("s must be real, finite and positive")
    return s


def _momentum(s, a, u):
    root_u = np.sqrt(u)
    return (np.sqrt(np.maximum((s - (a + root_u)**2) *
                              (s - (a - root_u)**2), 0.0)) /
            (2.0 * np.sqrt(s)) * (np.sqrt(s) >= a + root_u))


def _barrier(q, ell, radius):
    z2 = (q * radius)**2
    if ell == 0:
        return 1.0
    if ell == 1:
        return 1.0 / (1.0 + z2)
    return 1.0 / (9.0 + 3.0 * z2 + z2*z2)


class RunningWidth:
    """One BW state, with pole-matched partial widths and physical thresholds.

    5S is the primary model; 4S is the alternative. 3S is retained for
    diagnostics because it is too close to the two delta poles. 6S is blocked
    until channel partial widths are supplied by a separate input update.
    """

    def __init__(self, state="5S", *, radius=4.0, dstar_d1_l=0,
                 quadrature_nodes=1024):
        if not isinstance(state, str) or state not in STATES:
            raise ValueError("BW state must be one of 3S, 4S, 5S, 6S")
        if state == "6S":
            raise ValueError("6S running width unavailable: channel partial widths are missing")
        if np.ndim(radius) != 0 or np.asarray(radius).dtype.kind not in "iuf":
            raise ValueError("radius must be a finite nonnegative scalar in GeV^-1")
        radius = float(radius)
        if not np.isfinite(radius) or radius < 0:
            raise ValueError("radius must be a finite nonnegative scalar in GeV^-1")
        if (isinstance(dstar_d1_l, (bool, np.bool_)) or
                not isinstance(dstar_d1_l, (int, np.integer)) or dstar_d1_l not in (0, 2)):
            raise ValueError("dstar_d1_l must be 0 (S wave) or 2 (D wave)")
        if (isinstance(quadrature_nodes, (bool, np.bool_)) or
                not isinstance(quadrature_nodes, (int, np.integer)) or quadrature_nodes < 16):
            raise ValueError("quadrature_nodes must be an integer >= 16")
        self.state = state
        self.inputs = dict(STATES[state])
        self.mass = self.inputs["mass"]
        self.radius = radius
        self.dstar_d1_l = int(dstar_d1_l)
        self.quadrature_nodes = int(quadrature_nodes)
        self.pole_width = self.inputs["open_width"] + self.inputs["rest_width"]
        self.phi_pole = self.phase_space(self.mass**2)
        phi_ref = self.phase_space(self.inputs["reference_mass"]**2)
        if np.any(self.phi_pole <= 0) or np.any(phi_ref <= 0):
            raise ValueError("all supplied channels must be open at the reference and pole masses")
        shifted = np.asarray(self.inputs["partial_widths"]) * self.phi_pole / phi_ref
        self.pole_partial_widths = self.inputs["open_width"] * shifted / np.sum(shifted)

    def _stable(self, s, a, b, ell):
        q = _momentum(s, a, b*b)
        return q**(2*ell + 1) / np.sqrt(s) * _barrier(q, ell, self.radius)

    def _convolution(self, s, a, b, width, lower, ell):
        hi = np.maximum((np.sqrt(s) - a)**2, lower)
        active = np.sqrt(s) > a + np.sqrt(lower)
        sin2, jac, weights = _legendre(self.quadrature_nodes)
        u = lower + (hi[:, None] - lower) * sin2
        measure = (hi[:, None] - lower) * jac
        rho = b*width / (np.pi * ((u - b*b)**2 + (b*width)**2))
        q = _momentum(s[:, None], a, u)
        integrand = (rho * q**(2*ell + 1) / np.sqrt(s[:, None]) *
                     _barrier(q, ell, self.radius) * measure)
        return np.sum(weights * integrand, axis=1) * active

    def _broad(self, s, spectators, poles, width, daughters, ell):
        result = np.zeros_like(s)
        for k in range(2):
            result += 0.5 * (
                2.0/3.0 * self._convolution(
                    s, spectators[k], poles[k], width, (daughters[1-k] + PIC)**2, ell) +
                1.0/3.0 * self._convolution(
                    s, spectators[k], poles[k], width, (daughters[k] + PI0)**2, ell))
        return result

    def phase_space(self, s):
        """Return s.shape + (n_channels,), averaging charges before normalization."""
        s = _s_values(s)
        shape = s.shape
        flat = s.reshape(-1)
        result = np.empty((flat.size, len(self.inputs["channels"])))
        for start in range(0, flat.size, 128):
            ss = flat[start:start + 128]
            ddstar = sum(self._stable(ss, a, b, 1) for a, b in zip(D, DSTAR)) / 2.0
            if self.state == "3S":
                values = [ddstar]
            else:
                dsdsstar = self._stable(ss, DS, DSSTAR, 1)
                vv = sum(self._stable(ss, a, a, 1) for a in DSTAR) / 2.0
                vsvs = self._stable(ss, DSSTAR, DSSTAR, 1)
                if self.state == "4S":
                    values = [ddstar, dsdsstar, vv, vsvs]
                else:
                    values = [ddstar, dsdsstar,
                              self._broad(ss, D, D0STAR, D0_WIDTH, D, 0),
                              self._stable(ss, DS, DS0STAR, 0),
                              sum(self._stable(ss, a, b, 2) for a, b in zip(D, D2STAR)) / 2.0,
                              vv, vsvs,
                              self._broad(ss, DSTAR, D1, D1_WIDTH, DSTAR, self.dstar_d1_l)]
            result[start:start + len(ss)] = np.asarray(values).T
        return result.reshape(shape + (result.shape[-1],))

    def partial_widths(self, s):
        return self.phase_space(s) / self.phi_pole * self.pole_partial_widths

    def gamma(self, s):
        return self.inputs["rest_width"] + np.sum(self.partial_widths(s), axis=-1)

    def bw(self, s):
        """Raw BW(s) in GeV^-2; no finite-interval area renormalization."""
        s = _s_values(s)
        width = self.mass * self.gamma(s)
        return width / (np.pi * ((s - self.mass**2)**2 + width**2))

    def thresholds(self):
        points = [(a + b)**2 for a, b in zip(D, DSTAR)]
        if self.state != "3S":
            points += [(DS + DSSTAR)**2, *(4*a*a for a in DSTAR), (2*DSSTAR)**2]
        if self.state == "5S":
            points += [(DS + DS0STAR)**2, *((a + b)**2 for a, b in zip(D, D2STAR))]
            for spectators, daughters in ((D, D), (DSTAR, DSTAR)):
                for k in range(2):
                    points += [(spectators[k] + daughters[1-k] + PIC)**2,
                               (spectators[k] + daughters[k] + PI0)**2]
        return sorted(set(points))

    def borel_column(self, t, lower, upper, *, epsabs=2e-14, epsrel=2e-13):
        """Integrate BW independently of the continuum hat grid."""
        t = _s_values(t)
        if t.ndim != 1 or t.size == 0:
            raise ValueError("Borel grid must be a nonempty positive vector")
        lower = positive_scalar(lower, "lower")
        upper = positive_scalar(upper, "upper")
        if upper <= lower:
            raise ValueError("upper must exceed lower")
        epsabs = positive_scalar(epsabs, "epsabs")
        epsrel = positive_scalar(epsrel, "epsrel")
        points = sorted({v for v in self.thresholds() + [self.mass**2] if lower < v < upper})
        value, error, info = quad_vec(
            lambda s: np.exp(-s/t)/t * float(self.bw(s)), lower, upper,
            points=points, epsabs=epsabs, epsrel=epsrel, norm="max",
            full_output=True, limit=2000)
        if not info.success:
            raise RuntimeError("BW quadrature failed: " + str(info.message))
        if not np.all(np.isfinite(value)) or np.any(value <= 0):
            raise FloatingPointError("BW Borel column must be positive and finite")
        return value, dict(method="adaptive s quad_vec with threshold breakpoints",
                           absolute_error_estimate=float(error), evaluations=int(info.neval),
                           epsabs=epsabs, epsrel=epsrel, area_renormalized=False,
                           points=points)

    def metadata(self):
        return dict(state=self.state, status=self.inputs["status"], mass_GeV=self.mass,
                    reference_mass_GeV=self.inputs["reference_mass"],
                    pole_width_GeV=self.pole_width, open_width_GeV=self.inputs["open_width"],
                    rest_width_GeV=self.inputs["rest_width"],
                    rest_components_GeV=dict(gg=self.inputs["gg_width"],
                                             gamma_gamma=self.inputs["gamma_gamma_width"]),
                    channel_scope=self.inputs["scope"],
                    channels=list(self.inputs["channels"]),
                    reference_partial_widths_GeV=list(self.inputs["partial_widths"]),
                    pole_partial_widths_GeV=self.pole_partial_widths.tolist(),
                    radius_GeV_inverse=self.radius, DstarD1_L=self.dstar_d1_l,
                    quadrature_nodes=self.quadrature_nodes, area_renormalized=False,
                    partial_width_source="EPJC 75 (2015) 313, Table 3, model b",
                    width_policy="combined open-charm and constant residual inputs",
                    daughter_width_policy="constant neutral width also used for charged partner",
                    charge_policy="equal charge-group average; no extra conjugate factor",
                    sources=list(SOURCES))
