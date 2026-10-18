"""Conditional far-future scenarios.  Every scenario is a *model*, not a prediction about the real universe.

The expanding branch is obtained from t(a) quadrature; the recollapsing branch is mirrored about the turnaround.
An independent second-order ODE for a(t) is provided for cross-checking (``integrate_scale_factor_ode``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import integrate

from ..background import cosmic_time as ct
from ..background.dark_energy import de_density_factor, w_of_a
from ..background.friedmann import FLRWBackground
from ..config import CosmologyParams


@dataclass(frozen=True)
class Scenario:
    key: str
    label: str
    params: CosmologyParams
    assumption: str


def _base(**kw) -> CosmologyParams:
    return CosmologyParams(**kw)


def default_scenarios(h: float = 0.674, om: float = 0.315) -> list[Scenario]:
    r = CosmologyParams(h=h, omega_m=om).omega_r
    return [
        Scenario("lcdm", "LCDM (Lambda > 0)", _base(h=h, omega_m=om),
                 "Cosmological constant; flat; expands forever toward de Sitter."),
        Scenario("quint", "Constant w = -0.8", _base(h=h, omega_m=om, w0=-0.8),
                 "Flat; dark energy with constant w > -1 (not a cosmological constant)."),
        Scenario("phantom12", "Phantom w = -1.2", _base(h=h, omega_m=om, w0=-1.2),
                 "Flat; null-energy-condition-violating fluid; finite-time singularity in the classical model."),
        Scenario("phantom15", "Phantom w = -1.5", _base(h=h, omega_m=om, w0=-1.5),
                 "Flat; stronger phantom; earlier classical singularity."),
        Scenario("eds", "Einstein-de Sitter (Omega_m = 1 - Omega_r)", _base(h=h, omega_m=1.0 - r, omega_de=0.0, omega_k=0.0),
                 "Hypothetical flat universe without dark energy; decelerates forever."),
        Scenario("open", "Open matter-only (Omega_m = 0.315, no Lambda)", _base(h=h, omega_m=om, omega_de=0.0),
                 "Hypothetical negatively curved universe; approaches coasting expansion."),
        Scenario("closed", "Closed matter-only (Omega_m = 2)", _base(h=h, omega_m=2.0, omega_de=0.0),
                 "Hypothetical positively curved universe; recollapses."),
        Scenario("negL", "Negative Lambda (Omega_Lambda = -0.5)", _base(h=h, omega_m=om, omega_de=-0.5),
                 "Hypothetical anti-de Sitter-like Lambda < 0; recollapses regardless of curvature."),
    ]


@dataclass
class Trajectory:
    t_gyr: np.ndarray
    a: np.ndarray
    H_over_H0: np.ndarray
    t_end_gyr: float | None     # time of the big rip / big crunch (None if expansion is eternal)


def classify_fate(bg: FLRWBackground) -> dict:
    """Classify the asymptotic classical fate of the model by explicit computation."""
    amax = ct.a_turnaround(bg)
    if amax is not None:
        t_turn = ct.time_dimensionless(bg, amax * (1 - 1e-12))
        return {"fate": "recollapse (big crunch)", "a_max": amax,
                "t_end_gyr": 2.0 * t_turn * ct.hubble_time_gyr(bg)}
    from ..geometry.horizons import asymptotic_log_slope
    f = lambda x: x / np.sqrt(float(bg.Q(x)))
    # time integral int f da ~ a^p: finite (big rip) iff p < -1 (measured numerically, margin 0.01)
    finite = asymptotic_log_slope(f) < -1.01
    if finite:
        tail, _ = integrate.quad(f, 1.0, np.inf, limit=500, epsabs=1e-12, epsrel=1e-10)
        t_end = (ct.age_gyr(bg) + tail * ct.hubble_time_gyr(bg))
        return {"fate": "finite-time singularity (big rip)", "a_max": np.inf, "t_end_gyr": t_end}
    a_hi = 1e9
    h_ratio = float(bg.E(a_hi) / bg.E(1e6))
    q_inf = float(bg.deceleration(a_hi))
    if abs(h_ratio - 1.0) < 1e-3 and q_inf < -0.99:
        return {"fate": "asymptotic de Sitter", "a_max": np.inf, "t_end_gyr": None, "q_inf": q_inf}
    if q_inf < -0.05:
        return {"fate": "eternal accelerated power-law expansion", "a_max": np.inf, "t_end_gyr": None, "q_inf": q_inf}
    if abs(q_inf) <= 0.05:
        return {"fate": "eternal coasting expansion (q -> 0)", "a_max": np.inf, "t_end_gyr": None, "q_inf": q_inf}
    return {"fate": "eternal decelerating expansion", "a_max": np.inf, "t_end_gyr": None, "q_inf": q_inf}


def _cumulative_time_gyr(bg: FLRWBackground, a_grid: np.ndarray) -> np.ndarray:
    """t(a_i) accumulated segment by segment (robust when a spans many decades)."""
    th = ct.hubble_time_gyr(bg)
    f = lambda x: x / np.sqrt(float(bg.Q(x)))
    t = np.empty_like(a_grid)
    t[0] = ct.time_dimensionless(bg, float(a_grid[0]))
    for i in range(1, len(a_grid)):
        t[i] = t[i - 1] + integrate.quad(f, a_grid[i - 1], a_grid[i], epsabs=1e-14, epsrel=1e-11, limit=200)[0]
    return t * th


def trajectory(bg: FLRWBackground, a_start: float = 0.05, n: int = 500, a_stop: float | None = None) -> Trajectory:
    """a(t), H(t)/H0 from the present-day-normalised quadrature; handles rip, crunch and eternal expansion."""
    amax = ct.a_turnaround(bg)
    fate = classify_fate(bg)
    if amax is not None:
        a_top = amax * (1.0 - 1e-9)
        # cluster points near the turnaround where a(t) is flat
        u = np.linspace(0.0, 1.0, n // 2)
        a_up = a_start * (a_top / a_start) ** u
        t_up = _cumulative_time_gyr(bg, a_up)
        t_turn = float(t_up[-1])
        a_dn = a_up[::-1][1:]
        t_dn = 2.0 * t_turn - t_up[::-1][1:]
        a = np.concatenate([a_up, a_dn])
        t = np.concatenate([t_up, t_dn])
        H = np.concatenate([bg.E(a_up), -bg.E(a_dn)])
        return Trajectory(t, a, H, fate["t_end_gyr"])
    top = a_stop if a_stop is not None else (1e12 if fate["t_end_gyr"] is not None else 1e3)
    a = np.logspace(np.log10(a_start), np.log10(top), n)
    t = _cumulative_time_gyr(bg, a)
    return Trajectory(t, a, bg.E(a), fate["t_end_gyr"])


def integrate_scale_factor_ode(bg: FLRWBackground, t_span_dimless: float, rtol: float = 1e-11):
    """Independent check: integrate a'' = -(1/2) a sum_i Omega_i a^{-3(1+w_i)}(1+3w_i) forward from a=1, a'=1 (tau = H0 t).

    Uses the acceleration equation (second order) rather than the Friedmann constraint, so it verifies the quadrature.
    Returns scipy OdeResult with y = (a, da/dtau).
    """
    p = bg.params

    def rhs(tau, y):
        a, v = y
        w = float(w_of_a(a, p.w0, p.wa))
        acc = -0.5 * a * (2.0 * bg.om_r * a**-4 + bg.om_m * a**-3
                          + bg.om_de * float(de_density_factor(a, p.w0, p.wa)) * (1.0 + 3.0 * w))
        return [v, acc]

    return integrate.solve_ivp(rhs, [0.0, t_span_dimless], [1.0, 1.0], method="DOP853", rtol=rtol, atol=1e-14,
                               dense_output=True)
