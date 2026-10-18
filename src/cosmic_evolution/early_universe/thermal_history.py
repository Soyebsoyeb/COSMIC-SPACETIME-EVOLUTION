"""Thermal history of the radiation-dominated universe.

APPROXIMATION (labelled): the effective number of relativistic degrees of freedom g_*(T) is a step
function that switches off a species when T drops below its mass (T = m).  The real transitions
(electroweak crossover, QCD crossover) are smooth and require lattice-QCD equations of state.
Counting: 106.75 (full Standard Model) -> remove top (12*7/8=10.5) -> remove W+-,Z,H (10) -> remove b
(10.5) -> remove tau (3.5) -> remove c (10.5) -> 61.75 -> QCD confinement -> 17.25 (gamma, e, mu, nu, pi)
-> mu and pi annihilate -> 10.75 -> e+- annihilation (late values).
Time-temperature relation uses t = 1/(2H) (radiation domination); step thresholds of g_* are ignored when
tracking entropy exchange, giving O(10%) errors close to thresholds.
"""
from __future__ import annotations

import math

import numpy as np

from .. import constants as K
from ..config import CosmologyParams
from ..units import GeV_to_K, K_to_GeV
from .entropy import radiation_energy_density, G_STAR_S_TODAY

# (threshold temperature in GeV, g_* above... listed from high T to low T):  value applies for T >= threshold
_THRESH_GEV = np.array([173.0, 100.0, 4.2, 1.78, 1.3, 0.17, 0.1])
_G_BELOW = np.array([96.25, 86.25, 75.75, 72.25, 61.75, 17.25, 10.75])  # value for T just below that threshold
_G_TOP = 106.75
T_EE_GEV = 5.11e-4  # e+- annihilation at T ~ m_e (GeV)


def g_star_late(N_eff: float = K.N_EFF_STANDARD) -> float:
    """g_* after e+- annihilation: 2 + (7/4)(4/11)^{4/3} N_eff."""
    return 2.0 + 7.0 / 4.0 * (4.0 / 11.0) ** (4.0 / 3.0) * N_eff


def g_star(T_GeV, N_eff: float = K.N_EFF_STANDARD):
    """Energy-density relativistic degrees of freedom g_*(T) (step approximation)."""
    T = np.atleast_1d(np.asarray(T_GeV, dtype=float))
    out = np.empty_like(T)
    for i, t in enumerate(T):
        if t >= _THRESH_GEV[0]:
            out[i] = _G_TOP
        elif t >= T_EE_GEV:
            idx = np.searchsorted(-_THRESH_GEV, -t, side="right") - 1  # last threshold >= t
            out[i] = _G_BELOW[idx]
        else:
            out[i] = g_star_late(N_eff)
    return out if np.ndim(T_GeV) else float(out[0])


def g_star_s(T_GeV, N_eff: float = K.N_EFF_STANDARD):
    """Entropy degrees of freedom g_{*s}(T): equals g_* above e+- annihilation; 3.909 afterwards."""
    T = np.atleast_1d(np.asarray(T_GeV, dtype=float))
    g = np.atleast_1d(np.asarray(g_star(T, N_eff)))
    out = np.where(T < T_EE_GEV, G_STAR_S_TODAY, g)
    return out if np.ndim(T_GeV) else float(out[0])


def temperature_of_z(z, T0: float = K.T_CMB0):
    """Photon temperature T(z) = T0 (1+z) (adiabatic expansion of a blackbody)."""
    z = np.asarray(z, dtype=float)
    if np.any(z <= -1):
        raise ValueError("z must be > -1")
    return T0 * (1.0 + z)


def hubble_rate_radiation(T_GeV, N_eff: float = K.N_EFF_STANDARD):
    """H(T) in 1/s for a purely radiation-dominated universe with g_*(T) (flat)."""
    T_K = GeV_to_K(T_GeV)
    rho_c2 = radiation_energy_density(T_K, g_star(T_GeV, N_eff))  # J/m^3
    return np.sqrt(8.0 * math.pi * K.G / 3.0 * rho_c2 / K.c**2)


def time_of_temperature_s(T_GeV, N_eff: float = K.N_EFF_STANDARD):
    """Cosmic time t = 1/(2H) in seconds in the radiation era."""
    return 1.0 / (2.0 * hubble_rate_radiation(T_GeV, N_eff))


def temperature_of_time_gev(t_s, N_eff: float = K.N_EFF_STANDARD):
    """Invert t(T) by bisection in log T (monotonic)."""
    from scipy import optimize
    f = lambda lt, t: np.log(time_of_temperature_s(10 ** lt, N_eff)) - np.log(t)
    t_s = np.atleast_1d(np.asarray(t_s, dtype=float))
    out = [10 ** optimize.brentq(f, -12, 19, args=(t,)) for t in t_s]
    return np.array(out) if np.ndim(t_s) else out[0]


def epoch_table(params: CosmologyParams):
    """Characteristic epochs: (name, T_GeV approx, t_s). Temperatures at T~mass/transition scales."""
    rows = []
    for name, T in [("Electroweak scale (crossover; T ~ 100 GeV)", 100.0),
                    ("QCD confinement (crossover; T ~ 0.17 GeV)", 0.17),
                    ("Neutrino decoupling (T ~ 1 MeV)", 1e-3),
                    ("e+- annihilation (T ~ m_e)", T_EE_GEV),
                    ("Nucleosynthesis window (T ~ 0.07 MeV)", 7e-5)]:
        rows.append((name, T, float(time_of_temperature_s(T, params.N_eff))))
    return rows
