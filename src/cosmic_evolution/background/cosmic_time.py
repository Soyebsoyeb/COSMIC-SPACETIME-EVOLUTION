"""Cosmic time, age, lookback time, conformal time and the time-to-scale-factor map.

t(a) = (1/H0) * integral_0^a  a' da' / sqrt(Q(a')),    Q = a^4 E^2.
The integrand tends to zero as a -> 0 when radiation is present, and to sqrt(a/Om_m) otherwise,
so the lower limit 0 is exact (no finite cutoff).  Results are returned in Gyr.
"""
from __future__ import annotations

import numpy as np
from scipy import integrate, optimize

from .. import constants as K
from ..units import H0_to_inv_s
from .friedmann import FLRWBackground


def _t_integrand(bg: FLRWBackground):
    def f(a: float) -> float:
        q = float(bg.Q(a))
        if q <= 0:
            return np.inf
        return a / np.sqrt(q)
    return f


def hubble_time_gyr(bg: FLRWBackground) -> float:
    return 1.0 / H0_to_inv_s(bg.H0) / K.Gyr


def time_dimensionless(bg: FLRWBackground, a: float, epsabs: float = 1e-13, epsrel: float = 1e-11) -> float:
    """H0 * t(a) for a on the expanding branch."""
    if a <= 0:
        raise ValueError("a must be positive")
    val, _ = integrate.quad(_t_integrand(bg), 0.0, a, epsabs=epsabs, epsrel=epsrel, limit=400)
    return val


def cosmic_time_gyr(bg: FLRWBackground, a):
    a_arr = np.atleast_1d(np.asarray(a, dtype=float))
    th = hubble_time_gyr(bg)
    out = np.array([time_dimensionless(bg, float(x)) for x in a_arr]) * th
    return out if np.ndim(a) else float(out[0])


def age_gyr(bg: FLRWBackground, **kw) -> float:
    """Age of the universe within the model: t(a=1)."""
    return time_dimensionless(bg, 1.0, **kw) * hubble_time_gyr(bg)


def lookback_time_gyr(bg: FLRWBackground, z):
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    t0 = age_gyr(bg)
    out = np.array([t0 - float(cosmic_time_gyr(bg, 1.0 / (1.0 + zz))) for zz in z_arr])
    return out if np.ndim(z) else float(out[0])


def a_turnaround(bg: FLRWBackground, a_max: float = 1e4) -> float | None:
    """Scale factor where E^2 = 0 (expansion halts), or None if the model expands forever up to a_max."""
    grid = np.logspace(0, np.log10(a_max), 4000)
    q = bg.Q(grid)
    idx = np.where(q <= 0)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    if i == 0:
        return None
    return float(optimize.brentq(lambda x: float(bg.Q(x)), grid[i - 1], grid[i], xtol=1e-14))


def conformal_time_mpc(bg: FLRWBackground, a):
    """Conformal time eta(a) = c * integral_0^a da'/(a'^2 H) expressed as a comoving length in Mpc."""
    from .. import units
    dh = units.hubble_distance_Mpc(bg.H0)
    a_arr = np.atleast_1d(np.asarray(a, dtype=float))

    def f(x):
        return 1.0 / np.sqrt(float(bg.Q(x)))
    out = np.array([integrate.quad(f, 0.0, float(x), epsabs=1e-13, epsrel=1e-11, limit=400)[0] for x in a_arr]) * dh
    return out if np.ndim(a) else float(out[0])
