"""Particle, event and Hubble horizons.

Comoving particle horizon  chi_p(a) = (c/H0) int_0^a da'/sqrt(Q(a'))     (always finite with radiation)
Comoving event horizon     chi_e(a) = (c/H0) int_a^{a_end} da'/sqrt(Q(a'))  where a_end = inf (or the
turnaround scale factor for recollapsing models).  The event horizon is finite iff the integral converges;
we test this numerically by comparing tail integrals over successive decades (documented in docs/COSMIC_HORIZONS.md).
Hubble radius  R_H = c/H  (proper), comoving  c/(aH).
"""
from __future__ import annotations

import numpy as np
from scipy import integrate

from .. import units
from ..background.cosmic_time import a_turnaround
from ..background.friedmann import FLRWBackground


def _inv_sqrtQ(bg: FLRWBackground):
    def f(a: float) -> float:
        q = float(bg.Q(a))
        return 1.0 / np.sqrt(q) if q > 0 else np.inf
    return f


def particle_horizon_mpc(bg: FLRWBackground, a: float = 1.0) -> float:
    """Comoving particle horizon at scale factor a, in Mpc."""
    dh = units.hubble_distance_Mpc(bg.H0)
    val, _ = integrate.quad(_inv_sqrtQ(bg), 0.0, a, epsabs=1e-13, epsrel=1e-11, limit=400)
    return val * dh


def asymptotic_log_slope(f, a1: float = 1e6, a2: float = 1e9) -> float:
    """d ln f / d ln a between two large scale factors (the power-law index of the integrand)."""
    return float(np.log(f(a2) / f(a1)) / np.log(a2 / a1))


def event_horizon_exists(bg: FLRWBackground, a_start: float = 1.0) -> bool:
    """True if int_a^inf da/sqrt(Q) converges (or the model recollapses at finite a).

    The integrand behaves as a power law a^p at large a; the integral converges iff p < -1.  The index is measured
    numerically between a = 1e6 and 1e9 and a margin of 0.01 separates convergence from the logarithmic
    (marginal) case p = -1 realised by a pure cosmological constant in the *time* integral.
    """
    if a_turnaround(bg) is not None:
        return True
    return asymptotic_log_slope(_inv_sqrtQ(bg)) < -1.01


def event_horizon_mpc(bg: FLRWBackground, a: float = 1.0) -> float:
    """Comoving event horizon at scale factor a in Mpc; ``inf`` if it does not exist."""
    dh = units.hubble_distance_Mpc(bg.H0)
    amax = a_turnaround(bg)
    if amax is not None:
        val, _ = integrate.quad(_inv_sqrtQ(bg), a, amax, epsabs=1e-13, epsrel=1e-11, limit=400)
        return val * dh
    if not event_horizon_exists(bg):
        return float("inf")
    val, _ = integrate.quad(_inv_sqrtQ(bg), a, np.inf, epsabs=1e-13, epsrel=1e-11, limit=400)
    return val * dh


def hubble_radius_mpc(bg: FLRWBackground, a):
    """Proper Hubble radius c/H(a) in Mpc."""
    return units.hubble_distance_Mpc(bg.H0) / bg.E(a)


def comoving_hubble_radius_mpc(bg: FLRWBackground, a):
    """Comoving Hubble radius c/(aH) in Mpc."""
    a = np.asarray(a, dtype=float)
    return hubble_radius_mpc(bg, a) / a
