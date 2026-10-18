"""Convergence diagnostics for the background integrals."""
from __future__ import annotations

import numpy as np

from ..background import cosmic_time as ct
from ..background.friedmann import FLRWBackground
from ..geometry import horizons
from .integration import simpson_lna, trapezoid_lna


def age_vs_lower_cutoff(bg: FLRWBackground, a_min_values) -> np.ndarray:
    """Age (Gyr) when the integral starts at a_min instead of 0 (finite-cutoff error)."""
    from scipy import integrate
    th = ct.hubble_time_gyr(bg)
    f = lambda a: a / np.sqrt(float(bg.Q(a)))
    return np.array([integrate.quad(f, amin, 1.0, epsabs=1e-14, epsrel=1e-12, limit=500)[0] * th
                     for amin in a_min_values])


def age_vs_grid(bg: FLRWBackground, n_values, rule: str = "simpson", a_lo: float = 1e-10) -> np.ndarray:
    """Age (Gyr) from fixed log-spaced grids with n points (trapezoid or Simpson), integrand a/sqrt(Q) da."""
    th = ct.hubble_time_gyr(bg)
    fn = (lambda a: a / np.sqrt(bg.Q(a)))
    rule_fn = simpson_lna if rule == "simpson" else trapezoid_lna
    return np.array([rule_fn(fn, a_lo, 1.0, int(n)) * th for n in n_values])


def tolerance_sweep(bg: FLRWBackground, epsrel_values) -> np.ndarray:
    return np.array([ct.time_dimensionless(bg, 1.0, epsabs=0.0, epsrel=e) * ct.hubble_time_gyr(bg)
                     for e in epsrel_values])


def horizon_vs_upper_cutoff(bg: FLRWBackground, a_max_values) -> np.ndarray:
    """Event horizon (Mpc) when the upper limit is a finite a_max (cutoff error versus a_max -> infinity).

    Integrated in u = ln a (integrand a/sqrt(Q)) with breakpoints at every decade for numerical robustness."""
    from scipy import integrate
    from .. import units
    dh = units.hubble_distance_Mpc(bg.H0)
    f = lambda u: np.exp(u) / np.sqrt(float(bg.Q(np.exp(u))))
    out = []
    for amax in a_max_values:
        edges = np.linspace(0.0, np.log(amax), max(2, int(np.ceil(np.log10(amax))) + 1))
        out.append(sum(integrate.quad(f, lo, hi, epsabs=1e-15, epsrel=1e-13, limit=200)[0]
                       for lo, hi in zip(edges[:-1], edges[1:])) * dh)
    return np.array(out)
