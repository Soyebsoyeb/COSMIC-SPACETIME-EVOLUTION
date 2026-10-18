"""Linear matter power spectrum and sigma(R).

    Delta_m^2(k, a) = (4/25) [k^2 / (Omega_m H0^2/c^2)]^2 D_+(a)^2 T(k)^2 Delta_R^2(k),
    Delta_R^2(k) = A_s (k/k_pivot)^{n_s - 1},   P(k) = 2 pi^2 Delta_m^2 / k^3,

with D_+ normalised to a in the matter era.  k in 1/Mpc, P in Mpc^3.  Regime of validity: linear scales
(Delta^2 << 1), neglecting massive neutrinos, BAO (no-wiggle T), running, and radiation-era corrections to the
amplitude conversion at k >~ 0.5 1/Mpc.  This is NOT a non-linear or CMB calculation.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
from scipy import integrate

from .. import units
from ..background.friedmann import FLRWBackground
from ..perturbations.growth_factor import solve_growth
from . import transfer_functions as tf


def primordial_delta2(k, p):
    k = np.asarray(k, dtype=float)
    return p.A_s * (k / p.k_pivot) ** (p.n_s - 1.0)


def matter_delta2(k, bg: FLRWBackground, a: float = 1.0, transfer: Callable | None = None, growth=None):
    p = bg.params
    transfer = transfer or tf.eh_nowiggle
    growth = growth or solve_growth(bg, a_f=max(a, 1.0))
    Dp = float(growth.D(a))
    k = np.asarray(k, dtype=float)
    H0_over_c = bg.H0 / (units.hubble_distance_Mpc(bg.H0) * bg.H0)  # = 1/(c/H0) in 1/Mpc
    pref = (4.0 / 25.0) * (k**2 / (p.omega_m * H0_over_c**2)) ** 2
    return pref * Dp**2 * transfer(k, p) ** 2 * primordial_delta2(k, p)


def matter_power(k, bg: FLRWBackground, a: float = 1.0, transfer: Callable | None = None, growth=None):
    k = np.asarray(k, dtype=float)
    return 2.0 * np.pi**2 / k**3 * matter_delta2(k, bg, a, transfer, growth)


def tophat_window(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = x < 1e-3
    out[small] = 1.0 - x[small] ** 2 / 10.0
    xs = x[~small]
    out[~small] = 3.0 * (np.sin(xs) - xs * np.cos(xs)) / xs**3
    return out


def sigma_R(R_mpc: float, bg: FLRWBackground, a: float = 1.0, transfer: Callable | None = None, growth=None,
            lnk_min: float = np.log(1e-5), lnk_max: float = np.log(1e3), n: int = 20000) -> float:
    """RMS linear fluctuation in top-hat spheres of comoving radius R (Mpc)."""
    lnk = np.linspace(lnk_min, lnk_max, n)
    k = np.exp(lnk)
    d2 = matter_delta2(k, bg, a, transfer, growth)
    integrand = d2 * tophat_window(k * R_mpc) ** 2
    return float(np.sqrt(integrate.simpson(integrand, x=lnk)))


def sigma8(bg: FLRWBackground, transfer: Callable | None = None, growth=None) -> float:
    return sigma_R(8.0 / bg.params.h, bg, 1.0, transfer, growth)
