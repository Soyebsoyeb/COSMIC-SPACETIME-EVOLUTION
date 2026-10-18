"""Matter transfer functions (fitting formulae).  APPROXIMATE, no baryon acoustic oscillations.

  * bbks: Bardeen-Bond-Kaiser-Szalay (1986) shape with the Sugiyama (1995) baryon correction to the shape parameter.
  * eh_nowiggle: Eisenstein & Hu (1998) zero-oscillation ("no-wiggle") fitting formula.
Coefficients are transcribed from the literature; they are validated here only by internal consistency
(T -> 1 as k -> 0, mutual agreement of the two fits to ~10-20% over the range of interest, and the resulting sigma_8).
Neither captures massive neutrinos, BAO or non-standard radiation content.  For precision use CAMB or CLASS.
k is in 1/Mpc throughout.
"""
from __future__ import annotations

import numpy as np

from ..config import CosmologyParams


def bbks(k, p: CosmologyParams):
    k = np.asarray(k, dtype=float)
    om = p.omega_m
    ob = p.omega_b_h2 / p.h**2
    gamma = om * p.h * np.exp(-ob * (1.0 + np.sqrt(2.0 * p.h) / om))
    theta = p.T_cmb / 2.7
    q = k * theta**2 / (gamma * p.h)
    q = np.maximum(q, 1e-30)
    L = np.log(1.0 + 2.34 * q) / (2.34 * q)
    return L * (1.0 + 3.89 * q + (16.1 * q) ** 2 + (5.46 * q) ** 3 + (6.71 * q) ** 4) ** -0.25


def eh_nowiggle(k, p: CosmologyParams):
    k = np.asarray(k, dtype=float)
    wm = p.omega_m * p.h**2
    wb = p.omega_b_h2
    fb = wb / wm
    theta = p.T_cmb / 2.7
    s = 44.5 * np.log(9.83 / wm) / np.sqrt(1.0 + 10.0 * wb**0.75)           # Mpc
    a_g = 1.0 - 0.328 * np.log(431.0 * wm) * fb + 0.38 * np.log(22.3 * wm) * fb**2
    g_eff = p.omega_m * p.h * (a_g + (1.0 - a_g) / (1.0 + (0.43 * k * s) ** 4))
    q = np.maximum(k * theta**2 / (g_eff * p.h), 1e-30)
    L0 = np.log(2.0 * np.e + 1.8 * q)
    C0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
    return L0 / (L0 + C0 * q**2)
