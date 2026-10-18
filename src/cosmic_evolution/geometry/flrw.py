"""FLRW geometry: curvature-dependent transverse comoving distance and distance measures.

Metric: ds^2 = -c^2 dt^2 + a^2 [ dr^2/(1-k r^2) + r^2 dOmega^2 ].
Comoving line-of-sight distance D_C = (c/H0) int_0^z dz'/E(z');
transverse comoving distance D_M = D_H/sqrt(Om_k) sinh(sqrt(Om_k) D_C/D_H) (Om_k>0), D_C (Om_k=0),
D_H/sqrt(-Om_k) sin(sqrt(-Om_k) D_C/D_H) (Om_k<0); D_L = (1+z) D_M; D_A = D_M/(1+z).
"""
from __future__ import annotations

import numpy as np
from scipy import integrate

from .. import units
from ..background.friedmann import FLRWBackground


def comoving_distance_mpc(bg: FLRWBackground, z, epsabs: float = 1e-12, epsrel: float = 1e-11):
    """Line-of-sight comoving distance D_C(z) in Mpc: (c/H0) int_a^1 da / sqrt(Q(a))."""
    dh = units.hubble_distance_Mpc(bg.H0)
    z_arr = np.atleast_1d(np.asarray(z, dtype=float))
    if np.any(z_arr < 0):
        raise ValueError("distances are defined here for z >= 0")

    def f(a):
        q = float(bg.Q(a))
        if q <= 0:
            raise ValueError("Q(a) <= 0 along the light path")
        return 1.0 / np.sqrt(q)
    out = np.array([integrate.quad(f, 1.0 / (1.0 + zz), 1.0, epsabs=epsabs, epsrel=epsrel, limit=400)[0]
                    for zz in z_arr]) * dh
    return out if np.ndim(z) else float(out[0])


def transverse_comoving_distance_mpc(bg: FLRWBackground, z):
    dh = units.hubble_distance_Mpc(bg.H0)
    dc = np.asarray(comoving_distance_mpc(bg, z), dtype=float)
    ok = bg.om_k
    if abs(ok) < 1e-12:
        return dc
    s = np.sqrt(abs(ok))
    if ok > 0:
        return dh / s * np.sinh(s * dc / dh)
    return dh / s * np.sin(s * dc / dh)


def luminosity_distance_mpc(bg: FLRWBackground, z):
    return (1.0 + np.asarray(z, dtype=float)) * transverse_comoving_distance_mpc(bg, z)


def angular_diameter_distance_mpc(bg: FLRWBackground, z):
    return transverse_comoving_distance_mpc(bg, z) / (1.0 + np.asarray(z, dtype=float))


def distance_modulus(bg: FLRWBackground, z):
    return 5.0 * np.log10(luminosity_distance_mpc(bg, z)) + 25.0
