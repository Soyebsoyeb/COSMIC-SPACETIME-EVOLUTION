"""Spatial curvature of FLRW slices and the Ricci scalar of the background.

k = sign(-Omega_k) with physical curvature radius R_curv = (c/H0)/sqrt(|Omega_k|) today.
Ricci scalar for FLRW (signature -+++): R = 6 (a_ddot/a + H^2 + k c^2/a^2).
"""
from __future__ import annotations

import numpy as np

from .. import units
from ..background.friedmann import FLRWBackground


def curvature_sign(bg: FLRWBackground, tol: float = 1e-10) -> int:
    ok = bg.om_k
    if abs(ok) < tol:
        return 0
    return -1 if ok > 0 else 1


def curvature_radius_mpc(bg: FLRWBackground) -> float:
    """Present curvature radius a0 R_c = (c/H0)/sqrt(|Omega_k|); infinite for a flat model."""
    ok = abs(bg.om_k)
    if ok < 1e-12:
        return float("inf")
    return units.hubble_distance_Mpc(bg.H0) / np.sqrt(ok)


def ricci_scalar_over_H0sq(bg: FLRWBackground, a):
    """R c^2 / H0^2 (c = 1 units, signature -+++): R = 6 (a_ddot/a + H^2 + k/a^2).

    Using a_ddot/a = -q H^2 and k/a^2 = -Omega_k a^-2 H0^2:  R/H0^2 = 6 [ (1 - q) E^2 - Om_k a^-2 ].
    """
    a = np.asarray(a, dtype=float)
    q = bg.deceleration(a)
    return 6.0 * ((1.0 - q) * bg.E2(a) - bg.om_k * a**-2)
