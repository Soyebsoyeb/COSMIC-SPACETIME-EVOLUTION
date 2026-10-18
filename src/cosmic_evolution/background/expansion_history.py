"""Tabulated expansion history: one-stop table of background quantities on a redshift grid."""
from __future__ import annotations

import numpy as np

from ..geometry import flrw
from .cosmic_time import age_gyr, cosmic_time_gyr
from .friedmann import FLRWBackground


def expansion_table(bg: FLRWBackground, z) -> dict[str, np.ndarray]:
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)
    fr = bg.omega_fractions(a)
    return {
        "z": z,
        "a": a,
        "t_Gyr": np.asarray(cosmic_time_gyr(bg, a)),
        "H_km_s_Mpc": np.asarray(bg.H_of_z(z)),
        "E": np.asarray(bg.E_of_z(z)),
        "q": np.asarray(bg.deceleration(a)),
        "Omega_r": fr["radiation"],
        "Omega_m": fr["matter"],
        "Omega_k": fr["curvature"],
        "Omega_de": fr["dark_energy"],
        "D_C_Mpc": np.asarray(flrw.comoving_distance_mpc(bg, z)),
        "D_L_Mpc": np.asarray(flrw.luminosity_distance_mpc(bg, z)),
        "D_A_Mpc": np.asarray(flrw.angular_diameter_distance_mpc(bg, z)),
    }


def acceleration_onset_a(bg: FLRWBackground) -> float | None:
    """Largest scale factor in (1e-3, 1) region where q crosses zero from above (first onset of acceleration)."""
    from scipy import optimize
    grid = np.logspace(-3, 0, 3000)
    q = bg.deceleration(grid)
    s = np.where((q[:-1] > 0) & (q[1:] <= 0))[0]
    if len(s) == 0:
        return None
    i = s[-1]
    return float(optimize.brentq(lambda x: float(bg.deceleration(x)), grid[i], grid[i + 1], xtol=1e-14))
