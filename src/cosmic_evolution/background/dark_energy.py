"""Dark-energy equation-of-state models (CPL parametrisation).

w(a) = w0 + wa (1 - a).  The density evolves as

    rho_de(a) / rho_de(1) = a^{-3 (1 + w0 + wa)} exp[-3 wa (1 - a)]

which follows from the continuity equation d ln rho / d ln a = -3 (1 + w(a)).
w0 = -1, wa = 0 is a cosmological constant.  Constant-w models have wa = 0.
"""
from __future__ import annotations

import numpy as np


def w_of_a(a, w0: float = -1.0, wa: float = 0.0):
    a = np.asarray(a, dtype=float)
    return w0 + wa * (1.0 - a)


def de_density_factor(a, w0: float = -1.0, wa: float = 0.0):
    """rho_de(a)/rho_de(a=1) for the CPL model (exact solution of the continuity equation)."""
    a = np.asarray(a, dtype=float)
    return a ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a))
