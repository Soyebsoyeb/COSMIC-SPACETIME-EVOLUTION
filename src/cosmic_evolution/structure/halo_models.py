"""Halo mass functions from the Press-Schechter and Sheth-Tormen ansatz, built on linear theory.

APPROXIMATE / MODEL-DEPENDENT: these are analytic fits calibrated to N-body simulations in the literature; they
are not derived from the homogeneous background.  delta_c = 1.686 is the linear threshold of spherical collapse in
an Einstein-de Sitter universe (weakly cosmology dependent).  No halo profiles, baryonic physics or assembly bias.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import special

from .. import constants as K
from ..background.friedmann import FLRWBackground
from ..units import critical_density_kg_m3
from .power_spectrum import sigma_R

DELTA_C = 1.686
ST_A, ST_a, ST_p = 0.3222, 0.707, 0.3


def mean_matter_density_msun_mpc3(bg: FLRWBackground) -> float:
    rho_c = critical_density_kg_m3(bg.H0)  # kg/m^3 (present)
    return bg.om_m * rho_c * K.Mpc**3 / K.M_SUN


def radius_of_mass(M_msun, bg: FLRWBackground):
    """Lagrangian comoving radius (Mpc) enclosing mass M at the mean density."""
    rho = mean_matter_density_msun_mpc3(bg)
    return (3.0 * np.asarray(M_msun, dtype=float) / (4.0 * math.pi * rho)) ** (1.0 / 3.0)


def f_press_schechter(nu):
    nu = np.asarray(nu, dtype=float)
    return math.sqrt(2.0 / math.pi) * nu * np.exp(-0.5 * nu**2)


def f_sheth_tormen(nu):
    nu = np.asarray(nu, dtype=float)
    a = ST_a
    return ST_A * math.sqrt(2.0 * a / math.pi) * (1.0 + (a * nu**2) ** (-ST_p)) * nu * np.exp(-0.5 * a * nu**2)


def st_normalisation_analytic() -> float:
    """A such that integral f(nu) dln(nu) = 1:  A = [1 + 2^{-p} Gamma(1/2 - p)/sqrt(pi)]^{-1}."""
    return 1.0 / (1.0 + 2.0 ** (-ST_p) * special.gamma(0.5 - ST_p) / math.sqrt(math.pi))


def mass_function(M_msun, bg: FLRWBackground, kind: str = "ST", a: float = 1.0, transfer=None, growth=None):
    """dn/dlnM in Mpc^-3 (comoving) and sigma(M)."""
    M = np.atleast_1d(np.asarray(M_msun, dtype=float))
    rho = mean_matter_density_msun_mpc3(bg)
    sig = np.array([sigma_R(float(radius_of_mass(m, bg)), bg, a, transfer, growth) for m in M])
    lnM = np.log(M)
    dlns = np.gradient(np.log(1.0 / sig), lnM)
    nu = DELTA_C / sig
    f = f_sheth_tormen(nu) if kind.upper() == "ST" else f_press_schechter(nu)
    return rho / M * f * np.abs(dlns), sig
