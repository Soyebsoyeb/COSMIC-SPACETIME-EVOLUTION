"""Radiation thermodynamics in SI units: number, energy and entropy densities and the baryon-to-photon ratio."""
from __future__ import annotations

import math

import numpy as np

from .. import constants as K
from ..config import CosmologyParams
from ..units import critical_density_kg_m3


def _thermal_inv_length(T_K):
    """k_B T / (hbar c) in 1/m."""
    return K.k_B * np.asarray(T_K, dtype=float) / (K.hbar * K.c)


def photon_number_density(T_K):
    """n_gamma = (2 zeta(3)/pi^2) (kT/hbar c)^3 in m^-3."""
    return 2.0 * K.ZETA3 / math.pi**2 * _thermal_inv_length(T_K) ** 3


def radiation_energy_density(T_K, g_star=2.0):
    """rho c^2 = (pi^2/30) g_* (kT)^4/(hbar c)^3 in J/m^3."""
    T = np.asarray(T_K, dtype=float)
    return math.pi**2 / 30.0 * g_star * (K.k_B * T) * _thermal_inv_length(T) ** 3


def entropy_density(T_K, g_star_s=2.0):
    """s = (2 pi^2/45) g_{*s} k_B (kT/hbar c)^3 in J/(K m^3)."""
    return 2.0 * math.pi**2 / 45.0 * g_star_s * K.k_B * _thermal_inv_length(T_K) ** 3


def baryon_to_photon_ratio(params: CosmologyParams) -> float:
    """eta = n_b/n_gamma from Omega_b h^2 and T_cmb."""
    rho_b = params.omega_b_h2 * critical_density_kg_m3(100.0)  # kg/m^3 (h^2 absorbed)
    n_b = rho_b / K.m_p
    return float(n_b / photon_number_density(params.T_cmb))


G_STAR_S_TODAY = 2.0 + 7.0 / 4.0 * (4.0 / 11.0) * 3.0  # = 3.909 for three neutrino species
