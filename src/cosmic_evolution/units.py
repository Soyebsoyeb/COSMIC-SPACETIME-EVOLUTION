"""Unit conversion helpers. Internal convention: dimensionless time tau = H0 t, distances in c/H0."""
from __future__ import annotations

import numpy as np

from . import constants as K


def H0_to_inv_s(H0_km_s_Mpc: float) -> float:
    """Convert H0 from km/s/Mpc to 1/s."""
    if H0_km_s_Mpc <= 0:
        raise ValueError("H0 must be positive")
    return H0_km_s_Mpc * K.km / K.Mpc


def hubble_time_Gyr(H0_km_s_Mpc: float) -> float:
    """1/H0 in Gyr."""
    return 1.0 / H0_to_inv_s(H0_km_s_Mpc) / K.Gyr


def hubble_distance_Mpc(H0_km_s_Mpc: float) -> float:
    """c/H0 in Mpc."""
    if H0_km_s_Mpc <= 0:
        raise ValueError("H0 must be positive")
    return K.c / K.km / H0_km_s_Mpc


def critical_density_kg_m3(H0_km_s_Mpc: float) -> float:
    """rho_crit = 3 H0^2 / (8 pi G) in kg/m^3."""
    H = H0_to_inv_s(H0_km_s_Mpc)
    return 3.0 * H**2 / (8.0 * np.pi * K.G)


def Gyr_to_s(t_gyr):
    return np.asarray(t_gyr) * K.Gyr


def s_to_Gyr(t_s):
    return np.asarray(t_s) / K.Gyr


def Mpc_to_Gly(d_mpc):
    """Mpc to billion light-years."""
    return np.asarray(d_mpc) * K.Mpc / (K.c * K.yr) / 1e9


def GeV_to_K(E_GeV):
    """Temperature in kelvin corresponding to k_B T = E."""
    return np.asarray(E_GeV) * 1e9 * K.eV / K.k_B


def K_to_GeV(T_K):
    return np.asarray(T_K) * K.k_B / (1e9 * K.eV)


def K_to_MeV(T_K):
    return np.asarray(T_K) * K.k_B / (1e6 * K.eV)
