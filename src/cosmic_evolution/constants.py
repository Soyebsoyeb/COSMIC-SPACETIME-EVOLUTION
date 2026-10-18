"""Physical constants (SI) taken from ``scipy.constants`` (CODATA) plus derived quantities.

All quantities are SI unless the name states otherwise.
"""
from __future__ import annotations

import math

from scipy import constants as _sc
from scipy.special import zeta

c = _sc.c                      # m / s
G = _sc.G                      # m^3 / (kg s^2)
hbar = _sc.hbar                # J s
h_planck = _sc.h               # J s
k_B = _sc.k                    # J / K
eV = _sc.eV                    # J
m_e = _sc.m_e                  # kg
m_p = _sc.m_p                  # kg
m_H = m_p + m_e - 13.598 * eV / c**2  # kg, hydrogen atom mass (binding energy correction)
sigma_T = _sc.physical_constants["Thomson cross section"][0]  # m^2
sigma_SB = _sc.sigma           # W / (m^2 K^4)

Mpc = 3.0856775814913673e22    # m (IAU 2015: 648000/pi au)
Gpc = 1e3 * Mpc
yr = 365.25 * 86400.0          # s, Julian year
Gyr = 1e9 * yr
km = 1e3

# Planck units (non-reduced)
l_planck = math.sqrt(hbar * G / c**3)
t_planck = math.sqrt(hbar * G / c**5)
m_planck = math.sqrt(hbar * c / G)
E_planck_J = m_planck * c**2
E_planck_GeV = E_planck_J / (1e9 * eV)
# Reduced Planck mass (8 pi G = 1/M_P^2)
M_P_reduced_kg = math.sqrt(hbar * c / (8 * math.pi * G))
M_P_reduced_GeV = M_P_reduced_kg * c**2 / (1e9 * eV)

# CMB temperature: adopted fiducial value (Fixsen 2009 quotes 2.725 +/- 0.001 K)
T_CMB0 = 2.7255                # K
N_EFF_STANDARD = 3.046         # standard-model effective neutrino number (instantaneous-decoupling + QED corrections)
ZETA3 = float(zeta(3.0))

# Conversions
GeV_to_inv_s = 1e9 * eV / hbar  # 1 GeV expressed as angular frequency (1/s)
H0_100 = 100.0 * km / Mpc       # 100 km/s/Mpc in 1/s
c_over_H0_100_Mpc = c / (100.0 * km) # Mpc, c / (100 km/s/Mpc)

M_SUN = 1.98841e30             # kg, nominal solar mass (IAU 2015 nominal GM_sun / CODATA G, rounded)
