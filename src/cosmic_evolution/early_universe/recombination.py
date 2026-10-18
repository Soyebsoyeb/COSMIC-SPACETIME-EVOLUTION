"""Hydrogen recombination with the Peebles three-level-atom model (approximate).

Physics included: hydrogen only (helium treated as a neutral spectator providing Y_p of the baryon mass),
case-B recombination with the Pequignot-type fit used in RECFAST (including its empirical fudge factor F),
Ly-alpha escape (Sobolev) and the two-photon 2s-1s decay (Peebles C factor), T_matter = T_radiation.
NOT included: helium recombination, reionization, Compton-heating decoupling of T_m, radiative-transfer
corrections, higher-level resonances.  Accuracy is therefore of order a few per cent in x_e near
z ~ 1000 and worse in the freeze-out tail.  A validated recombination code (RECFAST/HyRec/CosmoRec) or a
Boltzmann solver (CAMB/CLASS) must be used for precision work.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import integrate

from .. import constants as K
from ..background.friedmann import FLRWBackground
from ..geometry import flrw
from ..units import H0_to_inv_s, critical_density_kg_m3

B1_EV = 13.598
B2_EV = B1_EV / 4.0
NU_ALPHA_EV = 0.75 * B1_EV
LAMBDA_ALPHA = K.h_planck * K.c / (NU_ALPHA_EV * K.eV)   # m
LAMBDA_2G = 8.227                                         # 1/s


def alpha_B(T_K, F: float = 1.14):
    """Case-B recombination coefficient in m^3/s (RECFAST-style fit with fudge F)."""
    t = np.asarray(T_K, dtype=float) / 1e4
    a, b, c, d = 4.309, -0.6166, 0.6703, 0.5300
    return F * 1e-19 * a * t**b / (1.0 + c * t**d)


def _n_H0(bg: FLRWBackground) -> float:
    p = bg.params
    rho_b = p.omega_b_h2 * critical_density_kg_m3(100.0)
    return (1.0 - p.Y_p) * rho_b / K.m_H


def saha_x(z, bg: FLRWBackground):
    """Saha-equilibrium ionisation fraction of hydrogen (an equilibrium approximation)."""
    z = np.asarray(z, dtype=float)
    T = bg.params.T_cmb * (1.0 + z)
    kT = K.k_B * T
    S = (K.m_e * kT / (2.0 * math.pi * K.hbar**2)) ** 1.5 * np.exp(-B1_EV * K.eV / kT)
    nH = _n_H0(bg) * (1.0 + z) ** 3
    A = S / nH
    return (-A + np.sqrt(A**2 + 4.0 * A)) / 2.0


@dataclass
class RecombinationResult:
    z: np.ndarray
    x_e: np.ndarray
    tau: np.ndarray
    visibility: np.ndarray    # g(z) = dtau/dz exp(-tau)   [per unit z]
    z_star: float             # peak of the visibility function
    z_tau1: float             # tau = 1 (hydrogen-only, no reionization)
    fwhm: float
    r_s_mpc: float            # comoving sound horizon at z_star
    D_M_mpc: float            # comoving angular-diameter distance to z_star
    theta_star: float         # r_s / D_M  (radians)


def solve_recombination(bg: FLRWBackground, z_start: float = 1800.0, z_end: float = 0.0, F: float = 1.14,
                        n_grid: int = 4000, rtol: float = 1e-9) -> RecombinationResult:
    p = bg.params
    nH0 = _n_H0(bg)
    H0s = H0_to_inv_s(bg.H0)

    def rhs(z, u):
        x = u[0]
        T = p.T_cmb * (1.0 + z)
        kT = K.k_B * T
        H = H0s * float(bg.E_of_z(z))
        nH = nH0 * (1.0 + z) ** 3
        al = float(alpha_B(T, F))
        be = al * (K.m_e * kT / (2.0 * math.pi * K.hbar**2)) ** 1.5 * math.exp(-B2_EV * K.eV / kT)
        Kf = LAMBDA_ALPHA**3 / (8.0 * math.pi * H)
        C = (1.0 + Kf * LAMBDA_2G * nH * (1.0 - x)) / (1.0 + Kf * (LAMBDA_2G + be) * nH * (1.0 - x))
        dxdt = C * (be * (1.0 - x) * math.exp(-NU_ALPHA_EV * K.eV / kT) - nH * al * x**2)
        return [-dxdt / ((1.0 + z) * H)]

    x0 = float(saha_x(z_start, bg))
    z_eval = np.linspace(z_start, z_end, n_grid)
    sol = integrate.solve_ivp(rhs, [z_start, z_end], [x0], method="Radau", rtol=rtol, atol=1e-14, t_eval=z_eval)
    if not sol.success:
        raise RuntimeError(sol.message)
    z = sol.t[::-1]
    x = sol.y[0][::-1]
    # optical depth tau(z) = int_0^z n_e sigma_T c dz'/((1+z')H)
    nH = nH0 * (1.0 + z) ** 3
    H = H0s * bg.E_of_z(z)
    dtau_dz = x * nH * K.sigma_T * K.c / ((1.0 + z) * H)
    tau = integrate.cumulative_trapezoid(dtau_dz, z, initial=0.0)
    vis = dtau_dz * np.exp(-tau)
    i = int(np.argmax(vis))
    z_star = float(z[i])
    half = vis[i] / 2.0
    lo = np.where(vis[:i] <= half)[0]
    hi = np.where(vis[i:] <= half)[0]
    fwhm = float(z[i + hi[0]] - z[lo[-1]]) if len(lo) and len(hi) else float("nan")
    z_tau1 = float(np.interp(1.0, tau, z))
    # sound horizon and angular scale
    om_g = p.omega_gamma
    om_b = p.omega_b_h2 / p.h**2

    def cs_over_H(zz):
        R = 0.75 * om_b / om_g / (1.0 + zz)
        return 1.0 / (math.sqrt(3.0 * (1.0 + R)) * float(bg.E_of_z(zz)))
    dh = K.c / K.km / bg.H0
    rs = integrate.quad(cs_over_H, z_star, np.inf, epsabs=1e-12, epsrel=1e-10, limit=500)[0] * dh
    DM = float(flrw.transverse_comoving_distance_mpc(bg, z_star))
    return RecombinationResult(z, x, tau, vis, z_star, z_tau1, fwhm, rs, DM, rs / DM)
