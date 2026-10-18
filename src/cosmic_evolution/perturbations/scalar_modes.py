"""Linear scalar (curvature) perturbations during inflation: Mukhanov-Sasaki equation in e-fold form.

For the comoving curvature perturbation R_k on a single-field background (M_P = 1):

    R'' + (3 - eps + eta) R' + (k / aH)^2 R = 0,     ' = d/dN, eta = d ln eps / dN.

R is gauge invariant (comoving gauge) and conserved on super-horizon scales for adiabatic single-field
perturbations.  The Bunch-Davies initial condition is imposed deep inside the horizon,
    R -> exp(-i k tau) / (z sqrt(2k)),  z = a sqrt(2 eps),
with the first-order WKB correction dR/dN = -(i k/aH + 1 + eta/2) R.  The dimensionless spectrum is
P_R = k^3 |R_k|^2 / (2 pi^2) evaluated well outside the horizon.
LIMITATIONS: linear theory, single field, no isocurvature, no non-Gaussianity, no backreaction.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import integrate
from scipy.interpolate import CubicSpline

from ..early_universe.inflation import InflationSolution, hubble_of_N


@dataclass
class ModeBackground:
    """Background functions of N: eps(N), eta(N), H(N) (M_P units), plus the integration range."""
    eps: Callable[[float], float]
    eta: Callable[[float], float]
    H: Callable[[float], float]
    N_min: float
    N_max: float

    @classmethod
    def from_inflation(cls, sol: InflationSolution, V0: float) -> "ModeBackground":
        eps = CubicSpline(sol.N, sol.eps_H)
        eta = CubicSpline(sol.N, sol.eta_H)
        lnH = CubicSpline(sol.N, np.log(hubble_of_N(sol, V0)))
        return cls(lambda n: float(eps(n)), lambda n: float(eta(n)), lambda n: math.exp(float(lnH(n))),
                   float(sol.N[0]), float(sol.N[-1]))

    @classmethod
    def quasi_de_sitter(cls, eps: float, H_init: float = 1e-5, N_total: float = 80.0) -> "ModeBackground":
        """Constant-eps background: H = H_init exp(-eps N) (exact power-law inflation)."""
        return cls(lambda n: eps, lambda n: 0.0, lambda n: H_init * math.exp(-eps * n), 0.0, N_total)


def k_of_exit(bg: ModeBackground, N_exit: float) -> float:
    """Comoving wavenumber (M_P units, a(N=0)=1) that crosses the horizon k = aH at N_exit."""
    return math.exp(N_exit) * bg.H(N_exit)


def _find_N_for_x(bg: ModeBackground, k: float, x_target: float, N_lo: float, N_hi: float) -> float:
    from scipy.optimize import brentq
    f = lambda n: k * math.exp(-n) / bg.H(n) - x_target
    return brentq(f, N_lo, N_hi, xtol=1e-12)


def evolve_curvature_mode(bg: ModeBackground, k: float, x_start: float = 400.0, x_end: float = 0.01,
                          rtol: float = 1e-10, atol: float = 1e-40) -> float:
    """Return P_R(k) from numerical integration of the Mukhanov-Sasaki equation."""
    N0 = _find_N_for_x(bg, k, x_start, bg.N_min, bg.N_max)
    N1 = _find_N_for_x(bg, k, x_end, N0, bg.N_max)
    a0 = math.exp(N0)
    eps0 = bg.eps(N0)
    amp = 1.0 / (a0 * math.sqrt(2.0 * eps0) * math.sqrt(2.0 * k))
    x0 = k / (a0 * bg.H(N0))
    c0 = -(1j * x0 + 1.0 + 0.5 * bg.eta(N0))
    R0 = amp + 0j
    dR0 = c0 * R0

    def rhs(N, y):
        R = y[0] + 1j * y[1]
        dR = y[2] + 1j * y[3]
        x = k * math.exp(-N) / bg.H(N)
        d2 = -(3.0 - bg.eps(N) + bg.eta(N)) * dR - x * x * R
        return [dR.real, dR.imag, d2.real, d2.imag]

    sol = integrate.solve_ivp(rhs, [N0, N1], [R0.real, R0.imag, dR0.real, dR0.imag], method="DOP853",
                              rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    R = sol.y[0, -1] + 1j * sol.y[1, -1]
    return float(k**3 / (2.0 * math.pi**2) * abs(R) ** 2)


def slow_roll_amplitude(bg: ModeBackground, N_exit: float) -> float:
    """Leading-order P_R = H^2/(8 pi^2 eps) at horizon exit (k = aH)."""
    return bg.H(N_exit) ** 2 / (8.0 * math.pi**2 * bg.eps(N_exit))
