"""Reheating: perturbative inflaton decay toy model and bookkeeping formulas.

Model (labelled MODEL-DEPENDENT): a homogeneous inflaton condensate oscillating about the minimum acts as a fluid with
constant equation of state w_phi (w = 0 for a quadratic minimum, <w> = 0 for the Starobinsky oscillations as a
common approximation), decaying with constant rate Gamma into radiation:

    d rho_phi / dt = -3 H (1 + w_phi) rho_phi - Gamma rho_phi,
    d rho_r   / dt = -4 H rho_r + Gamma rho_phi,        3 M_P^2 H^2 = rho_phi + rho_r.

With x = rho_phi/rho_end, y = rho_r/rho_end, G = Gamma/H_end and e-folds N:
    dx/dN = -3(1+w) x - (G/Htilde) x,  dy/dN = -4 y + (G/Htilde) x,  Htilde = sqrt(x + y).
The reheating temperature is conventionally defined by Gamma = H(T_reh):
    T_reh = (90 / (pi^2 g_*))^{1/4} sqrt(Gamma M_P).
The numerical value of Gamma is a free parameter; no observational value of T_reh is implied.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import integrate


def reheating_temperature_estimate(gamma_over_M_P: float, g_star: float = 106.75) -> float:
    """T_reh/M_P from Gamma = H: (90/(pi^2 g))^{1/4} sqrt(Gamma/M_P)."""
    if gamma_over_M_P <= 0:
        raise ValueError("Gamma must be positive")
    return (90.0 / (math.pi**2 * g_star)) ** 0.25 * math.sqrt(gamma_over_M_P)


def temperature_from_density(rho_over_MP4: float, g_star: float = 106.75) -> float:
    """T/M_P from rho = (pi^2/30) g T^4."""
    return (30.0 * rho_over_MP4 / (math.pi**2 * g_star)) ** 0.25


@dataclass
class ReheatingSolution:
    N: np.ndarray
    x: np.ndarray            # rho_phi / rho_end
    y: np.ndarray            # rho_r / rho_end
    gamma_over_Hend: float
    w: float

    @property
    def total(self):
        return self.x + self.y

    @property
    def radiation_fraction(self):
        return self.y / self.total


def solve_perturbative_reheating(gamma_over_Hend: float, w: float = 0.0, N_max: float | None = None,
                                 rtol: float = 1e-10) -> ReheatingSolution:
    """Integrate the decay equations using comoving variables X = x e^{3(1+w)N}, Y = y e^{4N} (well scaled)."""
    if gamma_over_Hend <= 0:
        raise ValueError("gamma_over_Hend must be positive")
    if gamma_over_Hend >= 1.0:
        raise ValueError("require Gamma < H_end for a reheating epoch to exist")
    if N_max is None:
        N_max = 2.0 / (3.0 * (1.0 + w)) * math.log(1.0 / gamma_over_Hend) + 10.0
    a3 = 3.0 * (1.0 + w)

    def rhs(N, u):
        X, Y = u
        Ht = math.sqrt(max(X * math.exp(-a3 * N) + Y * math.exp(-4.0 * N), 1e-300))
        g = gamma_over_Hend / Ht
        return [-g * X, g * X * math.exp((1.0 - 3.0 * w) * N)]

    N_eval = np.linspace(0.0, N_max, 4000)
    sol = integrate.solve_ivp(rhs, [0, N_max], [1.0, 0.0], method="Radau", rtol=rtol, atol=1e-30,
                              t_eval=N_eval)
    if not sol.success:
        raise RuntimeError(sol.message)
    x = sol.y[0] * np.exp(-a3 * sol.t)
    y = sol.y[1] * np.exp(-4.0 * sol.t)
    return ReheatingSolution(sol.t, x, y, gamma_over_Hend, w)


def equality_point(sol: ReheatingSolution) -> tuple[float, float]:
    """(N, rho_total/rho_end) where rho_r = rho_phi (end of inflaton domination)."""
    d = sol.y - sol.x
    i = int(np.argmax(d > 0))
    if i == 0:
        raise RuntimeError("no crossing found")
    f = d[i - 1] / (d[i - 1] - d[i])
    return float(sol.N[i - 1] + f * (sol.N[i] - sol.N[i - 1])), float(sol.total[i - 1] + f * (sol.total[i] - sol.total[i - 1]))
