"""Linear growth of matter perturbations.

Equation (sub-horizon, GR, scale-independent, pressureless matter perturbing in a smooth background):

    D'' + (1 - q) D' - (3/2) Omega_m(a) D = 0,     ' = d/d ln a.

Initial conditions deep in the radiation era use the Meszaros solution D = 1 + (3/2) a/a_eq (for matter
perturbations in a radiation background, ignoring Lambda and curvature at that epoch).  D is then normalised so that
D_+ -> a in the matter era.
ACCURACY: the initial conditions neglect curvature at a_i, so D_+ carries a relative normalisation error of order
a_i |Omega_k|/Omega_m (2e-4 for a_i = 1e-4 in a strongly curved model; zero for flat models).
LIMITATIONS: scale independent; no massive-neutrino free streaming; no baryon-photon coupling; dark energy is
assumed smooth (valid for Lambda and for fluids with sound speed c_s = c on sub-horizon scales); no modified gravity.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import integrate

from ..background.friedmann import FLRWBackground


@dataclass
class GrowthSolution:
    a: np.ndarray
    D_plus: np.ndarray     # normalised -> a in matter era
    f: np.ndarray          # d ln D / d ln a

    def D(self, a):
        """D_+(a) by cubic-spline interpolation in ln a (error ~1e-9 on the default grid)."""
        from scipy.interpolate import CubicSpline
        return CubicSpline(np.log(self.a), self.D_plus)(np.log(np.asarray(a, dtype=float)))

    def D_today_normalised(self) -> np.ndarray:
        """D normalised to unity at a = 1."""
        return self.D_plus / np.interp(1.0, self.a, self.D_plus)


def solve_growth(bg: FLRWBackground, a_i: float = 1e-4, a_f: float = 1.0, n: int = 800,
                 rtol: float = 1e-10, atol: float = 1e-14) -> GrowthSolution:
    if not (0 < a_i < a_f):
        raise ValueError("require 0 < a_i < a_f")
    a_eq = bg.a_equality_mr()
    y_i = a_i / a_eq
    D0, dD0 = 1.0 + 1.5 * y_i, 1.5 * y_i

    def rhs(x, u):
        a = np.exp(x)
        q = float(bg.deceleration(a))
        om = float(bg.omega_fractions(a)["matter"])
        return [u[1], -(1.0 - q) * u[1] + 1.5 * om * u[0]]

    x_eval = np.linspace(np.log(a_i), np.log(a_f), n)
    sol = integrate.solve_ivp(rhs, [x_eval[0], x_eval[-1]], [D0, dD0], method="DOP853", rtol=rtol, atol=atol,
                              t_eval=x_eval)
    if not sol.success:
        raise RuntimeError(sol.message)
    norm = a_eq / 1.5
    D = sol.y[0] * norm
    f = sol.y[1] / sol.y[0]
    return GrowthSolution(np.exp(sol.t), D, f)


def heath_integral_growth(bg: FLRWBackground, a: float) -> float:
    """Closed-form-by-quadrature growth for matter + Lambda + curvature (NO radiation): D = (5/2) Om E(a) int_0^a da'/(a' E)^3.

    Valid only when radiation is negligible and dark energy is a cosmological constant (Heath 1977)."""
    f = lambda x: 1.0 / (x * float(bg.E(x))) ** 3
    val, _ = integrate.quad(f, 0.0, a, epsabs=1e-14, epsrel=1e-12, limit=400)
    return 2.5 * bg.om_m * float(bg.E(a)) * val
