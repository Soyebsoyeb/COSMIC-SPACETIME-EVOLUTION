"""Single-field slow-roll inflation: exact background integration and slow-roll observables.

Units: reduced Planck mass M_P = 1 (8 pi G = 1), hbar = c = 1.  V(phi) = V0 * v(phi) with a dimensionless
shape v.  Using e-folds N = ln a as the clock the field equation becomes

    phi'' = -(3 - eps_H) (phi' + V_phi / V),   eps_H = phi'^2 / 2,   H^2 = V / (3 - eps_H),

where primes are d/dN.  This is independent of the amplitude V0, which is fixed afterwards by the scalar
amplitude A_s = H^2/(8 pi^2 eps_H) at horizon exit of the pivot scale.

Models (both are single-field, minimally coupled, canonical kinetic term):
  * QuadraticModel:   V = (1/2) m^2 phi^2         (V0 = m^2, v = phi^2/2).  Disfavoured by r limits.
  * StarobinskyModel: V = (3/4) M^2 (1 - exp(-sqrt(2/3) phi))^2   (V0 = 3 M^2/4).
LIMITATIONS: these are two representative potentials.  No claim is made that either is the true inflationary
model; the Planck-scale origin of the potential, initial conditions and eternal-inflation issues are not addressed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import integrate, optimize

from .. import constants as K
from ..units import K_to_GeV


class InflationModel:
    name = "abstract"
    default_phi_i = 1.0

    def v(self, phi):
        raise NotImplementedError

    def dv(self, phi):
        raise NotImplementedError

    def d2v(self, phi):
        raise NotImplementedError

    # --- potential slow-roll parameters (scale invariant under V0) -------------------------
    def eps_V(self, phi):
        return 0.5 * (self.dv(phi) / self.v(phi)) ** 2

    def eta_V(self, phi):
        return self.d2v(phi) / self.v(phi)

    def dlnv(self, phi):
        return self.dv(phi) / self.v(phi)


class QuadraticModel(InflationModel):
    name = "quadratic"
    default_phi_i = 19.0  # ~ 90 e-folds (N ~ phi^2/4)

    def v(self, phi):
        return 0.5 * np.asarray(phi, dtype=float) ** 2

    def dv(self, phi):
        return np.asarray(phi, dtype=float)

    def d2v(self, phi):
        return np.ones_like(np.asarray(phi, dtype=float))


class StarobinskyModel(InflationModel):
    name = "Starobinsky"
    default_phi_i = 6.2  # ~ 100 e-folds (N ~ (3/4) exp(sqrt(2/3) phi))
    _s = math.sqrt(2.0 / 3.0)

    def v(self, phi):
        return (1.0 - np.exp(-self._s * np.asarray(phi, dtype=float))) ** 2

    def dv(self, phi):
        e = np.exp(-self._s * np.asarray(phi, dtype=float))
        return 2.0 * (1.0 - e) * self._s * e

    def d2v(self, phi):
        e = np.exp(-self._s * np.asarray(phi, dtype=float))
        return 2.0 * self._s**2 * (2.0 * e**2 - e)


@dataclass
class InflationSolution:
    model: InflationModel
    N: np.ndarray          # e-folds from the start of the integration
    phi: np.ndarray
    dphi_dN: np.ndarray
    eps_H: np.ndarray
    eta_H: np.ndarray      # d ln eps_H / dN
    N_total: float         # e-folds from start to eps_H = 1
    phi_end: float

    def N_before_end(self):
        return self.N_total - self.N

    def phi_at(self, n_before_end: float) -> float:
        """Field value n_before_end e-folds before the end of inflation (cubic-spline interpolation)."""
        if not (0.0 < n_before_end < self.N_total):
            raise ValueError("n_before_end outside the integrated range")
        from scipy.interpolate import CubicSpline
        cs = CubicSpline(self.N, self.phi)
        return float(cs(self.N_total - n_before_end))


def solve_background(model: InflationModel, phi_i: float | None = None, rtol: float = 1e-11,
                     atol: float = 1e-13, N_max: float = 400.0) -> InflationSolution:
    """Integrate the exact e-fold equations from the slow-roll attractor until eps_H = 1."""
    phi_i = model.default_phi_i if phi_i is None else float(phi_i)
    if phi_i <= 0:
        raise ValueError("phi_i must be positive")
    y0 = [phi_i, -float(model.dlnv(phi_i))]  # slow-roll attractor phi' = -V_phi/V

    def rhs(N, y):
        phi, dphi = y
        eps = 0.5 * dphi**2
        return [dphi, -(3.0 - eps) * (dphi + float(model.dlnv(phi)))]

    def end(N, y):
        return 0.5 * y[1] ** 2 - 1.0
    end.terminal = True
    end.direction = 1

    sol = integrate.solve_ivp(rhs, [0.0, N_max], y0, method="DOP853", rtol=rtol, atol=atol,
                              events=end, dense_output=True, max_step=0.05)
    if sol.t_events[0].size == 0:
        raise RuntimeError("inflation did not end (eps_H < 1 throughout); increase N_max or change phi_i")
    N_tot = float(sol.t_events[0][0])
    N = np.linspace(0.0, N_tot, 6000)
    phi, dphi = sol.sol(N)
    eps = 0.5 * dphi**2
    ddphi = -(3.0 - eps) * (dphi + model.dlnv(phi))
    eta_H = 2.0 * ddphi / dphi
    return InflationSolution(model, N, phi, dphi, eps, eta_H, N_tot, float(sol.sol(N_tot)[0]))


def slow_roll_observables(model: InflationModel, phi_star: float) -> dict[str, float]:
    """Leading-order potential slow-roll results at horizon exit."""
    e = float(model.eps_V(phi_star))
    et = float(model.eta_V(phi_star))
    return {"eps_V": e, "eta_V": et, "n_s": 1.0 - 6.0 * e + 2.0 * et, "r": 16.0 * e}


def hubble_slow_roll_observables(sol: InflationSolution, n_before_end: float) -> dict[str, float]:
    """First-order Hubble-flow results: n_s - 1 = -2 eps_H - eta_H, r = 16 eps_H."""
    from scipy.interpolate import CubicSpline
    n = sol.N_total - n_before_end
    e = float(CubicSpline(sol.N, sol.eps_H)(n))
    et = float(CubicSpline(sol.N, sol.eta_H)(n))
    return {"eps_H": e, "eta_H": et, "n_s": 1.0 - 2.0 * e - et, "r": 16.0 * e}


def scalar_amplitude_normalisation(sol: InflationSolution, n_star: float, A_s: float) -> float:
    """V0 (units M_P^4) such that P_R = H^2/(8 pi^2 eps_H) = A_s at n_star e-folds before the end."""
    from scipy.interpolate import CubicSpline
    n = sol.N_total - n_star
    eps = float(CubicSpline(sol.N, sol.eps_H)(n))
    phi = float(CubicSpline(sol.N, sol.phi)(n))
    v = float(sol.model.v(phi))
    return A_s * 8.0 * math.pi**2 * eps * (3.0 - eps) / v


def hubble_of_N(sol: InflationSolution, V0: float):
    """H(N) in units of M_P along the solution."""
    return np.sqrt(V0 * sol.model.v(sol.phi) / (3.0 - sol.eps_H))


# ---- connecting the end of inflation to today ------------------------------------------------------
M_P_GEV = K.M_P_reduced_GeV
GEV_TO_INV_MPC = 1e9 * K.eV / (K.hbar * K.c) * K.Mpc   # 1 GeV in 1/Mpc  (natural units, hbar=c=1)
T0_GEV = float(K_to_GeV(K.T_CMB0))
G_S0 = 2.0 + 7.0 / 4.0 * (4.0 / 11.0) * 3.0


def pivot_efolds(model: InflationModel, A_s: float, k_pivot_mpc: float, w_reh: float = 1.0 / 3.0,
                 T_reh_GeV: float | None = None, g_reh: float = 106.75, sol: InflationSolution | None = None):
    """Self-consistent number of e-folds N* between horizon exit of k_pivot and the end of inflation.

    Assumptions (model dependent): reheating is an epoch of constant equation of state w_reh lasting from
    rho_end to rho_reh = (pi^2/30) g_reh T_reh^4; afterwards entropy is conserved (g_{*s} from the SM).
    T_reh_GeV=None means instantaneous reheating, rho_reh = rho_end.  The value of T_reh is an INPUT.
    """
    sol = sol or solve_background(model)

    def solve_for(nstar: float):
        V0 = scalar_amplitude_normalisation(sol, nstar, A_s)
        v_end = float(model.v(sol.phi_end))
        H_end = math.sqrt(V0 * v_end / (3.0 - 1.0))
        rho_end = 3.0 * H_end**2
        if T_reh_GeV is None:
            rho_reh = rho_end
        else:
            rho_reh = math.pi**2 / 30.0 * g_reh * (T_reh_GeV / M_P_GEV) ** 4
            if rho_reh > rho_end:
                raise ValueError("T_reh exceeds the instantaneous-reheating temperature")
        T_reh = (30.0 * rho_reh / (math.pi**2 * g_reh)) ** 0.25 * M_P_GEV
        N_reh = math.log(rho_end / rho_reh) / (3.0 * (1.0 + w_reh)) if rho_end > rho_reh else 0.0
        ln_a_reh = math.log((G_S0 / g_reh) ** (1.0 / 3.0) * T0_GEV / T_reh)
        ln_a_end = -N_reh + ln_a_reh
        from scipy.interpolate import CubicSpline
        n = sol.N_total - nstar
        eps = float(CubicSpline(sol.N, sol.eps_H)(n))
        phi = float(CubicSpline(sol.N, sol.phi)(n))
        Hk = math.sqrt(V0 * float(model.v(phi)) / (3.0 - eps))
        ln_k = ln_a_end - nstar + math.log(Hk * M_P_GEV * GEV_TO_INV_MPC)
        return ln_k, dict(V0=V0, H_end=H_end, rho_end=rho_end, T_reh_GeV=T_reh, N_reh=N_reh,
                          ln_a_end=ln_a_end, H_star=Hk)

    f = lambda n: solve_for(n)[0] - math.log(k_pivot_mpc)
    nstar = optimize.brentq(f, 20.0, min(95.0, sol.N_total - 1.0), xtol=1e-10)
    _, info = solve_for(nstar)
    info["N_star"] = nstar
    return info
