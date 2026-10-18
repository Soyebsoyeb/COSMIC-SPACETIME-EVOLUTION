"""Primordial tensor perturbations (gravitational waves) during inflation.

Each polarisation h_lambda obeys h'' + (3 - eps) h' + (k/aH)^2 h = 0 (' = d/dN) with Bunch-Davies initial
condition |h| = 2/(a M_P sqrt(2k)).  The total dimensionless spectrum is P_t = 2 k^3 |h|^2/(2 pi^2); the
de Sitter limit is P_t = 2 H^2/(pi^2 M_P^2).
"""
from __future__ import annotations

import math

from scipy import integrate

from .scalar_modes import ModeBackground, _find_N_for_x


def evolve_tensor_mode(bg: ModeBackground, k: float, x_start: float = 400.0, x_end: float = 0.01,
                       rtol: float = 1e-10, atol: float = 1e-40) -> float:
    N0 = _find_N_for_x(bg, k, x_start, bg.N_min, bg.N_max)
    N1 = _find_N_for_x(bg, k, x_end, N0, bg.N_max)
    a0 = math.exp(N0)
    amp = 2.0 / (a0 * math.sqrt(2.0 * k))
    x0 = k / (a0 * bg.H(N0))
    h0 = amp + 0j
    dh0 = -(1j * x0 + 1.0) * h0

    def rhs(N, y):
        h = y[0] + 1j * y[1]
        dh = y[2] + 1j * y[3]
        x = k * math.exp(-N) / bg.H(N)
        d2 = -(3.0 - bg.eps(N)) * dh - x * x * h
        return [dh.real, dh.imag, d2.real, d2.imag]

    sol = integrate.solve_ivp(rhs, [N0, N1], [h0.real, h0.imag, dh0.real, dh0.imag], method="DOP853",
                              rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    h = sol.y[0, -1] + 1j * sol.y[1, -1]
    return float(2.0 * k**3 / (2.0 * math.pi**2) * abs(h) ** 2)


def slow_roll_tensor_amplitude(bg: ModeBackground, N_exit: float) -> float:
    return 2.0 * bg.H(N_exit) ** 2 / math.pi**2
