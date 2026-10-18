"""Quadrature helpers with explicit error reporting."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import integrate


@dataclass(frozen=True)
class QuadResult:
    value: float
    abserr: float
    neval: int


def quad_checked(f: Callable[[float], float], a: float, b: float, epsabs: float = 1e-12,
                 epsrel: float = 1e-10, limit: int = 500) -> QuadResult:
    """scipy.quad with a hard failure if the estimated error exceeds the requested tolerance by 1e3."""
    val, err, info = integrate.quad(f, a, b, epsabs=epsabs, epsrel=epsrel, limit=limit, full_output=1)[:3]
    if not np.isfinite(val):
        raise FloatingPointError("integral is not finite")
    if err > 1e3 * max(epsabs, epsrel * abs(val)):
        raise RuntimeError(f"quadrature error estimate {err:.2e} too large")
    return QuadResult(float(val), float(err), int(info["neval"]))


def trapezoid_lna(f_over_a: Callable[[np.ndarray], np.ndarray], a_lo: float, a_hi: float, n: int) -> float:
    """integral_{a_lo}^{a_hi} g(a) da with g = f_over_a, evaluated on a uniform ln a grid (trapezoid rule, O(h^2))."""
    x = np.linspace(np.log(a_lo), np.log(a_hi), n)
    a = np.exp(x)
    return float(np.trapezoid(f_over_a(a) * a, x))


def simpson_lna(f_over_a: Callable[[np.ndarray], np.ndarray], a_lo: float, a_hi: float, n: int) -> float:
    """Same as ``trapezoid_lna`` with Simpson's rule (O(h^4)); n must be odd."""
    if n % 2 == 0:
        n += 1
    x = np.linspace(np.log(a_lo), np.log(a_hi), n)
    a = np.exp(x)
    return float(integrate.simpson(f_over_a(a) * a, x=x))
