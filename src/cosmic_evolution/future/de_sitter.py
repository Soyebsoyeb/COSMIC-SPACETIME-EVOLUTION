"""Asymptotic de Sitter expansion for a positive cosmological constant.

H_Lambda = H0 sqrt(Omega_Lambda).  In the pure de Sitter limit a(t) ~ exp(H_Lambda t), the comoving event
horizon is c/(a H_Lambda), the proper event horizon is c/H_Lambda, and a static observer sees a thermal horizon
temperature T = hbar H_Lambda / (2 pi k_B) (Gibbons-Hawking 1977, a semiclassical result for the pure de Sitter
spacetime; its relevance to the real universe is model dependent).
"""
from __future__ import annotations

import math

from .. import constants as K
from ..background.friedmann import FLRWBackground
from ..units import H0_to_inv_s


def hubble_lambda_inv_s(bg: FLRWBackground) -> float:
    if bg.om_de <= 0 or bg.params.w0 != -1.0 or bg.params.wa != 0.0:
        raise ValueError("de Sitter asymptote requires a positive cosmological constant")
    return H0_to_inv_s(bg.H0) * math.sqrt(bg.om_de)


def proper_event_horizon_gly(bg: FLRWBackground) -> float:
    """c/H_Lambda in billion light-years."""
    return K.c / hubble_lambda_inv_s(bg) / (K.c * K.yr) / 1e9


def gibbons_hawking_temperature_K(bg: FLRWBackground) -> float:
    return K.hbar * hubble_lambda_inv_s(bg) / (2.0 * math.pi * K.k_B)


def e_folding_time_gyr(bg: FLRWBackground) -> float:
    return 1.0 / hubble_lambda_inv_s(bg) / K.Gyr
