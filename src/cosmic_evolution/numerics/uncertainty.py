"""Illustrative parameter-uncertainty propagation.

ILLUSTRATIVE ONLY: parameters are sampled as *independent* Gaussians with the 1-sigma widths quoted in the
Planck 2018 VI abstract (H0 = 67.4 +- 0.5, Omega_m = 0.315 +- 0.007).  Real parameter posteriors are correlated and
non-Gaussian; use published chains for quantitative error budgets.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Callable

import numpy as np

from ..background import cosmic_time as ct
from ..background.friedmann import FLRWBackground
from ..config import CosmologyParams


def monte_carlo(fn: Callable[[FLRWBackground], float], base: CosmologyParams, n: int = 300, seed: int = 12345,
                sigma_h: float = 0.005, sigma_om: float = 0.007) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = []
    for h, om in zip(base.h + sigma_h * rng.standard_normal(n), base.omega_m + sigma_om * rng.standard_normal(n)):
        # keep Omega_b h^2 fixed (physical baryon density), flat closure
        p = replace(base, h=float(h), omega_m=float(om))
        out.append(fn(FLRWBackground(p)))
    return np.array(out)


def log_sensitivity(fn: Callable[[FLRWBackground], float], base: CosmologyParams, param: str, rel_step: float = 1e-3) -> float:
    """d ln f / d ln param by central differences."""
    v = getattr(base, param)
    up = fn(FLRWBackground(replace(base, **{param: v * (1 + rel_step)})))
    dn = fn(FLRWBackground(replace(base, **{param: v * (1 - rel_step)})))
    return float((np.log(up) - np.log(dn)) / (np.log(1 + rel_step) - np.log(1 - rel_step)))
