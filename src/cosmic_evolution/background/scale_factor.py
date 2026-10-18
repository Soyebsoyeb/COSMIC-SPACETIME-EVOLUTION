"""Redshift <-> scale factor transformations (a(t0) = 1)."""
from __future__ import annotations

import numpy as np


def a_of_z(z):
    z = np.asarray(z, dtype=float)
    if np.any(z <= -1.0):
        raise ValueError("redshift must satisfy z > -1")
    return 1.0 / (1.0 + z)


def z_of_a(a):
    a = np.asarray(a, dtype=float)
    if np.any(a <= 0.0):
        raise ValueError("scale factor must be positive")
    return 1.0 / a - 1.0
