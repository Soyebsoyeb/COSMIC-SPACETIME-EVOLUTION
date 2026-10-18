"""Cosmological parameter container with validation and YAML loading.

Fiducial values follow the base-LCDM numbers quoted in the abstract of Planck 2018 results VI
(H0 = 67.4, Omega_m = 0.315, n_s = 0.965, Omega_b h^2 = 0.0224, A_s from ln(10^10 A_s) ~ 3.04).
They are adopted as rounded fiducial inputs, not re-derived here.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import yaml

from . import constants as K


@dataclass(frozen=True)
class CosmologyParams:
    """Parameters of a homogeneous FLRW model with radiation, matter, curvature and dark energy.

    Omega_k is *derived* (closure: sum of density parameters = 1) unless ``flat`` is False and
    ``omega_k`` is given explicitly together with a dark-energy density that closes the budget.
    """

    h: float = 0.674
    omega_m: float = 0.315
    omega_b_h2: float = 0.0224
    omega_de: float | None = None      # if None: fixed by closure with omega_k
    omega_k: float | None = None       # if None and omega_de None -> flat (0)
    T_cmb: float = K.T_CMB0
    N_eff: float = K.N_EFF_STANDARD
    w0: float = -1.0
    wa: float = 0.0
    n_s: float = 0.965
    A_s: float = 2.1e-9
    k_pivot: float = 0.05              # 1/Mpc
    sigma8_planck: float = 0.811       # quoted late-time derived value (Planck 2018 VI abstract)
    Y_p: float = 0.245                 # primordial helium mass fraction INPUT (not computed here)
    name: str = "fiducial"
    # derived (filled in __post_init__ through object.__setattr__)
    omega_r: float = field(init=False, default=0.0)
    omega_de_eff: float = field(init=False, default=0.0)
    omega_k_eff: float = field(init=False, default=0.0)

    def __post_init__(self) -> None:
        if not (0.2 < self.h < 1.5):
            raise ValueError(f"h={self.h} outside plausible range (0.2, 1.5)")
        if self.omega_m < 0:
            raise ValueError("omega_m must be non-negative")
        if self.omega_b_h2 < 0 or self.omega_b_h2 / self.h**2 > self.omega_m + 1e-12:
            raise ValueError("require 0 <= Omega_b <= Omega_m")
        if self.T_cmb <= 0 or self.N_eff < 0:
            raise ValueError("T_cmb must be positive and N_eff non-negative")
        # Radiation: photons + massless neutrinos (documented approximation: massive neutrinos
        # are treated as radiation at all times).
        rho_gamma_h2 = _omega_gamma_h2(self.T_cmb)
        om_r = rho_gamma_h2 * (1.0 + 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0) * self.N_eff) / self.h**2
        object.__setattr__(self, "omega_r", om_r)
        if self.omega_de is None and self.omega_k is None:
            ok, ode = 0.0, 1.0 - self.omega_m - om_r
        elif self.omega_de is None:
            ok = float(self.omega_k)  # type: ignore[arg-type]
            ode = 1.0 - self.omega_m - om_r - ok
        elif self.omega_k is None:
            ode = float(self.omega_de)
            ok = 1.0 - self.omega_m - om_r - ode
        else:
            ode, ok = float(self.omega_de), float(self.omega_k)
            if abs(self.omega_m + om_r + ode + ok - 1.0) > 1e-9:
                raise ValueError("density parameters do not sum to unity")
        object.__setattr__(self, "omega_de_eff", ode)
        object.__setattr__(self, "omega_k_eff", ok)

    @property
    def omega_gamma(self) -> float:
        """Photon density parameter."""
        return _omega_gamma_h2(self.T_cmb) / self.h**2

    @property
    def H0(self) -> float:
        """Hubble constant in km/s/Mpc."""
        return 100.0 * self.h

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_yaml(cls, path: str | Path, section: str = "cosmology") -> "CosmologyParams":
        with open(path, "r", encoding="utf-8") as fh:
            cfg = yaml.safe_load(fh)
        data = dict(cfg[section])
        allowed = {f for f in cls.__dataclass_fields__ if cls.__dataclass_fields__[f].init}
        unknown = set(data) - allowed
        if unknown:
            raise ValueError(f"unknown cosmology keys: {sorted(unknown)}")
        return cls(**data)


def _omega_gamma_h2(T_cmb: float) -> float:
    """Photon density parameter times h^2 from the Stefan-Boltzmann law."""
    rho_gamma = 4.0 * K.sigma_SB * T_cmb**4 / K.c**3  # kg/m^3 (u/c^2 with u = 4 sigma T^4 / c)
    rho_crit_h2 = 3.0 * (K.H0_100) ** 2 / (8.0 * math.pi * K.G)
    return rho_gamma / rho_crit_h2


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)
