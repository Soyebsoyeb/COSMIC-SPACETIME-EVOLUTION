"""Friedmann background solver for radiation + matter + curvature + dark energy.

Conventions
-----------
* a(t0) = 1; E(a) = H(a)/H0.
* All Omega_i are *present-day* density parameters; closure sum_i Omega_i = 1 gives E(1) = 1.
* Work with  Q(a) = a^4 E^2(a) = Om_r + Om_m a + Om_k a^2 + Om_de a^4 f_de(a)
  because it stays finite as a -> 0 (radiation era) and all key integrals take the form
  integral of  a^p / sqrt(Q) da.
* Equation of state of curvature is w_k = -1/3 (it contributes no acceleration term).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CosmologyParams
from .dark_energy import de_density_factor, w_of_a


@dataclass(frozen=True)
class FLRWBackground:
    """Homogeneous and isotropic background model (Lambda or CPL dark energy)."""

    params: CosmologyParams

    # ---- basic densities -------------------------------------------------------------------
    @property
    def om_r(self) -> float:
        return self.params.omega_r

    @property
    def om_m(self) -> float:
        return self.params.omega_m

    @property
    def om_k(self) -> float:
        return self.params.omega_k_eff

    @property
    def om_de(self) -> float:
        return self.params.omega_de_eff

    @property
    def H0(self) -> float:
        return self.params.H0

    def Q(self, a):
        """a^4 E^2(a)."""
        a = np.asarray(a, dtype=float)
        with np.errstate(over="ignore"):
            return (self.om_r + self.om_m * a + self.om_k * a**2
                    + self.om_de * a**4 * de_density_factor(a, self.params.w0, self.params.wa))

    def E2(self, a):
        a = np.asarray(a, dtype=float)
        return self.Q(a) / a**4

    def E(self, a):
        e2 = self.E2(a)
        if np.any(e2 < 0):
            raise ValueError("E^2(a) < 0: the model has no real expansion rate at this scale factor")
        return np.sqrt(e2)

    def H(self, a):
        """Hubble rate in km/s/Mpc."""
        return self.H0 * self.E(a)

    def E_of_z(self, z):
        return self.E(1.0 / (1.0 + np.asarray(z, dtype=float)))

    def H_of_z(self, z):
        return self.H0 * self.E_of_z(z)

    # ---- density fractions and deceleration -----------------------------------------------
    def density_terms(self, a):
        """Return dict of Omega_i(a) E^2(a) contributions, i.e. rho_i(a)/rho_crit,0."""
        a = np.asarray(a, dtype=float)
        return {
            "radiation": self.om_r * a**-4,
            "matter": self.om_m * a**-3,
            "curvature": self.om_k * a**-2,
            "dark_energy": self.om_de * de_density_factor(a, self.params.w0, self.params.wa),
        }

    def omega_fractions(self, a):
        """Density parameters Omega_i(a) = rho_i(a)/rho_crit(a); they sum to one."""
        terms = self.density_terms(a)
        e2 = sum(terms.values())
        return {k: v / e2 for k, v in terms.items()}

    def deceleration(self, a):
        """q = -a a_ddot / a_dot^2 = (1/2) sum_i Omega_i(a) (1 + 3 w_i(a)); curvature drops out."""
        f = self.omega_fractions(a)
        w_de = w_of_a(a, self.params.w0, self.params.wa)
        return 0.5 * (2.0 * f["radiation"] + f["matter"] + (1.0 + 3.0 * w_de) * f["dark_energy"])

    # ---- special epochs --------------------------------------------------------------------
    def a_equality_mr(self) -> float:
        """Matter-radiation equality: Om_r a^-4 = Om_m a^-3 (exact for this two-fluid statement)."""
        if self.om_m <= 0:
            raise ValueError("no matter component")
        return self.om_r / self.om_m
