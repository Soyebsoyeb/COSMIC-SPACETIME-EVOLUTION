"""Big Bang nucleosynthesis: inputs, textbook equilibrium relations and an explicit external-solver interface.

SCOPE (read this): this module does NOT compute primordial abundances.  A quantitative prediction of
D/H, 3He/H, Y_p and 7Li/H requires a nuclear reaction network with weak-rate corrections (e.g. PRIMAT,
AlterBBN, PArthENoPE).  Implemented here are
  * the baryon-to-photon ratio eta computed from the background parameters,
  * the equilibrium neutron-to-proton ratio n/p = exp(-Q/T) (valid only while weak interactions are fast),
  * the identity Y_p = 2 x/(1+x) with x = n/p at the onset of nucleosynthesis (all surviving neutrons in 4He),
    where x is a USER-SUPPLIED input - it is not predicted here,
  * serialisation of the inputs for external solvers, and a backend hook that raises if no backend exists.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ..config import CosmologyParams
from .entropy import baryon_to_photon_ratio

Q_NP_MEV = 1.29333  # neutron-proton mass difference in MeV (PDG value, rounded)


@dataclass(frozen=True)
class BBNInputs:
    eta: float
    N_eff: float
    omega_b_h2: float
    Y_p_assumed_in_recombination: float


class BBNBackendUnavailable(RuntimeError):
    """Raised when an external nuclear-network backend is requested but not installed."""


def inputs_from_params(params: CosmologyParams) -> BBNInputs:
    return BBNInputs(eta=baryon_to_photon_ratio(params), N_eff=params.N_eff,
                     omega_b_h2=params.omega_b_h2, Y_p_assumed_in_recombination=params.Y_p)


def np_equilibrium_ratio(T_MeV):
    """Equilibrium n/p = exp(-Q/T)."""
    T = np.asarray(T_MeV, dtype=float)
    if np.any(T <= 0):
        raise ValueError("temperature must be positive")
    return np.exp(-Q_NP_MEV / T)


def helium_mass_fraction_from_np(n_over_p: float) -> float:
    """Y_p = 2 (n/p)/(1 + n/p): mass fraction if every neutron ends in 4He (an identity, not a prediction)."""
    if n_over_p < 0:
        raise ValueError("n/p must be non-negative")
    return 2.0 * n_over_p / (1.0 + n_over_p)


def export_inputs(params: CosmologyParams, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"bbn_inputs": asdict(inputs_from_params(params)),
                                "note": "Inputs for an external BBN network; no abundances computed."},
                               indent=2))
    return path


def run_backend(name: str, params: CosmologyParams):
    """Hook for external BBN solvers.  No backend is bundled; always raises BBNBackendUnavailable."""
    raise BBNBackendUnavailable(
        f"BBN backend '{name}' is not installed or wrapped. Use export_inputs() and run an external network.")
