"""Primordial scalar and tensor spectra for a given inflation model, computed numerically."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..early_universe import inflation as inf
from .scalar_modes import ModeBackground, evolve_curvature_mode, k_of_exit, slow_roll_amplitude
from .tensor_modes import evolve_tensor_mode, slow_roll_tensor_amplitude


@dataclass
class PrimordialSpectra:
    k_mpc: np.ndarray       # physical comoving wavenumbers in 1/Mpc (pivot mapped to k_pivot)
    P_s: np.ndarray         # numerical Mukhanov-Sasaki
    P_t: np.ndarray
    P_s_slowroll: np.ndarray
    P_t_slowroll: np.ndarray
    n_s_numeric: float      # log-derivative of numerical P_s at the pivot
    r_numeric: float
    N_star: float
    model_name: str
    A_s_numeric: float


def compute_spectra(model: inf.InflationModel, A_s: float = 2.1e-9, k_pivot_mpc: float = 0.05,
                    T_reh_GeV: float | None = None, w_reh: float = 1.0 / 3.0, n_k: int = 25,
                    delta_N: float = 8.0) -> PrimordialSpectra:
    sol = inf.solve_background(model)
    info = inf.pivot_efolds(model, A_s, k_pivot_mpc, w_reh=w_reh, T_reh_GeV=T_reh_GeV, sol=sol)
    V0 = info["V0"]
    N_star = info["N_star"]
    bg = ModeBackground.from_inflation(sol, V0)
    N_piv = sol.N_total - N_star
    N_exits = np.linspace(N_piv - delta_N / 2, N_piv + delta_N / 2, n_k)
    ks = np.array([k_of_exit(bg, n) for n in N_exits])
    Ps = np.array([evolve_curvature_mode(bg, k) for k in ks])
    Pt = np.array([evolve_tensor_mode(bg, k) for k in ks])
    Ps_sr = np.array([slow_roll_amplitude(bg, n) for n in N_exits])
    Pt_sr = np.array([slow_roll_tensor_amplitude(bg, n) for n in N_exits])
    k_pivot_internal = k_of_exit(bg, N_piv)
    k_phys = ks / k_pivot_internal * k_pivot_mpc
    lnk = np.log(k_phys)
    i = int(np.argmin(abs(N_exits - N_piv)))
    dlnP = np.gradient(np.log(Ps), lnk)
    n_s = 1.0 + float(dlnP[i])
    return PrimordialSpectra(k_phys, Ps, Pt, Ps_sr, Pt_sr, n_s, float(Pt[i] / Ps[i]), N_star, model.name,
                             float(np.interp(np.log(k_pivot_mpc), lnk, Ps)))
