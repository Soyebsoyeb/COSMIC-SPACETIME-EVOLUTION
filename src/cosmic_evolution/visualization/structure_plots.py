"""Figures for linear matter clustering and halo mass functions."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from ..background.friedmann import FLRWBackground
from ..config import CosmologyParams
from ..perturbations.growth_factor import solve_growth
from ..structure import halo_models as hm
from ..structure import power_spectrum as pk
from ..structure import transfer_functions as tf
from .publication_style import AQUA, BLUE, FigureResult, INK2, MUTED, ORANGE, VIOLET, apply_style, figure

apply_style()


@figure("fig18_matter_power_spectrum")
def fig18() -> FigureResult:
    bg = FLRWBackground(CosmologyParams())
    p = bg.params
    g = solve_growth(bg)
    k = np.logspace(-4, 1, 400)                       # 1/Mpc
    kh = k / p.h
    fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.5))
    axs[0].loglog(kh, tf.eh_nowiggle(k, p), color=BLUE, label="Eisenstein-Hu (no wiggle)")
    axs[0].loglog(kh, tf.bbks(k, p), color=ORANGE, ls="--", label="BBKS + Sugiyama shape")
    axs[0].set_xlabel("k [h Mpc$^{-1}$]")
    axs[0].set_ylabel("transfer function T(k)")
    axs[0].set_title("Matter transfer function fits (approximate)")
    axs[0].set_ylim(1e-4, 1.5)
    axs[0].legend(loc="lower left", fontsize=8.5)
    data = {"k_h_per_Mpc": kh}
    for a, c, ls, lab in ((1.0, BLUE, "-", "z = 0"), (0.5, ORANGE, "--", "z = 1"), (0.25, AQUA, "-.", "z = 3")):
        P = pk.matter_power(k, bg, a, growth=g) * p.h**3
        axs[1].loglog(kh, P, color=c, ls=ls, label=lab)
        data[f"P_{lab.replace(' ', '').replace('=', '')}_Mpc_over_h_cubed"] = P
    s8 = pk.sigma8(bg, growth=g)
    s8b = pk.sigma8(bg, tf.bbks, growth=g)
    axs[1].set_xlabel("k [h Mpc$^{-1}$]")
    axs[1].set_ylabel("linear P(k) [(Mpc/h)$^3$]")
    axs[1].set_title(f"Linear matter power spectrum; $\\sigma_8$ = {s8:.3f} (EH), {s8b:.3f} (BBKS)")
    axs[1].set_xlim(1e-4 / p.h, 10 / p.h)
    axs[1].legend(loc="lower left", fontsize=8.5)
    cap = (f"Left: two analytic fits to the matter transfer function (no baryon acoustic oscillations). Right: linear matter power spectrum built from "
           f"A_s = 2.1e-9, n_s = 0.965, the Eisenstein-Hu shape and the linear growth factor. Resulting sigma_8 = {s8:.3f} (EH) / {s8b:.3f} (BBKS) versus "
           f"the Planck-reported 0.811 +- 0.006; the few-percent excess reflects the no-wiggle fit, massless-neutrino treatment and fitting-formula errors. "
           f"Linear regime only; this does not include non-linear evolution or CMB physics.")
    return FigureResult("fig18_matter_power_spectrum", "Transfer function and linear matter power spectrum", cap, fig, data, "approximate")


@figure("fig19_halo_mass_function")
def fig19() -> FigureResult:
    bg = FLRWBackground(CosmologyParams())
    g = solve_growth(bg)
    M = np.logspace(9, 16, 36)
    fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.4))
    data = {"M_Msun": M}
    for a, ls, tag in ((1.0, "-", "z = 0"), (0.5, "--", "z = 1")):
        for kind, c in (("ST", BLUE), ("PS", ORANGE)):
            dn, sig = hm.mass_function(M, bg, kind, a, growth=g)
            axs[0].loglog(M, dn, color=c, ls=ls, label=f"{kind}, {tag}")
            data[f"dn_dlnM_{kind}_{tag.replace(' ', '').replace('=', '')}"] = dn
        data[f"sigma_{tag.replace(' ', '').replace('=', '')}"] = sig
    axs[0].set_ylim(1e-12, 1e0)
    axs[0].set_xlabel("halo mass M [M$_\\odot$]")
    axs[0].set_ylabel("dn/dlnM [Mpc$^{-3}$] (comoving)")
    axs[0].set_title("Halo mass functions from linear theory + collapse ansatz")
    axs[0].legend(fontsize=8.5)
    nu = np.logspace(-1, 0.7, 300)
    axs[1].loglog(nu, hm.f_sheth_tormen(nu), color=BLUE, label="Sheth-Tormen")
    axs[1].loglog(nu, hm.f_press_schechter(nu), color=ORANGE, ls="--", label="Press-Schechter")
    axs[1].set_ylim(1e-4, 1.2)
    axs[1].set_xlabel("$\\nu = \\delta_c/\\sigma(M)$")
    axs[1].set_ylabel("multiplicity function $f(\\nu)$")
    axs[1].set_title("Both normalised to unit mass fraction")
    axs[1].legend(fontsize=8.5)
    cap = ("Halo mass functions from the Press-Schechter and Sheth-Tormen analytic ansatz using sigma(M) from the linear power spectrum, delta_c = 1.686 and the fit "
           "parameters A = 0.3222, a = 0.707, p = 0.3 from the literature. These are phenomenological fits calibrated on N-body simulations, NOT derived from the "
           "homogeneous background; they illustrate the exponential cutoff of massive halos and the redshift evolution only. Accuracy of ~10-20% at best.")
    return FigureResult("fig19_halo_mass_function", "Halo mass functions", cap, fig, data, "model")
