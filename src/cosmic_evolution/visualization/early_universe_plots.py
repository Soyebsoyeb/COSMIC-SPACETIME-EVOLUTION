"""Figures for the early universe: inflation, reheating, thermal history, BBN relations, recombination."""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np

from .. import constants as K
from .. import units
from ..background.friedmann import FLRWBackground
from ..config import CosmologyParams
from ..early_universe import bbn_interface as bbn
from ..early_universe import inflation as inf
from ..early_universe import recombination as rc
from ..early_universe import reheating as rh
from ..early_universe import thermal_history as th
from ..early_universe.entropy import baryon_to_photon_ratio
from ..perturbations import primordial_spectrum as ps
from .publication_style import (AQUA, BLUE, FigureResult, INK2, LINESTYLES, MUTED, ORANGE, RED, VIOLET, apply_style,
                                figure)

apply_style()
_cache: dict = {}


def _sol(model):
    key = ("sol", model.name)
    if key not in _cache:
        _cache[key] = inf.solve_background(model)
    return _cache[key]


def _models():
    return [(inf.StarobinskyModel(), BLUE, "-"), (inf.QuadraticModel(), ORANGE, "--")]


@figure("fig09_inflation_background")
def fig09() -> FigureResult:
    fig, axs = plt.subplots(1, 3, figsize=(13.4, 4.4))
    data = {}
    for model, c, ls in _models():
        sol = _sol(model)
        info = inf.pivot_efolds(model, 2.1e-9, 0.05, sol=sol)
        nb = sol.N_before_end()
        phi = np.linspace(0.05, sol.phi[0], 400)
        axs[0].plot(phi, model.v(phi) / model.v(sol.phi[0]), color=c, ls=ls, label=model.name)
        axs[0].plot([sol.phi_end], [float(model.v(sol.phi_end)) / float(model.v(sol.phi[0]))], "o", color=c, ms=6)
        ph_star = sol.phi_at(info["N_star"])
        axs[0].plot([ph_star], [float(model.v(ph_star)) / float(model.v(sol.phi[0]))], "s", color=c, ms=6, mfc="white")
        axs[1].plot(nb, sol.phi, color=c, ls=ls, label=model.name)
        axs[2].semilogy(nb, sol.eps_H, color=c, ls=ls, label=f"{model.name}: $\\epsilon_H$")
        mm = nb < (nb.max() - 8.0)   # hide the initial-condition transient of eta_H
        axs[2].semilogy(nb[mm], np.abs(sol.eta_H)[mm], color=c, ls=ls, lw=1.0, alpha=0.6)
        data[f"N_before_end_{model.name}"] = nb
        data[f"phi_{model.name}"] = sol.phi
        data[f"epsH_{model.name}"] = sol.eps_H
    axs[0].set_xlabel("$\\phi$ [$M_P$]")
    axs[0].set_ylabel("$V(\\phi)/V(\\phi_{start})$")
    axs[0].set_title("Potentials (circle: end; square: pivot exit)")
    axs[0].legend()
    axs[1].set_xlim(100, 0)
    axs[1].set_xlabel("e-folds before end of inflation")
    axs[1].set_ylabel("$\\phi$ [$M_P$]")
    axs[1].set_title("Field evolution (exact background)")
    axs[1].legend()
    axs[2].set_xlim(100, 0)
    axs[2].axhline(1.0, color=MUTED, lw=0.8)
    axs[2].text(98, 1.15, "end of inflation: $\\epsilon_H = 1$", fontsize=8.5, color=INK2)
    axs[2].set_xlabel("e-folds before end of inflation")
    axs[2].set_ylabel("$\\epsilon_H = -\\dot H/H^2$ (thick); $|\\eta_H|$ (thin)")
    axs[2].set_title("Hubble slow-roll parameters")
    axs[2].set_ylim(1e-5, 3)
    cap = ("Exact numerical background for two representative single-field models (M_P = 1): quadratic V = m^2 phi^2/2 and Starobinsky-type "
           "V = (3/4) M^2 (1 - exp(-sqrt(2/3) phi))^2. The e-fold equation phi'' = -(3 - eps)(phi' + V_phi/V) is integrated from the slow-roll "
           "attractor until eps_H = 1. These potentials are model assumptions, not established physics, and initial conditions are not derived.")
    return FigureResult("fig09_inflation_background", "Inflationary background evolution", cap, fig, data, "model")


def comoving_hubble_composite(model, T_reh_GeV=None, w_reh=0.0, A_s=2.1e-9, k_p=0.05):
    sol = _sol(model)
    info = inf.pivot_efolds(model, A_s, k_p, w_reh=w_reh, T_reh_GeV=T_reh_GeV, sol=sol)
    bg = FLRWBackground(CosmologyParams())
    conv = inf.M_P_GEV * inf.GEV_TO_INV_MPC
    V0 = info["V0"]
    ln_a_end = info["ln_a_end"]
    # inflation
    lna_inf = sol.N - sol.N_total + ln_a_end
    H_inf = np.sqrt(V0 * model.v(sol.phi) / (3 - sol.eps_H)) * conv                     # 1/Mpc
    rows = [(lna_inf, 1.0 / (np.exp(lna_inf) * H_inf))]
    H_end = info["H_end"] * conv
    ln_a_reh = ln_a_end + info["N_reh"]
    if info["N_reh"] > 0:
        lna = np.linspace(ln_a_end, ln_a_reh, 200)
        H = H_end * np.exp(-1.5 * (1 + w_reh) * (lna - ln_a_end))
        rows.append((lna, 1.0 / (np.exp(lna) * H)))
    # radiation era with g_*(T), entropy conservation, then the full background
    lna = np.linspace(ln_a_reh, 0.0, 1500)
    a = np.exp(lna)
    H = np.empty_like(a)
    H0s = bg.H0 / (units.hubble_distance_Mpc(bg.H0) * bg.H0)
    for i, aa in enumerate(a):
        T = K.T_CMB0 / aa
        for _ in range(4):
            gs = th.g_star_s(float(units.K_to_GeV(T)))
            T = K.T_CMB0 / aa * (inf.G_S0 / gs) ** (1 / 3)
        g = th.g_star(float(units.K_to_GeV(T)))
        corr = (g / th.g_star_late()) * (inf.G_S0 / gs) ** (4 / 3)
        E2 = bg.om_r * corr * aa ** -4 + bg.om_m * aa ** -3 + bg.om_de * 1.0
        H[i] = H0s * math.sqrt(E2)
    rows.append((lna, 1.0 / (a * H)))
    return rows, info


@figure("fig10_comoving_hubble_radius")
def fig10() -> FigureResult:
    model = inf.StarobinskyModel()
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    data = {}
    for (T, lab, c, ls) in ((None, "instantaneous reheating", BLUE, "-"), (1e9, "T_reh = 1e9 GeV, matter-like reheating", ORANGE, "--")):
        rows, info = comoving_hubble_composite(model, T)
        x = np.concatenate([r[0] for r in rows])
        y = np.concatenate([r[1] for r in rows])
        ax.semilogy(x, y, color=c, ls=ls, label=f"{lab}  (N* = {info['N_star']:.1f})")
        data[f"lna_{'inst' if T is None else 'Treh1e9'}"] = x
        data[f"inv_aH_Mpc_{'inst' if T is None else 'Treh1e9'}"] = y
    kp = 0.05
    ax.axhline(1 / kp, color=MUTED, lw=0.9, ls="-.")
    ax.text(-122, 1 / kp * 2.5, "CMB pivot scale 1/k = 20 Mpc (comoving)", fontsize=8.5, color=INK2)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.text(-1.2, 3e-24, "today", fontsize=8.5, color=INK2, rotation=90, ha="right")
    ax.set_xlabel("ln a   [a = 1 today]")
    ax.set_ylabel("comoving Hubble radius (aH)$^{-1}$ [Mpc]")
    ax.set_title("Comoving Hubble radius: shrinks during inflation, grows afterwards")
    ax.legend(loc="upper right", fontsize=8.5)
    ax.set_xlim(-124, 3)
    ax.set_ylim(1e-26, 1e36)
    cap = ("Comoving Hubble radius (aH)^-1 built from the exact Starobinsky-type background during inflation (amplitude fixed by A_s = 2.1e-9 "
           "at the pivot), a reheating epoch of equation of state w = 0 with the stated T_reh (or instantaneous reheating), radiation domination "
           "with the SM g_*(T) step approximation and entropy conservation, then the matter and Lambda eras of the fiducial background. Scales exit "
           "the horizon during inflation (curve decreasing) and re-enter later; the horizon-problem solution corresponds to the 1/k line lying above "
           "the curve at early times. T_reh is an input, not a measurement.")
    return FigureResult("fig10_comoving_hubble_radius", "Comoving Hubble radius through cosmic history", cap, fig, data, "model")


@figure("fig13_inflation_ns_r")
def fig13() -> FigureResult:
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    data = {}
    for model, c, ls in _models():
        sol = _sol(model)
        N = np.linspace(40, 70, 61)
        ns, r = [], []
        for n in N:
            o = inf.slow_roll_observables(model, sol.phi_at(n))
            ns.append(o["n_s"])
            r.append(o["r"])
        ax.plot(ns, r, color=c, ls=ls, label=f"{model.name} (N = 40 to 70)")
        n50, n60 = inf.slow_roll_observables(model, sol.phi_at(50)), inf.slow_roll_observables(model, sol.phi_at(60))
        ax.plot([n50["n_s"], n60["n_s"]], [n50["r"], n60["r"]], "o", color=c, ms=6)
        info = inf.pivot_efolds(model, 2.1e-9, 0.05, sol=sol)
        o = inf.slow_roll_observables(model, sol.phi_at(info["N_star"]))
        ax.plot([o["n_s"]], [o["r"]], "D", color=c, ms=7, mfc="white")
        data[f"ns_{model.name}"] = np.array(ns)
        data[f"r_{model.name}"] = np.array(r)
    ax.axvspan(0.965 - 0.004, 0.965 + 0.004, color=AQUA, alpha=0.25, lw=0)
    ax.text(0.9655, 0.0007, "Planck 2018 $n_s$ = 0.965 +- 0.004 (68%)", rotation=90, fontsize=8, color=INK2, va="bottom")
    ax.axhline(0.06, color=RED, ls="-.", lw=1.2)
    ax.text(0.940, 0.063, "r < 0.06 (95%, Planck+BICEP2/Keck 2015)", color=INK2, fontsize=8.5)
    ax.set_yscale("log")
    ax.set_ylim(5e-4, 0.5)
    ax.set_xlim(0.93, 0.99)
    ax.set_xlabel("scalar spectral index $n_s$")
    ax.set_ylabel("tensor-to-scalar ratio r")
    ax.set_title("Slow-roll predictions of two potentials")
    ax.legend(loc="lower left", fontsize=8.5)
    cap = ("Leading-order potential slow-roll predictions n_s = 1 - 6 eps_V + 2 eta_V and r = 16 eps_V evaluated N e-folds before the end of the numerically "
           "integrated background. Filled circles: N = 50 and 60; open diamonds: the self-consistent N* for instantaneous reheating. Only the "
           "published constraint values (n_s band, r bound) are drawn; no likelihood contours are reproduced. The quadratic model lies above the bound, "
           "the Starobinsky-type model below it, subject to the model assumptions.")
    return FigureResult("fig13_inflation_ns_r", "Inflationary observables in the n_s-r plane", cap, fig, data, "model")


@figure("fig14_reheating")
def fig14() -> FigureResult:
    fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.4))
    data = {}
    for g, c, ls in ((1e-3, BLUE, "-"), (1e-6, ORANGE, "--"), (1e-9, AQUA, "-.")):
        s = rh.solve_perturbative_reheating(g)
        axs[0].semilogy(s.N, s.x, color=c, ls=ls, label=f"$\\rho_\\phi$, $\\Gamma/H_{{end}}$ = {g:.0e}")
        axs[0].semilogy(s.N, s.y, color=c, ls=ls, lw=0.9, alpha=0.65)
        data[f"N_{g:.0e}"] = s.N
        data[f"rho_phi_{g:.0e}"] = s.x
        data[f"rho_r_{g:.0e}"] = s.y
    axs[0].set_ylim(1e-30, 3)
    axs[0].set_xlabel("e-folds since end of inflation N")
    axs[0].set_ylabel("$\\rho / \\rho_{end}$   (thick: inflaton, thin: radiation)")
    axs[0].set_title("Perturbative decay: energy transfer")
    axs[0].legend(fontsize=8)
    gs = np.logspace(-10, -2, 17)
    ratios = []
    for g in gs:
        s = rh.solve_perturbative_reheating(g)
        Neq, rho_eq = rh.equality_point(s)
        Teq = rh.temperature_from_density(rho_eq * 3.0 / 2.0)           # Hend = 1: rho_end = 3
        ratios.append(Teq / rh.reheating_temperature_estimate(g))
    ratios = np.array(ratios)
    axs[1].semilogx(gs, ratios, "o-", color=BLUE)
    axs[1].set_ylim(0, 1)
    axs[1].set_xlabel("$\\Gamma / H_{end}$")
    axs[1].set_ylabel("$T(\\rho_r = \\rho_\\phi) / T_{reh}$ with $T_{reh}$ defined by $\\Gamma = H$")
    axs[1].set_title("Scaling check of the T$_{reh}$ definition")
    cap = ("Perturbative reheating toy model (constant decay rate, matter-like inflaton oscillations) with total density obeying the Friedmann constraint. "
           "Left: inflaton and radiation densities for three decay rates; reheating (radiation = inflaton) occurs near H ~ Gamma. Right: the temperature at "
           "radiation-inflaton equality relative to the conventional T_reh = (90/(pi^2 g))^(1/4) sqrt(Gamma M_P) is independent of Gamma (about 0.64), "
           "confirming the Gamma^(1/2) scaling. The decay rate is a free parameter; no reheating temperature is inferred from data.")
    return FigureResult("fig14_reheating", "Perturbative reheating", cap, fig, {**data, "gamma_over_Hend": gs, "Teq_over_Treh": ratios}, "model")


@figure("fig15_relativistic_degrees_of_freedom")
def fig15() -> FigureResult:
    T = np.logspace(-6, 3, 1500)
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    g = th.g_star(T)
    gs = th.g_star_s(T)
    ax.semilogx(T, g, color=BLUE, label="$g_*(T)$ (energy density)")
    ax.semilogx(T, gs, color=ORANGE, ls="--", label="$g_{*s}(T)$ (entropy)")
    for name, Tm in (("e+- annih.", 5.11e-4), ("QCD", 0.17), ("W, Z, H", 100.0), ("top", 173.0)):
        ax.axvline(Tm, color=MUTED, lw=0.7, ls=":")
        ax.text(Tm * 1.08, 6, name, rotation=90, fontsize=8, color=INK2)
    ax.set_xlabel("temperature T [GeV]")
    ax.set_ylabel("effective degrees of freedom")
    ax.set_title("Relativistic degrees of freedom (step approximation)")
    ax.legend(loc="upper left", fontsize=8.5)
    cap = ("Step-function approximation to g_*(T) and g_*s(T) from Standard-Model particle counting with thresholds at T = m. The real electroweak and QCD "
           "transitions are smooth crossovers, so the steps are an approximation (documented in docs/THERMAL_HISTORY.md). Late-time values 3.38 (energy) "
           "and 3.909 (entropy) follow from e+- annihilation heating the photons relative to neutrinos.")
    return FigureResult("fig15_relativistic_degrees_of_freedom", "Relativistic degrees of freedom", cap, fig,
                        {"T_GeV": T, "g_star": g, "g_star_s": gs}, "approximate")


@figure("fig16_recombination")
def fig16() -> FigureResult:
    bg = FLRWBackground(CosmologyParams())
    r = rc.solve_recombination(bg)
    m = (r.z > 200) & (r.z < 1800)
    fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.4))
    axs[0].semilogy(r.z[m], r.x_e[m], color=BLUE, label="Peebles three-level (hydrogen)")
    axs[0].semilogy(r.z[m], rc.saha_x(r.z[m], bg), color=ORANGE, ls="--", label="Saha equilibrium")
    axs[0].axvline(r.z_star, color=MUTED, lw=0.8)
    axs[0].set_xlabel("redshift z")
    axs[0].set_ylabel("free-electron fraction $x_e$ (hydrogen only)")
    axs[0].set_title("Ionisation history")
    axs[0].set_ylim(1e-4, 1.5)
    axs[0].legend(fontsize=8.5)
    mz = (r.z > 700) & (r.z < 1500)
    axs[1].plot(r.z[mz], r.visibility[mz] * 1e3, color=BLUE)
    axs[1].axvline(r.z_star, color=ORANGE, lw=1.0, ls="--")
    axs[1].annotate(f"peak z* = {r.z_star:.0f}\n(tau = 1 at z = {r.z_tau1:.0f})\nFWHM = {r.fwhm:.0f}", (r.z_star, r.visibility.max() * 1e3),
                    xytext=(-110, -50), textcoords="offset points", fontsize=8.5, arrowprops=dict(arrowstyle="-", color=INK2))
    axs[1].set_xlabel("redshift z")
    axs[1].set_ylabel("visibility $g(z) = d\\tau/dz\\, e^{-\\tau}$  [$10^{-3}$]")
    axs[1].set_title("Photon last scattering (no reionization)")
    cap = (f"Approximate recombination: Peebles three-level atom with case-B coefficients, Ly-alpha escape and two-photon decay, hydrogen only, "
           f"no reionization, T_m = T_gamma. The visibility function peaks at z = {r.z_star:.0f}; the comoving sound horizon at the peak is "
           f"{r.r_s_mpc:.1f} Mpc and 100 theta* = {100 * r.theta_star:.4f}, to be compared with the 100 theta_MC = 1.0409 +- 0.0003 reported by Planck 2018; "
           f"the ~0.6% offset reflects the approximations (helium, fudge factor, radiative transfer). This is NOT a CMB anisotropy calculation.")
    return FigureResult("fig16_recombination", "Recombination and last scattering", cap, fig,
                        {"z": r.z, "x_e": r.x_e, "tau": r.tau, "visibility": r.visibility}, "approximate")


@figure("fig17_primordial_spectra")
def fig17() -> FigureResult:
    fig, axs = plt.subplots(1, 3, figsize=(14.0, 4.4))
    data = {}
    for model, c, ls in _models():
        S = ps.compute_spectra(model, n_k=17)
        axs[0].loglog(S.k_mpc, S.P_s, color=c, ls=ls, label=f"{model.name}: numerical Mukhanov-Sasaki")
        axs[0].loglog(S.k_mpc, S.P_s_slowroll, color=c, ls=":", lw=1.0, label="leading-order slow roll")
        axs[1].semilogx(S.k_mpc, S.P_s / S.P_s_slowroll, color=c, ls=ls, label=f"scalar, {model.name}")
        axs[1].semilogx(S.k_mpc, S.P_t / S.P_t_slowroll, color=c, ls=ls, lw=0.9, alpha=0.6)
        axs[2].loglog(S.k_mpc, S.P_t, color=c, ls=ls, label=f"{model.name}: r = {S.r_numeric:.4f}, $n_s$ = {S.n_s_numeric:.4f}")
        data[f"k_{model.name}"] = S.k_mpc
        data[f"Ps_{model.name}"] = S.P_s
        data[f"Pt_{model.name}"] = S.P_t
    axs[0].set_xlabel("k [Mpc$^{-1}$]")
    axs[0].set_ylabel("$\\mathcal{P}_\\mathcal{R}(k)$")
    axs[0].set_title("Primordial scalar spectrum")
    axs[0].legend(fontsize=7.5)
    axs[1].set_xlabel("k [Mpc$^{-1}$]")
    axs[1].set_ylabel("numerical / leading-order slow roll")
    axs[1].set_title("Next-to-leading-order correction (thin: tensors)")
    axs[1].legend(fontsize=7.5)
    axs[2].set_xlabel("k [Mpc$^{-1}$]")
    axs[2].set_ylabel("$\\mathcal{P}_t(k)$")
    axs[2].set_title("Primordial tensor spectrum")
    axs[2].legend(fontsize=7.5)
    cap = ("Primordial spectra from direct numerical integration of the Mukhanov-Sasaki (scalar) and tensor mode equations on the exact inflationary background, "
           "with Bunch-Davies initial conditions and the pivot k = 0.05 /Mpc fixed by the self-consistent N* for instantaneous reheating. The leading-order "
           "slow-roll amplitude is accurate to ~2-3% for scalars and <0.01% for tensors in these models; n_s and r quoted are numerical log-derivative and ratio "
           "at the pivot. Amplitude normalisation A_s = 2.1e-9 was imposed with the leading-order formula, so the numerical A_s differs by the same 2-3%.")
    return FigureResult("fig17_primordial_spectra", "Primordial scalar and tensor spectra", cap, fig, data, "model")


@figure("fig23_bbn_relations")
def fig23() -> FigureResult:
    p = CosmologyParams()
    T = np.logspace(-2, 1.5, 300)
    npr = bbn.np_equilibrium_ratio(T)
    x = np.logspace(-3, 0, 200)
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.3))
    axs[0].loglog(T, npr, color=BLUE)
    axs[0].axvline(bbn.Q_NP_MEV, color=MUTED, lw=0.8, ls=":")
    axs[0].text(bbn.Q_NP_MEV * 1.05, 2e-4, "Q = 1.293 MeV", fontsize=8.5, color=INK2)
    axs[0].axhline(1.0, color=MUTED, lw=0.6)
    axs[0].set_xlabel("temperature T [MeV]")
    axs[0].set_ylabel("equilibrium n/p = exp(-Q/T)")
    axs[0].set_title("Equilibrium neutron-to-proton ratio")
    axs[0].set_ylim(1e-4, 2)
    axs[1].semilogx(x, [bbn.helium_mass_fraction_from_np(v) for v in x], color=ORANGE)
    axs[1].set_xlabel("n/p at the onset of nucleosynthesis (INPUT)")
    axs[1].set_ylabel("$Y_p = 2x/(1+x)$ if all neutrons enter $^4$He")
    axs[1].set_title("Identity, not a prediction")
    axs[1].text(0.003, 0.80, f"computed baryon-to-photon ratio\n$\\eta$ = {baryon_to_photon_ratio(p):.3e}\nfor $\\Omega_b h^2$ = {p.omega_b_h2}", fontsize=9, color=INK2)
    cap = ("Textbook relations used in the BBN interface. Left: equilibrium n/p, valid only while weak interactions keep neutrons and protons in "
           "equilibrium; weak freeze-out occurs at roughly MeV temperatures and is NOT computed here. Right: the mass fraction of helium-4 if every neutron "
           "were bound into 4He, as a function of an externally supplied n/p. No primordial abundance (D/H, Y_p, 7Li/H) is predicted by this repository; "
           "a nuclear reaction network (PRIMAT, AlterBBN, PArthENoPE) is required.")
    return FigureResult("fig23_bbn_relations", "BBN interface relations", cap, fig,
                        {"T_MeV": T, "np_equilibrium": npr}, "established")
