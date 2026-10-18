"""Background, distance, growth, future-scenario, convergence and sensitivity figures.

Every function is registered in ``FIGURES`` and returns a ``FigureResult`` with the arrays that were plotted,
so that the CSV datasets and the figures come from the same calculation.
"""
from __future__ import annotations

from dataclasses import replace

import matplotlib.pyplot as plt
import numpy as np

from .. import constants as K
from .. import units
from ..background import cosmic_time as ct
from ..background import expansion_history as eh
from ..config import CosmologyParams
from ..background.friedmann import FLRWBackground
from ..future import dark_energy_scenarios as des
from ..geometry import flrw, horizons
from ..numerics import convergence as cv
from ..numerics import uncertainty as unc
from ..perturbations import growth_factor as gf
from .publication_style import (AQUA, BLUE, FigureResult, GRID, INK2, LINESTYLES, MAGENTA, MUTED, ORANGE, RED,
                                SERIES, VIOLET, apply_style, figure)

apply_style()


def fiducial() -> FLRWBackground:
    return FLRWBackground(CosmologyParams())


def _model(**kw) -> FLRWBackground:
    return FLRWBackground(CosmologyParams(**kw))


def _eds() -> FLRWBackground:
    r = CosmologyParams().omega_r
    return FLRWBackground(CosmologyParams(omega_m=1.0 - r, omega_de=0.0, omega_k=0.0))


def _epoch_scale_factors(bg):
    a_eq = bg.a_equality_mr()
    a_lm = (bg.om_m / bg.om_de) ** (1.0 / 3.0)
    return a_eq, a_lm


# ----------------------------------------------------------------------------------------------------------
@figure("fig01_scale_factor_vs_time")
def fig01() -> FigureResult:
    bg = fiducial()
    a = np.logspace(-8, 3, 400)
    t = np.asarray(ct.cosmic_time_gyr(bg, a))
    a_eq, a_lm = _epoch_scale_factors(bg)
    t_eq = float(ct.cosmic_time_gyr(bg, a_eq))
    t_lm = float(ct.cosmic_time_gyr(bg, a_lm))
    t0 = ct.age_gyr(bg)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.loglog(t, a, color=BLUE, label="Fiducial flat LCDM (background)")
    # asymptotic guides
    tt = np.logspace(-6, -3, 3)
    ax.loglog(tt, 0.0 + (tt / t_eq) ** 0.5 * a_eq * 0.55, color=MUTED, ls=":", lw=1)
    ax.axvspan(t.min(), t_eq, color=ORANGE, alpha=0.08, lw=0)
    ax.axvspan(t_eq, t_lm, color=AQUA, alpha=0.08, lw=0)
    ax.axvspan(t_lm, t.max(), color=VIOLET, alpha=0.08, lw=0)
    ax.plot([t0], [1.0], "o", color=RED, ms=6, zorder=5)
    ax.annotate(f"today: t0 = {t0:.2f} Gyr", (t0, 1.0), xytext=(-110, 22), textcoords="offset points",
                arrowprops=dict(arrowstyle="-", color=INK2), fontsize=9)
    ax.text(3e-9, 5e-6, "radiation\ndominated\n$a\\propto t^{1/2}$", color=INK2, fontsize=9)
    ax.text(1e-3, 3e-2, "matter\ndominated\n$a\\propto t^{2/3}$", color=INK2, fontsize=9)
    ax.text(2e1, 4e-1, "Lambda\ndominated\n$a\\propto e^{H_\\Lambda t}$", color=INK2, fontsize=9)
    ax.set_xlim(1e-9, t.max())
    ax.set_ylim(1e-8, 1e3)
    ax.set_xlabel("cosmic time t [Gyr]")
    ax.set_ylabel("scale factor a(t)   [a(t0) = 1]")
    ax.set_title("Scale factor versus cosmic time (homogeneous background only)")
    ax.legend(loc="lower right")
    cap = ("Scale factor of the fiducial flat LCDM background from the exact quadrature t(a). Shaded bands mark "
           "radiation, matter and cosmological-constant domination (boundaries at a_eq = Om_r/Om_m and a = (Om_m/Om_L)^(1/3)). "
           "The far-future branch assumes a cosmological constant persists forever (conditional on the model).")
    return FigureResult("fig01_scale_factor_vs_time", "Scale factor versus cosmic time", cap, fig,
                        {"t_Gyr": t, "a": a}, "numerical")


@figure("fig02_hubble_parameter_vs_redshift")
def fig02() -> FigureResult:
    bg = fiducial()
    eds = _eds()
    z = np.logspace(-3, 4, 400)
    zl = np.linspace(0, 3, 300)
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.3))
    ax = axs[0]
    ax.loglog(z, bg.H_of_z(z), color=BLUE, label="LCDM (fiducial)")
    ax.loglog(z, eds.H_of_z(z), color=ORANGE, ls="--", label="Einstein-de Sitter (hypothetical)")
    ax.loglog(z, bg.H0 * np.sqrt(bg.om_m) * (1 + z) ** 1.5, color=MUTED, ls=":", label="matter only: $H_0\\sqrt{\\Omega_m}(1+z)^{3/2}$")
    ax.set_xlabel("redshift z")
    ax.set_ylabel("H(z) [km s$^{-1}$ Mpc$^{-1}$]")
    ax.set_title("Expansion rate H(z)")
    ax.legend(loc="upper left", fontsize=8)
    ax = axs[1]
    adot = bg.H_of_z(zl) / (1 + zl)
    ax.plot(zl, adot, color=BLUE, label="LCDM: $\\dot a\\,H_0^{-1}$... in km/s/Mpc")
    ax.plot(zl, eds.H_of_z(zl) / (1 + zl), color=ORANGE, ls="--", label="Einstein-de Sitter")
    zmin = zl[int(np.argmin(adot))]
    ax.plot([zmin], [adot.min()], "o", color=RED, ms=6)
    ax.annotate(f"minimum of $\\dot a$ at z = {zmin:.2f}\n(onset of acceleration)", (zmin, adot.min()),
                xytext=(25, 40), textcoords="offset points", arrowprops=dict(arrowstyle="-", color=INK2), fontsize=8.5)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("$\\dot a$ = H(z)/(1+z)   [km s$^{-1}$ Mpc$^{-1}$]")
    ax.set_title("Cosmic expansion speed: $\\dot a$ decreases, then increases")
    ax.legend(loc="upper right", fontsize=8)
    ax.get_legend().get_texts()[0].set_text("LCDM (fiducial)")
    cap = ("Left: H(z) for the fiducial LCDM model, an Einstein-de Sitter comparison model and the matter-only power law. "
           "Right: a_dot = H/(1+z) (in units of the Hubble constant) whose minimum marks the transition from deceleration "
           "to acceleration, here z = %.2f." % zmin)
    return FigureResult("fig02_hubble_parameter_vs_redshift", "Hubble parameter versus redshift", cap, fig,
                        {"z": z, "H_LCDM": bg.H_of_z(z), "H_EdS": eds.H_of_z(z)}, "numerical")


@figure("fig03_density_fractions")
def fig03() -> FigureResult:
    bg1 = fiducial()
    bg2 = FLRWBackground(CosmologyParams(omega_m=0.30, omega_k=0.05))
    a = np.logspace(-6, 1, 600)
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.3), sharey=True)
    data = {"a": a}
    for ax, bg, ttl in zip(axs, (bg1, bg2), ("Fiducial: flat LCDM", "Illustrative curved model: $\\Omega_k = 0.05$ (not a fit)")):
        fr = bg.omega_fractions(a)
        for (key, lab), c, ls in zip((("radiation", "radiation"), ("matter", "matter"), ("curvature", "curvature"),
                                      ("dark_energy", "Lambda")), (ORANGE, BLUE, VIOLET, AQUA), LINESTYLES):
            ax.semilogx(a, fr[key], color=c, ls=ls, label=lab)
            data[f"{ttl[:4]}_{key}"] = fr[key]
        ax.axvline(1.0, color=MUTED, lw=0.8)
        ax.text(1.1, 0.5, "today", color=INK2, rotation=90, va="center", fontsize=8.5)
        ax.set_xlabel("scale factor a")
        ax.set_title(ttl)
        ax.set_ylim(-0.08, 1.08)
    axs[0].set_ylabel("density parameter $\\Omega_i(a)$")
    axs[0].legend(loc="center left", fontsize=8.5)
    cap = ("Density parameters Omega_i(a) = rho_i/rho_crit(a) for radiation, matter, spatial curvature and a cosmological constant. "
           "The fractions sum to one at every a (checked in tests). The right panel uses an illustrative non-flat closure to show the "
           "curvature term; it is not a fit to data. Neutrinos are treated as massless radiation.")
    return FigureResult("fig03_density_fractions", "Density fractions versus scale factor", cap, fig, data, "numerical")


@figure("fig04_deceleration_parameter")
def fig04() -> FigureResult:
    bg = fiducial()
    a = np.logspace(-5, 1, 500)
    z = np.linspace(0, 2.5, 300)
    a_acc = eh.acceleration_onset_a(bg)
    z_acc = 1 / a_acc - 1
    a_acc_an = (bg.om_m / (2 * bg.om_de)) ** (1 / 3)
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.3))
    ax = axs[0]
    ax.semilogx(a, bg.deceleration(a), color=BLUE)
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.axhline(1.0, color=MUTED, lw=0.6, ls=":")
    ax.axhline(0.5, color=MUTED, lw=0.6, ls=":")
    ax.axhline(-1.0, color=MUTED, lw=0.6, ls=":")
    ax.text(1.5e-5, 1.03, "radiation: q = 1", fontsize=8.5, color=INK2)
    ax.text(1.5e-5, 0.53, "matter: q = 1/2", fontsize=8.5, color=INK2)
    ax.text(1.5e-5, -0.95, "de Sitter: q = -1", fontsize=8.5, color=INK2)
    ax.plot([a_acc], [0.0], "o", color=RED, ms=6)
    ax.set_xlabel("scale factor a")
    ax.set_ylabel("deceleration parameter q")
    ax.set_title("q(a) across all epochs, fiducial LCDM")
    ax = axs[1]
    for (w, c, ls) in ((-1.0, BLUE, "-"), (-0.8, ORANGE, "--"), (-1.2, AQUA, "-.")):
        b = _model(w0=w)
        ax.plot(z, b.deceleration(1 / (1 + z)), color=c, ls=ls, label=f"w = {w}")
    ax.plot(z, _eds().deceleration(1 / (1 + z)), color=VIOLET, ls=":", label="no dark energy")
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.annotate(f"onset of acceleration (LCDM)\nz = {z_acc:.3f}  (analytic, radiation neglected: "
                f"z = {1 / a_acc_an - 1:.3f})", (z_acc, 0), xytext=(15, 30), textcoords="offset points",
                arrowprops=dict(arrowstyle="-", color=INK2), fontsize=8.5)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("q(z)")
    ax.set_title("Acceleration transition for constant w")
    ax.legend(loc="lower right", fontsize=8.5)
    cap = ("Left: deceleration parameter q = (1/2) sum Omega_i (1+3w_i) over all epochs, tending to 1, 1/2 and -1 in radiation, matter and "
           "Lambda domination. Right: q(z) for constant dark-energy equation of state w and for a universe without dark energy. The "
           "numerical zero crossing for LCDM agrees with the analytic estimate a = (Om_m/2Om_L)^(1/3).")
    return FigureResult("fig04_deceleration_parameter", "Deceleration parameter and acceleration transition", cap, fig,
                        {"a": a, "q_LCDM": bg.deceleration(a)}, "numerical")


@figure("fig05_photon_temperature_scaling")
def fig05() -> FigureResult:
    from ..early_universe import thermal_history as th
    bg = fiducial()
    z = np.logspace(-2, 16, 400)
    T = th.temperature_of_z(z, bg.params.T_cmb)
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.loglog(1 + z, T, color=BLUE)
    marks = [("recombination (visibility peak)", 1079.0), ("matter-radiation equality", 1 / bg.a_equality_mr() - 1)]
    for name, T_GeV in (("e+- annihilation (T ~ m_e)", th.T_EE_GEV), ("QCD crossover (~0.17 GeV)", 0.17),
                        ("electroweak scale (~100 GeV)", 100.0)):
        marks.append((name, units.GeV_to_K(T_GeV) / bg.params.T_cmb - 1))
    for i, (name, zz) in enumerate(marks):
        ax.plot([1 + zz], [bg.params.T_cmb * (1 + zz)], "o", color=ORANGE, ms=5)
        ax.annotate(name, (1 + zz, bg.params.T_cmb * (1 + zz)), xytext=(8, -14 - 2 * i), textcoords="offset points",
                    fontsize=8.5, color=INK2)
    ax.set_xlabel("1 + z")
    ax.set_ylabel("photon temperature T [K]")
    ax.set_title("Photon temperature scaling T = T0 (1 + z)")
    sec = ax.secondary_yaxis("right", functions=(lambda x: x * K.k_B / K.eV, lambda x: x * K.eV / K.k_B))
    sec.set_ylabel("k$_B$T [eV]")
    sec.spines["right"].set_visible(True)
    ax.set_ylim(1, 1e30)
    cap = ("Adiabatic cooling of the CMB blackbody, T(z) = T0 (1+z) with T0 = 2.7255 K (adopted fiducial; Fixsen 2009 quotes 2.725 +- 0.001 K). "
           "Markers locate characteristic epochs computed from the model (recombination from the Peebles calculation; other temperatures are "
           "textbook scales, and the QCD/electroweak transitions are smooth crossovers, not sharp events). No temperature measurements are plotted.")
    return FigureResult("fig05_photon_temperature_scaling", "Photon temperature scaling", cap, fig,
                        {"z": z, "T_K": T}, "established")


@figure("fig07_comoving_distance")
def fig07() -> FigureResult:
    bg = fiducial()
    eds = _eds()
    z = np.logspace(-2, np.log10(1100), 160)
    dc = flrw.comoving_distance_mpc(bg, z) / 1e3
    dce = flrw.comoving_distance_mpc(eds, z) / 1e3
    b2 = _model(w0=-0.8)
    dcw = flrw.comoving_distance_mpc(b2, z) / 1e3
    ph = horizons.particle_horizon_mpc(bg) / 1e3
    eh_ = horizons.event_horizon_mpc(bg) / 1e3
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.semilogx(z, dc, color=BLUE, label="LCDM (fiducial)")
    ax.semilogx(z, dcw, color=ORANGE, ls="--", label="constant w = -0.8")
    ax.semilogx(z, dce, color=VIOLET, ls=":", label="Einstein-de Sitter")
    ax.axhline(ph, color=MUTED, lw=0.8, ls="-.")
    ax.text(0.012, ph + 0.25, f"particle horizon (LCDM) {ph:.1f} Gpc = {units.Mpc_to_Gly(ph * 1e3):.0f} Gly", fontsize=8.5, color=INK2)
    ax.axhline(eh_, color=MUTED, lw=0.8, ls="-.")
    ax.text(0.012, eh_ + 0.25, f"event horizon (LCDM) {eh_:.1f} Gpc = {units.Mpc_to_Gly(eh_ * 1e3):.0f} Gly", fontsize=8.5, color=INK2)
    ax.set_xlabel("redshift z")
    ax.set_ylabel("comoving distance D$_C$ [Gpc]")
    ax.set_title("Comoving distance to redshift z")
    ax.legend(loc="center right")
    cap = ("Line-of-sight comoving distance D_C(z) = (c/H0) int dz/E(z) for three background models (quadrature, tolerance 1e-11 relative). "
           "Horizontal lines: comoving particle horizon today (radiation included) and comoving event horizon today (exists only for the "
           "accelerating LCDM case). Distances are model-conditional.")
    return FigureResult("fig07_comoving_distance", "Comoving distance versus redshift", cap, fig,
                        {"z": z, "D_C_LCDM_Gpc": dc, "D_C_w-0.8_Gpc": dcw, "D_C_EdS_Gpc": dce}, "numerical")


@figure("fig08_luminosity_angular_distance")
def fig08() -> FigureResult:
    bg = fiducial()
    z = np.logspace(-2, np.log10(1100), 160)
    dl = flrw.luminosity_distance_mpc(bg, z) / 1e3
    da = flrw.angular_diameter_distance_mpc(bg, z) / 1e3
    eds = _eds()
    dl2 = flrw.luminosity_distance_mpc(eds, z) / 1e3
    da2 = flrw.angular_diameter_distance_mpc(eds, z) / 1e3
    i = int(np.argmax(da))
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.3))
    axs[0].loglog(z, dl, color=BLUE, label="LCDM")
    axs[0].loglog(z, dl2, color=VIOLET, ls=":", label="Einstein-de Sitter")
    axs[0].set_xlabel("redshift z")
    axs[0].set_ylabel("luminosity distance D$_L$ [Gpc]")
    axs[0].set_title("D$_L$ = (1+z) D$_M$")
    axs[0].legend(loc="upper left")
    axs[1].semilogx(z, da, color=BLUE, label="LCDM")
    axs[1].semilogx(z, da2, color=VIOLET, ls=":", label="Einstein-de Sitter")
    axs[1].plot([z[i]], [da[i]], "o", color=RED, ms=6)
    axs[1].annotate(f"maximum D$_A$ = {da[i]:.2f} Gpc at z = {z[i]:.2f}", (z[i], da[i]), xytext=(10, -40),
                    textcoords="offset points", arrowprops=dict(arrowstyle="-", color=INK2), fontsize=8.5)
    axs[1].set_xlabel("redshift z")
    axs[1].set_ylabel("angular-diameter distance D$_A$ [Gpc]")
    axs[1].set_title("D$_A$ = D$_M$/(1+z): turns over")
    axs[1].legend(loc="lower right")
    cap = ("Luminosity distance (left) and angular-diameter distance (right) for the flat fiducial model and an Einstein-de Sitter comparison. "
           "D_A has a maximum, so distant objects can subtend larger angles; the location of the maximum depends on the model. "
           "Distance duality D_L = (1+z)^2 D_A holds by construction in the implemented metric geometry.")
    return FigureResult("fig08_luminosity_angular_distance", "Luminosity and angular-diameter distances", cap, fig,
                        {"z": z, "D_L_Gpc": dl, "D_A_Gpc": da}, "numerical")


@figure("fig11_linear_growth_factor")
def fig11() -> FigureResult:
    bg = fiducial()
    models = [("LCDM (fiducial)", bg, BLUE, "-"), ("constant w = -0.8", _model(w0=-0.8), ORANGE, "--"),
              ("open, no dark energy", _model(omega_m=0.315, omega_de=0.0), VIOLET, ":"), ("Einstein-de Sitter", _eds(), AQUA, "-.")]
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.3))
    data = {}
    for lab, b, c, ls in models:
        g = gf.solve_growth(b)
        axs[0].plot(g.a, g.D_plus / g.a, color=c, ls=ls, label=lab)
        axs[1].plot(g.a, g.f, color=c, ls=ls, label=lab)
        if lab.startswith("LCDM"):
            data = {"a": g.a, "D_plus": g.D_plus, "f": g.f}
            om = b.omega_fractions(g.a)["matter"]
            axs[1].plot(g.a, om**0.55, color=MUTED, lw=1.0, ls=(0, (1, 2)), label="$\\Omega_m(a)^{0.55}$ (LCDM)")
    axs[0].set_xscale("log")
    axs[1].set_xscale("log")
    axs[0].set_xlim(1e-2, 1)
    axs[1].set_xlim(1e-2, 1)
    axs[0].set_xlabel("scale factor a")
    axs[0].set_ylabel("$D_+(a)/a$   [$D_+\\to a$ in matter era]")
    axs[0].set_title("Growth suppression relative to matter domination")
    axs[0].legend(loc="lower left", fontsize=8.5)
    axs[1].set_xlabel("scale factor a")
    axs[1].set_ylabel("growth rate f = d ln D / d ln a")
    axs[1].set_title("Linear growth rate")
    axs[1].legend(loc="lower left", fontsize=8.5)
    cap = ("Linear, scale-independent growth of pressureless matter perturbations from the equation D'' + (1-q) D' = (3/2) Om_m(a) D with "
           "Meszaros initial conditions. Valid on sub-horizon, linear scales with smooth dark energy; no neutrino free streaming, no modified "
           "gravity. This is NOT a nonlinear structure-formation calculation and does not predict galaxies or halos by itself.")
    return FigureResult("fig11_linear_growth_factor", "Linear growth factor", cap, fig, data, "approximate")


def _fig12_data():
    sc = des.default_scenarios()
    out = []
    for s in sc:
        bg = FLRWBackground(s.params)
        tr = des.trajectory(bg, a_start=0.2, n=420)
        out.append((s, bg, tr, des.classify_fate(bg)))
    return out


@figure("fig12_future_scenarios")
def fig12() -> FigureResult:
    res = _fig12_data()
    colors = [BLUE, ORANGE, AQUA, VIOLET, MAGENTA, RED, "#6b6b6b", "#b08900"]
    fig, axs = plt.subplots(1, 2, figsize=(11.2, 4.6))
    data = {}
    for (s, bg, tr, fate), c, ls in zip(res, colors, [LINESTYLES[i % 4] for i in range(len(res))]):
        m = (tr.t_gyr < 120) & np.isfinite(tr.t_gyr)
        axs[0].plot(tr.t_gyr[m], tr.a[m], color=c, ls=ls, label=s.label)
        axs[1].plot(tr.t_gyr[m], tr.H_over_H0[m], color=c, ls=ls)
        if fate["t_end_gyr"] is not None and fate["t_end_gyr"] < 120:
            axs[0].plot([fate["t_end_gyr"]], [0.25 if np.isfinite(fate["a_max"]) else 6e3], marker="x", color=c, ms=8, mew=2, clip_on=False)
        data[f"t_Gyr_{s.key}"] = tr.t_gyr
        data[f"a_{s.key}"] = tr.a
    axs[0].set_yscale("log")
    axs[0].set_ylim(0.2, 1e4)
    axs[0].set_xlim(0, 120)
    axs[0].axvline(ct.age_gyr(fiducial()), color=MUTED, lw=0.8)
    axs[0].text(ct.age_gyr(fiducial()) + 1, 0.3, "today", fontsize=8.5, color=INK2)
    axs[0].set_xlabel("cosmic time t [Gyr]")
    axs[0].set_ylabel("scale factor a(t)")
    axs[0].set_title("Scale factor in conditional future scenarios")
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="upper center", ncol=4, fontsize=8.2, bbox_to_anchor=(0.5, 0.0))
    axs[1].set_yscale("symlog", linthresh=1.0)
    axs[1].set_ylim(-20, 3e2)
    axs[1].set_xlim(0, 120)
    axs[1].axhline(0, color=MUTED, lw=0.8)
    axs[1].set_xlabel("cosmic time t [Gyr]")
    axs[1].set_ylabel("H(t) / H0")
    axs[1].set_title("Expansion rate (negative = contraction)")
    cap = ("Classical future expansion of eight scenarios, each defined by explicit density parameters (see docs/FUTURE_COSMOLOGY.md). "
           "Crosses mark finite-time singularities (big rip for w < -1, big crunch for recollapsing models). These are conditional outcomes of "
           "assumed models and classical General Relativity, NOT predictions of the actual final state of the universe; vacuum decay, "
           "quantum gravity and unknown physics are excluded.")
    return FigureResult("fig12_future_scenarios", "Comparison of conditional future scenarios", cap, fig, data, "model")


@figure("fig20_numerical_convergence")
def fig20() -> FigureResult:
    bg = fiducial()
    ref = ct.age_gyr(bg)
    amin = np.logspace(-14, -3, 12)
    e_cut = np.abs(cv.age_vs_lower_cutoff(bg, amin) - ref) / ref
    ns = np.array([21, 41, 81, 161, 321, 641, 1281, 2561])
    e_tr = np.abs(cv.age_vs_grid(bg, ns, "trapezoid") - ref) / ref
    e_si = np.abs(cv.age_vs_grid(bg, ns, "simpson") - ref) / ref
    eps = np.logspace(-4, -11, 8)
    e_tol = np.abs(cv.tolerance_sweep(bg, eps) - ref) / ref
    amx = np.logspace(1, 7, 7)
    ehr = horizons.event_horizon_mpc(bg)
    e_eh = np.abs(cv.horizon_vs_upper_cutoff(bg, amx) - ehr) / ehr
    floor = 1e-17
    fig, axs = plt.subplots(2, 2, figsize=(10.4, 8.0))
    fig.subplots_adjust(hspace=0.42, wspace=0.28)
    ax = axs[0, 0]
    ax.loglog(amin, np.maximum(e_cut, floor), "o-", color=BLUE)
    pred = (amin ** 2 / (2 * np.sqrt(bg.om_r))) / (ref / ct.hubble_time_gyr(bg))
    ax.loglog(amin, pred, "--", color=ORANGE, label="analytic: $a_{min}^2/(2\\sqrt{\\Omega_r})$ / age")
    ax.set_xlabel("lower cutoff a$_{min}$")
    ax.set_ylabel("relative error in age")
    ax.set_title("Finite lower cutoff (radiation era)")
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    ax.loglog(ns, np.maximum(e_tr, floor), "o-", color=BLUE, label="trapezoid in ln a")
    ax.loglog(ns, np.maximum(e_si, floor), "s--", color=ORANGE, label="Simpson in ln a")
    ax.loglog(ns, e_tr[0] * (ns / ns[0]) ** -2.0, ":", color=MUTED, label="slope -2")
    ax.loglog(ns, e_si[0] * (ns / ns[0]) ** -4.0, ":", color=VIOLET, label="slope -4")
    ax.set_xlabel("grid points N (log-spaced a from 1e-10 to 1)")
    ax.set_ylabel("relative error in age")
    ax.set_title("Fixed-grid quadrature order of convergence")
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    ax.loglog(eps, np.maximum(e_tol, floor), "o-", color=BLUE)
    ax.set_xlabel("requested relative tolerance (adaptive quad)")
    ax.set_ylabel("relative error in age")
    ax.set_title("Adaptive-quadrature tolerance sweep")
    ax.invert_xaxis()
    ax = axs[1, 1]
    ax.loglog(amx, np.maximum(e_eh, floor), "o-", color=BLUE, label="numerical")
    ax.loglog(amx, 1.0 / (np.sqrt(bg.om_de) * amx) / (ehr / (units.hubble_distance_Mpc(bg.H0))), "--", color=ORANGE,
              label="analytic tail $1/(\\sqrt{\\Omega_\\Lambda} a_{max})$")
    ax.set_xlabel("upper cutoff a$_{max}$")
    ax.set_ylabel("relative error in event horizon")
    ax.set_title("Finite upper cutoff of the event-horizon integral")
    ax.legend(fontsize=8)
    cap = ("Numerical-convergence diagnostics for the background integrals of the fiducial model; the reference is the adaptive quadrature at "
           "tolerance 1e-11. Errors below 1e-16 are floored for display. (a) Starting the age integral at a_min instead of 0 removes ~a_min^2/(2 sqrt(Om_r)). "
           "(b) Fixed-grid rules converge with the expected orders. (c) Adaptive tolerance sweep (the achieved error is far below loose requests because the adaptive rule overshoots). (d) Truncating the event-horizon integral at a_max "
           "leaves the analytic tail 1/(sqrt(Om_L) a_max).")
    return FigureResult("fig20_numerical_convergence", "Numerical convergence diagnostics", cap, fig,
                        {"a_min": amin, "err_cutoff": e_cut, "N": ns.astype(float), "err_trap": e_tr, "err_simpson": e_si}, "numerical")


@figure("fig21_parameter_sensitivity")
def fig21() -> FigureResult:
    base = CosmologyParams()
    from ..geometry import flrw as fl

    quantities = {
        "age t0": lambda b: ct.age_gyr(b),
        "D_C(z=1)": lambda b: fl.comoving_distance_mpc(b, 1.0),
        "D_C(z=1100)": lambda b: fl.comoving_distance_mpc(b, 1100.0),
        "particle horizon": lambda b: horizons.particle_horizon_mpc(b),
        "1/a_eq": lambda b: 1.0 / b.a_equality_mr(),
        "event horizon": lambda b: horizons.event_horizon_mpc(b),
    }
    params = ["h", "omega_m"]
    sens = np.array([[unc.log_sensitivity(fn, base, p) for p in params] for fn in quantities.values()])
    ages = unc.monte_carlo(lambda b: ct.age_gyr(b), base, n=300)
    fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.4))
    y = np.arange(len(quantities))
    axs[0].barh(y - 0.19, sens[:, 0], height=0.36, color=BLUE, label="h")
    axs[0].barh(y + 0.19, sens[:, 1], height=0.36, color=ORANGE, label="$\\Omega_m$ (flat closure, $\\Omega_b h^2$ fixed)", hatch="//", edgecolor="white")
    axs[0].set_yticks(y)
    axs[0].set_yticklabels(list(quantities))
    axs[0].invert_yaxis()
    axs[0].axvline(0, color=INK2, lw=0.8)
    axs[0].set_xlabel("logarithmic sensitivity  d ln X / d ln parameter")
    axs[0].set_title("Parameter sensitivities of background quantities")
    axs[0].legend(loc="lower right", fontsize=8)
    axs[1].hist(ages, bins=24, color=BLUE, alpha=0.85, edgecolor="white")
    axs[1].axvline(np.mean(ages), color=ORANGE, lw=1.6)
    axs[1].set_xlabel("age of the universe [Gyr]")
    axs[1].set_ylabel("Monte Carlo draws (n = 300)")
    axs[1].set_title(f"Illustrative age spread: {np.mean(ages):.3f} +- {np.std(ages):.3f} Gyr")
    cap = ("Left: central-difference logarithmic sensitivities of background quantities to h and Om_m (the density-parameter derivative includes the "
           "flat closure Om_L = 1 - Om_m - Om_r with Om_b h^2 held fixed). Right: age distribution from INDEPENDENT Gaussian draws of h and Om_m with "
           "widths 0.5 km/s/Mpc and 0.007 quoted in the Planck 2018 VI abstract. Real posteriors are correlated; this is illustrative, not an error budget.")
    return FigureResult("fig21_parameter_sensitivity", "Parameter sensitivity and illustrative uncertainty", cap, fig,
                        {"age_Gyr_mc": ages}, "approximate")


@figure("fig22_conformal_diagrams")
def fig22() -> FigureResult:
    cases = [("Flat LCDM (fiducial)", fiducial(), 1e3), ("Einstein-de Sitter", _eds(), 3e1),
             ("Closed matter-only ($\\Omega_m = 2$)", _model(omega_m=2.0, omega_de=0.0), None)]
    fig, axs = plt.subplots(1, 3, figsize=(13.4, 4.6))
    data = {}
    for ax, (ttl, bg, amax_cut) in zip(axs, cases):
        top = ct.a_turnaround(bg)
        a_hi = top * (1 - 1e-9) if top is not None else amax_cut
        a = np.logspace(-9, np.log10(a_hi), 500)
        eta = np.asarray(ct.conformal_time_mpc(bg, a)) / 1e3        # Gpc
        eta0 = float(ct.conformal_time_mpc(bg, 1.0)) / 1e3
        etaf = eta[-1]
        if top is not None:
            eta_inf = np.inf          # recollapsing: event horizon not drawn (expansion phase only)
        else:
            eh_tail = horizons.event_horizon_mpc(bg, a_hi) / 1e3
            eta_inf = float(ct.conformal_time_mpc(bg, a_hi)) / 1e3 + eh_tail if np.isfinite(eh_tail) else np.inf
        chi_hub = np.asarray(horizons.comoving_hubble_radius_mpc(bg, a)) / 1e3
        ax.plot(eta, eta, color=BLUE, label="particle horizon of an observer at $\\chi = 0$")
        m = eta <= eta0
        ax.plot(eta[m], eta0 - eta[m], color=ORANGE, ls="--", label="our past light cone")
        m2 = eta >= eta0
        if np.isfinite(eta_inf):
            ax.plot(eta[m2], np.maximum(eta_inf - eta[m2], 0), color=AQUA, ls="-.", label="event horizon (if it exists)")
        ax.plot(eta, chi_hub, color=VIOLET, ls=":", label="comoving Hubble radius c/(aH)")
        ax.axvline(eta0, color=MUTED, lw=0.8)
        ax.set_xlim(0, min(1.05 * etaf, 60))
        ax.set_ylim(0, max(1.1 * eta0, 5))
        ax.set_xlabel("conformal time $\\eta$ [comoving Gpc]")
        ax.set_title(ttl, fontsize=10)
        data[f"eta_Gpc_{ttl[:4]}"] = eta
    axs[0].set_ylabel("comoving distance $\\chi$ [Gpc]")
    axs[0].legend(loc="upper left", fontsize=7.4)
    cap = ("Conformal-time diagrams (light rays are 45 degree lines, chi = +-eta). Curves: particle horizon, the past light cone of an observer "
           "today (vertical line marks conformal time today), the event horizon when finite, and the comoving Hubble radius c/(aH). For the accelerating "
           "LCDM case the comoving Hubble radius grows again after matter domination while the event horizon shrinks; "
           "the Einstein-de Sitter case has no finite event horizon. The closed model is shown up to maximum expansion only. Spatial curvature is not drawn.")
    return FigureResult("fig22_conformal_diagrams", "Conformal diagrams and horizons", cap, fig, data, "numerical")
