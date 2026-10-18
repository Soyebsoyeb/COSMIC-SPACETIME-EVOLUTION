"""Logarithmic cosmic timeline with epistemic-status labelling."""
from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np

from .. import constants as K
from ..background import cosmic_time as ct
from ..background.friedmann import FLRWBackground
from ..config import CosmologyParams
from ..early_universe import inflation as inf
from ..early_universe import thermal_history as th
from ..future.de_sitter import e_folding_time_gyr
from .publication_style import (AQUA, BLUE, FigureResult, INK2, MAGENTA, MUTED, ORANGE, RED, VIOLET, apply_style, figure)

apply_style()

CATEGORIES = {
    "established": ("Established theory", BLUE, None),
    "numerical": ("Computed by this code", AQUA, None),
    "approximate": ("Approximate, stated validity", VIOLET, "//"),
    "model": ("Model-dependent", ORANGE, ".."),
    "speculative": ("Speculative / unresolved", RED, "xx"),
}


def epoch_list():
    p = CosmologyParams()
    bg = FLRWBackground(p)
    gyr = K.Gyr
    a_eq = bg.a_equality_mr()
    a_lm = (bg.om_m / bg.om_de) ** (1 / 3)
    a_acc = (bg.om_m / (2 * bg.om_de)) ** (1 / 3)
    t = lambda a: float(ct.cosmic_time_gyr(bg, a)) * gyr
    t0 = ct.age_gyr(bg) * gyr
    sol = inf.solve_background(inf.StarobinskyModel())
    info = inf.pivot_efolds(inf.StarobinskyModel(), p.A_s, p.k_pivot, sol=sol)
    t_end_inf = 1.0 / (2.0 * info["H_end"] * inf.M_P_GEV * K.GeV_to_inv_s)   # radiation-era-like estimate t ~ 1/2H
    dur = float(np.trapezoid(1.0 / (np.sqrt(info["V0"] * inf.StarobinskyModel().v(sol.phi) / (3 - sol.eps_H)) * inf.M_P_GEV * K.GeV_to_inv_s), sol.N))
    return [
        ("Planck era: classical GR not applicable", 1e-45, K.t_planck, "speculative"),
        ("Inflation (Starobinsky-type model; duration only)", t_end_inf - min(dur, 0.9 * t_end_inf) , t_end_inf, "model"),
        ("Reheating (T_reh unconstrained between here and ~MeV)", t_end_inf, float(th.time_of_temperature_s(1e-3)), "model"),
        ("Electroweak and QCD crossovers", float(th.time_of_temperature_s(100.0)), float(th.time_of_temperature_s(0.1)), "approximate"),
        ("Neutrino decoupling, e+- annihilation", float(th.time_of_temperature_s(1e-3)), float(th.time_of_temperature_s(5e-4)), "approximate"),
        ("Nucleosynthesis window (network not run)", float(th.time_of_temperature_s(1e-4)), float(th.time_of_temperature_s(3e-5)), "established"),
        ("Radiation domination", float(th.time_of_temperature_s(1e-3)), t(a_eq), "numerical"),
        ("Matter domination", t(a_eq), t(a_lm), "numerical"),
        ("Recombination / last scattering", t(1 / 1200.0), t(1 / 1000.0), "approximate"),
        ("Structure formation (qualitative, not computed)", t(1 / 31.0), t0, "speculative"),
        ("Lambda domination and acceleration", t(a_acc), t0 * 1e3, "numerical"),
        ("Far future: conditional on constant Lambda", t0, 1e19, "model"),
    ], t0


@figure("fig06_cosmic_timeline")
def fig06() -> FigureResult:
    rows, t0 = epoch_list()
    fig, ax = plt.subplots(figsize=(11.0, 6.2))
    data = {}
    for i, (name, a, b, cat) in enumerate(rows):
        lab, color, hatch = CATEGORIES[cat]
        ax.barh(i, math.log10(b) - math.log10(a), left=math.log10(a), height=0.62, color=color, alpha=0.85,
                hatch=hatch, edgecolor="white", linewidth=0.6)
        la, lb = math.log10(a), math.log10(b)
        if lb < 8.0:
            ax.text(lb + 0.4, i, name, va="center", ha="left", fontsize=8.6)
        else:
            ax.text(la - 0.4, i, name, va="center", ha="right", fontsize=8.6)
        data[f"t_start_s_{i}"] = np.array([a])
        data[f"t_end_s_{i}"] = np.array([b])
    ax.axvline(math.log10(t0), color="black", lw=1.0)
    ax.text(math.log10(t0) + 0.25, -0.75, "today (13.8 Gyr)", fontsize=8.5)
    ax.set_yticks([])
    ax.set_xlim(-46, 21)
    ax.set_ylim(len(rows) - 0.4, -1.0)
    ax.set_xlabel("cosmic time since the model's t = 0  [log$_{10}$ seconds]")
    ax.set_title("Cosmic timeline with epistemic status of each epoch")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, hatch=h, ec="white", alpha=0.85) for (_, c, h) in CATEGORIES.values()]
    ax.legend(handles, [v[0] for v in CATEGORIES.values()], loc="lower left", fontsize=8.5, title="status")
    ax.grid(axis="y", visible=False)
    cap = ("Logarithmic timeline of the universe's epochs. Colour AND hatch encode the epistemic category: established theory; computed by this code; "
           "approximate with stated validity; model dependent; speculative. Times of post-inflationary epochs come from the radiation-era t = 1/(2H) relation with the "
           "g_* step approximation (about 10% uncertain near thresholds) and from the fiducial background; inflation timing is shown for one illustrative model and its "
           "time origin is not physically meaningful; the Planck-era bar spans the regime where classical GR cannot be trusted. Structure-formation epoch is qualitative "
           "(starting at z ~ 30 as an indicative marker) and not computed by the background solver.")
    return FigureResult("fig06_cosmic_timeline", "Cosmic timeline with model dependence", cap, fig, data, "model")
