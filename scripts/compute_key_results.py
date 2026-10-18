#!/usr/bin/env python
"""Compute headline numbers for the fiducial model and write results/tables/key_results.{csv,md}."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

from cosmic_evolution import FLRWBackground  # noqa: E402
from cosmic_evolution.background import cosmic_time as ct  # noqa: E402
from cosmic_evolution.background import expansion_history as eh  # noqa: E402
from cosmic_evolution.config import CosmologyParams  # noqa: E402
from cosmic_evolution.early_universe import entropy, inflation as inf, recombination as rc  # noqa: E402
from cosmic_evolution.future import de_sitter as ds  # noqa: E402
from cosmic_evolution.geometry import flrw, horizons  # noqa: E402
from cosmic_evolution.perturbations import growth_factor as gf  # noqa: E402
from cosmic_evolution.structure import power_spectrum as pk  # noqa: E402
from cosmic_evolution import units  # noqa: E402


def compute(config: Path | None = None):
    p = CosmologyParams.from_yaml(config) if config else CosmologyParams()
    bg = FLRWBackground(p)
    g = gf.solve_growth(bg)
    rec = rc.solve_recombination(bg)
    a_acc = eh.acceleration_onset_a(bg)
    rows = [
        # (quantity, value, unit, category, method)
        ("Age of the universe t0", ct.age_gyr(bg), "Gyr", "numerical", "exact quadrature of the Friedmann equation"),
        ("Hubble time 1/H0", ct.hubble_time_gyr(bg), "Gyr", "numerical", "unit conversion"),
        ("Omega_r (photons + massless neutrinos)", p.omega_r, "", "numerical", "Stefan-Boltzmann with N_eff"),
        ("Matter-radiation equality z_eq", 1 / bg.a_equality_mr() - 1, "", "numerical", "Omega_r a^-4 = Omega_m a^-3 (massless-nu approximation)"),
        ("Onset of acceleration z_acc (q = 0)", 1 / a_acc - 1, "", "numerical", "root of q(a)"),
        ("Cosmic time at onset of acceleration", float(ct.cosmic_time_gyr(bg, a_acc)), "Gyr", "numerical", "t(a_acc)"),
        ("Comoving particle horizon today", horizons.particle_horizon_mpc(bg) / 1e3, "Gpc", "numerical", "int_0^1 da/sqrt(Q)"),
        ("Comoving event horizon today", horizons.event_horizon_mpc(bg) / 1e3, "Gpc", "numerical", "int_1^inf da/sqrt(Q)"),
        ("Comoving distance to z = 1100", flrw.comoving_distance_mpc(bg, 1100.0) / 1e3, "Gpc", "numerical", "quadrature"),
        ("Angular-diameter distance to z = 1100", flrw.angular_diameter_distance_mpc(bg, 1100.0), "Mpc", "numerical", "quadrature"),
        ("Linear growth suppression D+(1)", float(g.D(1.0)), "", "approximate", "linear growth ODE, scale independent"),
        ("sigma_8 (EH no-wiggle, A_s = 2.1e-9)", pk.sigma8(bg, growth=g), "", "approximate", "fitting-formula transfer function; Planck-quoted 0.811 +- 0.006"),
        ("Baryon-to-photon ratio eta", entropy.baryon_to_photon_ratio(p), "", "numerical", "Omega_b h^2 and T_cmb"),
        ("Recombination visibility peak z*", rec.z_star, "", "approximate", "Peebles three-level hydrogen model"),
        ("Comoving sound horizon at z*", rec.r_s_mpc, "Mpc", "approximate", "Peebles model + background"),
        ("100 theta* = 100 r_s/D_M", 100 * rec.theta_star, "", "approximate", "Planck 2018 reports 100 theta_MC = 1.0409 +- 0.0003"),
        ("de Sitter horizon radius c/H_Lambda", ds.proper_event_horizon_gly(bg), "Gly", "numerical", "asymptotic future (conditional on Lambda)"),
        ("Gibbons-Hawking temperature", ds.gibbons_hawking_temperature_K(bg), "K", "established", "hbar H_Lambda/(2 pi k_B), pure de Sitter"),
    ]
    for model in (inf.StarobinskyModel(), inf.QuadraticModel()):
        sol = inf.solve_background(model)
        info = inf.pivot_efolds(model, p.A_s, p.k_pivot, sol=sol)
        o = inf.slow_roll_observables(model, sol.phi_at(info["N_star"]))
        rows += [
            (f"[{model.name}] N* (instantaneous reheating)", info["N_star"], "e-folds", "model", "self-consistent pivot e-folds"),
            (f"[{model.name}] n_s at N*", o["n_s"], "", "model", "leading-order slow roll, exact background"),
            (f"[{model.name}] r at N*", o["r"], "", "model", "16 eps_V"),
            (f"[{model.name}] H_end", info["H_end"] * inf.M_P_GEV, "GeV", "model", "amplitude fixed by A_s = 2.1e-9"),
        ]
    return rows


def main(argv=None) -> int:
    out = ROOT / "results" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    cfg = ROOT / "configs" / "fiducial_lcdm.yaml"
    rows = compute(cfg)
    with open(out / "key_results.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("# Key results for configs/fiducial_lcdm.yaml; category = epistemic status (see docs/PHYSICS.md)\n")
        w = csv.writer(fh)
        w.writerow(["quantity", "value", "unit", "category", "method"])
        for r in rows:
            w.writerow([r[0], f"{r[1]:.8g}", r[2], r[3], r[4]])
    lines = ["| Quantity | Value | Unit | Status | Method |", "|---|---:|---|---|---|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]:.6g} | {r[2]} | {r[3]} | {r[4]} |")
    (out / "key_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
