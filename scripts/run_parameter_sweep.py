#!/usr/bin/env python
"""Configuration-driven parameter sweep over (h, Omega_m) and a dark-energy model comparison."""
from __future__ import annotations

import csv
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import yaml  # noqa: E402

from cosmic_evolution import FLRWBackground  # noqa: E402
from cosmic_evolution.background import cosmic_time as ct  # noqa: E402
from cosmic_evolution.background import expansion_history as eh  # noqa: E402
from cosmic_evolution.config import CosmologyParams  # noqa: E402
from cosmic_evolution.geometry import flrw, horizons  # noqa: E402


def row(bg):
    a_acc = eh.acceleration_onset_a(bg)
    return {"age_Gyr": ct.age_gyr(bg), "z_eq": 1 / bg.a_equality_mr() - 1,
            "z_acc": (1 / a_acc - 1) if a_acc else float("nan"),
            "particle_horizon_Gpc": horizons.particle_horizon_mpc(bg) / 1e3,
            "D_A_z1100_Mpc": float(flrw.angular_diameter_distance_mpc(bg, 1100.0))}


def main() -> int:
    out = ROOT / "results" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    sweep_cfg = yaml.safe_load((ROOT / "configs" / "parameter_sweep.yaml").read_text())
    base = CosmologyParams.from_yaml(ROOT / "configs" / sweep_cfg["base"])
    with open(out / "parameter_sweep.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("# Sweep of flat LCDM over (h, Omega_m); Omega_b h^2 held fixed; configs/parameter_sweep.yaml\n")
        w = None
        for h in sweep_cfg["sweep"]["h"]:
            for om in sweep_cfg["sweep"]["omega_m"]:
                r = row(FLRWBackground(replace(base, h=h, omega_m=om)))
                if w is None:
                    w = csv.DictWriter(fh, fieldnames=["h", "omega_m", *r.keys()])
                    w.writeheader()
                w.writerow({"h": h, "omega_m": om, **{k: f"{v:.8g}" for k, v in r.items()}})
    mc = yaml.safe_load((ROOT / "configs" / "model_comparison.yaml").read_text())
    base = CosmologyParams.from_yaml(ROOT / "configs" / mc["base"])
    with open(out / "model_comparison.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write("# Dark-energy equation-of-state comparison (flat closure); configs/model_comparison.yaml\n")
        w = None
        for m in mc["models"]:
            r = row(FLRWBackground(replace(base, w0=m["w0"], wa=m["wa"])))
            if w is None:
                w = csv.DictWriter(fh, fieldnames=["model", "w0", "wa", *r.keys()])
                w.writeheader()
            w.writerow({"model": m["name"], "w0": m["w0"], "wa": m["wa"], **{k: f"{v:.8g}" for k, v in r.items()}})
    print("wrote parameter_sweep.csv and model_comparison.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
