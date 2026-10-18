#!/usr/bin/env python
"""Generate every figure, dataset and metadata file for the repository.

Usage:  python scripts/generate_all_figures.py [--only fig01 fig02 ...] [--outdir results]
Deterministic: fixed seeds, fixed tolerances; figures embed no timestamps.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib  # noqa: E402
import numpy as np  # noqa: E402
import scipy  # noqa: E402

from cosmic_evolution import __version__  # noqa: E402
from cosmic_evolution.config import CosmologyParams  # noqa: E402
from cosmic_evolution.visualization import publication_style as ps  # noqa: E402
from cosmic_evolution.visualization import (cosmology_plots, early_universe_plots, structure_plots, timeline)  # noqa: E402,F401


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--outdir", default=str(ROOT / "results"))
    args = ap.parse_args(argv)
    out = Path(args.outdir)
    registry = {}
    p = CosmologyParams()
    for name, fn in sorted(ps.FIGURES.items()):
        if args.only and not any(name.startswith(o) for o in args.only):
            continue
        t0 = time.time()
        res = fn()
        paths = ps.save_figure(res, out / "figures")
        hdr = [f"dataset for {name}: {res.title}", f"fiducial parameters: h={p.h}, Omega_m={p.omega_m}, Omega_b h^2={p.omega_b_h2}, "
               f"T_cmb={p.T_cmb} K, N_eff={p.N_eff}, A_s={p.A_s}, n_s={p.n_s}", "units: see column names; time in Gyr or s, distances in Mpc or Gpc, k in 1/Mpc",
               f"epistemic category: {res.category}", f"code version: {__version__}"]
        csv = ps.write_dataset_csv(name, res.data, out / "datasets", hdr)
        registry[name] = {"title": res.title, "caption": res.caption, "category": res.category,
                          "png": str(Path(paths["png"]).relative_to(out.parent)), "pdf": str(Path(paths["pdf"]).relative_to(out.parent)),
                          "dataset": str(csv.relative_to(out.parent)) if csv else None, "seconds": round(time.time() - t0, 2)}
        print(f"[ok] {name}  ({time.time() - t0:.1f} s)")
    meta = out / "metadata"
    meta.mkdir(parents=True, exist_ok=True)
    existing = {}
    if args.only and (meta / "figure_registry.json").exists():
        existing = json.loads((meta / "figure_registry.json").read_text())
    existing.update(registry)
    (meta / "figure_registry.json").write_text(json.dumps(existing, indent=2, sort_keys=True))
    (meta / "environment.json").write_text(json.dumps({
        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__, "package_version": __version__,
        "fiducial_parameters": p.to_dict()}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
