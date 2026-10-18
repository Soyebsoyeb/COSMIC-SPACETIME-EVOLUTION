"""Command-line interface:  cosmic-evolution {summary,figures,sweep,validate,reproduce,fate}."""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from . import FLRWBackground
from .background import cosmic_time as ct
from .config import CosmologyParams

_ROOT = Path(__file__).resolve().parents[2]


def _script(name: str, *args: str) -> int:
    path = _ROOT / "scripts" / name
    if not path.exists():
        print(f"script {name} not found next to the package (run from a source checkout)", file=sys.stderr)
        return 2
    return subprocess.call([sys.executable, str(path), *args])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="cosmic-evolution", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summary", help="print background quantities for a YAML config")
    s.add_argument("--config", default=None)
    sub.add_parser("figures", help="regenerate all figures and datasets")
    sub.add_parser("sweep", help="run the parameter sweep and model comparison")
    sub.add_parser("validate", help="run the numerical validation report")
    sub.add_parser("reproduce", help="run everything")
    sub.add_parser("fate", help="print the conditional-fate table")
    args = ap.parse_args(argv)
    if args.cmd == "summary":
        p = CosmologyParams.from_yaml(args.config) if args.config else CosmologyParams()
        bg = FLRWBackground(p)
        print(f"H0 = {p.H0:.3f} km/s/Mpc; Omega_m = {p.omega_m}; Omega_r = {p.omega_r:.4e}; Omega_de = {p.omega_de_eff:.5f}")
        print(f"age = {ct.age_gyr(bg):.4f} Gyr; z_eq = {1 / bg.a_equality_mr() - 1:.1f}")
        return 0
    if args.cmd == "fate":
        from .future.fate_classification import fate_table
        for r in fate_table():
            print(f"{r['label']:55s} {r['fate']}")
        return 0
    return _script({"figures": "generate_all_figures.py", "sweep": "run_parameter_sweep.py",
                    "validate": "validate_numerics.py", "reproduce": "reproduce_results.py"}[args.cmd])


if __name__ == "__main__":
    raise SystemExit(main())
