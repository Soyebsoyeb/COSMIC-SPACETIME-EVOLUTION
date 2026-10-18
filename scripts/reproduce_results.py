#!/usr/bin/env python
"""Reproduce every result: validation, key results, sweeps, figures; save a config snapshot and checksums."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(script: str, *args: str) -> None:
    print(f"==> {script} {' '.join(args)}", flush=True)
    subprocess.run([PY, str(ROOT / "scripts" / script), *args], check=True)


def main() -> int:
    run("validate_numerics.py")
    run("compute_key_results.py")
    run("run_parameter_sweep.py")
    run("generate_all_figures.py")
    meta = ROOT / "results" / "metadata"
    meta.mkdir(parents=True, exist_ok=True)
    shutil.copy(ROOT / "configs" / "fiducial_lcdm.yaml", meta / "config_snapshot_fiducial_lcdm.yaml")
    sums = {str(p.relative_to(ROOT)): sha256(p) for p in sorted((ROOT / "results").rglob("*"))
            if p.is_file() and p.suffix in {".csv", ".md"} and "checksums" not in p.name}
    (meta / "checksums_tables_datasets.json").write_text(json.dumps(sums, indent=2, sort_keys=True))
    print("reproduction complete; checksums written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
