"""Build (if missing) and execute the notebooks in notebooks/ in place."""
from __future__ import annotations

import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    for path in sorted((ROOT / "notebooks").glob("*.ipynb")):
        nb = nbformat.read(path, as_version=4)
        NotebookClient(nb, timeout=600, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}}).execute()
        nbformat.write(nb, path)
        print("executed", path.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
