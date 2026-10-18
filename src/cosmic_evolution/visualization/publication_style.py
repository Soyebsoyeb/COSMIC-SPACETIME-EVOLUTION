"""Deterministic publication style, validated categorical palette and figure registry/saving helpers.

Palette: first six slots of a validated categorical palette (blue, orange, aqua, violet, magenta, red).  Identity is
never colour-only: every multi-series plot also uses distinct line styles and/or direct labels and legends.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BLUE, ORANGE, AQUA, VIOLET, MAGENTA, RED = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4", "#e34948"
YELLOW, GREEN = "#eda100", "#008300"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e3e2dd"
SERIES = [BLUE, ORANGE, AQUA, VIOLET, MAGENTA, RED]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (1, 1))]


def apply_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 100, "savefig.dpi": 200, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
        "axes.titleweight": "regular", "axes.labelsize": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.spines.top": False,
        "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.axisbelow": True, "lines.linewidth": 1.8, "legend.frameon": False, "legend.fontsize": 9,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.hashsalt": "cosmic-evolution", "axes.prop_cycle":
        matplotlib.cycler(color=SERIES), "mathtext.default": "regular",
    })


@dataclass
class FigureResult:
    name: str
    title: str
    caption: str
    fig: "plt.Figure"
    data: dict = field(default_factory=dict)          # column name -> 1D array (written to CSV)
    category: str = "numerical"                        # see docs: established / numerical / approximate / model / speculative
    script: str = "scripts/generate_all_figures.py"


FIGURES: dict[str, Callable[[], FigureResult]] = {}


def figure(name: str):
    def deco(fn):
        FIGURES[name] = fn
        return fn
    return deco


def save_figure(res: FigureResult, outdir: str | Path) -> dict[str, str]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / f"{res.name}.png"
    pdf = outdir / f"{res.name}.pdf"
    res.fig.savefig(png, bbox_inches="tight")
    res.fig.savefig(pdf, bbox_inches="tight", metadata={"CreationDate": None, "Creator": "cosmic_evolution"})
    plt.close(res.fig)
    return {"png": str(png), "pdf": str(pdf)}


def write_dataset_csv(name: str, data: dict, outdir: str | Path, header_lines: list[str]) -> Path | None:
    import numpy as np
    if not data:
        return None
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cols = {k: np.atleast_1d(np.asarray(v, dtype=float)) for k, v in data.items()}
    # group columns by length: write one CSV per distinct length to keep tables rectangular
    groups: dict[int, dict] = {}
    for k, v in cols.items():
        groups.setdefault(len(v), {})[k] = v
    last = None
    for i, (n, g) in enumerate(sorted(groups.items(), key=lambda kv: -kv[0])):
        path = outdir / (f"{name}.csv" if i == 0 else f"{name}__{i}.csv")
        with open(path, "w", encoding="utf-8") as fh:
            for line in header_lines:
                fh.write(f"# {line}\n")
            fh.write(",".join(g.keys()) + "\n")
            for row in zip(*g.values()):
                fh.write(",".join(f"{x:.10e}" for x in row) + "\n")
        last = path
    return last


def label_line(ax, x, y, text, color=INK2, dx=0, dy=0, **kw):
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", color=color, fontsize=8.5, **kw)
