"""Create the four tutorial notebooks (unexecuted). Run scripts/execute_notebooks.py afterwards."""
from __future__ import annotations

from pathlib import Path

import nbformat as nbf

OUT = Path(__file__).resolve().parents[1] / "notebooks"


def md(s):
    return nbf.v4.new_markdown_cell(s)


def code(s):
    return nbf.v4.new_code_cell(s)


NOTEBOOKS = {
"01_flrw_background.ipynb": [
md("# FLRW background\n\nAge, expansion history and horizons of the fiducial flat LCDM model. Label: **numerical** (exact quadrature of the Friedmann equation)."),
code("import numpy as np\nimport matplotlib.pyplot as plt\nfrom cosmic_evolution.config import CosmologyParams\nfrom cosmic_evolution.background.friedmann import FLRWBackground\nfrom cosmic_evolution.background.cosmic_time import age_gyr\nfrom cosmic_evolution.background.expansion_history import expansion_table\nfrom cosmic_evolution.geometry import horizons\n\nbg = FLRWBackground(CosmologyParams())\nprint('age [Gyr]:', age_gyr(bg))\nprint('particle horizon today [Mpc]:', horizons.particle_horizon_mpc(bg))\nprint('event horizon today [Mpc]:', horizons.event_horizon_mpc(bg))"),
md("## Exact checks\n\nIn an Einstein-de Sitter universe $H_0 t_0 = 2/3$ exactly."),
code("eds = FLRWBackground(CosmologyParams(omega_m=1.0 - CosmologyParams().omega_r, h=0.7))\nfrom cosmic_evolution.background.cosmic_time import time_dimensionless\nprint(time_dimensionless(eds, 1.0), 2/3, '(small difference is the radiation density)')"),
md("## Expansion history"),
code("z = np.logspace(-2, 3, 200)\nT = expansion_table(bg, z)\nfig, ax = plt.subplots(1, 2, figsize=(10, 3.5))\nax[0].loglog(1 + z, T['H_km_s_Mpc']); ax[0].set_xlabel('1+z'); ax[0].set_ylabel('H [km/s/Mpc]')\nax[1].semilogx(1 + z, T['q']); ax[1].axhline(0, color='k', lw=0.5); ax[1].set_xlabel('1+z'); ax[1].set_ylabel('q')\nplt.tight_layout(); plt.show()"),
],
"02_thermal_history.ipynb": [
md("# Thermal history\n\nRelativistic degrees of freedom and the time-temperature relation. Label: **approximate** (step function at T = m, instantaneous decoupling)."),
code("import numpy as np\nimport matplotlib.pyplot as plt\nfrom cosmic_evolution.config import CosmologyParams\nfrom cosmic_evolution.early_universe import thermal_history as th\n\nT = np.logspace(-6, 3, 400)\nplt.semilogx(T, th.g_star(T), label='g*')\nplt.semilogx(T, th.g_star_s(T), '--', label='g*s')\nplt.xlabel('T [GeV]'); plt.ylabel('degrees of freedom'); plt.legend(); plt.show()"),
code("for name, Tg, t in th.epoch_table(CosmologyParams()):\n    print(f'{name:50s} T = {Tg:9.3g} GeV   t = {t:9.3g} s')"),
md("BBN is **not computed** in this repository; only the textbook relations $n/p=e^{-Q/T}$ and $Y_p=2x/(1+x)$ are provided, with $x$ an input."),
code("from cosmic_evolution.early_universe import bbn_interface as b\nprint([n for n in dir(b) if not n.startswith('_')])"),
],
"03_inflation_models.ipynb": [
md("# Inflation models\n\nExact e-fold background evolution for two potentials. Label: **model dependent**."),
code("import numpy as np\nimport matplotlib.pyplot as plt\nfrom cosmic_evolution.early_universe import inflation as inf\n\nfor model in (inf.StarobinskyModel(), inf.QuadraticModel()):\n    sol = inf.solve_background(model)\n    n = 55.0\n    phi = sol.phi_at(n)\n    sr = inf.slow_roll_observables(model, phi)\n    hs = inf.hubble_slow_roll_observables(sol, n)\n    print(model.name, 'N_total = %.1f' % sol.N_total, 'n_s = %.4f' % hs['n_s'], 'r = %.4f' % hs['r'], '| potential slow roll:', round(sr['n_s'], 4), round(sr['r'], 4))"),
code("sol = inf.solve_background(inf.StarobinskyModel())\nplt.semilogy(sol.N_before_end(), sol.eps_H); plt.gca().invert_xaxis()\nplt.xlabel('e-folds before end'); plt.ylabel('epsilon_H'); plt.show()"),
md("The quadratic potential gives $r\\approx0.14$, above the observational upper limit of order $0.06$: the code reproduces this known tension, it does not hide it."),
],
"04_future_scenarios.ipynb": [
md("# Future scenarios\n\nEvery result is **conditional on the assumed dark-energy model**."),
code("from cosmic_evolution.future.fate_classification import fate_table\nfor row in fate_table():\n    t = row['t_end_gyr']\n    print(f\"{row['label']:50s} {row['fate']:45s} t_end = {'none' if t is None else round(t, 1)}\")"),
code("import numpy as np\nimport matplotlib.pyplot as plt\nfrom cosmic_evolution.config import CosmologyParams\nfrom cosmic_evolution.background.friedmann import FLRWBackground\nfrom cosmic_evolution.future.fate_classification import big_rip_time_estimate_gyr\nbg = FLRWBackground(CosmologyParams(w0=-1.2))\nprint('big-rip time estimate [Gyr] (approximate):', big_rip_time_estimate_gyr(bg))"),
],
}

if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for name, cells in NOTEBOOKS.items():
        nb = nbf.v4.new_notebook()
        nb.cells = cells
        nb.metadata["kernelspec"] = {"display_name": "Python 3", "language": "python", "name": "python3"}
        nbf.write(nb, OUT / name)
        print("wrote", name)
