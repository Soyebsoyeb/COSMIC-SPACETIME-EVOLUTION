"""Summary table of conditional fates for the default scenario set."""
from __future__ import annotations

from ..background.friedmann import FLRWBackground
from .dark_energy_scenarios import Scenario, classify_fate, default_scenarios


def fate_table(scenarios: list[Scenario] | None = None) -> list[dict]:
    rows = []
    for sc in scenarios or default_scenarios():
        bg = FLRWBackground(sc.params)
        res = classify_fate(bg)
        rows.append({"key": sc.key, "label": sc.label, "fate": res["fate"],
                     "t_end_gyr": res.get("t_end_gyr"), "a_max": res.get("a_max"),
                     "Omega_k": bg.om_k, "Omega_de": bg.om_de, "w0": sc.params.w0})
    return rows


def big_rip_time_estimate_gyr(bg: FLRWBackground) -> float:
    """Caldwell-Kamionkowski-Weinberg estimate t_rip - t0 ~ (2/3) / (|1+w| H0 sqrt(1 - Om_m)), neglecting matter
    at late times; approximate and valid for w < -1 with Om_m treated as negligible near the rip."""
    from ..background.cosmic_time import hubble_time_gyr
    w = bg.params.w0
    if w >= -1.0:
        raise ValueError("big rip estimate requires w < -1")
    return (2.0 / 3.0) / (abs(1.0 + w) * (1.0 - bg.om_m) ** 0.5) * hubble_time_gyr(bg)
