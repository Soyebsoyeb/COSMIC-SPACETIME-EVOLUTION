import math

import numpy as np
import pytest
from scipy import integrate

from cosmic_evolution import constants as K
from cosmic_evolution.background import cosmic_time as ct
from cosmic_evolution.background.friedmann import FLRWBackground
from cosmic_evolution.config import CosmologyParams
from cosmic_evolution.future import de_sitter as ds
from cosmic_evolution.future import dark_energy_scenarios as des
from cosmic_evolution.future import fate_classification as fc

from conftest import no_rad


def test_fate_table_expected_outcomes():
    expected = {"lcdm": "asymptotic de Sitter", "quint": "eternal accelerated power-law expansion",
                "phantom12": "finite-time singularity (big rip)", "phantom15": "finite-time singularity (big rip)",
                "eds": "eternal decelerating expansion", "open": "eternal coasting expansion (q -> 0)",
                "closed": "recollapse (big crunch)", "negL": "recollapse (big crunch)"}
    rows = {r["key"]: r for r in fc.fate_table()}
    assert set(rows) == set(expected)
    for k, v in expected.items():
        assert rows[k]["fate"] == v, k


def test_pure_phantom_rip_time_integral_exact():
    # dark energy only, w = -3/2: int_1^inf da/(a H/H0) = 2/(3|1+w|) = 4/3  (Caldwell et al. 2003, Hubble units)
    bg = no_rad(omega_m=0.0, omega_de=1.0, omega_k=0.0, w0=-1.5)
    val, _ = integrate.quad(lambda a: a / math.sqrt(float(bg.Q(a))), 1.0, np.inf, epsabs=1e-13, epsrel=1e-11)
    assert val == pytest.approx(4.0 / 3.0, rel=1e-7)


def test_big_rip_estimate_close_to_numeric():
    sc = {s.key: s for s in des.default_scenarios()}
    bg = FLRWBackground(sc["phantom15"].params)
    rows = {r["key"]: r for r in fc.fate_table()}
    numeric = rows["phantom15"]["t_end_gyr"] - ct.age_gyr(bg)
    est = fc.big_rip_time_estimate_gyr(bg)
    # the estimate neglects matter near the rip: documented ~3% overestimate for w = -1.5, Omega_m = 0.315
    assert est == pytest.approx(numeric, rel=0.05)
    assert est > numeric
    with pytest.raises(ValueError):
        fc.big_rip_time_estimate_gyr(FLRWBackground(CosmologyParams()))


def test_closed_recollapse_analytic_with_negligible_radiation():
    bg = no_rad(omega_m=2.0, omega_de=0.0)
    res = des.classify_fate(bg)
    assert res["a_max"] == pytest.approx(2.0, rel=1e-9)
    assert res["t_end_gyr"] == pytest.approx(2 * math.pi * ct.hubble_time_gyr(bg), rel=1e-5)


def test_negative_lambda_always_recollapses():
    bg = FLRWBackground(CosmologyParams(omega_m=0.315, omega_de=-0.5))
    assert ct.a_turnaround(bg) is not None


def test_ode_matches_quadrature_for_all_expanding_scenarios():
    for s in des.default_scenarios():
        bg = FLRWBackground(s.params)
        amax = ct.a_turnaround(bg)
        if amax is not None:
            continue
        fate = des.classify_fate(bg)
        if fate["t_end_gyr"] is not None:
            continue
        sol = des.integrate_scale_factor_ode(bg, 0.5)
        a = sol.sol(0.5)[0]
        dt = ct.time_dimensionless(bg, a) - ct.time_dimensionless(bg, 1.0)
        assert dt == pytest.approx(0.5, rel=1e-7), s.key


def test_trajectory_monotonic_and_endpoints():
    sc = {s.key: s for s in des.default_scenarios()}
    tr = des.trajectory(FLRWBackground(sc["lcdm"].params), a_start=0.2, n=100)
    assert np.all(np.diff(tr.t_gyr) > 0) and np.all(np.diff(tr.a) > 0)
    tr = des.trajectory(FLRWBackground(sc["closed"].params), a_start=0.2, n=100)
    assert tr.a.max() == pytest.approx(2.0, rel=1e-3)
    assert tr.H_over_H0[-1] < 0 < tr.H_over_H0[0]


def test_de_sitter_quantities(bg):
    hl = ds.hubble_lambda_inv_s(bg)
    assert hl == pytest.approx(bg.H0 * 1e3 / K.Mpc * math.sqrt(bg.om_de), rel=1e-12)
    assert ds.proper_event_horizon_gly(bg) == pytest.approx(17.53, abs=0.02)
    assert ds.gibbons_hawking_temperature_K(bg) == pytest.approx(K.hbar * hl / (2 * math.pi * K.k_B), rel=1e-14)
    assert ds.gibbons_hawking_temperature_K(bg) < 1e-29
    with pytest.raises(ValueError):
        ds.hubble_lambda_inv_s(FLRWBackground(CosmologyParams(w0=-0.8)))


def test_lcdm_hubble_rate_approaches_h_lambda(bg):
    assert bg.E(1e4) == pytest.approx(math.sqrt(bg.om_de), rel=1e-6)


def test_eds_asymptotic_power_law():
    bg = no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)
    t = [ct.time_dimensionless(bg, a) for a in (10.0, 20.0)]
    assert math.log(20.0 / 10.0) / math.log(t[1] / t[0]) == pytest.approx(2.0 / 3.0, rel=1e-6)   # a ~ t^{2/3}
