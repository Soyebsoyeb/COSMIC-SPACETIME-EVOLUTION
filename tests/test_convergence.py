import math

import numpy as np
import pytest

from cosmic_evolution import units
from cosmic_evolution.background import cosmic_time as ct
from cosmic_evolution.geometry import horizons
from cosmic_evolution.numerics import convergence as cv
from cosmic_evolution.numerics import integration as ni
from cosmic_evolution.numerics import uncertainty as unc


def test_finite_lower_cutoff_matches_analytic_prediction(bg):
    ref = ct.age_gyr(bg)
    a_min = np.array([1e-6, 1e-5])      # deep in the radiation era (a_eq = 2.9e-4), where the pure-radiation estimate holds
    err = ref - cv.age_vs_lower_cutoff(bg, a_min)
    pred = a_min**2 / (2 * math.sqrt(bg.om_r)) * ct.hubble_time_gyr(bg)
    assert np.allclose(err, pred, rtol=3e-2)


def test_grid_rules_have_expected_orders(bg):
    ref = ct.age_gyr(bg)
    n = np.array([161, 321, 641, 1281])
    e_tr = np.abs(cv.age_vs_grid(bg, n, "trapezoid") - ref)
    e_si = np.abs(cv.age_vs_grid(bg, n, "simpson") - ref)
    s_tr = -np.polyfit(np.log(n), np.log(e_tr), 1)[0]
    s_si = -np.polyfit(np.log(n), np.log(e_si), 1)[0]
    assert s_tr == pytest.approx(2.0, abs=0.3)
    assert s_si > 3.3
    assert e_si[-1] < 1e-8 * ref


def test_tolerance_sweep_error_below_request(bg):
    ref = ct.age_gyr(bg)
    eps = np.array([1e-4, 1e-6, 1e-8, 1e-10])
    err = np.abs(cv.tolerance_sweep(bg, eps) - ref) / ref
    assert np.all(err <= eps)


def test_event_horizon_cutoff_tail_is_analytic(bg):
    ehr = horizons.event_horizon_mpc(bg)
    amax = np.array([1e3, 1e5])
    err = ehr - cv.horizon_vs_upper_cutoff(bg, amax)
    dh = units.hubble_distance_Mpc(bg.H0)
    assert np.allclose(err, dh / (math.sqrt(bg.om_de) * amax), rtol=1e-3)


def test_quad_checked_detects_nonfinite():
    with pytest.raises((FloatingPointError, RuntimeError, ZeroDivisionError)):
        ni.quad_checked(lambda x: 1.0 / x, 0.0, 1.0)


def test_quad_checked_value():
    r = ni.quad_checked(lambda x: x**2, 0.0, 3.0)
    assert r.value == pytest.approx(9.0, rel=1e-12)


def test_simpson_even_n_is_promoted_to_odd():
    v = ni.simpson_lna(lambda a: 1.0 / a**2, 1.0, math.e, 100)
    assert v == pytest.approx(1.0 - 1.0 / math.e, rel=1e-6)


def test_monte_carlo_is_reproducible_and_centred(params):
    f = lambda b: ct.age_gyr(b)
    a = unc.monte_carlo(f, params, n=40, seed=7)
    b = unc.monte_carlo(f, params, n=40, seed=7)
    assert np.array_equal(a, b)
    assert abs(a.mean() - ct.age_gyr(__import__("cosmic_evolution").FLRWBackground(params))) < 0.1


def test_log_sensitivity_of_age_to_h(params):
    s = unc.log_sensitivity(lambda b: ct.age_gyr(b), params, "h")
    assert s == pytest.approx(-1.0, abs=0.01)           # t0 ~ 1/H0 at fixed density parameters
