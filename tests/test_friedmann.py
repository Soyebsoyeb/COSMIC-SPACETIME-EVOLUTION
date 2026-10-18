import numpy as np
import pytest

from cosmic_evolution.background import dark_energy as de
from cosmic_evolution.background import expansion_history as eh
from cosmic_evolution.background import scale_factor as sf
from cosmic_evolution.background.friedmann import FLRWBackground
from cosmic_evolution.config import CosmologyParams

from conftest import no_rad


def test_present_day_normalisation(bg):
    assert bg.E(1.0) == pytest.approx(1.0, abs=1e-14)
    assert bg.H(1.0) == pytest.approx(bg.H0, rel=1e-14)


def test_fractions_sum_to_one(bg):
    a = np.logspace(-8, 3, 200)
    s = sum(bg.omega_fractions(a).values())
    assert np.allclose(s, 1.0, atol=1e-12)


def test_radiation_and_matter_scaling(bg):
    a = np.array([1e-8 * 0.999, 1e-8 * 1.001])
    slope = np.diff(np.log(bg.E2(a))) / np.diff(np.log(a))
    assert slope[0] == pytest.approx(-4.0, abs=1e-3)            # radiation era
    a = np.array([0.03 * 0.999, 0.03 * 1.001])
    slope = np.diff(np.log(bg.E2(a))) / np.diff(np.log(a))
    assert slope[0] == pytest.approx(-3.0, abs=0.05)            # matter era (radiation ~1%)


def test_lambda_density_constant():
    bg = FLRWBackground(CosmologyParams())
    d = bg.density_terms(np.array([1e-3, 1.0, 1e3]))["dark_energy"]
    assert np.allclose(d, d[0], rtol=1e-14)


def test_scaling_laws_exact_for_single_components():
    r = no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)
    assert r.E2(0.5) == pytest.approx(2.0 ** 3, rel=1e-12)       # matter: a^-3
    k = no_rad(omega_m=0.0, omega_de=0.0)
    assert k.E2(0.5) == pytest.approx(2.0 ** 2, rel=1e-12)       # curvature: a^-2


def test_deceleration_limits(bg):
    assert bg.deceleration(1e-9) == pytest.approx(1.0, abs=1e-3)
    assert bg.deceleration(0.05) == pytest.approx(0.5, abs=0.03)
    assert bg.deceleration(1e4) == pytest.approx(-1.0, abs=1e-9)


def test_deceleration_matches_numerical_second_derivative(bg):
    """q = -1 - dlnH/dlna  via finite differences of the independently computed H."""
    for a in (0.01, 0.3, 0.8, 2.0):
        h = 1e-5
        dlnH = (np.log(bg.E(a * (1 + h))) - np.log(bg.E(a * (1 - h)))) / (np.log(1 + h) - np.log(1 - h))
        assert bg.deceleration(a) == pytest.approx(-1.0 - dlnH, abs=1e-8)


def test_matter_radiation_equality(bg):
    a_eq = bg.a_equality_mr()
    assert bg.density_terms(a_eq)["radiation"] == pytest.approx(bg.density_terms(a_eq)["matter"], rel=1e-12)
    assert 1 / a_eq - 1 == pytest.approx(3419.0, rel=2e-3)       # regression: massless-nu approximation


def test_acceleration_onset_matches_analytic(bg):
    a_num = eh.acceleration_onset_a(bg)
    a_an = (bg.om_m / (2 * bg.om_de)) ** (1 / 3)               # radiation neglected
    assert a_num == pytest.approx(a_an, rel=1e-3)
    assert bg.deceleration(a_num) == pytest.approx(0.0, abs=1e-10)


def test_dark_energy_cpl_density_solves_continuity():
    w0, wa = -0.9, 0.3
    a = np.array([0.4, 0.5])
    h = 1e-6
    rho = lambda x: de.de_density_factor(x, w0, wa)
    dln = (np.log(rho(0.5 * (1 + h))) - np.log(rho(0.5 * (1 - h)))) / (np.log(1 + h) - np.log(1 - h))
    assert dln == pytest.approx(-3 * (1 + de.w_of_a(0.5, w0, wa)), rel=1e-8)
    assert de.de_density_factor(1.0, w0, wa) == pytest.approx(1.0)


def test_constant_w_density_power_law():
    assert de.de_density_factor(0.5, -0.8, 0.0) == pytest.approx(0.5 ** (-3 * 0.2), rel=1e-14)


def test_negative_e2_raises():
    bg = no_rad(omega_m=2.0, omega_de=0.0)       # closed: E^2 < 0 beyond a = 2
    with pytest.raises(ValueError):
        bg.E(3.0)


def test_redshift_scale_factor_roundtrip():
    z = np.array([0.0, 0.5, 10.0, 1100.0])
    assert np.allclose(sf.z_of_a(sf.a_of_z(z)), z, rtol=1e-14, atol=1e-14)
    with pytest.raises(ValueError):
        sf.a_of_z(-1.0)
    with pytest.raises(ValueError):
        sf.z_of_a(0.0)


def test_h_of_z_equals_h_of_a(bg):
    z = np.array([0.0, 1.0, 100.0])
    assert np.allclose(bg.H_of_z(z), bg.H(1 / (1 + z)), rtol=1e-14)
