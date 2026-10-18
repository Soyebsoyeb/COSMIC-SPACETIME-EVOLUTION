import math

import numpy as np
import pytest
from scipy import integrate

from cosmic_evolution import units
from cosmic_evolution.background import cosmic_time as ct
from cosmic_evolution.future import dark_energy_scenarios as des

from conftest import no_rad


def test_eds_age():
    bg = no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)
    assert ct.time_dimensionless(bg, 1.0) == pytest.approx(2.0 / 3.0, rel=1e-9)


def test_flat_lcdm_age_analytic():
    bg = no_rad(omega_m=0.3, omega_de=0.7, omega_k=0.0)
    ana = 2.0 / (3.0 * math.sqrt(0.7)) * math.asinh(math.sqrt(0.7 / 0.3))
    assert ct.time_dimensionless(bg, 1.0) == pytest.approx(ana, rel=1e-9)


def test_milne_age():
    bg = no_rad(omega_m=0.0, omega_de=0.0)
    assert ct.time_dimensionless(bg, 1.0) == pytest.approx(1.0, rel=1e-9)


def test_hubble_time_conversion():
    # 1/H0 = 9.778 Gyr / h  (derived from the Julian year and IAU Mpc)
    assert units.hubble_time_Gyr(100.0) == pytest.approx(9.7779, rel=1e-4)


def test_fiducial_age_regression(bg):
    assert ct.age_gyr(bg) == pytest.approx(13.79, abs=0.02)


def test_age_integration_converges(bg):
    ref = ct.age_gyr(bg)
    for eps in (1e-6, 1e-8, 1e-10):
        val = ct.time_dimensionless(bg, 1.0, epsabs=0.0, epsrel=eps) * ct.hubble_time_gyr(bg)
        assert val == pytest.approx(ref, rel=max(10 * eps, 1e-12))


def test_cosmic_time_ode_cross_check(bg):
    """Second-order acceleration equation (independent of the Friedmann constraint) reproduces t(a)."""
    sol = des.integrate_scale_factor_ode(bg, 0.5)
    a_half = sol.sol(0.5)[0]
    dt = ct.time_dimensionless(bg, a_half) - ct.time_dimensionless(bg, 1.0)
    assert dt == pytest.approx(0.5, rel=1e-8)


def test_lookback_plus_cosmic_time_equals_age(bg):
    z = 2.0
    assert ct.lookback_time_gyr(bg, z) + ct.cosmic_time_gyr(bg, 1 / (1 + z)) == pytest.approx(ct.age_gyr(bg), rel=1e-10)


def test_time_monotonic_in_a(bg):
    t = ct.cosmic_time_gyr(bg, np.logspace(-6, 1, 50))
    assert np.all(np.diff(t) > 0)


def test_radiation_era_time_formula(bg):
    a = 1e-7
    assert ct.time_dimensionless(bg, a) == pytest.approx(a**2 / (2 * math.sqrt(bg.om_r)), rel=2e-3)


def test_turnaround_closed_matter_universe():
    bg = no_rad(omega_m=2.0, omega_de=0.0)
    assert ct.a_turnaround(bg) == pytest.approx(2.0, rel=1e-9)
    fate = des.classify_fate(bg)
    assert fate["t_end_gyr"] / ct.hubble_time_gyr(bg) == pytest.approx(2 * math.pi, rel=1e-5)  # crunch at 2 pi / H0


def test_no_turnaround_for_lcdm(bg):
    assert ct.a_turnaround(bg) is None


def test_invalid_scale_factor_rejected(bg):
    with pytest.raises(ValueError):
        ct.time_dimensionless(bg, -1.0)


def test_conformal_time_equals_particle_horizon(bg):
    from cosmic_evolution.geometry import horizons
    assert ct.conformal_time_mpc(bg, 1.0) == pytest.approx(horizons.particle_horizon_mpc(bg), rel=1e-12)
