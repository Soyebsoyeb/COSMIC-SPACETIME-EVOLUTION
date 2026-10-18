import math

import numpy as np
import pytest

from cosmic_evolution import units
from cosmic_evolution.background.friedmann import FLRWBackground
from cosmic_evolution.config import CosmologyParams
from cosmic_evolution.geometry import curvature, flrw, horizons

from conftest import no_rad


def test_density_closure(params):
    s = params.omega_r + params.omega_m + params.omega_k_eff + params.omega_de_eff
    assert s == pytest.approx(1.0, abs=1e-14)


def test_invalid_parameters_raise():
    with pytest.raises(ValueError):
        CosmologyParams(h=-0.5)
    with pytest.raises(ValueError):
        CosmologyParams(omega_m=-0.1)
    with pytest.raises(ValueError):
        CosmologyParams(T_cmb=0.0)
    with pytest.raises(ValueError):
        CosmologyParams(omega_m=0.3, omega_de=0.9, omega_k=0.5)   # does not close
    with pytest.raises(ValueError):
        CosmologyParams(omega_m=0.04)                             # Omega_b = 0.049 > Omega_m


def test_eds_comoving_distance_analytic():
    bg = no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)
    dh = units.hubble_distance_Mpc(bg.H0)
    for z in (0.1, 1.0, 5.0):
        assert flrw.comoving_distance_mpc(bg, z) == pytest.approx(2 * dh * (1 - 1 / math.sqrt(1 + z)), rel=1e-8)


def test_milne_transverse_distance_analytic():
    bg = no_rad(omega_m=0.0, omega_de=0.0)          # empty open universe: Omega_k = 1
    dh = units.hubble_distance_Mpc(bg.H0)
    for z in (0.5, 2.0):
        assert flrw.comoving_distance_mpc(bg, z) == pytest.approx(dh * math.log(1 + z), rel=1e-8)
        assert flrw.transverse_comoving_distance_mpc(bg, z) == pytest.approx(dh * math.sinh(math.log(1 + z)), rel=1e-8)


def test_closed_universe_uses_sine():
    bg = no_rad(omega_m=2.0, omega_de=0.0)           # Omega_k = -1
    dh = units.hubble_distance_Mpc(bg.H0)
    dc = flrw.comoving_distance_mpc(bg, 1.0)
    assert flrw.transverse_comoving_distance_mpc(bg, 1.0) == pytest.approx(dh * math.sin(dc / dh), rel=1e-10)


def test_distance_duality(bg):
    z = np.array([0.1, 0.5, 2.0, 10.0])
    assert np.allclose(flrw.luminosity_distance_mpc(bg, z), (1 + z) ** 2 * flrw.angular_diameter_distance_mpc(bg, z), rtol=1e-12)


def test_distance_integral_consistency_with_independent_rule(bg):
    """Compare adaptive quadrature with a high-resolution Simpson rule in ln(1+z)."""
    from scipy.integrate import simpson
    z = 3.0
    x = np.linspace(0.0, math.log(1 + z), 4001)
    f = (1 + np.exp(x) - 1) ** 0 * np.exp(x) / bg.E_of_z(np.exp(x) - 1)
    ref = simpson(f, x=x) * units.hubble_distance_Mpc(bg.H0)
    assert flrw.comoving_distance_mpc(bg, z) == pytest.approx(ref, rel=1e-9)


def test_negative_redshift_rejected(bg):
    with pytest.raises(ValueError):
        flrw.comoving_distance_mpc(bg, -0.1)


def test_curvature_sign_and_radius():
    assert curvature.curvature_sign(no_rad(omega_m=0.3, omega_de=0.0)) == -1   # open
    assert curvature.curvature_sign(no_rad(omega_m=2.0, omega_de=0.0)) == 1    # closed
    assert curvature.curvature_sign(no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)) == 0
    assert math.isinf(curvature.curvature_radius_mpc(no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)))


def test_ricci_scalar_de_sitter():
    bg = no_rad(omega_m=0.0, omega_de=1.0, omega_k=0.0)     # pure de Sitter: R = 12 H^2
    assert curvature.ricci_scalar_over_H0sq(bg, 0.7) == pytest.approx(12.0, rel=1e-12)


def test_particle_horizon_eds_and_event_horizon_de_sitter():
    eds = no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)
    dh = units.hubble_distance_Mpc(eds.H0)
    assert horizons.particle_horizon_mpc(eds) == pytest.approx(2 * dh, rel=1e-8)
    assert math.isinf(horizons.event_horizon_mpc(eds))
    ds = no_rad(omega_m=0.0, omega_de=1.0, omega_k=0.0)
    assert horizons.event_horizon_mpc(ds) == pytest.approx(dh, rel=1e-7)


def test_fiducial_horizons_regression(bg):
    assert horizons.particle_horizon_mpc(bg) / 1e3 == pytest.approx(14.14, abs=0.05)   # Gpc, comoving
    assert horizons.event_horizon_mpc(bg) / 1e3 == pytest.approx(5.11, abs=0.05)
    assert horizons.event_horizon_exists(bg)


def test_hubble_radius_equals_c_over_h(bg):
    assert horizons.hubble_radius_mpc(bg, 1.0) == pytest.approx(units.hubble_distance_Mpc(bg.H0), rel=1e-12)
