import math

import numpy as np
import pytest

from cosmic_evolution.background.friedmann import FLRWBackground
from cosmic_evolution.config import CosmologyParams
from cosmic_evolution.early_universe import inflation as inf
from cosmic_evolution.perturbations import growth_factor as gf
from cosmic_evolution.perturbations import primordial_spectrum as ps
from cosmic_evolution.perturbations import scalar_modes as sm
from cosmic_evolution.perturbations import tensor_modes as tm
from cosmic_evolution.structure import halo_models as hm
from cosmic_evolution.structure import power_spectrum as pk
from cosmic_evolution.structure import transfer_functions as tf

from conftest import no_rad


# ---------------------------------------------------------------- mode equations
@pytest.mark.parametrize("N_exit", [10.0, 30.0])
def test_scalar_de_sitter_limit(N_exit):
    b = sm.ModeBackground.quasi_de_sitter(1e-4)
    k = sm.k_of_exit(b, N_exit)
    assert sm.evolve_curvature_mode(b, k) / sm.slow_roll_amplitude(b, N_exit) == pytest.approx(1.0, abs=1e-3)


def test_tensor_de_sitter_limit():
    b = sm.ModeBackground.quasi_de_sitter(1e-4)
    k = sm.k_of_exit(b, 20.0)
    assert tm.evolve_tensor_mode(b, k) / tm.slow_roll_tensor_amplitude(b, 20.0) == pytest.approx(1.0, abs=1e-3)


def test_mode_integration_converges_with_initial_depth():
    b = sm.ModeBackground.quasi_de_sitter(1e-3)
    k = sm.k_of_exit(b, 20.0)
    p1 = sm.evolve_curvature_mode(b, k, x_start=150.0)
    p2 = sm.evolve_curvature_mode(b, k, x_start=600.0)
    assert p1 == pytest.approx(p2, rel=3e-3)


def test_power_law_inflation_spectral_index():
    """Constant eps: P ~ H^2/eps exact slope n_s - 1 = -2 eps/(1-eps) at leading order."""
    eps = 0.01
    b = sm.ModeBackground.quasi_de_sitter(eps)
    ks = [sm.k_of_exit(b, n) for n in (20.0, 21.0)]
    P = [sm.evolve_curvature_mode(b, k) for k in ks]
    ns_m1 = math.log(P[1] / P[0]) / math.log(ks[1] / ks[0])
    assert ns_m1 == pytest.approx(-2 * eps / (1 - eps), rel=0.02)


def test_scalar_spectrum_full_model_consistency():
    S = ps.compute_spectra(inf.StarobinskyModel(), n_k=5, delta_N=4.0)
    assert S.n_s_numeric == pytest.approx(0.9645, abs=0.004)
    assert 0.003 < S.r_numeric < 0.004
    assert np.all(np.abs(S.P_s / S.P_s_slowroll - 1) < 0.04)             # next-to-leading-order correction
    assert np.all(np.abs(S.P_t / S.P_t_slowroll - 1) < 1e-3)
    assert S.A_s_numeric == pytest.approx(2.1e-9, rel=0.05)


def test_bad_wavenumber_rejected():
    b = sm.ModeBackground.quasi_de_sitter(1e-3, N_total=10.0)
    with pytest.raises(ValueError):
        sm.evolve_curvature_mode(b, 1e300)


# ---------------------------------------------------------------- growth
def test_growth_eds_exact():
    bg = no_rad(omega_m=1.0, omega_de=0.0, omega_k=0.0)
    g = gf.solve_growth(bg, a_i=1e-3)
    assert np.allclose(g.D_plus, g.a, rtol=1e-6)
    assert np.allclose(g.f, 1.0, atol=1e-6)


def test_growth_matches_heath_integral():
    bg = no_rad(omega_m=0.3, omega_de=0.7, omega_k=0.0)
    g = gf.solve_growth(bg, a_i=1e-4, n=200)
    for a in (0.2, 0.5, 1.0):
        assert g.D(a) == pytest.approx(gf.heath_integral_growth(bg, a), rel=1e-7)


def test_growth_matches_heath_for_open_universe():
    bg = no_rad(omega_m=0.3, omega_de=0.0)
    # initial conditions neglect curvature at a_i: error ~ a_i Omega_k/Omega_m, so start earlier
    g = gf.solve_growth(bg, a_i=1e-6, n=200)
    assert g.D(1.0) == pytest.approx(gf.heath_integral_growth(bg, 1.0), rel=1e-5)
    g4 = gf.solve_growth(bg, a_i=1e-4, n=200)
    assert abs(g4.D(1.0) / gf.heath_integral_growth(bg, 1.0) - 1) == pytest.approx(1.9e-4, rel=0.3)   # documented IC error


def test_growth_suppressed_by_dark_energy(bg):
    g = gf.solve_growth(bg)
    assert g.D_plus[-1] / g.a[-1] == pytest.approx(0.788, abs=0.01)         # regression
    assert np.all(np.diff(g.D_plus) > 0)
    om = bg.omega_fractions(g.a)["matter"]
    sel = g.a > 0.3
    assert np.allclose(g.f[sel], om[sel] ** 0.55, rtol=0.02)                # Linder gamma ~ 0.55


def test_growth_matter_era_amplitude_meszaros():
    """In the radiation era the growing mode is D ~ ln a-like then 1 + 3y/2; check the normalised solution -> a."""
    bg = FLRWBackground(CosmologyParams())
    g = gf.solve_growth(bg, a_i=1e-4)
    i = np.argmin(abs(g.a - 0.02))
    assert g.D_plus[i] / g.a[i] == pytest.approx(1.0, rel=0.02)


def test_growth_invalid_range(bg):
    with pytest.raises(ValueError):
        gf.solve_growth(bg, a_i=2.0, a_f=1.0)


# ---------------------------------------------------------------- transfer functions, power spectrum, halos
def test_transfer_functions_limits(params):
    k = np.array([1e-5, 1e-4])
    assert np.allclose(tf.eh_nowiggle(k, params), 1.0, atol=2e-3)
    assert np.allclose(tf.bbks(k, params), 1.0, atol=2e-3)
    kk = np.logspace(-3, 0, 20)
    ratio = tf.bbks(kk, params) / tf.eh_nowiggle(kk, params)
    assert np.all(np.abs(ratio - 1) < 0.06)
    assert np.all(np.diff(tf.eh_nowiggle(np.logspace(-3, 1, 50), params)) < 0)


def test_sigma8_in_expected_range(bg):
    g = gf.solve_growth(bg)
    s8 = pk.sigma8(bg, growth=g)
    # approximate fitting-formula pipeline: must be within ~5% of the Planck-quoted 0.811 (documented in VALIDATION.md)
    assert s8 == pytest.approx(0.811, rel=0.05)


def test_power_spectrum_growth_scaling(bg):
    g = gf.solve_growth(bg)
    k = np.array([0.01, 0.1])
    p1 = pk.matter_power(k, bg, 1.0, growth=g)
    p05 = pk.matter_power(k, bg, 0.5, growth=g)
    d = float(g.D(0.5) / g.D(1.0))
    assert np.allclose(p05 / p1, d**2, rtol=1e-12)


def test_large_scale_power_spectrum_slope(bg):
    k = np.array([1e-4, 2e-4])
    P = pk.matter_power(k, bg)
    assert math.log(P[1] / P[0]) / math.log(2.0) == pytest.approx(bg.params.n_s, abs=0.01)


def test_tophat_window_limits():
    assert pk.tophat_window(np.array([1e-6]))[0] == pytest.approx(1.0)
    assert pk.tophat_window(np.array([50.0]))[0] == pytest.approx(3 * (math.sin(50) - 50 * math.cos(50)) / 50**3)


def test_sheth_tormen_normalisation():
    assert hm.ST_A == pytest.approx(hm.st_normalisation_analytic(), abs=1e-4)
    nu = np.logspace(-14, 1.5, 400000)     # small-nu tail ~ nu^0.4 converges slowly, hence the deep lower limit
    from scipy.integrate import trapezoid
    assert trapezoid(hm.f_sheth_tormen(nu), np.log(nu)) == pytest.approx(1.0, abs=2e-3)
    assert trapezoid(hm.f_press_schechter(nu), np.log(nu)) == pytest.approx(1.0, abs=1e-6)


def test_mass_function_conserves_mass(bg):
    g = gf.solve_growth(bg)
    from scipy.special import erfc
    M = np.logspace(2, 17, 301)
    dn, sig = hm.mass_function(M, bg, "PS", 1.0, growth=g)
    rho = hm.mean_matter_density_msun_mpc3(bg)
    frac = np.trapezoid(dn * M / rho, np.log(M))
    nu = hm.DELTA_C / sig
    # mass fraction in the finite mass range equals the analytic PS integral over the corresponding nu range
    expected = erfc(nu.min() / math.sqrt(2)) - erfc(nu.max() / math.sqrt(2))
    assert frac == pytest.approx(expected, abs=3e-3)
    assert 0.8 < frac < 1.0                              # CDM sigma(M) saturates at small M, so < 1 is expected


def test_mean_matter_density(bg):
    # rho_m = 2.775e11 h^2 Omega_m Msun/Mpc^3
    assert hm.mean_matter_density_msun_mpc3(bg) == pytest.approx(2.775e11 * bg.params.h**2 * bg.om_m, rel=2e-3)
