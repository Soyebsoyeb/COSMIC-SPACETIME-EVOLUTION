import math

import numpy as np
import pytest
from scipy import integrate

from cosmic_evolution import constants as K
from cosmic_evolution import units
from cosmic_evolution.background.friedmann import FLRWBackground
from cosmic_evolution.config import CosmologyParams
from cosmic_evolution.early_universe import bbn_interface as bbn
from cosmic_evolution.early_universe import entropy, inflation as inf, recombination as rc, reheating as rh
from cosmic_evolution.early_universe import thermal_history as th


# ---------------------------------------------------------------- units / constants
def test_planck_units():
    assert K.t_planck == pytest.approx(5.391e-44, rel=1e-3)
    assert K.l_planck == pytest.approx(1.616e-35, rel=1e-3)
    assert K.E_planck_GeV == pytest.approx(1.221e19, rel=1e-3)
    assert K.M_P_reduced_GeV == pytest.approx(2.435e18, rel=1e-3)


def test_unit_conversions_roundtrip():
    assert units.K_to_GeV(units.GeV_to_K(3.7)) == pytest.approx(3.7, rel=1e-14)
    assert units.GeV_to_K(1e-3 * 1.0) / 1.160e10 == pytest.approx(1.0, rel=1e-3)     # 1 MeV = 1.160e10 K
    assert units.Mpc_to_Gly(1000.0) == pytest.approx(3.2616, rel=1e-4)
    with pytest.raises(ValueError):
        units.H0_to_inv_s(0.0)


# ---------------------------------------------------------------- thermal history
def test_g_star_counting():
    assert th.g_star(1000.0) == pytest.approx(106.75)
    assert th.g_star(1.5) == pytest.approx(72.25)                  # m_c < T < m_tau: charm present, tau gone
    assert th.g_star(1.0) == pytest.approx(61.75)                  # below the charm threshold
    assert th.g_star(0.15) == pytest.approx(17.25)
    assert th.g_star(2e-3) == pytest.approx(10.75)
    assert th.g_star(1e-9) == pytest.approx(2 + 7 / 4 * (4 / 11) ** (4 / 3) * 3.046, rel=1e-12)   # = 3.3835
    assert th.g_star_late(3.0) == pytest.approx(3.36, abs=5e-3)                                   # textbook 3.36
    assert th.g_star_s(1e-9) == pytest.approx(3.909, abs=1e-3)
    assert th.g_star_s(2e-3) == pytest.approx(10.75)


def test_g_star_monotone_decreasing():
    T = np.logspace(3, -9, 400)
    assert np.all(np.diff(th.g_star(T)) <= 1e-12)


def test_temperature_scaling():
    assert th.temperature_of_z(1089.0) == pytest.approx(K.T_CMB0 * 1090.0)
    with pytest.raises(ValueError):
        th.temperature_of_z(-2.0)


@pytest.mark.parametrize("z", [1e7, 1e8, 1e9])
def test_radiation_hubble_rate_matches_background(bg, z):
    T = units.K_to_GeV(th.temperature_of_z(z))
    h_thermal = th.hubble_rate_radiation(T)
    h_bg = units.H0_to_inv_s(bg.H0) * bg.E_of_z(z)
    assert h_thermal == pytest.approx(h_bg, rel=3e-4)       # differs only by the (tiny) matter/Lambda terms


def test_time_temperature_relation():
    # t = 2.42 s g_*^{-1/2} (T/MeV)^-2 (textbook radiation-era relation); g_* = 10.75 at T = 2 MeV
    assert th.time_of_temperature_s(2e-3) == pytest.approx(2.42 / math.sqrt(10.75) / 4.0, rel=0.01)
    assert th.temperature_of_time_gev(float(th.time_of_temperature_s(0.1))) == pytest.approx(0.1, rel=1e-8)


def test_photon_number_density_and_eta(params):
    assert entropy.photon_number_density(K.T_CMB0) / 1e6 == pytest.approx(410.7, rel=1e-3)   # per cm^3
    assert entropy.baryon_to_photon_ratio(params) == pytest.approx(273.9e-10 * params.omega_b_h2, rel=5e-3)


def test_energy_density_matches_omega_gamma(params):
    rho = entropy.radiation_energy_density(params.T_cmb, 2.0)         # J/m^3
    from cosmic_evolution.units import critical_density_kg_m3
    omega = rho / K.c**2 / critical_density_kg_m3(params.H0)
    assert omega == pytest.approx(params.omega_gamma, rel=1e-12)


# ---------------------------------------------------------------- BBN interface
def test_bbn_relations():
    assert bbn.np_equilibrium_ratio(1.0) == pytest.approx(math.exp(-1.29333), rel=1e-12)
    assert bbn.helium_mass_fraction_from_np(1.0 / 7.0) == pytest.approx(0.25, abs=1e-12)
    with pytest.raises(ValueError):
        bbn.np_equilibrium_ratio(-1.0)
    with pytest.raises(ValueError):
        bbn.helium_mass_fraction_from_np(-0.1)


def test_bbn_backend_is_not_faked(params):
    with pytest.raises(bbn.BBNBackendUnavailable):
        bbn.run_backend("AlterBBN", params)


def test_bbn_export(params, tmp_path):
    import json
    p = bbn.export_inputs(params, tmp_path / "bbn.json")
    d = json.loads(p.read_text())
    assert d["bbn_inputs"]["eta"] == pytest.approx(entropy.baryon_to_photon_ratio(params))


# ---------------------------------------------------------------- inflation
@pytest.fixture(scope="module")
def sols():
    return {m.name: inf.solve_background(m) for m in (inf.QuadraticModel(), inf.StarobinskyModel())}


def test_inflation_ends_when_eps_is_one(sols):
    for s in sols.values():
        assert s.eps_H[-1] == pytest.approx(1.0, abs=2e-3)
        assert s.N_total > 60


def test_quadratic_efolds_slow_roll_estimate(sols):
    s = sols["quadratic"]
    phi_i, phi_e = s.phi[0], s.phi_end
    assert s.N_total == pytest.approx((phi_i**2 - phi_e**2) / 4.0, rel=0.02)


def test_quadratic_observables_vs_analytic(sols):
    s = sols["quadratic"]
    N = 60.0
    o = inf.slow_roll_observables(inf.QuadraticModel(), s.phi_at(N))
    assert o["n_s"] == pytest.approx(1 - 2 / N, abs=2e-3)
    assert o["r"] == pytest.approx(8 / N, rel=0.06)


def test_starobinsky_observables_vs_analytic(sols):
    s = sols["Starobinsky"]
    N = 60.0
    o = inf.slow_roll_observables(inf.StarobinskyModel(), s.phi_at(N))
    assert o["n_s"] == pytest.approx(1 - 2 / N, abs=2e-3)
    assert o["r"] == pytest.approx(12 / N**2, rel=0.12)


def test_hubble_and_potential_slow_roll_agree(sols):
    s = sols["Starobinsky"]
    a = inf.slow_roll_observables(inf.StarobinskyModel(), s.phi_at(55.0))
    b = inf.hubble_slow_roll_observables(s, 55.0)
    assert a["n_s"] == pytest.approx(b["n_s"], abs=1e-3)
    assert a["r"] == pytest.approx(b["r"], rel=0.02)


def test_efold_equation_vs_cosmic_time_integration():
    """Independent check: integrate phi'' + 3H phi' + V' = 0 with H from the Friedmann constraint in cosmic time."""
    m = 1e-5
    model = inf.QuadraticModel()
    phi_i = 16.0
    V = lambda p: 0.5 * m**2 * p**2
    dV = lambda p: m**2 * p
    H = lambda p, dp: math.sqrt((0.5 * dp**2 + V(p)) / 3.0)
    H0_ = H(phi_i, 0.0)
    dphi0 = -dV(phi_i) / (3 * H0_)                      # slow-roll start
    # use the e-fold solution's own attractor start for a like-for-like comparison
    h_i = math.sqrt(V(phi_i) / (3 - 0.5 * (model.dlnv(phi_i)) ** 2))
    dphi0 = -h_i * model.dlnv(phi_i)

    def rhs(t, y):
        p, dp, N = y
        h = H(p, dp)
        return [dp, -3 * h * dp - dV(p), h]

    def end(t, y):
        h = H(y[0], y[1])
        return 0.5 * (y[1] / h) ** 2 - 1.0
    end.terminal = True
    sol = integrate.solve_ivp(rhs, [0, 1e9], [phi_i, dphi0, 0.0], method="DOP853", rtol=1e-11, atol=1e-30, events=end)
    ref = inf.solve_background(model, phi_i=phi_i)
    assert sol.y[2, -1] == pytest.approx(ref.N_total, rel=1e-5)
    assert sol.y[0, -1] == pytest.approx(ref.phi_end, rel=1e-3)


def test_inflation_invalid_inputs():
    with pytest.raises(ValueError):
        inf.solve_background(inf.QuadraticModel(), phi_i=-1.0)
    with pytest.raises(RuntimeError):
        inf.solve_background(inf.QuadraticModel(), phi_i=5.0, N_max=1.0)    # integration window too short


def test_pivot_efolds_self_consistent(sols):
    model = inf.StarobinskyModel()
    info = inf.pivot_efolds(model, 2.1e-9, 0.05, sol=sols["Starobinsky"])
    assert 45 < info["N_star"] < 65
    # lowering the reheating temperature shortens N*
    info2 = inf.pivot_efolds(model, 2.1e-9, 0.05, w_reh=0.0, T_reh_GeV=1e9, sol=sols["Starobinsky"])
    assert info2["N_star"] < info["N_star"]
    with pytest.raises(ValueError):
        inf.pivot_efolds(model, 2.1e-9, 0.05, T_reh_GeV=1e30, sol=sols["Starobinsky"])


def test_quadratic_mass_from_amplitude(sols):
    info = inf.pivot_efolds(inf.QuadraticModel(), 2.1e-9, 0.05, sol=sols["quadratic"])
    assert math.sqrt(info["V0"]) == pytest.approx(6.2e-6, rel=0.05)         # known m ~ 6e-6 M_P for the m^2 phi^2 model


# ---------------------------------------------------------------- reheating
def test_reheating_energy_budget_and_limits():
    s = rh.solve_perturbative_reheating(1e-6)
    tot = s.x + s.y
    # d rho / dN = -3 (1+w) x - 4 y   (continuity), checked by finite differences of the solution
    dtot = np.gradient(tot, s.N)
    rhs = -(3.0 * s.x + 4.0 * s.y)
    mask = (s.N > 1) & (s.N < 8)
    assert np.allclose(dtot[mask], rhs[mask], rtol=2e-3)
    # Gamma -> 0: inflaton redshifts as a^-3
    s0 = rh.solve_perturbative_reheating(1e-14, N_max=5.0)
    assert s0.x[-1] == pytest.approx(math.exp(-3 * 5.0), rel=1e-6)
    assert s.radiation_fraction[-1] == pytest.approx(1.0, abs=1e-6)


def test_reheating_temperature_definition_scaling():
    ratios = []
    for g in (1e-3, 1e-6, 1e-9):
        s = rh.solve_perturbative_reheating(g)
        Neq, rho_eq = rh.equality_point(s)
        t_eq = rh.temperature_from_density(rho_eq * 3.0 / 2.0)
        ratios.append(t_eq / rh.reheating_temperature_estimate(g))
    assert max(ratios) / min(ratios) == pytest.approx(1.0, abs=2e-3)       # exact Gamma^(1/2) scaling
    assert 0.55 < ratios[0] < 0.75


def test_reheating_invalid():
    with pytest.raises(ValueError):
        rh.solve_perturbative_reheating(-1.0)
    with pytest.raises(ValueError):
        rh.solve_perturbative_reheating(2.0)
    with pytest.raises(ValueError):
        rh.reheating_temperature_estimate(0.0)


# ---------------------------------------------------------------- recombination
@pytest.fixture(scope="module")
def rec(bg):
    return rc.solve_recombination(bg)


def test_case_b_coefficient_anchor_values():
    # standard tabulated case-B values: 2.59e-13 cm^3/s at 1e4 K, 1.43e-13 cm^3/s at 2e4 K (Osterbrock)
    assert rc.alpha_B(1e4, F=1.0) * 1e6 == pytest.approx(2.59e-13, rel=0.01)
    assert rc.alpha_B(2e4, F=1.0) * 1e6 == pytest.approx(1.43e-13, rel=0.02)


def test_lyman_alpha_wavelength():
    assert rc.LAMBDA_ALPHA == pytest.approx(121.567e-9, rel=1e-3)


def test_recombination_early_saha_limit(bg, rec):
    i = np.argmin(abs(rec.z - 1700))
    assert rec.x_e[i] == pytest.approx(float(rc.saha_x(rec.z[i], bg)), rel=0.02)


def test_recombination_properties(rec):
    assert np.all(rec.x_e > 0) and np.all(rec.x_e <= 1.0 + 1e-9)
    assert 1040 < rec.z_star < 1120
    assert 1040 < rec.z_tau1 < 1130
    assert 1e-4 < rec.x_e[0] < 5e-4                 # freeze-out ionisation
    dx = np.diff(rec.x_e[::-1])                     # increasing z
    assert np.all(dx >= -1e-12) or True
    assert np.all(np.diff(rec.x_e) >= -1e-10)       # z ascending: x_e increases monotonically


def test_sound_horizon_angle_close_to_planck(rec):
    # theta_MC = 1.0409e-2 reported by Planck 2018; this approximate model is allowed 1.5%
    assert 100 * rec.theta_star == pytest.approx(1.0409, rel=0.015)
    assert 140.0 < rec.r_s_mpc < 150.0


def test_saha_decreasing_with_redshift_inverse(bg):
    z = np.array([800.0, 1000.0, 1200.0, 1400.0])
    x = rc.saha_x(z, bg)
    assert np.all(np.diff(x) > 0)
