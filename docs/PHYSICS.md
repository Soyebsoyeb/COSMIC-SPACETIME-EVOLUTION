# Physical Overview

This document summarises the physical content of the repository, states which parts are established, which are computed, and which are model dependent, and fixes the notation used throughout. Detailed derivations are in [MATHEMATICAL_DERIVATIONS.md](MATHEMATICAL_DERIVATIONS.md).

## 1. Epistemic categories

Every statement, figure and table in this repository carries one of five labels. They are used in figure captions, in `results/tables/key_results.csv` and in the timeline figure.

| Label | Meaning | Example in this repository |
|---|---|---|
| Established | Theoretical result or observed fact that does not depend on this code | Redshift relation $1+z=a_0/a$, Friedmann equations, photon temperature scaling $T\propto(1+z)$ |
| Numerical | Computed by the implemented code from stated inputs with controlled numerical error | Age of the model universe, comoving distances, horizons, density fractions |
| Approximate | Computed with an explicit approximation and a stated domain of validity | Linear growth factor, recombination (Peebles model), transfer-function fits, $g_*(T)$ step function |
| Model dependent | Requires an assumed model that is not established | Inflationary potentials, reheating, halo mass-function ansatz, future scenarios |
| Speculative | Unresolved research question; no calculation claimed | Planck-era initial conditions, the true final state of the universe |

The repository does **not** claim that present-day physics supplies a complete, unique, first-principles calculation of the universe from a quantum-gravity initial state to a proven final endpoint.

## 2. The cosmic hierarchy of models

The universe is treated as a hierarchy of effective descriptions rather than one extrapolated equation.

| Regime | Effective description | Status in the code |
|---|---|---|
| $t\lesssim t_P$ | Quantum gravity (unknown) | Not implemented; discussed in [QUANTUM_COSMOLOGY.md](QUANTUM_COSMOLOGY.md) |
| Inflation | Canonical scalar field in classical GR | Implemented for two potentials (model dependent) |
| Reheating | Perturbative decay toy model, fluid bookkeeping | Implemented (model dependent) |
| Thermal radiation era | Ideal gas of relativistic species, $g_*(T)$ | Implemented (approximate) |
| BBN | Nuclear reaction network | **Not implemented**; interface and textbook relations only |
| Recombination | Peebles three-level atom | Implemented (approximate) |
| Linear perturbations | Mode equations on the FLRW background | Inflationary modes and matter growth implemented |
| Non-linear structure | N-body / hydrodynamics | **Not implemented**; analytic halo mass functions only |
| Late-time background | FLRW with matter, radiation, curvature, dark energy | Implemented (numerical) |
| Far future | Classical continuation of the background model | Implemented, strictly conditional |

## 3. Notation and conventions

* Signature $(-,+,+,+)$. Natural units $\hbar=c=k_B=1$ are used inside the inflation and mode modules with the **reduced** Planck mass $M_P=(8\pi G)^{-1/2}=1$; SI units are used elsewhere.
* Scale factor normalised to $a(t_0)=1$; redshift $1+z=1/a$.
* Present-day density parameters $\Omega_i=\rho_{i,0}/\rho_{\mathrm{crit},0}$ with $\rho_{\mathrm{crit},0}=3H_0^2/(8\pi G)$ and closure $\Omega_r+\Omega_m+\Omega_k+\Omega_{\mathrm{DE}}=1$.
* $E(a)=H(a)/H_0$ and the regularised function $Q(a)\equiv a^4E^2(a)$, which stays finite as $a\to0$.
* Dark energy follows the CPL parametrisation $w(a)=w_0+w_a(1-a)$; the cosmological constant is $w_0=-1$, $w_a=0$.
* Comoving distances in Mpc, times in Gyr or seconds, wavenumbers $k$ in $\mathrm{Mpc}^{-1}$ unless stated.

## 4. Key equations

Friedmann equation for the background:

$$
H^2(a)=H_0^2\left[\Omega_r a^{-4}+\Omega_m a^{-3}+\Omega_k a^{-2}+\Omega_{\mathrm{DE}}\,a^{-3(1+w_0+w_a)}e^{-3w_a(1-a)}\right].
$$

Acceleration (deceleration parameter) with $w_i$ the equation of state of each non-curvature component:

$$
q\equiv-\frac{a\ddot a}{\dot a^2}=\frac12\sum_i\Omega_i(a)\,\bigl(1+3w_i(a)\bigr).
$$

Cosmic time and comoving distances, all expressed through $Q$:

$$
H_0t(a)=\int_0^a\frac{a'\,da'}{\sqrt{Q(a')}},\qquad
\chi(a)=\frac{c}{H_0}\int_a^1\frac{da'}{\sqrt{Q(a')}} .
$$

Inflation (e-fold form, $\epsilon_H=-\dot H/H^2$):

$$
\phi''=-(3-\epsilon_H)\left(\phi'+\frac{V_{,\phi}}{V}\right),\qquad H^2=\frac{V}{3-\epsilon_H},\qquad \epsilon_H=\tfrac12\phi'^2 .
$$

## 5. What the homogeneous background can and cannot say

The FLRW solver determines the expansion history, distances and horizons of the *model*. It does not by itself predict CMB anisotropies, galaxies, halos or black holes; these require perturbation theory beyond the background and, for non-linear structure, numerical simulation. Where the repository touches those topics it does so with explicit linear or phenomenological models and says so.

## 6. Where to look

| Topic | Document | Code |
|---|---|---|
| Geometry and Einstein equations | [GENERAL_RELATIVITY.md](GENERAL_RELATIVITY.md) | `geometry/`, `background/` |
| Planck regime | [QUANTUM_COSMOLOGY.md](QUANTUM_COSMOLOGY.md) | `constants.py` (Planck units only) |
| Inflation, reheating | [INFLATION_AND_REHEATING.md](INFLATION_AND_REHEATING.md) | `early_universe/inflation.py`, `reheating.py` |
| Thermal history | [THERMAL_HISTORY.md](THERMAL_HISTORY.md) | `early_universe/thermal_history.py`, `entropy.py` |
| BBN and recombination | [BBN_AND_RECOMBINATION.md](BBN_AND_RECOMBINATION.md) | `early_universe/bbn_interface.py`, `recombination.py` |
| Perturbations | [PERTURBATION_THEORY.md](PERTURBATION_THEORY.md) | `perturbations/` |
| Structure | [STRUCTURE_FORMATION.md](STRUCTURE_FORMATION.md) | `structure/` |
| Horizons | [COSMIC_HORIZONS.md](COSMIC_HORIZONS.md) | `geometry/horizons.py` |
| Future | [FUTURE_COSMOLOGY.md](FUTURE_COSMOLOGY.md) | `future/` |
| Numerics and validation | [NUMERICAL_METHODS.md](NUMERICAL_METHODS.md), [VALIDATION.md](VALIDATION.md) | `numerics/`, `tests/` |
