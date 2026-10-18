# Thermal History of the Early Universe

Status: **approximate**. The thermal module counts relativistic degrees of freedom with a step function at $T=m$ and assumes instantaneous decoupling. It supports the background calculation and the BBN and recombination interfaces; it is not a Boltzmann solver.

## 1. Radiation-era relations (established)

For relativistic species in equilibrium with temperature $T$ (units $\hbar=c=k_B=1$):

$$
\rho=\frac{\pi^2}{30}g_*T^4,\qquad s=\frac{2\pi^2}{45}g_{*s}T^3,\qquad
H=\sqrt{\frac{8\pi G}{3}\rho},\qquad t=\frac1{2H}.
$$

The relation $t=1/(2H)$ is valid only while radiation dominates and $g_*$ is constant. The numerical constant in $t\approx2.42\,g_*^{-1/2}(T/\mathrm{MeV})^{-2}$ s is checked in `tests/test_early_universe.py`.

## 2. Degrees of freedom

Counting bosons with weight $g$ and fermions with weight $\tfrac78g$ for each species that is relativistic ($T\gtrsim m$):

| Temperature range | Species added or removed | $g_*$ |
|---|---|---:|
| above 175 GeV | all Standard Model particles | 106.75 |
| below top quark mass | top quark annihilates | 96.25 |
| below $W$, $Z$, Higgs | massive bosons removed | 86.25 |
| below bottom, charm, $\tau$ | heavy fermions removed | 75.75 to 61.75 |
| near the QCD transition | quarks and gluons to hadrons (pions, muons) | 17.25 |
| after muon annihilation | muons removed | 10.75 |
| after $e^\pm$ annihilation | only photons and neutrinos | 3.38 |

The code uses a step approximation (`early_universe/thermal_history.py`). The QCD transition is a smooth crossover in lattice calculations; the step there is a coarse approximation and the resulting error in $g_*$ during that interval is not quantified here. $g_{*s}$ after $e^\pm$ annihilation is $3.909$, while $g_*=3.38$ because neutrinos are colder than photons by $(4/11)^{1/3}$. Figure: `fig15_relativistic_degrees_of_freedom`.

## 3. Consistency with the background solver

The thermal Hubble rate $H(T)$ was compared with $H(z)$ from the FLRW module at $z=10^7$ to $10^9$ ($T\approx$ keV to MeV): agreement is at the $3\times10^{-4}$ level in the test, as expected because both use the same massless-neutrino radiation density and differ only through the $g_*$ step approximation. The temperature scaling $T=T_0(1+z)$ is shown in `fig05_photon_temperature_scaling`.

## 4. Neutrinos

$N_{\mathrm{eff}}=3.046$ is used as the standard-model value and neutrinos are treated as massless radiation at all times. **Neutrino masses are not modelled**; massive neutrinos change the matter-radiation equality and late-time growth at the sub-percent to percent level. This is a documented simplification, and the effect is listed in [RESEARCH_ROADMAP.md](RESEARCH_ROADMAP.md).

## 5. Baryon asymmetry

The baryon-to-photon ratio $\eta=n_b/n_\gamma=6.12\times10^{-10}$ is computed from $\Omega_bh^2$ and $T_{\mathrm{CMB}}$. The origin of the asymmetry (baryogenesis) is not modelled.

## 6. Not modelled

Electroweak and QCD phase-transition dynamics, dark-matter freeze-out, axions, gravitational-wave backgrounds from transitions, and neutrino decoupling with spectral distortions.
