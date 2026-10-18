# Perturbation Theory

## 1. Primordial perturbations from inflation (model dependent)

The comoving curvature perturbation $\mathcal R$ and the two tensor polarisations $h$ are integrated mode by mode in e-fold time (equations in [MATHEMATICAL_DERIVATIONS.md](MATHEMATICAL_DERIVATIONS.md), section 8):

$$
\mathcal R''+(3-\epsilon_H+\eta_H)\mathcal R'+\left(\frac{k}{aH}\right)^2\mathcal R=0,\qquad
h''+(3-\epsilon_H)h'+\left(\frac{k}{aH}\right)^2h=0.
$$

Initial conditions are Bunch-Davies, imposed at $k/(aH)=400$ (default `x_start`), and the spectra are read off after freezing. Outputs: $\mathcal P_{\mathcal R}(k)$, $\mathcal P_t(k)$, $n_s$ and $r$ from the numerical solutions (`perturbations/primordial_spectrum.py`; `fig17_primordial_spectra`).

**Checks.** In exact de Sitter the solver returns $\mathcal P_{\mathcal R}$ within $4\times10^{-5}$ of the analytic $H^2/(8\pi^2\epsilon M_P^2)$ (validation suite), and the tensor amplitude within the $10^{-3}$ tolerance of the slow-roll value (test suite).

## 2. Linear matter growth (approximate)

$$
D''+(1-q)D'-\tfrac32\Omega_m(a)D=0,\qquad D=D(\ln a).
$$

* Valid for **sub-horizon, scale-independent growth of pressureless matter**. Neutrino free-streaming, baryon-photon coupling before decoupling and super-horizon gauge effects are not modelled.
* Initial conditions from the Meszaros solution in the radiation era. For curved models the initial condition has a documented limited accuracy ($10^{-4}$ level at $a_i=10^{-4}$, improved at smaller $a_i$); see the docstring of `solve_growth`.
* Checks: Einstein-de Sitter gives $D=a$; flat $\Lambda$CDM matches the Heath integral to $4\times10^{-11}$; the growth index is close to the standard $\gamma\approx0.55$.

Fiducial result: $D_+(a=1)=0.788$ when normalised to $D\to a$ at early matter domination (`fig11_linear_growth_factor`).

## 3. What is not here

No Boltzmann hierarchy, no photon, baryon or neutrino perturbations, no CMB temperature or polarisation spectrum, no gravitational lensing, no second-order or non-Gaussian effects, no isocurvature modes. Users needing $C_\ell$ should use CAMB or CLASS; an interface to them is on the roadmap.
