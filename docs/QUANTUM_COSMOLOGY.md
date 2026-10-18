# Quantum Cosmology and the Planck Regime

**Module status: not implemented.** The repository contains no quantum-gravity or toy quantum-cosmology solver. This document states what is known, what is conjectured, and what the code does and does not do. Every claim below about the Planck epoch is labelled *speculative* unless marked otherwise.

## 1. Planck units (established definitions)

$$
\ell_P=\sqrt{\frac{\hbar G}{c^3}},\qquad t_P=\frac{\ell_P}{c},\qquad m_P=\sqrt{\frac{\hbar c}{G}},\qquad T_P=\frac{m_Pc^2}{k_B}.
$$

These are computed from CODATA constants in `constants.py` and listed in `results/tables/`. They mark the scale at which a classical metric description is expected to fail, because the quantum fluctuation of the geometry becomes comparable to the geometry itself. They are **dimensional-analysis scales, not a derived boundary**; the real transition scale is unknown.

## 2. The classical singularity

Extrapolating the radiation-dominated FLRW solution backwards gives $a\to0$, $\rho\to\infty$, $R\to\infty$ in finite proper time. The Hawking-Penrose theorems show that under the stated energy conditions and in classical GR this is generic geodesic incompleteness. The correct reading is that **classical GR ceases to be reliable** there, not that physical spacetime has a singularity. Nothing in this repository integrates toward $a=0$ beyond the radiation-era extrapolation, and the regularised function $Q(a)$ is used only to make the quadratures finite, which is a mathematical convenience and not a physical statement.

## 3. Candidate frameworks (literature overview, no implementation)

| Approach | Core idea | Status |
|---|---|---|
| Wheeler-DeWitt minisuperspace | Quantise $a$ and a matter field with $\hat H\Psi=0$; boundary conditions select a wave function (Hartle-Hawking, Vilenkin tunnelling) | Model dependent; factor ordering and boundary condition are choices |
| Loop quantum cosmology | Discrete geometry replaces the singularity by a bounce at $\rho\sim\rho_c\approx0.41\rho_P$ in the simplest model | Model dependent; relation to full loop quantum gravity unresolved |
| String-inspired scenarios | Pre-big-bang, ekpyrotic or brane constructions | Speculative |
| Eternal inflation / multiverse | Quantum fluctuations make inflation globally unending | Speculative; measure problem unresolved |
| Emergent / bouncing universes | Non-singular earlier phase | Speculative |

None of these is selected or tested in the code. Presenting any one of them as the true initial condition would not be honest.

## 4. What the code takes as the starting point

The earliest physical regime computed is slow-roll inflation of a canonical scalar field in classical GR ([INFLATION_AND_REHEATING.md](INFLATION_AND_REHEATING.md)), initialised at a field value chosen so that the pivot scale exits the horizon about 55 to 58 e-folds before the end of inflation. This is a model-dependent choice of initial condition, not a derivation from quantum gravity. The Hubble rate at the end of inflation implied by the observed scalar amplitude ($H_{\mathrm{end}}\sim7\times10^{12}$ GeV in the two models) is well below the reduced Planck mass $2.4\times10^{18}$ GeV, which is the consistency requirement for using classical GR.

## 5. Possible future extension

A minisuperspace toy model (for example a loop-quantum-cosmology effective Friedmann equation $H^2=\tfrac{8\pi G}{3}\rho(1-\rho/\rho_c)$) could be added as a clearly labelled model-dependent module. It is listed in [RESEARCH_ROADMAP.md](RESEARCH_ROADMAP.md) and is not part of version 0.1.0.
