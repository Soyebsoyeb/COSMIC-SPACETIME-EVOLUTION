# Validation Report

This file separates (A) checks against exact analytic results, (B) convergence checks, and (C) known approximation-level discrepancies. The machine-generated version is `results/logs/validation_report.md` (regenerate with `python scripts/validate_numerics.py`). No continuous-integration run has been performed on GitHub; results below were obtained in the development environment only.

## A and B. Analytic and convergence checks (13 of 13 pass)

| Check | Relative error | Tolerance |
|---|---:|---:|
| Einstein-de Sitter age $H_0t_0=2/3$ | $1.4\times10^{-14}$ | $10^{-9}$ |
| Flat $\Lambda$CDM age vs closed form | $4.9\times10^{-14}$ | $10^{-9}$ |
| Milne universe age | $1.0\times10^{-14}$ | $10^{-9}$ |
| EdS comoving distance, $z=3$ | $9.8\times10^{-15}$ | $10^{-9}$ |
| EdS particle horizon | $2.4\times10^{-14}$ | $10^{-9}$ |
| de Sitter event horizon | $1.0\times10^{-15}$ | $10^{-7}$ |
| Closed matter universe crunch time $2\pi/H_0$ | $2.2\times10^{-13}$ | $10^{-5}$ |
| Second-order acceleration ODE vs quadrature | $7.2\times10^{-13}$ | $10^{-8}$ |
| Growth ODE vs Heath integral | $4.4\times10^{-11}$ | $10^{-7}$ |
| Mukhanov-Sasaki vs de Sitter $\mathcal P_{\mathcal R}$ | $4.3\times10^{-5}$ | $10^{-3}$ |
| Simpson age, $N=321$ | $2.9\times10^{-8}$ | $10^{-6}$ |
| Simpson age, $N=1281$ | $1.2\times10^{-10}$ | $10^{-9}$ |
| Recombination $100\theta_*$ vs Planck | $6.3\times10^{-3}$ | $1.5\times10^{-2}$ |

The last row is a model-approximation comparison, not a numerical-accuracy test, and its tolerance was set from the expected accuracy of the Peebles model before the comparison was examined.

## Unit and regression tests

The `tests/` directory contains 113 tests (geometry, Friedmann, cosmic time, early universe, perturbations, future models, convergence) which all passed in the final development run. During development, nine test expectations were found to be wrong (for example the value of $g_*$ at a given temperature, the radiation-era time-temperature constant and the closed-universe growth initial-condition accuracy) while the code was right; in each case the expectation was corrected with a documented derivation, and no tolerance was loosened merely to obtain a pass. Reviewers are encouraged to read `tests/` directly.

## C. Known discrepancies with observation or literature

| Quantity | This code | Reference | Cause |
|---|---:|---:|---|
| $\sigma_8$ | 0.836 (EH), 0.847 (BBKS) | 0.811 $\pm$ 0.006 (Planck 2018) | fitted transfer function, massless neutrinos, rounded $A_s$ |
| $100\theta_*$ | 1.0475 | 1.0409 $\pm$ 0.0003 | Peebles model approximations |
| $z_*$ | 1078.9 | about 1090 | same |
| Numerical $A_s$ vs imposed | +2.6 percent | | leading-order normalisation formula |
| Big-rip estimate | +2.9 percent vs numeric | | approximation neglects matter |
| Quadratic inflation $r$ | 0.140 | observational limit of order $r<0.06$ | model is disfavoured; reproduced, not tuned |

## D. What has not been validated

BBN (not implemented), full CMB spectra (not implemented), non-linear structure (not implemented), massive neutrinos (not modelled), and transfer-function coefficients and RECFAST-style fit constants against their original papers (recalled from the literature; see [REFERENCES.md](REFERENCES.md)). Anyone using these for quantitative science should compare against CAMB or CLASS first.
