# Numerical Methods

## 1. Quadrature

All background integrals are written in terms of $Q(a)=a^4E^2$ so that integrands are finite at $a\to0$ and no singularities need to be removed by hand.

| Task | Method | Module |
|---|---|---|
| Time, distances, horizons | adaptive Gauss-Kronrod (`scipy.integrate.quad`) with the error estimate checked; default tolerances $\mathrm{epsabs}=10^{-12}$, $\mathrm{epsrel}=10^{-10}$ (the fiducial config requests $10^{-13}$ and $10^{-11}$) | `numerics/integration.py` |
| Cross-check | composite trapezoid and Simpson in $\ln a$ | `numerics/integration.py` |
| Infinite upper limit | split in $u=\ln a$ at decade breakpoints; finite-cutoff tail added analytically where applicable | `numerics/convergence.py` |
| Cumulative time for plots | segment-wise quadrature (avoids non-monotonic spline artefacts) | `future/fate_classification.py` |

## 2. ODE integration

| System | Integrator | Tolerances |
|---|---|---|
| Inflation background | `solve_ivp`, DOP853, in $N$ | rtol $10^{-11}$, atol $10^{-13}$ |
| Scalar and tensor modes | complex mode split into real parts, DOP853, in $N$ | rtol $10^{-10}$, atol $10^{-40}$ |
| Reheating | Radau (stiff), comoving variables | rtol $10^{-10}$, atol $10^{-30}$ |
| Recombination | Radau (stiff) in redshift | rtol $10^{-9}$, atol $10^{-14}$ |
| Growth factor | DOP853 in $\ln a$; cubic spline in $\ln a$ for evaluation | rtol $10^{-10}$, atol $10^{-14}$ |

## 3. Convergence studies

`scripts/validate_numerics.py` and `fig20_numerical_convergence` report:

* Simpson age error against the number of log-spaced points: $2.9\times10^{-8}$ at $N=321$ and $1.2\times10^{-10}$ at $N=1281$, consistent with fourth-order convergence.
* Age error against the lower integration cutoff $a_{\min}$ (within the stated validity region).
* Quadrature tolerance sweep.
* Event-horizon error against the upper cutoff $a_{\max}$ (integrated in $\ln a$).

## 4. Uncertainty and sensitivity

`numerics/uncertainty.py` provides log-sensitivities $\partial\ln X/\partial\ln p$ by central differences and a Monte Carlo propagation with **independent Gaussian** parameter errors. The Monte Carlo is **illustrative** only: it ignores parameter correlations and non-Gaussianity of the real posterior, so the resulting spreads are not confidence intervals. Seed 12345 makes runs reproducible. See `fig21_parameter_sensitivity`.

## 5. Reproducibility

Deterministic Matplotlib style and fixed seeds; `results/metadata/environment.json` records versions; `scripts/reproduce_results.py` writes SHA-256 checksums of all outputs. PNG files can differ at the byte level across Matplotlib or platform versions, so checksums are guidance and not a guarantee of bit-identical images; the CSV datasets are the numerical record.

## 6. Known limits

Floating-point double precision; no arbitrary-precision cross-check; no adaptive mesh in $k$ for mode integration (a fixed log grid is used).
