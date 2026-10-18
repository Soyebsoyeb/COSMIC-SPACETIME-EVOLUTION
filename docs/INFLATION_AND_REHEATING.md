# Inflation and Reheating

Status: **model dependent**. Slow-roll inflation is the leading but unconfirmed paradigm for generating the primordial perturbations; the potentials below are examples, not detections. Reheating is a toy fluid model.

## 1. Setup and units

A canonical scalar field in flat FLRW with reduced Planck mass $M_P=1$. The derivation of the e-fold equation is in [MATHEMATICAL_DERIVATIONS.md](MATHEMATICAL_DERIVATIONS.md), section 7:

$$
\phi''=-(3-\epsilon_H)\left(\phi'+\frac{V_{,\phi}}{V}\right),\qquad \epsilon_H=\tfrac12\phi'^2,\qquad H^2=\frac{V}{3-\epsilon_H}.
$$

The solution is integrated in $N=\ln a$ from slow-roll initial conditions until $\epsilon_H=1$. Because the equation depends on $V$ only through $V_{,\phi}/V$, the energy scale is fixed afterwards by matching $\mathcal P_{\mathcal R}(k_*)=A_s$.

## 2. Potentials implemented

| Model | $V(\phi)$ | Leading-order prediction |
|---|---|---|
| Quadratic (chaotic) | $\tfrac12m^2\phi^2$ | $n_s\simeq1-2/N$, $r\simeq8/N$ |
| Starobinsky | $V_0\left(1-e^{-\sqrt{2/3}\,\phi}\right)^2$ | $n_s\simeq1-2/N$, $r\simeq12/N^2$ |

## 3. Results at the pivot scale (instantaneous reheating, $A_s=2.1\times10^{-9}$)

| Model | $N_*$ | $n_s$ | $r$ | $H_{\mathrm{end}}$ [GeV] |
|---|---:|---:|---:|---:|
| Starobinsky | 55.72 | 0.9645 | 0.0036 | $7.2\times10^{12}$ |
| Quadratic | 57.54 | 0.9649 | 0.140 | $7.6\times10^{12}$ |

Source: `results/tables/key_results.csv`. The tensor-to-scalar ratio of the quadratic model exceeds the published observational upper limit (of order $r<0.06$ at 95 percent confidence from combined CMB data; see the cited Planck 2018 inflation paper and the later BICEP/Keck analysis for the current bound). The code therefore **reproduces the known tension** of this potential; it does not rank models statistically, since no likelihood is computed. Figure: `fig13_inflation_ns_r`.

## 4. Pivot e-folds from the thermal history

The number of e-folds between horizon exit of the pivot scale $k_*$ and the end of inflation is fixed by requiring that $k_*=a_*H_*$ equals today's $k_*=a_0 k_*$ after the full expansion history:

$$
\ln\frac{k_*}{a_0H_0}=-N_*+\ln\frac{a_{\mathrm{end}}}{a_{\mathrm{reh}}}+\ln\frac{a_{\mathrm{reh}}}{a_0}+\ln\frac{H_*}{H_0}.
$$

`pivot_efolds` solves this with a root find, including a reheating epoch of equation of state $w_{\mathrm{reh}}$ and temperature $T_{\mathrm{reh}}$, and conserves comoving entropy with $g_{*s}$. The result depends only logarithmically on $T_{\mathrm{reh}}$, which is why $N_*\approx55$ to $58$ across a wide range. The dependence on the reheating assumptions is scanned in `fig14_reheating`.

## 5. Reheating model

The condensate energy $\rho_\phi$ decays into radiation $\rho_r$ with constant rate $\Gamma$ (equations in [MATHEMATICAL_DERIVATIONS.md](MATHEMATICAL_DERIVATIONS.md), section 9). Comoving variables remove the large dynamic range, and a stiff Radau integrator is used. Outputs: the equality point $\rho_\phi=\rho_r$, the instantaneous-decay temperature estimate from $\Gamma=H$, and the actual temperature from $\rho_r$. The choice of constant $\Gamma$, the absence of preheating (parametric resonance) and the absence of thermalisation physics are **model dependent simplifications**.

## 6. Perturbations

Mode equations for $\mathcal R$ and tensor modes are integrated from Bunch-Davies initial conditions deep inside the horizon until super-horizon freezing; see [PERTURBATION_THEORY.md](PERTURBATION_THEORY.md). The numerical scalar amplitude is 2.6 percent higher than the imposed $A_s$ because the normalisation uses the leading-order slow-roll formula, while the mode integration includes next-order corrections. This is reported, not hidden, in [VALIDATION.md](VALIDATION.md).

## 7. Limits

* Single canonical field only; no multi-field, non-canonical, non-Gaussianity or isocurvature physics.
* No quantum corrections to the potential or the initial state beyond Bunch-Davies.
* The observed spectral index and amplitude constrain, but do not uniquely select, an inflaton potential.
