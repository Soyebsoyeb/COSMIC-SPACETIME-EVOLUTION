# Future Cosmology

**Every result in this document is conditional on the assumed model.** The code continues the classical FLRW background into the future under stated dark-energy assumptions. It does not predict which scenario is realised, and the physical nature of dark energy is unknown.

## 1. Scenarios and computed fates

Computed with `future/fate_classification.py` using the fiducial matter and radiation content and the dark-energy or curvature settings listed (`results/tables/`, `fig12_future_scenarios`).

| Scenario | Fate classification | Time of end (cosmic time) | Category |
|---|---|---:|---|
| $\Lambda$CDM ($\Lambda>0$) | asymptotic de Sitter | none | model dependent |
| Constant $w=-0.8$ | eternal accelerated power-law expansion | none | model dependent |
| Phantom $w=-1.2$ | finite-time singularity (big rip) | 71.7 Gyr | model dependent |
| Phantom $w=-1.5$ | finite-time singularity (big rip) | 37.3 Gyr | model dependent |
| Einstein-de Sitter ($\Omega_m\approx1$) | eternal decelerating expansion | none | model dependent, not our universe |
| Open matter-only ($\Omega_m=0.315$) | eternal coasting expansion ($q\to0$) | none | model dependent, not our universe |
| Closed matter-only ($\Omega_m=2$) | recollapse (big crunch) | 91.1 Gyr | model dependent, not our universe |
| Negative $\Lambda$ ($\Omega_\Lambda=-0.5$) | recollapse (big crunch) | 57.6 Gyr | model dependent |

The EdS, open and closed rows are comparison universes that disagree with observations; they illustrate the classification, not forecasts.

## 2. Classification method

The classification uses the asymptotic logarithmic slope of the integrand $1/\sqrt Q$ of the time integral at large $a$. If the slope is below $-1.01$ the integral $\int^\infty da/\sqrt Q$ for the event horizon converges (eternal accelerated expansion), and for $w<-1$ the time integral $\int^\infty a\,da/\sqrt Q$ converges, signalling a finite-time singularity. A turning point where $Q=0$ at finite $a$ indicates recollapse. This is a numerical criterion; a slope within a few percent of the threshold should be treated as inconclusive. For $\Lambda$ domination the slope is $-2$ for the event-horizon integrand and the fate is de Sitter.

## 3. Established and semiclassical results used

* Asymptotic de Sitter: $H\to H_\Lambda=H_0\sqrt{\Omega_\Lambda}$, horizon radius $c/H_\Lambda=17.5$ Gly (computed from the fiducial parameters).
* Gibbons-Hawking temperature $T=\hbar H_\Lambda/(2\pi k_B)=2.2\times10^{-30}$ K: established for exact de Sitter space in semiclassical gravity (Gibbons and Hawking 1977). Its role in the real universe is model dependent.
* Big-rip time estimate $t_{\mathrm{rip}}-t_0\simeq2/(3\lvert1+w\rvert H_0\sqrt{1-\Omega_m})$ (Caldwell, Kamionkowski and Weinberg 2003) is an approximation; the estimate exceeds the numerical value by roughly 1 to 3 percent in the cases tested (0.9 percent for the fiducial $w=-1.2$ model).

## 4. Not covered

Vacuum decay, proton decay and black-hole evaporation timescales, Poincare recurrence, Boltzmann brains, cyclic models, interacting dark energy and modified gravity are **not** modelled and no claim is made about them. The statement "the universe will expand forever" holds only within the listed models.
