# Cosmic Horizons and Causal Structure

## 1. Definitions

Using $\chi$ for comoving distance and $\tau$ for conformal time ($d\tau=dt/a$):

| Concept | Definition | Formula in the code |
|---|---|---|
| Particle horizon | Largest comoving distance from which signals could have reached us since the beginning | $\chi_p=D_H\int_0^a da'/\sqrt{Q}$ |
| Event horizon | Largest comoving distance from which signals sent now can ever reach us | $\chi_e=D_H\int_a^{a_*}da'/\sqrt{Q}$ |
| Hubble radius | $c/H$; sets whether a mode is within the horizon; **not** a causal boundary | `hubble_radius` |
| Comoving Hubble radius | $c/(aH)$; shrinks during inflation | `fig10_comoving_hubble_radius` |

The event horizon exists only when the integral to $a_*=\infty$ converges, which in a model with an asymptotically power-law $a(t)\propto t^p$ requires $p>1$ (accelerated expansion). `event_horizon_exists` tests this numerically via the asymptotic logarithmic slope of the integrand.

## 2. Fiducial values

| Quantity | Value | Status |
|---|---:|---|
| Particle horizon today | 14.14 Gpc | numerical |
| Event horizon today | 5.11 Gpc (about 16.7 Gly) | numerical, conditional on a cosmological constant |
| de Sitter Hubble radius $c/H_\Lambda$ | 17.5 Gly | numerical, conditional |

The particle horizon includes the radiation epoch, extrapolated to $a=0$ without inflation. With a preceding inflationary phase, the particle horizon of the *spacetime* is much larger than the value computed from the hot-big-bang solution alone, which is the usual horizon-problem resolution. The code reports the horizon of the FLRW model from $a=0$ and does not add the inflationary contribution.

## 3. Horizon problem

At recombination the comoving particle horizon is a small fraction of the distance to the last-scattering surface (about 280 Mpc against $D_C\approx13.9$ Gpc in the fiducial model, i.e. an angular scale of about one degree). Without an earlier period in which the comoving Hubble radius $1/(aH)$ decreases, widely separated CMB patches would never have been in causal contact. `fig10_comoving_hubble_radius` and `fig22_conformal_diagrams` illustrate this.

## 4. Conformal diagrams

`fig22_conformal_diagrams` plots light cones in $(\eta,\chi)$ for Einstein-de Sitter, flat $\Lambda$CDM, the closed matter model and the big-rip model. Event horizons are drawn only where they exist.

## 5. Cautions

Horizon radii depend on the assumed future expansion. The "event horizon today" is conditional on $\Lambda$ staying constant; with dynamical dark energy it can change or vanish.
