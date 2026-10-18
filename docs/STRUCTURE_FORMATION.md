# Structure Formation

Status: **approximate**. This module links the primordial spectrum to the linear matter power spectrum and to analytic halo abundance. It does not evolve non-linear gravity.

## 1. Linear matter power spectrum

$$
\Delta^2_m(k,a)=\frac{4}{25}\left(\frac{k^2}{\Omega_mH_0^2}\right)^2D_+^2(a)\,T^2(k)\,\Delta^2_{\mathcal R}(k),\qquad
P(k)=\frac{2\pi^2}{k^3}\Delta^2_m,
$$

where $k$ is in $\mathrm{Mpc}^{-1}$, $H_0$ is expressed as $H_0/c$ in $\mathrm{Mpc}^{-1}$, $D_+$ is normalised to $a$ during matter domination, and $\Delta^2_{\mathcal R}=A_s(k/k_*)^{n_s-1}$. Radiation-era corrections to the amplitude conversion at $k\gtrsim0.5\,\mathrm{Mpc}^{-1}$ are neglected.

### Transfer functions

| Fit | Reference | Verification status |
|---|---|---|
| BBKS with Sugiyama shape parameter | Bardeen et al. 1986; Sugiyama 1995 | coefficients recalled, not re-checked |
| Eisenstein-Hu no-wiggle | Eisenstein and Hu 1998 | coefficients recalled, not re-checked |

Both are fitting formulae with known limited accuracy; neither contains baryon acoustic oscillations in the form used (no-wiggle). They are plotted together in `fig18_matter_power_spectrum`.

## 2. Variance and sigma_8

$$
\sigma^2(R)=\int\frac{dk}{k}\Delta^2_m(k)\,W^2(kR),\qquad W(x)=\frac{3(\sin x-x\cos x)}{x^3}.
$$

Result for the fiducial parameters: $\sigma_8=0.836$ (Eisenstein-Hu) and $0.847$ (BBKS). Planck 2018 quotes $0.811\pm0.006$, so the fit-based value is roughly **3 percent high**. The offset has identifiable sources: massless neutrinos, an approximate transfer function, and the fixed amplitude $A_s=2.1\times10^{-9}$ (rounded). This is reported as a limitation, and no attempt was made to tune the result to $0.811$.

## 3. Halo mass function

Press-Schechter and Sheth-Tormen:

$$
\frac{dn}{d\ln M}=\frac{\bar\rho_m}{M}\,f(\nu)\left|\frac{d\ln\sigma}{d\ln M}\right|,\qquad \nu=\frac{\delta_c}{\sigma(M)},\ \delta_c=1.686,
$$

with $f_{\mathrm{PS}}=\sqrt{2/\pi}\,\nu e^{-\nu^2/2}$ and $f_{\mathrm{ST}}=A\sqrt{2a/\pi}\,[1+(a\nu^2)^{-p}]\,\nu e^{-a\nu^2/2}$, $A=0.3222$, $a=0.707$, $p=0.3$ (Sheth and Tormen 1999, arXiv:astro-ph/9901122). The mass-integral normalisation of Sheth-Tormen is verified in the tests (the fraction of mass in halos integrates to one over the full range). Both formulae are **calibrations and ansaetze, accurate to tens of percent at best compared with N-body**; see `fig19_halo_mass_function`.

## 4. Not included

N-body or hydrodynamic evolution, baryonic feedback, non-linear power spectra (halofit), galaxy bias, redshift-space distortions, reionisation, black-hole formation and stellar evolution. These are not claimed.
