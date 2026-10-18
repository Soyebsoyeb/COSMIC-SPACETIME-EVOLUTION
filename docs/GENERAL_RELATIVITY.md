# General Relativity and Spacetime Geometry

## 1. The FLRW line element

Homogeneity and isotropy of space single out the Friedmann-Lemaitre-Robertson-Walker metric

$$
ds^2=-c^2dt^2+a^2(t)\left[\frac{dr^2}{1-kr^2}+r^2\,d\Omega^2\right],\qquad d\Omega^2=d\theta^2+\sin^2\theta\,d\varphi^2 ,
$$

where $t$ is the proper time of comoving observers (cosmic time), $r$ is a comoving radial coordinate, and $k\in\{-1,0,+1\}$ labels the spatial curvature after rescaling $r$ so that $a$ carries the length scale.

| Quantity | Definition | Meaning |
|---|---|---|
| Scale factor $a(t)$ | Ratio of proper separations of comoving observers to their value today | Overall expansion of space |
| Comoving distance $\chi$ | $\chi=\int dr/\sqrt{1-kr^2}$ | Separation that is constant for comoving observers |
| Proper distance | $d_p=a(t)\chi$ at fixed $t$ | Not directly observable for distant objects |
| Hubble rate | $H=\dot a/a$ | Expansion rate; local recession law $v=Hd_p$ |
| Curvature radius | $a\,R_c$ with $R_c=1/\sqrt{\lvert k\rvert}$ | Defined for $k\neq0$; today $=(c/H_0)/\sqrt{\lvert\Omega_k\rvert}$ |

Convention: $a(t_0)=1$ (so comoving and present proper distances coincide), and $\Omega_k\equiv-kc^2/(a_0^2H_0^2)$, so $\Omega_k>0$ is an open (negatively curved) universe.

## 2. Redshift

A photon follows a null geodesic, $c\,dt=a\,d\chi$. Two successive wave crests emitted at $t_e$ and $t_e+\delta t_e$ and received at $t_0$ and $t_0+\delta t_0$ satisfy $\int_{t_e}^{t_0}c\,dt/a=\int_{t_e+\delta t_e}^{t_0+\delta t_0}c\,dt/a$, hence $\delta t_0/a_0=\delta t_e/a_e$ and

$$
1+z\equiv\frac{\lambda_0}{\lambda_e}=\frac{\nu_e}{\nu_0}=\frac{a(t_0)}{a(t_e)} .
$$

Setting $a(t_0)=1$ gives $1+z=1/a$. The choice $a_0=1$ is a convention: only ratios of scale factors are observable (in a flat universe the overall scale of $a$ and of comoving coordinates can be rescaled together).

## 3. Einstein field equations

$$
G_{\mu\nu}+\Lambda g_{\mu\nu}=\frac{8\pi G}{c^4}T_{\mu\nu}.
$$

* The **Einstein tensor** $G_{\mu\nu}=R_{\mu\nu}-\tfrac12Rg_{\mu\nu}$ encodes curvature; $R_{\mu\nu}$ is the Ricci tensor and $R$ the Ricci scalar.
* The **stress-energy tensor** $T_{\mu\nu}$ describes matter and radiation. For a perfect fluid with four-velocity $u^\mu$: $T_{\mu\nu}=(\rho+p/c^2)u_\mu u_\nu c^2+p\,g_{\mu\nu}$.
* The **cosmological constant** $\Lambda$ is equivalent to a perfect fluid with $\rho_\Lambda c^2=\Lambda c^4/(8\pi G)$ and $p_\Lambda=-\rho_\Lambda c^2$ ($w=-1$).
* The **contracted Bianchi identity** $\nabla^\mu G_{\mu\nu}=0$, together with $\nabla^\mu g_{\mu\nu}=0$, forces $\nabla^\mu T_{\mu\nu}=0$: local energy-momentum conservation is a consequence of the field equations, not an extra assumption, and $\Lambda$ must be constant.

## 4. Curvature of the FLRW metric

With $c=1$ the non-vanishing Christoffel symbols are

$$
\Gamma^0_{ij}=a\dot a\,\gamma_{ij},\qquad \Gamma^i_{0j}=\frac{\dot a}{a}\delta^i_j,\qquad \Gamma^i_{jk}=\text{(Christoffel symbols of }\gamma_{ij}),
$$

with $\gamma_{ij}$ the metric of the unit-curvature spatial slice. The Ricci tensor components are

$$
R_{00}=-3\frac{\ddot a}{a},\qquad R_{ij}=\bigl(a\ddot a+2\dot a^2+2k\bigr)\gamma_{ij},
$$

and the Ricci scalar

$$
R=6\left(\frac{\ddot a}{a}+\frac{\dot a^2}{a^2}+\frac{k}{a^2}\right).
$$

This expression is implemented in `geometry/curvature.py::ricci_scalar_over_H0sq` and tested against the exact de Sitter value $R=12H^2$.

The Einstein tensor is

$$
G_{00}=3\left(\frac{\dot a^2}{a^2}+\frac{k}{a^2}\right),\qquad G_{ij}=-\bigl(2a\ddot a+\dot a^2+k\bigr)\gamma_{ij}.
$$

## 5. From the field equations to the Friedmann equations

For comoving observers $u^\mu=(1,0,0,0)$, so $T_{00}=\rho$ and $T_{ij}=p\,a^2\gamma_{ij}$. The $00$ component of the field equations gives the first Friedmann equation, and the $ij$ components combined with it give the second:

$$
H^2+\frac{kc^2}{a^2}=\frac{8\pi G}{3}\rho+\frac{\Lambda c^2}{3},\qquad
\frac{\ddot a}{a}=-\frac{4\pi G}{3}\left(\rho+\frac{3p}{c^2}\right)+\frac{\Lambda c^2}{3}.
$$

The full derivation, including the continuity equation from $\nabla_\mu T^{\mu\nu}=0$, is in [MATHEMATICAL_DERIVATIONS.md](MATHEMATICAL_DERIVATIONS.md).

## 6. Singularities and their limits

The Hawking-Penrose singularity theorems show that, under energy conditions and with classical GR, generic expanding cosmologies have an initial singularity (geodesic incompleteness). The theorems say that **classical** GR cannot be continued past that point; they do not say that spacetime has a singularity in the full quantum theory. See [QUANTUM_COSMOLOGY.md](QUANTUM_COSMOLOGY.md). Phantom dark energy ($w<-1$) violates the null energy condition and admits a finite-time future singularity in the classical model ([FUTURE_COSMOLOGY.md](FUTURE_COSMOLOGY.md)).

## 7. Implementation map

| Concept | Module |
|---|---|
| $E(a)$, $H(a)$, density fractions, $q(a)$ | `background/friedmann.py` |
| CPL dark energy | `background/dark_energy.py` |
| Cosmic and conformal time | `background/cosmic_time.py` |
| $D_C$, $D_M$, $D_L$, $D_A$ for any sign of $\Omega_k$ | `geometry/flrw.py` |
| Curvature radius, Ricci scalar | `geometry/curvature.py` |
| Horizons, Hubble radius | `geometry/horizons.py` |
