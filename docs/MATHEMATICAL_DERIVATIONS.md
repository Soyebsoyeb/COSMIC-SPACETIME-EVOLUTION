# Mathematical Derivations

Throughout, $c=1$ unless stated and overdots denote $d/dt$. Each section ends with the repository module or test that depends on the result.

## 1. Friedmann equations from the Einstein equations

Take the FLRW metric with $a(t)$, curvature index $k$, and the perfect-fluid stress tensor $T^\mu{}_\nu=\mathrm{diag}(-\rho,p,p,p)$. The curvature components (see [GENERAL_RELATIVITY.md](GENERAL_RELATIVITY.md)) are $G_{00}=3(\dot a^2+k)/a^2$ and $G_{ij}=-(2a\ddot a+\dot a^2+k)\gamma_{ij}$.

**00 component.** $G_{00}+\Lambda g_{00}=8\pi G\,T_{00}$ with $g_{00}=-1$, $T_{00}=\rho$:

$$
\frac{3(\dot a^2+k)}{a^2}-\Lambda=8\pi G\rho\quad\Longrightarrow\quad H^2+\frac{k}{a^2}=\frac{8\pi G}{3}\rho+\frac{\Lambda}{3}.\tag{F1}
$$

**ij components.** $-(2a\ddot a+\dot a^2+k)+\Lambda a^2=8\pi G\,p\,a^2$. Substituting (F1) to eliminate $\dot a^2+k$,

$$
-2a\ddot a-\left(\frac{8\pi G}{3}\rho+\frac{\Lambda}{3}\right)a^2+\Lambda a^2=8\pi Gpa^2
\;\Longrightarrow\;
\frac{\ddot a}{a}=-\frac{4\pi G}{3}(\rho+3p)+\frac{\Lambda}{3}.\tag{F2}
$$

Defining the critical density $\rho_c=3H^2/(8\pi G)$ and the density parameters $\Omega_i=\rho_i/\rho_c$, $\Omega_k=-k/(aH)^2$ and $\Omega_\Lambda=\Lambda/(3H^2)$, (F1) reads $\sum_i\Omega_i+\Omega_k+\Omega_\Lambda=1$ at every time.

## 2. Energy conservation and the density scalings

$\nabla_\mu T^{\mu\nu}=0$ with $\nu=0$ gives, for the perfect fluid,

$$
\partial_t\rho+3\frac{\dot a}{a}(\rho+p)=0 .\tag{C}
$$

Derivation: $\nabla_\mu T^{\mu0}=\partial_\mu T^{\mu0}+\Gamma^\mu_{\mu\lambda}T^{\lambda0}+\Gamma^0_{\mu\lambda}T^{\mu\lambda}$. With $T^{00}=\rho$, $T^{ij}=p\gamma^{ij}/a^2$, $\Gamma^\mu_{\mu0}=3\dot a/a$ and $\Gamma^0_{ij}=a\dot a\gamma_{ij}$ one finds $\dot\rho+3(\dot a/a)\rho+a\dot a\gamma_{ij}\,p\gamma^{ij}/a^2=\dot\rho+3H(\rho+p)=0$.

For a constant equation of state $p=w\rho$, (C) integrates to $\rho\propto a^{-3(1+w)}$:

| Component | $w$ | $\rho\propto$ | Physical origin |
|---|---|---|---|
| Radiation | $1/3$ | $a^{-4}$ | Number density dilutes as $a^{-3}$ and each photon is redshifted, $E\propto a^{-1}$ |
| Matter | $0$ | $a^{-3}$ | Number dilutes as $a^{-3}$; rest energy is constant |
| Curvature | $-1/3$ | $a^{-2}$ | Effective density $-3k/(8\pi Ga^2)$ |
| Cosmological constant | $-1$ | const | Vacuum energy density does not dilute; work $-p\,dV=\rho\,dV$ exactly compensates volume growth |

For the CPL model $w(a)=w_0+w_a(1-a)$: $d\ln\rho/d\ln a=-3(1+w_0+w_a-w_aa)$, whose solution is
$\rho(a)/\rho(1)=a^{-3(1+w_0+w_a)}e^{-3w_a(1-a)}$. (`background/dark_energy.py`; test `test_dark_energy_cpl_density_solves_continuity`.)

## 3. The regularised expansion function $Q(a)$

Using (F1) and the scalings,

$$
E^2(a)=\Omega_ra^{-4}+\Omega_ma^{-3}+\Omega_ka^{-2}+\Omega_{\mathrm{DE}}f(a),\qquad Q(a)\equiv a^4E^2=\Omega_r+\Omega_ma+\Omega_ka^2+\Omega_{\mathrm{DE}}a^4f(a).
$$

$Q$ is finite and (with radiation) positive as $a\to0$, which removes the integrable singularities that appear in $1/E$.

## 4. Cosmic time, distances and horizons

From $H=\dot a/a$: $dt=da/(aH)$, hence

$$
H_0t(a)=\int_0^a\frac{da'}{a'E(a')}=\int_0^a\frac{a'\,da'}{\sqrt{Q(a')}} .
$$

Light satisfies $d\chi=-dt/a$, so the comoving distance to a source emitted at scale factor $a$ is

$$
D_C=\frac{c}{H_0}\int_a^1\frac{da'}{a'^2E(a')}=\frac{c}{H_0}\int_a^1\frac{da'}{\sqrt{Q(a')}} .
$$

Transverse comoving distance $D_M=\mathrm{sinn}(D_C)$: $D_C$ for $\Omega_k=0$; $(D_H/\sqrt{\Omega_k})\sinh(\sqrt{\Omega_k}D_C/D_H)$ for $\Omega_k>0$; $(D_H/\sqrt{-\Omega_k})\sin(\sqrt{-\Omega_k}D_C/D_H)$ for $\Omega_k<0$, with $D_H=c/H_0$. The angular-diameter and luminosity distances are $D_A=D_M/(1+z)$ and $D_L=(1+z)D_M$, so $D_L=(1+z)^2D_A$ (distance duality).

The particle horizon is the same integral from $0$: $\chi_p(a)=D_H\int_0^a da'/\sqrt{Q}$. The comoving event horizon is $\chi_e(a)=D_H\int_a^{a_\ast}da'/\sqrt{Q}$ with $a_\ast=\infty$ for eternal expansion or $a_{\max}$ for a recollapsing model.

*Analytic checks used in the tests:* for $\Omega_m=1$, $H_0t_0=2/3$, $D_C=2D_H(1-(1+z)^{-1/2})$, $\chi_p=2D_H$; for a flat $\Lambda$CDM without radiation $H_0t_0=\frac{2}{3\sqrt{\Omega_\Lambda}}\,\mathrm{arsinh}\sqrt{\Omega_\Lambda/\Omega_m}$; for the empty (Milne) universe $H_0t_0=1$ and $D_M=D_H\sinh\ln(1+z)$; for a closed matter-only universe $a_{\max}=\Omega_m/(\Omega_m-1)$ and the crunch occurs at $H_0t=\pi\Omega_m/(\Omega_m-1)^{3/2}$ (for $\Omega_m=2$: $2\pi$).

## 5. Deceleration parameter

Dividing (F2) by $H^2$ and using $\Omega_i=8\pi G\rho_i/(3H^2)$,

$$
q\equiv-\frac{\ddot a}{aH^2}=\frac{4\pi G}{3H^2}\sum_i\rho_i(1+3w_i)=\frac12\sum_i\Omega_i(1+3w_i).
$$

Curvature has $w_k=-1/3$ and drops out. Limits: $q=1$ (radiation), $1/2$ (matter), $-1$ (de Sitter). For flat $\Lambda$CDM neglecting radiation, $q=0$ when $\Omega_m a^{-3}=2\Omega_\Lambda$, i.e. $a_{\mathrm{acc}}=(\Omega_m/2\Omega_\Lambda)^{1/3}$.

## 6. Linear matter growth

For non-relativistic matter on sub-horizon scales in a smooth background, the continuity, Euler and Poisson equations give $\ddot\delta+2H\dot\delta=4\pi G\bar\rho_m\delta$. With $D=D(\ln a)$ and $\dot{}=H\,d/d\ln a$:

$$
\ddot\delta+2H\dot\delta=H^2\left[D''+\left(2+\frac{d\ln H}{d\ln a}\right)D'\right],\qquad
\frac{d\ln H}{d\ln a}=-(1+q),
$$

so, using $4\pi G\bar\rho_m=\tfrac32H^2\Omega_m(a)$,

$$
D''+(1-q)D'-\tfrac32\Omega_m(a)D=0 .
$$

In an Einstein-de Sitter universe ($q=1/2$, $\Omega_m=1$) the growing solution is $D=a$. In the radiation era the Meszaros solution $D\propto1+\tfrac32 a/a_{\mathrm{eq}}$ applies. For $\Lambda$CDM without radiation, $D\propto E(a)\int_0^a da'/(a'E)^3$ (Heath 1977); the code is tested against this integral. (`perturbations/growth_factor.py`)

## 7. Scalar-field inflation

For a canonical field $\phi$ with potential $V$ in a flat FLRW universe ($M_P=1$):
$H^2=\tfrac13(\tfrac12\dot\phi^2+V)$, $\ddot\phi+3H\dot\phi+V_{,\phi}=0$, and $\dot H=-\tfrac12\dot\phi^2$, so $\epsilon_H\equiv-\dot H/H^2=\dot\phi^2/(2H^2)$.

**E-fold form.** With $N=\ln a$, $d/dt=H\,d/dN$ and primes denoting $d/dN$: $\dot\phi=H\phi'$, $\ddot\phi=H^2(\phi''-\epsilon_H\phi')$. The field equation becomes $\phi''+(3-\epsilon_H)\phi'+V_{,\phi}/H^2=0$, and since $H^2=V/(3-\epsilon_H)$,

$$
\phi''=-(3-\epsilon_H)\left(\phi'+\frac{V_{,\phi}}{V}\right),\qquad \epsilon_H=\tfrac12\phi'^2 .
$$

The shape of the trajectory is independent of the amplitude of $V$ (the amplitude is fixed afterwards from $A_s$). Inflation ends when $\epsilon_H=1$ ($\ddot a=0$).

**Slow-roll parameters.** Potential parameters: $\epsilon_V=\tfrac12(V_{,\phi}/V)^2$, $\eta_V=V_{,\phi\phi}/V$. Hubble parameters: $\epsilon_H$, $\eta_H\equiv d\ln\epsilon_H/dN$. To leading order $\epsilon_H\simeq\epsilon_V$ and $\eta_H\simeq4\epsilon_V-2\eta_V$.

**Observables.** $\mathcal P_{\mathcal R}=H^2/(8\pi^2\epsilon_HM_P^2)$ at $k=aH$, $\mathcal P_t=2H^2/(\pi^2M_P^2)$, so $r=16\epsilon_H$. Since $d\ln k=(1-\epsilon_H)dN$, $n_s-1=d\ln\mathcal P_{\mathcal R}/d\ln k\simeq-2\epsilon_H-\eta_H=-6\epsilon_V+2\eta_V$.

**Quadratic potential.** $V=\tfrac12m^2\phi^2$: $\epsilon_V=\eta_V=2/\phi^2$ and $N\simeq(\phi^2-\phi_e^2)/4$, giving $n_s\simeq1-2/N$, $r\simeq8/N$. **Starobinsky-type potential:** $n_s\simeq1-2/N$, $r\simeq12/N^2$. The code compares the exact background results to these limits in `tests/test_early_universe.py`.

## 8. Mukhanov-Sasaki equation in e-fold time

The comoving curvature perturbation $\mathcal R$ obeys, with $v=z\mathcal R$ and $z=a\sqrt{2\epsilon_H}M_P$, the equation $v''+(k^2-z''/z)v=0$ in conformal time $\tau$ ($d\tau=dt/a$). Equivalently $\mathcal R_{\tau\tau}+2(z_\tau/z)\mathcal R_\tau+k^2\mathcal R=0$ with $z_\tau/z=aH(1+\eta_H/2)$. Converting with $d/d\tau=aH\,d/dN$ and $d(aH)/dN=(1-\epsilon_H)aH$:

$$
\mathcal R''+(3-\epsilon_H+\eta_H)\mathcal R'+\left(\frac{k}{aH}\right)^2\mathcal R=0 .
$$

Tensor modes ($z\to a$): $h''+(3-\epsilon_H)h'+(k/aH)^2h=0$. Bunch-Davies initial conditions: $\mathcal R\to e^{-ik\tau}/(z\sqrt{2k})$ and $h\to2e^{-ik\tau}/(aM_P\sqrt{2k})$, the latter for each of the two polarisations; $\mathcal P_t=2\cdot k^3\lvert h\rvert^2/(2\pi^2)$. In exact de Sitter the result is $\mathcal P_{\mathcal R}=H^2/(8\pi^2\epsilon M_P^2)$ and $\mathcal P_t=2H^2/(\pi^2M_P^2)$; the numerical solver reproduces the scalar value to $4\times10^{-5}$ and the tensor value to within the $10^{-3}$ test tolerance (`tests/test_perturbations.py`).

## 9. Reheating bookkeeping

A homogeneous oscillating condensate with equation of state $w$ decaying at constant rate $\Gamma$:

$$
\dot\rho_\phi=-3H(1+w)\rho_\phi-\Gamma\rho_\phi,\qquad \dot\rho_r=-4H\rho_r+\Gamma\rho_\phi,\qquad 3M_P^2H^2=\rho_\phi+\rho_r .
$$

Adding the equations gives $\dot\rho=-3H(\rho+p)$ with $p=w\rho_\phi+\rho_r/3$, consistent with (C). The reheating temperature is defined by $\Gamma=H$: $\rho=3M_P^2\Gamma^2=\frac{\pi^2}{30}g_*T_{\mathrm{reh}}^4$, hence
$T_{\mathrm{reh}}=(90/\pi^2g_*)^{1/4}\sqrt{\Gamma M_P}$.

## 10. Thermal relations

For relativistic species in equilibrium: $\rho=\frac{\pi^2}{30}g_*T^4$, $s=\frac{2\pi^2}{45}g_{*s}T^3$. In radiation domination $H=\sqrt{8\pi G\rho/3}$ and $a\propto t^{1/2}$ so $t=1/(2H)=0.301\,g_*^{-1/2}M_{\mathrm{Pl}}/T^2$ ($M_{\mathrm{Pl}}=1.22\times10^{19}$ GeV). After $e^\pm$ annihilation $T_\nu=(4/11)^{1/3}T_\gamma$ and $g_*=2+\frac74(4/11)^{4/3}N_{\mathrm{eff}}$, $g_{*s}=2+\frac74\cdot\frac4{11}\cdot3=3.909$. Entropy conservation $g_{*s}T^3a^3=$ const links the reheating temperature to the scale factor at the end of inflation (`early_universe/inflation.py::pivot_efolds`).

## 11. Recombination equations

Saha equilibrium for hydrogen: $\dfrac{x_e^2}{1-x_e}=\dfrac1{n_H}\left(\dfrac{m_ekT}{2\pi\hbar^2}\right)^{3/2}e^{-B_1/kT}$. Peebles' three-level model replaces equilibrium by

$$
\frac{dx_e}{dt}=\mathcal C\left[\beta(1-x_e)e^{-h\nu_\alpha/kT}-n_H\alpha_Bx_e^2\right],\quad
\mathcal C=\frac{1+K\Lambda_{2\gamma}n_H(1-x_e)}{1+K(\Lambda_{2\gamma}+\beta)n_H(1-x_e)},\quad K=\frac{\lambda_\alpha^3}{8\pi H},
$$

with $\beta=\alpha_B(2\pi m_ekT/h^2)^{3/2}e^{-B_2/kT}$ by detailed balance, so that the equilibrium limit reproduces Saha. The optical depth is $\tau(z)=\int_0^z\sigma_Tn_ec\,dz'/[(1+z')H]$ and the visibility function $g(z)=(d\tau/dz)e^{-\tau}$.

## 12. De Sitter horizon temperature

For pure de Sitter space with Hubble rate $H_\Lambda=H_0\sqrt{\Omega_\Lambda}$ the cosmological horizon has radius $c/H_\Lambda$ and the Gibbons-Hawking temperature $T=\hbar H_\Lambda/(2\pi k_B)$ (an established semiclassical result for exact de Sitter; its role in the real universe is model dependent).

## 13. Phantom singularity

For $w<-1$ with dark energy dominating, $H^2=H_0^2\Omega_{\mathrm{DE}}a^{-3(1+w)}$ grows with $a$. Then

$$
t_{\mathrm{rip}}-t_0=\int_1^\infty\frac{da}{aH}=\frac{1}{H_0\sqrt{\Omega_{\mathrm{DE}}}}\int_1^\infty a^{p-1}\,da,\qquad p=\tfrac32(1+w)<0,
$$

which converges to $1/(\lvert p\rvert H_0\sqrt{\Omega_{\mathrm{DE}}})=2/(3\lvert1+w\rvert H_0\sqrt{\Omega_{\mathrm{DE}}})$. For $\Omega_{\mathrm{DE}}=1$ and $w=-3/2$ this equals $4/(3H_0)$ (tested). When matter is not negligible the estimate of Caldwell, Kamionkowski and Weinberg (2003), $t_{\mathrm{rip}}-t_0\simeq2/(3\lvert1+w\rvert H_0\sqrt{1-\Omega_m})$, is approximate; the repository finds it overshoots the numerical value by roughly 1 to 3 percent in the cases tested.
