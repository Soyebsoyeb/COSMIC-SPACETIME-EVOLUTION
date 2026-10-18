# Big Bang Nucleosynthesis and Recombination

## Part A: Big Bang Nucleosynthesis

**Module status: interface only. The repository does not compute nuclear abundances.** Predicting the abundances of D, ${}^3$He, ${}^4$He and ${}^7$Li requires a coupled network of weak and nuclear rates (as in PArthENoPE, AlterBBN or PRIMAT), which is not implemented. Any BBN abundance figure in this repository would be invented, so none is shown.

### A.1 What is implemented (established relations)

1. The neutron-to-proton ratio in weak equilibrium, $n/p=\exp(-Q/T)$ with $Q=m_n-m_p=1.29333$ MeV.
2. The helium mass fraction as an identity when essentially all neutrons are bound into ${}^4$He: with $x=n/p$ at the onset of nucleosynthesis,

$$
Y_p=\frac{2x}{1+x}.
$$

Here $x$ is an **input** (typically about $1/7$ after neutron decay, giving $Y_p\approx0.25$), not a result of the code. `fig23_bbn_relations` plots these relations and is labelled accordingly.

3. A backend hook, `early_universe/bbn_interface.py`, raises `BBNBackendUnavailable` when a network backend is requested. A user can plug in an external code; none is shipped.

The observationally inferred primordial helium fraction and deuterium abundance (for example the review by Cyburt et al., arXiv:1505.01076) are cited as external constraints but are not reproduced here.

## Part B: Recombination

Status: **approximate**. Implemented in `early_universe/recombination.py`.

### B.1 Model

Peebles' three-level hydrogen atom (equations in [MATHEMATICAL_DERIVATIONS.md](MATHEMATICAL_DERIVATIONS.md), section 11) with a case-B recombination coefficient from a RECFAST-style fit including the fudge factor $F=1.14$. The Saha solution is provided for comparison and is accurate early on but wrongly predicts the freeze-out. The coefficient values used are *recalled* from the literature and have not been re-checked against the source papers here; the code is tested against two anchor values of $\alpha_B$ ($2.59\times10^{-13}$ cm$^3$/s at $10^4$ K and $1.43\times10^{-13}$ cm$^3$/s at $2\times10^4$ K).

### B.2 Outputs

| Quantity | This code | Reference |
|---|---:|---|
| Visibility peak $z_*$ | 1078.9 | Planck 2018 reports $z_*\approx1090$ (from a full recombination code) |
| Comoving sound horizon at $z_*$ | 145.2 Mpc | |
| $100\,\theta_*=100\,r_s/D_M$ | 1.0475 | $1.0409\pm0.0003$ ($100\theta_{MC}$, Planck 2018 VI) |

The 0.63 percent difference in $\theta_*$ is a **model approximation**: the Peebles model omits helium recombination, additional two-photon and Lyman-alpha radiative-transfer corrections, and uses $\Omega_r$ with massless neutrinos. It should not be mistaken for a numerical error, and it is more than ten times the Planck measurement uncertainty. Figure: `fig16_recombination`.

### B.3 Limits

No helium recombination, no reionisation, no dark-matter annihilation or other energy injection, no Compton-drag or Thomson-scattering polarisation physics, and no CMB anisotropy spectrum. The code produces the *ingredients* (visibility function and sound horizon), not $C_\ell$.
