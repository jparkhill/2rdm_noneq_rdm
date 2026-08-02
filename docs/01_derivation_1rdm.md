# Reproduction of the dissipative 1-RDM EOM from Mukherjee–Kutzelnigg normal ordering

**Status:** derived from first principles and verified numerically (`verify_1rdm.py`).
Independent of the published version — to be reconciled against the literature dossier.

## 1. Model

$$H = H_S + H_B + H_{SB},\qquad
H_S=\sum_{pq}h_{pq}a^\dagger_p a_q,\qquad
H_B=\sum_k \omega_k b^\dagger_k b_k,$$

$$H_{SB}=\sum_\alpha A_\alpha\otimes B_\alpha,\qquad
A_\alpha=\sum_{pq}A^\alpha_{pq}a^\dagger_p a_q,\qquad
B_\alpha=\sum_k c_{\alpha k}(b_k+b^\dagger_k).$$

The system side of the coupling is one-body and number-conserving, so
$N=\operatorname{tr}\gamma$ is conserved and the expected stationary state is
Fermi–Dirac at the bath temperature with $\mu$ fixed by $N$.

Bath correlation spectrum
$$\Gamma_\alpha(\omega)=\int dt\,e^{i\omega t}\langle B_\alpha(t)B_\alpha(0)\rangle
=2\pi J_\alpha(\omega)\,[\,\bar n(\omega)+1\,],$$
which obeys the KMS / detailed-balance identity
$$\boxed{\;\Gamma(-\omega)=e^{-\beta\omega}\Gamma(\omega)\;}$$
This single identity is the *only* bath property the Fermi–Dirac fixed point needs.

## 2. Where the hierarchy fails to close

With $\gamma_{pq}=\langle a^\dagger_q a_p\rangle$, the coherent part is exact and trivial:
$$i\,\dot\gamma = [h,\gamma].$$

The TCL2 (Born–Markov, second-order time-convolutionless) dissipator built from a
**quadratic** $A_\alpha$ is *quartic* in fermion operators. Contracting it down to the
1-RDM therefore produces $\langle a^\dagger_p a^\dagger_q a_s a_r\rangle$: the equation of
motion for $\gamma$ is the first member of a BBGKY-like hierarchy and does **not** close.
This is the structural fact that makes the whole programme necessary.

## 3. The MK closure

Generalized (Mukherjee–Kutzelnigg) normal ordering is taken with respect to the
*instantaneous, correlated, non-idempotent* $\gamma$, so the elementary contractions are

$$\langle a^\dagger_p a_q\rangle=\gamma_{qp},\qquad
\langle a_p a^\dagger_q\rangle=\delta_{pq}-\gamma_{pq},$$

and the extended Wick theorem gives the cumulant decomposition

$$^2\Gamma_{pq,rs}=\gamma_{pr}\gamma_{qs}-\gamma_{ps}\gamma_{qr}+\lambda_{pq,rs}.$$

**Closure: $\lambda_2 = 0$.** With an idempotent reference this collapses to ordinary
Wick; the MK version is what permits *fractional* occupations, and that is exactly what
lets the fixed point be Fermi–Dirac instead of a determinant.

## 4. Result

In the eigenbasis of $h$ ($h|p\rangle=\varepsilon_p|p\rangle$, $n_p=\gamma_{pp}$), with
the secular approximation:

$$\boxed{\;\dot n_p=\sum_q\Big[\,W_{p\leftarrow q}\,n_q(1-n_p)\;-\;W_{q\leftarrow p}\,n_p(1-n_q)\,\Big]\;}$$
$$W_{p\leftarrow q}=\sum_\alpha\big|A^\alpha_{pq}\big|^2\,\Gamma_\alpha(\varepsilon_q-\varepsilon_p),$$
$$\dot\gamma_{pq}\big|_{p\neq q}=-i(\varepsilon_p-\varepsilon_q)\gamma_{pq}
-\tfrac12(\Gamma_p+\Gamma_q)\gamma_{pq}.$$

**The factors $(1-n_p)$ are literally the MK hole contractions** $\langle a_p a^\dagger_p\rangle
= 1-\gamma_{pp}$. They are not inserted by hand; they fall out of normal ordering about
the correlated reference. Everything below follows from them.

### 4.1 Fermi–Dirac stationarity
Setting each bracket to zero and using KMS with $\Delta=\varepsilon_q-\varepsilon_p$:
$$\frac{W_{q\leftarrow p}}{W_{p\leftarrow q}}=\frac{\Gamma(-\Delta)}{\Gamma(\Delta)}=e^{-\beta\Delta}
\;\Longrightarrow\;
\frac{n_q/(1-n_q)}{n_p/(1-n_p)}=e^{-\beta(\varepsilon_q-\varepsilon_p)}
\;\Longrightarrow\;
\frac{n_p}{1-n_p}=e^{-\beta(\varepsilon_p-\mu)},$$
$$n_p=\frac{1}{1+e^{\beta(\varepsilon_p-\mu)}},\qquad \textstyle\sum_p n_p=N .$$

### 4.2 N-representability is automatic at rank 1 — *and here is why*
The 1-RDM N-representability conditions are just $0\le n_p\le 1$. The flow makes that box
**invariant**, because the gain term carries a factor $n$ and the loss term a factor $(1-n)$:

$$n_p=0\;\Rightarrow\;\dot n_p=\sum_q W_{p\leftarrow q}n_q\;\ge 0,\qquad
n_p=1\;\Rightarrow\;\dot n_p=-\sum_q W_{q\leftarrow p}(1-n_q)\;\le 0.$$

**The flux through the boundary of the representable set vanishes identically on that
boundary.** This is the property we must generalise to rank 2. It is a purely algebraic,
checkable statement — a good Lean target.

### 4.3 Lyapunov function
$F[n]=\sum_p\varepsilon_p n_p+\beta^{-1}\sum_p[n_p\ln n_p+(1-n_p)\ln(1-n_p)]$ decreases
monotonically and is minimised, at fixed $N$, exactly by the Fermi–Dirac $n$.

## 5. Numerical verification

`verify_1rdm.py`, $M=10$ spin-orbitals, random $h$ and two random Hermitian couplings,
ohmic bath ($\eta=0.05$, $\omega_c=6$), $\beta=1.4$, $N=4$, started from a
far-from-equilibrium Slater determinant (top four orbitals filled).

| diagnostic | MK closure kept | hole contractions dropped |
|---|---|---|
| drift in $N$ | $2.2\times10^{-14}$ | $3.0\times10^{-13}$ |
| min eigenvalue of $\gamma$ | $0.000000$ | $0.000000$ |
| **max eigenvalue of $\gamma$** | $1.000000$ | $\mathbf{2.8135}$ ← **N-representability violated** |
| $F$ monotone decreasing | **yes** | no |
| $\|n_\infty-n_{FD}(\beta,\mu)\|_\infty$ | $\mathbf{1.1\times10^{-10}}$ | $1.86$ |

$\mu=-0.5510$ recovered from $N$ alone; the converged occupations match Fermi–Dirac to
$10^{-10}$ across the whole spectrum, i.e. to integrator tolerance.

Dropping the MK hole contractions gives a Boltzmann-like steady state and drives the
lowest occupation to **2.81 electrons in one spin-orbital**. The same mechanism that
delivers Fermi–Dirac is the one that enforces representability. That is the thesis of
this project, verified at rank 1.
