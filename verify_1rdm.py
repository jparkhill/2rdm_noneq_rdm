"""Verification: does the MK-closed 1-RDM EOM thermalize to Fermi-Dirac,
conserve N, and preserve 0 <= gamma <= 1?  And what does dropping the MK hole
contractions do?"""
import numpy as np
from qdiss.rdm1 import DissipativeRDM1, fermi_dirac, solve_mu

rng = np.random.default_rng(20260802)
M, BETA = 10, 1.4

# random one-body hamiltonian
X = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
h = (X + X.conj().T) / 2 / np.sqrt(M) + np.diag(np.linspace(-2.0, 2.0, M))
# two independent bath coupling operators
cps = []
for _ in range(2):
    Y = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
    cps.append((Y + Y.conj().T) / 2)

def run(blocking, label):
    m = DissipativeRDM1(h, cps, BETA, eta=0.05, wc=6.0, blocking=blocking)
    # initial state: a far-from-equilibrium Slater determinant (top 4 orbitals filled)
    n0 = np.zeros(M); n0[[6, 7, 8, 9]] = 1.0
    g0 = np.diag(n0.astype(complex))
    g0[0, 1] = g0[1, 0] = 0.0
    t, gs = m.propagate(g0, 4000.0, n_out=400)

    N   = np.array([np.real(np.trace(g)) for g in gs])
    ev  = np.array([np.real(np.linalg.eigvalsh(g)) for g in gs])
    nd  = np.array([np.real(np.diag(g)) for g in gs])
    # exact Lyapunov functional of the fermionic Pauli master equation
    def F(n):
        n = np.clip(n, 1e-15, 1-1e-15)
        return m.eps@n + (n*np.log(n)+(1-n)*np.log(1-n)).sum()/BETA
    Fs = np.array([F(x) for x in nd])

    n_fin = nd[-1]
    mu = solve_mu(m.eps, BETA, N[-1])
    n_fd = fermi_dirac(m.eps, BETA, mu)
    coh = np.max(np.abs(gs[-1] - np.diag(np.diag(gs[-1]))))

    print(f"\n=== {label} ===")
    print(f"  N drift              : {abs(N-N[0]).max():.3e}   (N0 = {N[0]:.6f})")
    print(f"  min eigenvalue seen  : {ev.min():+.6e}")
    print(f"  max eigenvalue seen  : {ev.max():+.6e}")
    print(f"  residual coherence   : {coh:.3e}")
    print(f"  free energy monotone : {'YES' if np.all(np.diff(Fs) < 1e-9) else 'NO'}"
          f"   (dF total = {Fs[-1]-Fs[0]:+.6f})")
    print(f"  mu from N            : {mu:+.6f}")
    print(f"  ||n_final - n_FD||_inf = {np.abs(n_fin-n_fd).max():.3e}")
    print("   eps_p      n_final      n_FD(beta,mu)     diff")
    for e, a, b in zip(m.eps, n_fin, n_fd):
        print(f"  {e:+8.4f}  {a:12.8f}  {b:12.8f}  {a-b:+.2e}")
    return m, n_fin, n_fd

run(True,  "MK-closed (hole contractions kept)  -- reproduction target")
run(False, "hole contractions DROPPED (linear Redfield)")
