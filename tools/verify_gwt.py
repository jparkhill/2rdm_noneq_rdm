"""
Independent numerical verification of the Mukherjee-Kutzelnigg generalized
normal ordering (GNO) and generalized Wick theorem (GWT) contraction rules,
in an explicit Fock space over M spin orbitals with a random correlated
N-electron reference.

Nothing here is taken on faith from the literature: every identity quoted in
docs/07_mkno_primary_sources.md that could be checked numerically is checked
here, in the stated index/sign convention.

Conventions (verified, not assumed):
    a^{p q ...}_{r s ...}  =  a+_p a+_q ... a_s a_r          (upper=creation,
                                                              lower=annihilation,
                                                              annihilators reversed)
    gamma^{p...}_{q...}    =  <Psi| a^{p...}_{q...} |Psi>
    eta^{p...}_{q...}      =  <Psi| a_{q1} a_{q2} ... a+_{p2} a+_{p1} |Psi>   (hole RDM)
Run: python3 tools/verify_gwt.py
"""
import itertools, numpy as np

rng = np.random.default_rng(20260802)
M = 6           # spin orbitals
N = 3           # electrons

# ---------- Fock space ----------
dim = 1 << M
def occ(s): return [i for i in range(M) if (s >> i) & 1]
def cre(p):
    """matrix of a+_p"""
    A = np.zeros((dim, dim))
    for s in range(dim):
        if (s >> p) & 1: continue
        t = s | (1 << p)
        sign = (-1) ** bin(s & ((1 << p) - 1)).count("1")
        A[t, s] = sign
    return A
Ad = [cre(p) for p in range(M)]
Aa = [Ad[p].T.copy() for p in range(M)]

# sanity: anticommutators
for p in range(M):
    for q in range(M):
        assert np.allclose(Aa[p] @ Ad[q] + Ad[q] @ Aa[p], (p == q) * np.eye(dim))
        assert np.allclose(Ad[p] @ Ad[q] + Ad[q] @ Ad[p], 0)

# ---------- random correlated N-electron reference ----------
idx = [s for s in range(dim) if bin(s).count("1") == N]
psi = np.zeros(dim); v = rng.normal(size=len(idx))
for k, s in enumerate(idx): psi[s] = v[k]
psi /= np.linalg.norm(psi)

def ev(Op): return float(psi @ (Op @ psi))

def astr(upper, lower):
    """a^{upper}_{lower} = a+_{u1}..a+_{uk} a_{lk}..a_{l1}"""
    Op = np.eye(dim)
    for p in upper: Op = Op @ Ad[p]
    for q in reversed(lower): Op = Op @ Aa[q]
    return Op

def hstr(upper, lower):
    """hole string: a_{l1} a_{l2} ... a+_{uk} ... a+_{u1}"""
    Op = np.eye(dim)
    for q in lower: Op = Op @ Aa[q]
    for p in reversed(upper): Op = Op @ Ad[p]
    return Op

R = range(M)
g1 = np.array([[ev(astr([p], [q])) for q in R] for p in R])
g2 = np.zeros((M,)*4); g3 = np.zeros((M,)*6)
h1 = np.array([[ev(hstr([p], [q])) for q in R] for p in R])
h2 = np.zeros((M,)*4); h3 = np.zeros((M,)*6)
for p, q, r, s in itertools.product(R, repeat=4):
    g2[p, q, r, s] = ev(astr([p, q], [r, s]))
    h2[p, q, r, s] = ev(hstr([p, q], [r, s]))
for t in itertools.product(R, repeat=6):
    g3[t] = ev(astr(list(t[:3]), list(t[3:])))
    h3[t] = ev(hstr(list(t[:3]), list(t[3:])))

ok = lambda name, c: print(f"  [{'OK ' if c else 'FAIL'}] {name}")

print("== 1-RDM basics ==")
ok("trace gamma = N", abs(np.trace(g1) - N) < 1e-10)
ok("eta = 1 - gamma (index-transposed)", np.allclose(h1, np.eye(M) - g1.T))
w = np.linalg.eigvalsh(g1)
ok("gamma non-idempotent (correlated)", not np.allclose(g1 @ g1, g1))
print(f"        natural occupations: {np.round(np.sort(w)[::-1], 4)}")

# ---------- (c) cumulant decomposition ----------
print("== (c) cumulant decomposition, MK convention ==")
l2 = g2 - (np.einsum('pr,qs->pqrs', g1, g1) - np.einsum('ps,qr->pqrs', g1, g1))
ok("gamma2 = g^p_r g^q_s - g^p_s g^q_r + lambda2 (defines l2)", True)

def wedge3_from(g, l2):
    """Return (A,B) with gamma3 = A + B + lambda3, where
       A = sum_{perm sigma of the 3 lower indices} sgn * g g g   [= (1/6) g^g^g]
       B = sum_{i,j} (-1)^(i+j) g^{p_i}_{q_j} lambda2^{rest}     [= lambda2 ^ g]
       Written as explicit loops on purpose: no einsum index gymnastics to trust."""
    import itertools as it
    A = np.zeros((M,) * 6)
    B = np.zeros((M,) * 6)
    perms3 = [(p, int(np.linalg.det(np.eye(3)[list(p)]).round()))
              for p in it.permutations(range(3))]
    for up in it.product(R, repeat=3):
        for lo in it.product(R, repeat=3):
            a = 0.0
            for pi, sg in perms3:
                a += sg * g[up[0], lo[pi[0]]] * g[up[1], lo[pi[1]]] * g[up[2], lo[pi[2]]]
            A[up + lo] = a
            b = 0.0
            for i in range(3):
                for j in range(3):
                    ru = [up[k] for k in range(3) if k != i]
                    rl = [lo[k] for k in range(3) if k != j]
                    b += (-1) ** (i + j) * g[up[i], lo[j]] * l2[ru[0], ru[1], rl[0], rl[1]]
            B[up + lo] = b
    return A, B

A3, B3 = wedge3_from(g1, l2)
l3 = g3 - A3 - B3
print(f"        |lambda2| = {np.linalg.norm(l2):.6f}   |lambda3| = {np.linalg.norm(l3):.6f}")
ok("lambda2 nonzero for a correlated reference", np.linalg.norm(l2) > 1e-6)
ok("lambda3 nonzero for a correlated reference", np.linalg.norm(l3) > 1e-6)
# partial trace test lambda^{pr}_{qr} summed = delta^p_q (n_q^2 - n_q) in NO basis
wnat, U = np.linalg.eigh(g1)
l2n = np.einsum('pi,qj,rk,sl,pqrs->ijkl', U, U, U, U, l2)
lhs = np.einsum('prqr->pq', l2n)
rhs = np.diag(wnat**2 - wnat)
ok("cumulant partial trace  sum_r l^{pr}_{qr} = d^p_q (n^2-n)", np.allclose(lhs, rhs, atol=1e-9))

# ---------- (e) particle-hole conjugation ----------
print("== (e) particle-hole conjugation of the cumulants ==")
l2h = h2 - (np.einsum('pr,qs->pqrs', h1, h1) - np.einsum('ps,qr->pqrs', h1, h1))
ok("hole 2-cumulant == particle 2-cumulant  (eta2 = eta^eta + lambda2)",
   np.allclose(l2h, l2, atol=1e-9))
A3h, B3h = wedge3_from(h1, l2)          # note: same lambda2, hole gamma1
l3h = h3 - A3h - B3h
ok("hole 3-cumulant == MINUS particle 3-cumulant (eta3 = -l3 + l2^eta + eta^3/6)",
   np.allclose(l3h, -l3, atol=1e-9))

# ---------- (a)/(b) GNO and GWT ----------
print("== (a) GNO defining property, (b) GWT contraction rules ==")
def gno1(p, q):  return astr([p], [q]) - g1[p, q] * np.eye(dim)
def gno2(p, q, r, s):
    O = astr([p, q], [r, s]) - g2[p, q, r, s] * np.eye(dim)
    O -= g1[p, r] * gno1(q, s); O += g1[p, s] * gno1(q, r)
    O += g1[q, r] * gno1(p, s); O -= g1[q, s] * gno1(p, r)
    return O
bad1 = max(abs(ev(gno1(p, q))) for p in R for q in R)
bad2 = max(abs(ev(gno2(p, q, r, s))) for p, q, r, s in itertools.product(R, repeat=4))
ok(f"<{{a^p_q}}> = 0        (max |.| = {bad1:.2e})", bad1 < 1e-10)
ok(f"<{{a^pq_rs}}> = 0      (max |.| = {bad2:.2e})", bad2 < 1e-10)

# GWT for {a^p_q}{a^r_s}
maxerr_full = 0.0; maxerr_nolam = 0.0
for p, q, r, s in itertools.product(R, repeat=4):
    lhs = gno1(p, q) @ gno1(r, s)
    rhs = (gno2(p, r, q, s)
           + h1[r, q] * gno1(p, s)
           - g1[p, s] * gno1(r, q)
           + (g1[p, s] * h1[r, q] + l2[p, r, q, s]) * np.eye(dim))
    maxerr_full = max(maxerr_full, np.abs(lhs - rhs).max())
    rhs_nl = rhs - l2[p, r, q, s] * np.eye(dim)      # engine with gamma/eta only
    maxerr_nolam = max(maxerr_nolam, np.abs(lhs - rhs_nl).max())
ok(f"GWT  {{a^p_q}}{{a^r_s}} = {{a^pr_qs}} + eta^r_q{{a^p_s}} - g^p_s{{a^r_q}} + g^p_s eta^r_q + lambda^pr_qs"
   f"   (max err {maxerr_full:.2e})", maxerr_full < 1e-10)
ok(f"SAME rule with the lambda2 term DROPPED is WRONG (max err {maxerr_nolam:.2e})",
   maxerr_nolam > 1e-6)

# ---------- (d) positivity ----------
print("== (d) positivity: is lambda2 the obstruction? ==")
def pairmat(T):
    """map antisymmetric T^{pq}_{rs} to a matrix on ordered pairs p<q"""
    prs = [(p, q) for p in R for q in R if p < q]
    return np.array([[T[p, q, r, s] for (r, s) in prs] for (p, q) in prs])
D  = pairmat(g2); Q = pairmat(h2)
Dq = pairmat(np.einsum('pr,qs->pqrs', g1, g1) - np.einsum('ps,qr->pqrs', g1, g1))
Qq = pairmat(np.einsum('pr,qs->pqrs', h1, h1) - np.einsum('ps,qr->pqrs', h1, h1))
ok(f"exact D = 2-RDM is PSD (min eig {np.linalg.eigvalsh(D).min():+.3e})",
   np.linalg.eigvalsh(D).min() > -1e-10)
ok(f"exact Q = 2-hole RDM is PSD (min eig {np.linalg.eigvalsh(Q).min():+.3e})",
   np.linalg.eigvalsh(Q).min() > -1e-10)
ok(f"lambda2=0 surrogate  D_qf = g^g  is PSD (min eig {np.linalg.eigvalsh(Dq).min():+.3e})",
   np.linalg.eigvalsh(Dq).min() > -1e-10)
ok(f"lambda2=0 surrogate  Q_qf = e^e  is PSD (min eig {np.linalg.eigvalsh(Qq).min():+.3e})",
   np.linalg.eigvalsh(Qq).min() > -1e-10)
# random-lambda stress test: is g^g + L ever indefinite for admissible-looking L?
worst = 0.0
for _ in range(400):
    Lr = rng.normal(size=(M,)*4) * 0.15
    Lr = Lr - Lr.transpose(1, 0, 2, 3); Lr = Lr - Lr.transpose(0, 1, 3, 2)
    Lr = 0.5 * (Lr + Lr.transpose(2, 3, 0, 1))
    worst = min(worst, np.linalg.eigvalsh(Dq + pairmat(Lr)).min())
ok(f"adding an arbitrary antisym/hermitian lambda2 CAN break D>=0 (worst min eig {worst:+.3e})",
   worst < -1e-8)

# theorem check: 0 <= gamma <= 1 => gamma ^ gamma PSD, for many random gammas
bad = 0
for _ in range(300):
    X = rng.normal(size=(M, M)); S = X @ X.T
    G = S / (np.linalg.eigvalsh(S).max() * (1 + rng.random()))   # 0 <= G <= 1
    W = pairmat(np.einsum('pr,qs->pqrs', G, G) - np.einsum('ps,qr->pqrs', G, G))
    if np.linalg.eigvalsh(W).min() < -1e-10: bad += 1
ok("for every 0<=gamma<=1 tested, gamma^gamma is PSD (300 random draws)", bad == 0)
# ---------- quasi-normal ordering (Kong-Nooijen-Mukherjee) vs. true GNO ----------
print("== gamma/eta-only engine == quasi-normal ordering: does <{..}> still vanish? ==")
def qno2(p, q, r, s):
    """'{a^pq_rs}' as produced by an engine whose only contractions are gamma and eta:
       the fully-contracted piece is the gamma-wedge-gamma product, NOT the true 2-RDM."""
    O = astr([p, q], [r, s]) - (g1[p, r] * g1[q, s] - g1[p, s] * g1[q, r]) * np.eye(dim)
    O -= g1[p, r] * gno1(q, s); O += g1[p, s] * gno1(q, r)
    O += g1[q, r] * gno1(p, s); O -= g1[q, s] * gno1(p, r)
    return O
resid = np.array([[[[ev(qno2(p, q, r, s)) for s in R] for r in R] for q in R] for p in R])
ok(f"quasi-normal-ordered <{{a^pq_rs}}> = lambda2 exactly (max dev {np.abs(resid - l2).max():.2e})",
   np.allclose(resid, l2, atol=1e-10))
print(f"        i.e. the GNO defining property <{{...}}>=0 FAILS by exactly lambda2;")
print(f"        max |residual| here = {np.abs(resid).max():.4f}")

# ---------- the exactness boundary: quasi-free (Gaussian) references ----------
print("== exactness boundary: is gamma/eta-only Wick exact for a NON-IDEMPOTENT gamma? ==")
eps = rng.normal(size=M)                       # mean-field orbital energies
nocc = 1.0 / (1.0 + np.exp(eps))               # Fermi-Dirac occupations, beta=1
w = np.array([np.prod([nocc[i] if (s >> i) & 1 else 1 - nocc[i] for i in range(M)])
              for s in range(dim)])
rho = np.diag(w)
tev = lambda Op: float(np.trace(rho @ Op))
gG = np.array([[tev(astr([p], [q])) for q in R] for p in R])
gG2 = np.zeros((M,)*4); gG3 = np.zeros((M,)*6)
for p, q, r, s in itertools.product(R, repeat=4):
    gG2[p, q, r, s] = tev(astr([p, q], [r, s]))
for t in itertools.product(R, repeat=6):
    gG3[t] = tev(astr(list(t[:3]), list(t[3:])))
lG2 = gG2 - (np.einsum('pr,qs->pqrs', gG, gG) - np.einsum('ps,qr->pqrs', gG, gG))
AG, BG = wedge3_from(gG, lG2)
lG3 = gG3 - AG - BG
ok("mean-field thermal state: gamma is NON-idempotent", not np.allclose(gG @ gG, gG))
ok(f"mean-field thermal state: lambda2 = 0  (|l2| = {np.linalg.norm(lG2):.2e})",
   np.linalg.norm(lG2) < 1e-10)
ok(f"mean-field thermal state: lambda3 = 0  (|l3| = {np.linalg.norm(lG3):.2e})",
   np.linalg.norm(lG3) < 1e-10)
print("        => gamma/eta-only Wick is EXACT for quasi-free states, idempotent or not;")
print("           it is the lambda_k = 0 approximation for any correlated reference.")

print("\nDone.")
