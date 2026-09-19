#!/usr/bin/env python3
r"""Route F, episode 8: poles of the extended (X, delta) system with a dynamical delta.

Same system/recipe as tools/f0delta_poles.py (H4 jitter geometry, sto-3g, M = 8, N = 4, six lowest FCI singlets,
biorthogonal matching on two-body transition densities). The extended Jacobian is over (X, delta_c), 2*n2 = 1568
real variables; delta is first relaxed to its stationary value at the ground state (fixed point of the delta equation
with X frozen), so the linearisation is at a genuine stationary point of the extended flow. Targets are matched on their
X components only. Also reported: the spectrum of the delta sector (fictitious modes) -- real parts ~ -1/tau -- and the
largest |Re w| among modes with appreciable X weight.
"""
import sys, time
import numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "tools")
import f0delta_poles as fp
from f0delta import F0Delta, DynamicDelta
from commutator_k3 import from_mat2
from qdiss.eom import fullci
from scipy.optimize import linear_sum_assignment

tau, kappa, sigma = [float(x) for x in sys.argv[1:4]] if len(sys.argv) > 3 else (0.5, 10.0, 1e-3)
Hs = fp.op_matrix(fp.Hfn); Hs = 0.5 * (Hs + Hs.conj().T)
S2s = fp.op_matrix(fp.S2fn); S2s = 0.5 * (S2s + S2s.conj().T)
es, vec = np.linalg.eigh(Hs)
s2 = np.real(np.einsum("ik,ij,jk->k", vec.conj(), S2s, vec)); sing = np.where(np.abs(s2) < 1e-6)[0]
exc = (es - es[sing[0]]).real; sing_exc = np.array([exc[k] for k in sing[1:]])
keep = np.argsort(sing_exc)[:6]; targets = sing_exc[keep]; tstate = [sing[1:][k] for k in keep]
H2 = []
for k in tstate:
    T = fp.trans2(vec[:, sing[0]], vec[:, k]); H2.append(0.5 * (T + np.transpose(T.conj(), (2, 3, 0, 1))))
M, N = fp.M, fp.N
base = F0Delta(fp.h, fp.V, M, N); dd = DynamicDelta(base, tau=tau, kappa=kappa, sigma=sigma)
g2_0 = fullci.rdm(fp.embed(vec[:, sing[0]]), M, 2)
# relax delta to its fixed point at frozen X (Picard on the linear-in-delta force; the force is affine in delta_c
# since grad Phi is linear in the flow): solve by iterating tau*dc/dt = 0
# the force is affine in delta_c: force(dc) = f0 + A dc. Build A by finite differences (n2 evaluations) and solve.
dc0 = np.zeros(base.n2); _, fz = dd.rhs(g2_0, dc0)
A = np.zeros((base.n2, base.n2)); t0 = time.time()
for k in range(base.n2):
    e = np.zeros(base.n2); e[k] = 1e-4
    A[:, k] = (dd.rhs(g2_0, e)[1] - fz) / 1e-4
dc = np.linalg.solve(A, -fz); _, f = dd.rhs(g2_0, dc)
print("delta fixed point (linear solve, %.0fs): |force| %.1e, |delta| %.4f, Phi %.2e" % (time.time() - t0, np.linalg.norm(f), np.linalg.norm(dc), dd.last_phi))
n2 = fp.n2
def rhs_vec(y):
    g2 = y[:n2] @ np.array(fp.B2).reshape(n2, -1); g2 = g2.reshape((M,) * 4) if False else None
def unpack(y):
    return sum(c * b for c, b in zip(y[:n2], fp.B2)) if False else None
# coordinates: X in the B2 basis (fp.B2 are 4-index tensors), delta in base.B coefficients
B2 = np.array(fp.B2); nrm = fp.nrm
def to_coords(g2):
    return np.array([float(np.real(np.tensordot(B2[m].conj(), g2, axes=4))) for m in range(n2)]) / nrm
def from_coords(c):
    return np.einsum("k,kpvqw->pvqw", c, B2)
x0 = to_coords(g2_0); y0 = np.concatenate([x0, dc])
def F(y):
    g2 = from_coords(y[:n2]); gd, dcd = dd.rhs(g2, y[n2:])
    return np.concatenate([to_coords(gd), dcd])
f0 = F(y0); print("|F(y0)| = %.2e (X part %.2e, delta part %.2e)" % (np.linalg.norm(f0), np.linalg.norm(f0[:n2]), np.linalg.norm(f0[n2:])))
J = np.zeros((2 * n2, 2 * n2)); t0 = time.time()
for k in range(2 * n2):
    e = np.zeros(2 * n2); e[k] = fp.EPS
    J[:, k] = (F(y0 + e) - f0) / fp.EPS
print("Jacobian %.0fs" % (time.time() - t0))
w, R = np.linalg.eig(J); Rinv = np.linalg.inv(R)
xw = np.linalg.norm(R[:n2], axis=0) / np.linalg.norm(R, axis=0)          # X-weight of each mode
print("modes: %d ; with X-weight > 0.5: %d ; Re w range of those [%.2e, %.2e] ; delta-sector Re w median %.3f (expect ~ -1/tau = %.3f)"
      % (len(w), (xw > 0.5).sum(), w.real[xw > 0.5].min(), w.real[xw > 0.5].max(), np.median(w.real[xw <= 0.5]), -1 / tau))
keepm = np.where(w.imag > 1e-6)[0]
S = np.zeros((6, len(keepm)))
for a in range(6):
    x = np.concatenate([to_coords(H2[a]), np.zeros(n2)]); c = Rinv @ x
    for b, i in enumerate(keepm): S[a, b] = abs(c[i]) * np.linalg.norm(R[:, i]) / np.linalg.norm(x)
r, cc = linear_sum_assignment(-S)
row = np.zeros(6); sc = np.zeros(6); re = np.zeros(6)
for a, b in zip(r, cc): row[a] = w[keepm[b]].imag; sc[a] = S[a, b]; re[a] = w[keepm[b]].real
print("FCI          | %s" % "  ".join("%8.4f" % x for x in targets))
print("dyn-delta    | %s  mean|err| %.4f" % ("  ".join("%8.4f" % x for x in row), np.mean(np.abs(row - targets))))
print("             | %s  <- weight" % "  ".join("%8.3f" % x for x in sc))
print("             | %s  <- Re w" % "  ".join("%+8.1e" % x for x in re))
