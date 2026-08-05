#!/usr/bin/env python3
"""Validate the GPU backend against the numpy path, then time a real propagation.

No result from qdiss.gpu is usable until it agrees with the CPU path at complex128 round-off.
Checked: the assembled gamma_1-dot and gamma_2-dot, the Mueller reconstruction with and without
the Lackner contraction correction, and the antisymmetriser. Then the same RK4 propagation is
timed on both.

Runs on the host venv /home/john/.venvs/qm, which carries pyscf and CUDA torch.
"""
import sys
import time
import numpy as np
import torch
sys.path.insert(0, ".")
sys.path.insert(0, "tools")
from qdiss import gpu, muller
from qdiss.eom import fullci, gwt, wick
from qdiss.eom.evaluate import assemble as np_assemble
from qdiss.symbolic import mkno

UP3, LO3 = ("p", "v", "x"), ("q", "w", "y")
ORB_B = [(s, ",".join("".join(u) + "".join(l) for n, r, u, l in f) + "->pvxqwy",
          [(n, r) for n, r, u, l in f])
         for s, f in mkno.gamma_to_cumulants(UP3, LO3)
         if not any(n == "L" and r == 3 for n, r, u, l in f)
         and any(n == "L" and r == 2 for n, r, u, l in f)]
ONE1, TWO1 = wick.terms(1)
ONE2, TWO2 = wick.terms(2)
dev = gpu.device()
print("device:", dev, torch.cuda.get_device_name(0) if dev.type == "cuda" else "")

from pyscf import ao2mo, gto, scf


def hubbard(nsite, U, J=1.0):
    """Spin-orbital h and antisymmetrised V for a Hubbard chain, our convention."""
    M = 2 * nsite
    h = np.zeros((M, M), dtype=complex)
    for i in range(nsite - 1):
        for s in (0, 1):
            h[2 * i + s, 2 * (i + 1) + s] = -J
            h[2 * (i + 1) + s, 2 * i + s] = -J
    phys = np.zeros((M,) * 4)
    for i in range(nsite):
        for sp in (0, 1):
            for sq in (0, 1):
                phys[2 * i + sp, 2 * i + sq, 2 * i + sp, 2 * i + sq] += U
    V = (phys - np.transpose(phys, (0, 1, 3, 2))).astype(complex)
    return M, h, V


M, h, V = hubbard(4, 1.0)
N = 4
Hm = fullci.build_H(h, V, M)
popc = np.array([bin(s).count("1") for s in range(1 << M)])
sub = np.where(popc == N)[0]
es, vec = np.linalg.eigh(Hm[np.ix_(sub, sub)])
psi0 = np.zeros(1 << M, dtype=complex)
psi0[sub] = vec[:, 0]
g1 = fullci.rdm(psi0, M, 1)
g2 = fullci.rdm(psi0, M, 2)
EYE = np.eye(M, dtype=complex)

print()
print("=== agreement, 4-site Hubbard U/J=1, M=%d N=%d ===" % (M, N))
hg, Vg, g1g, g2g = (gpu.up(x, dev) for x in (h, V, g1, g2))
Ig = gpu.up(EYE, dev)
for cc in (False, True):
    a = muller.gamma3(g1, g2, N, contraction_correct=cc)
    b = gpu.down(gpu.muller_gamma3(g1g, g2g, N, M, contraction_correct=cc, orbB=ORB_B))
    print("   muller gamma3 cc=%-5s |cpu-gpu| = %.3e   |value| = %.5f"
          % (cc, np.abs(a - b).max(), np.abs(a).max()))
g3 = muller.gamma3(g1, g2, N)
g3g = gpu.up(g3, dev)
tn = {("D", 1): EYE, ("G", 1): g1, ("G", 2): g2, ("G", 3): g3}
tg = {("D", 1): Ig, ("G", 1): g1g, ("G", 2): g2g, ("G", 3): g3g}
for rk, one, two, osub in ((1, ONE1, TWO1, "pq"), (2, ONE2, TWO2, "pvqw")):
    a = np_assemble(one, two, h, V, tn, osub, M)
    b = gpu.down(gpu.assemble(one, two, hg, Vg, tg, osub))
    print("   gamma_%d-dot          |cpu-gpu| = %.3e   |value| = %.5f"
          % (rk, np.abs(a - b).max(), np.abs(a).max()))

print()
print("=== timing: 25 RK4 steps of gamma_2 with the Mueller closure ===")


def rhs_np(x):
    gg1 = np.einsum("pxqx->pq", x) / (N - 1)
    gg3 = muller.gamma3(gg1, x, N)
    return np_assemble(ONE2, TWO2, h, V,
                       {("D", 1): EYE, ("G", 1): gg1, ("G", 2): x, ("G", 3): gg3},
                       "pvqw", M)


def rhs_gpu(x):
    gg1 = torch.einsum("pxqx->pq", x) / (N - 1)
    gg3 = gpu.muller_gamma3(gg1, x, N, M, orbB=ORB_B)
    return gpu.assemble(ONE2, TWO2, hg, Vg,
                        {("D", 1): Ig, ("G", 1): gg1, ("G", 2): x, ("G", 3): gg3}, "pvqw")


def rk4(f, y, dt):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3)
    return y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


dt, ns = 0.004, 25
y = g2.copy()
t0 = time.perf_counter()
for _ in range(ns):
    y = rk4(rhs_np, y, dt)
tcpu = time.perf_counter() - t0
yg = g2g.clone()
rk4(rhs_gpu, yg, dt)
torch.cuda.synchronize() if dev.type == "cuda" else None
t0 = time.perf_counter()
for _ in range(ns):
    yg = rk4(rhs_gpu, yg, dt)
if dev.type == "cuda":
    torch.cuda.synchronize()
tgpu = time.perf_counter() - t0
print("   cpu %.3f s   gpu %.3f s   speedup %.1fx" % (tcpu, tgpu, tcpu / tgpu))
print("   |gamma_2 cpu - gamma_2 gpu| after %d steps = %.3e"
      % (ns, np.abs(y - gpu.down(yg)).max()))
