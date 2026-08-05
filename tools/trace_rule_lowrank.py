import sys
import numpy as np
sys.path.insert(0, ".")
from qdiss.eom import fullci
from qdiss.symbolic import mkno

def wedge(a, b):
    return np.einsum("pq,vw->pvqw", a, b) - np.einsum("pw,vq->pvqw", a, b)

def anti1(A, X):
    return A @ X + X @ A

def anti2(A, X):
    return (np.einsum("py,yvqw->pvqw", A, X) + np.einsum("vy,pyqw->pvqw", A, X)
            + np.einsum("pvyw,yq->pvqw", X, A) + np.einsum("pvqy,yw->pvqw", X, A))

T3 = [(s, [("".join(u), "".join(l)) for n, r, u, l in f], [(n, r) for n, r, u, l in f])
      for s, f in mkno.gamma_to_cumulants(("p","v","x"), ("q","w","y"))
      if not any(n == "L" and r == 3 for n, r, u, l in f)]

print("Is the rank-2 identity also anti(E, lam_{k-1}) plus a lower-cumulant CHAIN?")
print("%8s | %13s %13s %13s" % ("(M,N)", "k=2 resid", "k=3 resid", "|target|"))
for M, N in ((6, 3), (8, 4), (6, 4)):
    rng = np.random.default_rng(9)
    psi = fullci.random_state(M, N, rng)
    g1 = fullci.rdm(psi, M, 1)
    g2 = fullci.rdm(psi, M, 2)
    g3 = fullci.rdm(psi, M, 3)
    l2 = g2 - wedge(g1, g1)
    tot = np.zeros((M,)*6, dtype=complex)
    for s, pairs, keys in T3:
        tot = tot + s * np.einsum(",".join(u+l for u, l in pairs) + "->pvxqwy",
                                  *[{("G",1): g1, ("L",2): l2}[k] for k in keys], optimize=True)
    l3 = g3 - tot
    I = np.eye(M, dtype=complex)
    E = g1 - 0.5 * I
    # k=2 : sum_x lam2 = anti1(E, lam1) - lam1 chain lam1,  lam1 = g1
    k2 = np.abs(np.einsum("pxqx->pq", l2) - (anti1(E, g1) - g1 @ g1)).max()
    # k=3 : sum_y lam3 = anti2(E, lam2), no chain possible
    k3 = np.abs(np.einsum("pvyqwy->pvqw", l3) - anti2(E, l2)).max()
    print("%8s | %13.3e %13.3e %13.4f"
          % ("(%d,%d)" % (M, N), k2, k3, np.abs(np.einsum("pvyqwy->pvqw", l3)).max()))
print()
print("k=2  sum_x lam2 = anti(E, lam1) - (lam1 * lam1)          chain coeff -1")
print("k=3  sum_y lam3 = anti(E, lam2)                          no chain exists")
print("k=4  sum_z lam4 = anti(E, lam3) + (9/2)(chainA - chainB) chain of lam2 * lam2")
print()
print("The chain at rank k is lam_{k/2} * lam_{k/2}: it exists for k=2 and k=4 and there is")
print("no lam_{3/2} for k=3, which is exactly where the fit found no chain.")
