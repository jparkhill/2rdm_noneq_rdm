import numpy as np, itertools
def build(M):
    def ad(p,v):
        out=np.zeros_like(v); bit=1<<p
        for s in np.nonzero(v)[0]:
            if s&bit: continue
            out[s|bit]+=((-1.0)**bin(s&(bit-1)).count("1"))*v[s]
        return out
    def aa(p,v):
        out=np.zeros_like(v); bit=1<<p
        for s in np.nonzero(v)[0]:
            if not(s&bit): continue
            out[s^bit]+=((-1.0)**bin(s&(bit-1)).count("1"))*v[s]
        return out
    return ad,aa

print("%7s %10s %10s %10s | %9s %9s %9s | %9s"%("(M,N)","tr g2 raw","N(N-1)","trD pair","ours","comb","pairnorm","N/2"))
for M,N in ((4,2),(6,2),(8,4),(10,4),(12,6)):
    dim=1<<M; ad,aa=build(M); K,n=M//2,N//2
    vac=np.zeros(dim,dtype=complex); vac[0]=1.0
    psi=vac
    for _ in range(n):
        acc=np.zeros_like(psi)
        for k in range(K): acc=acc+ad(2*k,ad(2*k+1,psi))
        psi=acc
    psi=psi/np.linalg.norm(psi)
    # raw gamma_2 exactly as tools/*.py rdm(psi,2) builds it
    g2=np.zeros((M,)*4,dtype=complex)
    for p,v,q,w in itertools.product(*([range(M)]*4)):
        g2[p,v,q,w]=np.vdot(psi, ad(p,ad(v,aa(w,aa(q,psi)))))
    tr_raw=float(np.real(np.einsum("pvpv->",g2)))
    PAIRS=[(p,v) for p in range(M) for v in range(p+1,M)]
    D=np.array([[2.0*g2[p,v,q,w] for (q,w) in PAIRS] for (p,v) in PAIRS])
    D=0.5*(D+D.conj().T)
    ev=float(np.linalg.eigvalsh(D).max())
    print("%7s %10.4f %10d %10.4f | %9.4f %9.4f %9.4f | %9.1f"
          %("(%d,%d)"%(M,N),tr_raw,N*(N-1),float(np.real(np.trace(D))),
            ev, ev/2, ev/(2*(N-1)), N/2))
