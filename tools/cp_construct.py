import sys, numpy as np
from scipy.optimize import minimize
sys.path.insert(0,"."); sys.path.insert(0,"tools")
src=open("tools/sdp_cpfit.py").read().split("nt = np.linalg.norm(target)")[0]
exec(src)
nt=np.linalg.norm(target)
H=0.5*(target+target.conj().T)
w,U=np.linalg.eigh(0.5*(rho+rho.conj().T))
Kb=U[:,w<1e-9]; Rb=U[:,w>=1e-9]; nr=Rb.shape[1]
PR=Rb@Rb.conj().T
NJ=nr*nr
def unpack(x):
    i=0
    kr=x[i:i+NP*NP].reshape(NP,NP); i+=NP*NP
    ki=x[i:i+NP*NP].reshape(NP,NP); i+=NP*NP
    KH=kr+kr.T+1j*(ki-ki.T)
    Ls=[]
    for _ in range(NJ):
        ar=x[i:i+nr*nr].reshape(nr,nr); i+=nr*nr
        ai=x[i:i+nr*nr].reshape(nr,nr); i+=nr*nr
        Ls.append(Rb@(ar+1j*ai)@Rb.conj().T)   # supported on range(rho)
    return KH,Ls
def gen(KH,Ls):
    d=-1j*(KH@rho-rho@KH)
    for L in Ls:
        LdL=L.conj().T@L
        d=d+L@rho@L.conj().T-0.5*(LdL@rho+rho@LdL)
    return d
def obj(x):
    KH,Ls=unpack(x); r=gen(KH,Ls)-target
    return float(np.real(np.vdot(r,r)))
npar=2*NP*NP+2*NJ*nr*nr
best=None
for t in range(3):
    x0=np.random.default_rng(200+t).normal(size=npar)*0.05
    res=minimize(obj,x0,method="L-BFGS-B",options={"maxiter":8000,"maxfun":120000})
    if best is None or res.fun<best.fun: best=res
    print("  trial %d residual %.6e (%.3f%%)" % (t, np.sqrt(max(res.fun,0)), 100*np.sqrt(max(res.fun,0))/nt))
KH,Ls=unpack(best.x); fit=gen(KH,Ls)
print()
print("CONSTRAINED CP fit (jumps supported on range(rho)): residual %.4e  = %.4f%% of |exact|" % (np.linalg.norm(fit-target), 100*np.linalg.norm(fit-target)/nt))
print("  d(tr rho^2)  exact %+.4e   fit %+.4e" % (2*np.trace(rho@target).real, 2*np.trace(rho@fit).real))
