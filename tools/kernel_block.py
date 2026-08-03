import sys, numpy as np
sys.path.insert(0,"."); sys.path.insert(0,"tools")
src=open("tools/sdp_cpfit.py").read().split("nt = np.linalg.norm(target)")[0]
exec(src)
nt=np.linalg.norm(target)
H=0.5*(target+target.conj().T)
w,U=np.linalg.eigh(0.5*(rho+rho.conj().T))
K=U[:,w<1e-9]; R=U[:,w>=1e-9]
XKK = K.conj().T @ H @ K
XKR = K.conj().T @ H @ R
ev=np.linalg.eigvalsh(0.5*(XKK+XKK.conj().T))
print("pair dim %d   rank rho %d   dim ker %d" % (NP, R.shape[1], K.shape[1]))
print("kernel block X_KK: min eig %+.3e  max eig %+.3e  |X_KK| %.3e" % (ev.min(), ev.max(), abs(XKK).max()))
print("off-diagonal block |X_KR| = %.3e" % abs(XKR).max())
print("|exact rho-dot| = %.4f" % nt)
print()
print("GKSL image condition at rho is <v|X|v> >= 0 on ker rho, i.e. X_KK PSD.")
print("verdict:", "SATISFIED (X_KK is PSD)" if ev.min() > -1e-9 else "VIOLATED")
