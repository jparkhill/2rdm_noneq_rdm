#!/usr/bin/env python3
r"""Route F, episode 9: dipole spectra from the 4.8 fs propagations (nonlinear, self-consistent poles).

Parses the text logs of tools/long_propagate.py / tools/dyndelta_gpu.py (dipole <mu_x>(t) every 0.5 a.u., and the exact
dipole alongside), removes the mean, applies a Hann window, zero-pads, and plots |FT| against the FCI singlet excitation
energies of this geometry (vertical lines; the kick is x-polarised so only x-dipole-allowed singlets appear).
Resolution: 200 a.u. window -> Hann FWHM ~ 2*2pi/200 = 0.06 Ha; peaks closer than that merge.
"""
import sys, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "."); sys.path.insert(0, "tools")

RUNS = [("exact (FCI)", "results/long_muller.txt", 8, "k"),
        ("muller2405", "results/long_muller.txt", 7, "tab:orange"),
        ("F0  (lambda_3 = 0), diverges 2.7 fs", "results/long_f0.txt", 7, "tab:gray"),
        ("dyn-delta + barrier k_neg=1e3", "results/long_gpu_dyndelta_k100n1e3.txt", 7, "tab:blue"),
        ("dyn-delta + barrier k_neg=1e4", "results/long_gpu_dyndelta_k100n1e4.txt", 7, "tab:green"),
        ("dyn-delta + softplus barrier", "results/long_gpu_dyndelta_k100n1e4sp.txt", 7, "tab:red")]


def parse(path, col):
    t, mu = [], []
    for line in open(path):
        if re.match(r"\s+[0-9.]+ \|", line):
            f = line.replace("|", " ").split()
            t.append(float(f[0])); mu.append(float(f[col]))
    return np.array(t), np.array(mu)


def spectrum(t, mu, tmax=200.0):
    m = t <= tmax
    t, mu = t[m], mu[m] - mu[m].mean()
    w = np.hanning(len(t)); dt = t[1] - t[0]
    n = 16 * len(t)
    F = np.fft.rfft(mu * w, n=n) * dt
    om = 2 * np.pi * np.fft.rfftfreq(n, d=dt)
    return om, np.abs(F)


def main():
    # FCI singlet excitations for this geometry (H4 trapezoid as in commutator_structure.H4)
    import commutator_structure as cs
    from molecular_2rdm import build_H_sector, molecule, sector
    from qdiss.eom import fullci
    atom, N, spin = cs.H4
    mol, mf, M, h, V, _ = molecule(atom, "sto-3g", 0, spin); dets = sector(M, N)
    es, vec = np.linalg.eigh(build_H_sector(h, V, M, dets))
    mu1 = cs.dipole(mol, mf, M); Mu = cs.sector_operator(mu1, M, dets)
    tdm = np.abs(vec[:, 0].conj() @ Mu @ vec)          # |<0|mu_x|k>|
    exc = es - es[0]
    allowed = [(e, f) for e, f in zip(exc, tdm) if f > 1e-3 and e > 1e-6]
    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=False)
    ax = axes[0]
    for label, path, col, color in RUNS:
        try:
            t, mu = parse(path, col)
        except FileNotFoundError:
            continue
        tmax = 200.0 if "F0" not in label else 90.0
        om, S = spectrum(t, mu, tmax)
        ax.plot(om, S / S.max() if "exact" not in label else S / S.max(), color=color, lw=2.2 if "exact" in label else 1.2, label=label + ("" if tmax == 200 else "  [0-90 a.u. only]"))
    for e, f in allowed:
        ax.axvline(e, color="k", ls=":", lw=0.8, alpha=0.6)
    ax.set_xlim(0, 1.8); ax.set_xlabel("frequency (Ha)"); ax.set_ylabel("|FT of <mu_x>(t)|, normalised")
    ax.set_title("H4 trapezoid, M=8 N=4, x-kicked ground state, 4.8 fs: dipole spectrum vs FCI dipole-allowed singlets (dotted)")
    ax.legend(fontsize=8)
    ax = axes[1]
    for label, path, col, color in RUNS:
        try:
            t, mu = parse(path, col)
        except FileNotFoundError:
            continue
        ax.plot(t * 0.02419, mu, color=color, lw=2.2 if "exact" in label else 1.0, label=label)
    ax.set_xlabel("t (fs)"); ax.set_ylabel("<mu_x>(t) (a.u.)"); ax.set_ylim(2.6, 4.2); ax.legend(fontsize=8)
    ax.set_title("dipole moment along the propagation")
    plt.tight_layout(); plt.savefig("results/dipole_spectra.png", dpi=130)
    # peak table
    print("FCI x-allowed singlet excitations (Ha) and |<0|mu_x|k>|: " + "  ".join("%.4f(%.3f)" % (e, f) for e, f in allowed if e < 1.8))
    for label, path, col, color in RUNS:
        try:
            t, mu = parse(path, col)
        except FileNotFoundError:
            continue
        om, S = spectrum(t, mu, 200.0 if "F0" not in label else 90.0)
        pk = [om[i] for i in range(1, len(S) - 1) if S[i] > S[i - 1] and S[i] > S[i + 1] and S[i] > 0.05 * S.max() and om[i] < 1.8]
        print("%-40s peaks: %s" % (label, "  ".join("%.4f" % p for p in pk)))


if __name__ == "__main__":
    main()
