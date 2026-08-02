"""Regression tests for the rank-1 dissipative RDM EOM.

These pin the result that motivates the whole project (docs/01_derivation_1rdm.md):
keeping the Mukherjee-Kutzelnigg hole contractions makes 0 <= gamma <= 1 an
invariant region and sends the system to Fermi-Dirac; dropping them breaks both.

If these fail, something foundational has regressed -- fix before proceeding.
"""
import numpy as np
import pytest

from qdiss.rdm1 import DissipativeRDM1, fermi_dirac, solve_mu


M = 10          # spin-orbitals
N = 4           # electrons
BETA = 1.4
SEED = 0


def _model(blocking=True):
    rng = np.random.default_rng(SEED)
    a = rng.normal(size=(M, M))
    h = 0.5 * (a + a.T)
    couplings = []
    for _ in range(2):
        c = rng.normal(size=(M, M))
        couplings.append(0.5 * (c + c.T))
    return DissipativeRDM1(h, couplings, beta=BETA, eta=0.05, wc=6.0,
                           blocking=blocking)


def _initial_gamma(model):
    """Far-from-equilibrium: the TOP N orbitals filled, in the h eigenbasis."""
    g = np.zeros((M, M), dtype=complex)
    for i in range(M - N, M):
        g[i, i] = 1.0
    return g


def _final(model, t_end=4000.0):
    g0 = _initial_gamma(model)
    ts, gs = model.propagate(g0, t_end=t_end, n_out=50)
    return gs[-1]


@pytest.fixture(scope="module")
def evolved_mk():
    m = _model(blocking=True)
    return m, _final(m)


def test_particle_number_conserved(evolved_mk):
    model, g = evolved_mk
    assert abs(np.trace(g).real - N) < 1e-8


def test_n_representability_invariant(evolved_mk):
    """0 <= gamma <= 1 must hold -- the whole point of the MK closure."""
    _, g = evolved_mk
    w = np.linalg.eigvalsh(0.5 * (g + g.conj().T))
    assert w.min() > -1e-8, f"negative occupation {w.min()}"
    assert w.max() < 1.0 + 1e-8, f"occupation above one: {w.max()}"


def test_converges_to_fermi_dirac(evolved_mk):
    """Occupations must match FD at the bath beta, with mu fixed by N alone."""
    model, g = evolved_mk
    n = np.diag(g).real
    mu = solve_mu(model.eps, BETA, N)
    n_fd = fermi_dirac(model.eps, BETA, mu)
    assert np.max(np.abs(np.sort(n) - np.sort(n_fd))) < 1e-6


def test_free_energy_is_monotone(evolved_mk):
    """The Lyapunov function must decrease along the flow."""
    model, _ = evolved_mk
    ts, gs = model.propagate(_initial_gamma(model), t_end=2000.0, n_out=60)
    F = np.array([model.free_energy(g) for g in gs])
    dF = np.diff(F)
    assert dF.max() < 1e-7, f"free energy increased by {dF.max()}"


def test_dropping_hole_contractions_breaks_representability():
    """The control: without the (1-gamma) blocking factors, occupations blow past 1.

    This is the negative result from docs/01 table -- it is what makes the
    positive result meaningful, so it is pinned too.
    """
    model = _model(blocking=False)
    g = _final(model)
    w = np.linalg.eigvalsh(0.5 * (g + g.conj().T))
    assert w.max() > 1.0 + 1e-3, (
        "expected N-representability violation without hole contractions, "
        f"but max occupation was {w.max()}"
    )
