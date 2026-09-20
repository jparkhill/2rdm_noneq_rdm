# Prior art: what must be cited, and which relationships must be checked before any claim of novelty

Status: OPEN. Written 2026-08-13. Nothing here has been read/verified today; these are the obligations.
Each item: what we use, the likely source, and the specific question to settle by reading.

## 1. RDMs as Gram matrices; the D, Q, G conditions
- Coleman (1963) Rev. Mod. Phys. 35, 668; Garrod & Percus (1964) J. Math. Phys. 5, 1756 — D, Q, G as Gram
  matrices of `a a|Ψ>`, `a†a†|Ψ>`, `a†a|Ψ>`. Cite for §3 of DYNDELTA-SPEC and tools/channel_cones.py.
- Mazziotti (ed.) Adv. Chem. Phys. 134 (2007) for the modern N-representability picture.

## 2. The "XR theorem" (docs/XR-POSITIVITY.md, commit bd24aee): exact flow of Mat_k is i(S†X − XS)
Ingredients are standard: Gram-matrix tangency (elementary); the form i(S†X − XS) = i[K,X] − ½{Γ,X} is the
non-jump part of GKSL (Gorini–Kossakowski–Sudarshan 1976; Lindblad 1976). **To check:** whether the specific
refactoring — the (k+1)-annihilator term of the two-body commutator factoring through the k-annihilated span,
so the exact BBGKY flow has this form at every rank — appears in
  - Mazziotti, "Anti-Hermitian contracted Schrödinger equation" PRL 97, 143002 (2006) and follow-ups on
    positivity-preserving 2-RDM propagation (Mazziotti, Chem. Phys. Lett. / PRA ~2007–2012);
  - Lackner, Březinová, Sato, Ishikawa, Burgdörfer, PRA 91, 023412 (2015) and PRA 95, 033414 (2017)
    (TD-2RDM with purification; "contraction-consistent" reconstruction);
  - Joost, Schlünzen, Bonitz, PRB 101, 245101 (2020) and Joost et al. PRB 105, 165155 (2022) (G1–G2 scheme,
    purification, "positivity" discussions);
  - Akbari, Hashemi, Rubio, Nieminen, van Leeuwen, PRB 85, 235121 (2012) (TD-2RDM instability analysis) —
    this one is the closest to our "F0 diverges" result and MUST be cited for episode 9.
Verdict needed: is the factored form stated there? If yes, cite and drop the name "XR theorem".

## 3. Cumulant reconstruction / λ₃ = 0 flow
- Mazziotti PRA 57, 4219 (1998); Kutzelnigg & Mukherjee JCP 110, 2800 (1999) (MKNO cumulants; the separable
  γ₃ we use IS their λ₃=0 reconstruction).
- Nakatsuji & Yasuda PRL 76, 1039 (1996); Valdemoro's reconstructions (Adv. Quantum Chem.).
- Our closed channel form F0 (episode 5) is the SAME flow as "propagate γ₂ with γ₃ → Mazziotti/Kutzelnigg
  separable part", written as pair/ph matrix products. **To check:** whether anyone has written the λ₃=0 rank-2
  EOM in the pp-ladder/ph-ring sandwich form (the (1−n_p−n_v)V and (n_p−n_q)V structures are textbook in
  Green's-function language: Bethe–Goldstone / pp-RPA and ph-RPA / GW; cite Fetter–Walecka or Ring–Schuck,
  "The Nuclear Many-Body Problem", for those). The identification is likely novel as an RDM statement; the
  structures are not.
- Mueller-type reconstruction we compare against: arXiv:2405.05093 (with the Lackner correction, above).

## 4. Extensivity = connectedness (docs/CONSTRAINT-LEDGER.md)
- Standard cumulant/linked-cluster theorem (Kubo 1962, J. Phys. Soc. Jpn. 17, 1100; coupled-cluster
  literature, Bartlett & Musiał RMP 2007 for size-extensivity). Cite; not novel.

## 5. Lie-algebraic framing (route F, episodes 1–3)
- Universality theorem used for "no polynomial Lie algebra contains u(M) and an interaction": Bravyi & Kitaev,
  Ann. Phys. 298, 210 (2002); Jozsa & Miyake, Proc. R. Soc. A 464, 3089 (2008); Brod & Galvão PRA 84, 022310.
- Lie–Poisson / Hamiltonian closures of moment hierarchies: Morrison, Rev. Mod. Phys. 70, 467 (1998);
  Perin, Chandre, Morrison, Tassi, Ann. Phys. 348, 50 (2014); metriplectic: Morrison Physica D 18, 410 (1986);
  GENERIC: Grmela & Öttinger PRE 56, 6620 (1997).
- Marcinkiewicz theorem (no finite cumulant truncation beyond 2 is a distribution): Marcinkiewicz, Math. Z. 44,
  612 (1939); quantum analogue used loosely — **check** the precise fermionic statement (Rajagopal & Sudarshan
  PRA 10, 1852 (1974) is the usual reference).

## 6. Dynamical constraint enforcement (episodes 8–9)
- Car & Parrinello PRL 55, 2471 (1985) — the fictitious-dynamics idea we borrowed.
- Penalty/barrier and softplus smoothing: standard optimisation (Nocedal & Wright); cite generically.
- Purification-by-projection alternatives we did NOT use but compared conceptually: Mazziotti PRE 65, 026704
  (2002) (purification of 2-RDMs); Lackner et al. 2015 above; Joost et al. 2022 above.
  **To check:** whether any of them use a smooth, dynamical (non-projective) enforcement. Our claim would be
  "first smooth dynamical positivity enforcement with an intact linear-response spectrum"; verify nobody has.

## 7. Bath / dissipation (docs/BATH-DERIVATION.md)
- TCL2: Breuer & Petruccione, "The Theory of Open Quantum Systems" (2002); Davies, Commun. Math. Phys. 39, 91
  (1974) for the secular limit and Gibbs fixed point.
- Rank-1 precedent: Parkhill et al. JCP 142, 134113 (2015) — our own; the rank-2 contraction with γ₃ is the
  new step, and "EE2" (unpublished) the failed prior attempt.
- Related: Bonitz group's dissipative G1–G2; Schlünzen, Joost, Bonitz PRL 124, 076601 (2020).

## 8. Pole/linear-response methodology
- Biorthogonal projection for non-normal Jacobians: standard (Trefethen & Embree, Spectra and Pseudospectra).
- ERPA / extended RPA from correlated 2-RDMs (comparison methods for poles): Chatterjee & Pernal JCP 137,
  204109 (2012); Rowe RMP 40, 153 (1968).

## Procedure
Before any external write-up: read items 2, 3, 6 in full; record verdicts here with page references; adjust
naming ("XR theorem" → cite or rename) accordingly. Until then, docs in this repo may use the internal names but
must not claim novelty.
