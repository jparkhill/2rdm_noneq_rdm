# 07 — Mukherjee–Kutzelnigg generalized normal ordering and the generalized Wick
# theorem: primary sources, exact algebraic statements, and an audit of `MKNO.nb`

Compiled 2026-08-02. Mirrors live in `refs/papers/mkno/`; manifest entries appended to
`refs/papers/lit/MANIFEST.md`. Every numerical claim tagged **[NUM]** is reproduced by
`tools/verify_gwt.py` (run it; it prints pass/fail per identity).

## 0. Evidence tags used throughout

| tag | meaning |
|---|---|
| **[P-FULL]** | read in the primary paper's full text |
| **[P-ABS]** | read only the primary paper's *publisher abstract* (full text paywalled) |
| **[S-FULL]** | read in a secondary but authoritative full text that restates the primary result and cites it |
| **[NUM]** | independently verified numerically here, in explicit Fock space |
| **[INF]** | my inference — not stated in any source I read |

## 1. Access ledger — be clear about what could NOT be obtained

**All four requested primaries are closed access.** Checked against Unpaywall, OpenAlex,
Semantic Scholar, arXiv and Crossref on 2026-08-02; no open-access location exists for any
of them, and no legal copy was fetched. **I read none of these in full text.**

| # | citation (Crossref-verified) | DOI | status |
|---|---|---|---|
| 1 | W. Kutzelnigg, D. Mukherjee, *Normal order and extended Wick theorem for a multiconfiguration reference wave function*, **J. Chem. Phys. 107, 432–449 (1997)** | 10.1063/1.474405 | paywalled; **abstract obtained** |
| 2 | D. Mukherjee, *Normal ordering and a Wick-like reduction theorem for fermions with respect to a multi-determinantal reference state*, **Chem. Phys. Lett. 274, 561–566 (1997)** | 10.1016/S0009-2614(97)00714-8 | paywalled; **no abstract obtainable** |
| 3 | W. Kutzelnigg, D. Mukherjee, *Cumulant expansion of the reduced density matrices*, **J. Chem. Phys. 110, 2800–2809 (1999)** | 10.1063/1.478189 | paywalled; **abstract obtained** |
| 4 | W. Kutzelnigg, D. Mukherjee, *Irreducible Brillouin conditions and contracted Schrödinger equations for n-electron systems. **II. Spin-free formulation***, **J. Chem. Phys. 116, 4787–4801 (2002)** | 10.1063/1.1448827 | paywalled; **abstract obtained** |

> **Bibliographic correction to the brief.** JCP **116**, 4787 (2002) is *Part II* and is
> spin-free-specific. *Part I* — the one that actually derives the IBC_k / ICSE_k in terms
> of cumulants — is **D. Mukherjee, W. Kutzelnigg, J. Chem. Phys. 114, 2047–2061 (2001)**,
> doi:10.1063/1.1337058. There is an earlier companion: W. Kutzelnigg, D. Mukherjee,
> *Chem. Phys. Lett.* **317**, 567 (2000). If we cite "the IBC/ICSE paper" we should cite I,
> not II. **[P-ABS]**

Also paywalled, abstract only: Kong/Nooijen/Mukherjee JCP **132**, 234107 (2010)
(10.1063/1.3439395); Shamasundar JCP **131**, 174109 (2009) (10.1063/1.3256237);
Kutzelnigg/Shamasundar/Mukherjee *Mol. Phys.* **108**, 433 (2010) (10.1080/00268970903547926);
Hanauer & Köhn *Chem. Phys.* **401**, 50 (2012) (10.1016/j.chemphys.2011.09.024);
Neuscamman/Yanai/Chan JCP **130**, 124102 (2009) (10.1063/1.3086932).

**Could not obtain at all:** Misiewicz, Turney, Schaefer, *Reduced Density Matrix Cumulants:
The Combinatorics of Size-Consistency and Generalized Normal Ordering*, JCTC **16**, 6150
(2020). The journal version is paywalled and the ChemRxiv preprint
(10.26434/chemrxiv.12197886) sits behind a Cloudflare challenge that refuses non-browser
clients. Nothing below depends on it. Also not obtained: Evangelista, *Annu. Rev. Phys.
Chem.* **70**, 349 (2019) (Annual Reviews blocks scripted download).

**What I did obtain in full text** (all mirrored, all open access):

| key | citation | why it matters |
|---|---|---|
| `2205.01178v1` | F. A. Evangelista, *Automatic derivation of many-body theories based on general Fermi vacua* (Wick&d), JCP **157**, 062502 (2022) | complete modern restatement of MK GNO + GWT, explicitly citing KM-1997 |
| `margocsy_thesis_elte2022` | Á. Margócsy, PhD thesis, ELTE (2022), §2.1.3–2.1.6 | **the most complete open derivation**: GNO, GWT, worked contraction algebra, particle–hole conjugation of cumulants |
| `hanauer_thesis_mainz2013` | M. Hanauer, PhD thesis, Mainz (2013), §3.2.2 and ch. 6 | independent statement of the same rules; the extended version of Hanauer & Köhn 2012; numerical consequences of cumulant truncation |
| `1602.05667v1` | Li & Evangelista, *MR-LDSRG(2)*, JCP **144**, 164114 (2016) | third independent statement of the GWT rules |
| `1512.06956`, `hergert_frontiers2020_tour`, `1408.6555` | Hergert et al., IMSRG / MR-IMSRG | nuclear-physics restatement — independent convention cross-check |
| `1911.06836` | Lischka et al. / Evangelista et al., *Chem. Rev.* **120**, 5878 (2020) | context |

---

## 2. (a) The definition of generalized normal ordering

**Defining property.** For an arbitrary (correlated, non-idempotent, possibly
multideterminantal) reference $\Psi_0$, MK define the normal-ordered product
$\{\hat a^{pq\cdots}_{rs\cdots}\}$ by the requirement that its expectation value in that
same reference vanishes:

$$\boxed{\;\langle \Psi_0 | \{\hat a^{pq\cdots}_{rs\cdots}\} | \Psi_0\rangle = 0\;}$$

with the operator-string convention
$\hat a^{pq\cdots}_{rs\cdots} = \hat a^\dagger_p \hat a^\dagger_q \cdots \hat a_s \hat a_r$
(upper = creation, lower = annihilation, annihilators read right-to-left).
**[S-FULL]** Evangelista 2022, Eq. (2) and Eq. (1), citing KM-1997.
Same statement **[S-FULL]** Margócsy §2.1.5 Eq. (2.1.37); Hanauer §3.2.2.

The KM-1997 abstract states the three guiding principles verbatim: *"normal order operators
with respect to any reference state must be expressible as linear combinations of those with
respect to the genuine vacuum, that the vacuum expectation value of a normal order operator
must vanish (with respect to the vacuum to which it is in normal order), and that the
well-known formalism for a single Slater determinant as physical vacuum must be contained as
a special case."* **[P-ABS]**

**The definition is recursive, and this is where the engine's assumption enters.** Applying
the vanishing condition from rank 1 upward gives, in MK convention:

$$\hat a^{p}_{q} = \{\hat a^{p}_{q}\} + \gamma^{p}_{q}$$
$$\hat a^{pq}_{rs} = \{\hat a^{pq}_{rs}\}
+ \gamma^{p}_{r}\{\hat a^{q}_{s}\} - \gamma^{p}_{s}\{\hat a^{q}_{r}\}
- \gamma^{q}_{r}\{\hat a^{p}_{s}\} + \gamma^{q}_{s}\{\hat a^{p}_{r}\}
+ \gamma^{pq}_{rs}$$

**[S-FULL]** Margócsy Eqs. (2.1.40)–(2.1.41); condensed form
$\hat a_K=\sum_{K'=0}^{K}\gamma_{K'}\wedge\{\hat a_{K-K'}\}$, Eq. (2.1.42). **[NUM]** verified.

> **Typo in the source, flagged.** The extracted text of Margócsy Eq. (2.1.41) reads
> `... + γ^p_s :a^p_q:` for the fourth γ term. That is a typo (repeated index `p`, and it
> duplicates the second term). The correct term is $+\gamma^{r}_{s}\{\hat a^{p}_{q}\}$;
> that is what makes $\langle\{\hat a^{pq}_{rs}\}\rangle=0$ hold. **[NUM]** confirmed
> (max $|\langle\{\hat a^{pq}_{rs}\}\rangle| = 1.6\times10^{-16}$ over all indices).

**Crucially: the fully-contracted term in the rank-2 relation is the *full* 2-RDM
$\gamma^{pq}_{rs}$, not $\gamma\wedge\gamma$.** An engine whose only contractions are
$\gamma$ and $\eta$ produces $\gamma^p_r\gamma^q_s-\gamma^p_s\gamma^q_r$ there instead. This
single substitution is the whole audit, and it is quantified in §6.

**Elementary (pairwise) contractions.** These *are* what the brief describes:
$$\contraction{}{\hat a^\dagger_p}{}{\hat a_q}\hat a^\dagger_p \hat a_q = \gamma^p_q
\equiv \langle\Psi_0|\hat a^\dagger_p\hat a_q|\Psi_0\rangle,\qquad
\hat a_q \hat a^\dagger_p = \eta^p_q = \delta^p_q - \gamma^p_q .$$
**[S-FULL]** Evangelista 2022 Eqs. (4)–(5); Li & Evangelista 2016 Eqs. (20)–(21);
Margócsy Eq. (2.1.38); Hanauer Eqs. (3.50)–(3.51). Four independent statements, identical.
**[NUM]** verified ($\eta = 1-\gamma^{\mathsf T}$ in the index convention above).

---

## 3. (b) The generalized Wick theorem, in full — **the load-bearing answer**

### 3.1 The answer

**GWT is NOT limited to pairwise $\gamma/\eta$ contractions.** It generates multi-leg
contractions, and a $2k$-leg contraction ($k\ge 2$) evaluates to an element of the $k$-body
density cumulant $\lambda_k$. This is stated, in identical form, by every full-text source I
read.

$$\hat q_1\hat q_2\cdots = \{\hat q_1\hat q_2\cdots\}
+ \sum_{\text{single pairs}} + \sum_{\text{double pairs}}
+ \sum_{\text{single 4-leg}} + \sum_{\text{triple pairs}}
+ \sum_{\text{single pairs}}\sum_{\text{single 4-leg}} + \cdots$$

**[S-FULL]** Evangelista 2022, Eq. (3). Note the term structure: pair contractions and
4-leg contractions appear *at the same order* and are summed together.

The 4-leg contraction is:
$$\contraction{}{\hat a^\dagger_p}{\hat a^\dagger_q \hat a_s}{\hat a_r}
\hat a^\dagger_p\hat a^\dagger_q\hat a_s\hat a_r \;=\; \lambda^{pq}_{rs}
\;=\; \gamma^{pq}_{rs} - \gamma^{p}_{r}\gamma^{q}_{s} + \gamma^{p}_{s}\gamma^{q}_{r}$$
and multi-leg contractions are antisymmetric under permutation of the operators involved,
$\hat a^\dagger_p\hat a^\dagger_q\hat a_s\hat a_r = -\,\hat a^\dagger_p\hat a_s\hat a^\dagger_q\hat a_r$.
**[S-FULL]** Evangelista 2022, Eqs. (6)–(7).

Hanauer states the same as a general rule: *"a joint contraction of $2n$ operators with
$n\ge2$ always yields an $n$-body cumulant"*, $\lambda^{rs}_{pq}$, $\lambda^{rst}_{opq}$, …
**[S-FULL]** Hanauer §3.2.2, Eqs. (3.52)–(3.53).

Margócsy states it as $\overbrace{\alpha_{p_1}\cdots\alpha_{p_{2K}}} = \langle
\alpha_{p_1}\cdots\alpha_{p_{2K}}\rangle_c$ — "the MRGNO expansion is obtained by considering
all possible $K$-particle contractions" — **[S-FULL]** Eq. (2.1.43).

Li & Evangelista: *"in addition to pairwise contractions, new multi-legged contractions
appear … these new contractions correspond to elements of the $k$-body density cumulant
$\lambda_k$."* **[S-FULL]** arXiv:1602.05667, text after Eq. (21).

Hergert (nuclear): the MR-IMSRG uses *"a generalized normal ordering that is extended with
contractions of higher rank, namely the irreducible $k$-body density matrices
$\lambda^{(k)}$"*, with $\lambda^{pq}_{rs}\equiv\rho^{pq}_{rs}-\rho^p_r\rho^q_s+\rho^q_r\rho^p_s$.
**[S-FULL]** Frontiers 2020 Eqs. (25)–(26) — same sign convention.

### 3.2 The explicit rank-1 × rank-1 product rule (the one to test `MKNO.nb` against)

$$\boxed{\;\{\hat a^{p}_{q}\}\{\hat a^{r}_{s}\}
= \{\hat a^{pr}_{qs}\}
+ \eta^{r}_{q}\{\hat a^{p}_{s}\}
- \gamma^{p}_{s}\{\hat a^{r}_{q}\}
+ \gamma^{p}_{s}\eta^{r}_{q}
+ \lambda^{pr}_{qs}\;}$$

**[S-FULL]** Margócsy Eq. (2.1.52). **[NUM]** verified as a full operator identity in Fock
space (max element-wise error $4.4\times10^{-16}$); with the $\lambda^{pr}_{qs}$ term
deleted the identity fails by $1.4\times10^{-1}$ on the same random reference.

The rank-2 × rank-1 product $\{\hat a^{pr}_{qs}\}\{\hat a^{t}_{u}\}$ is written out in full
in Margócsy Eq. (2.1.53) — 19 terms, of which **8 carry a $\lambda_2$ factor and 1 carries
$\lambda_3$**. That is the single most useful line in the open literature for regression-testing
a hand-built engine, and it is in `refs/papers/mkno/margocsy_thesis_elte2022.pdf`.

### 3.3 The general product rule and its one exception

> *"Product of two MRGNO operators $\{S_1\}$ and $\{S_2\}$ is evaluated by considering all
> possible contractions between them, each contraction producing a $K$-body cumulant (with
> $\gamma_1=\lambda_1$ understood). **The one exception under this rule is the
> pair-contraction between an annihilation operator of $\{S_1\}$ and a creation operator of
> $\{S_2\}$, which yields $\eta_1$.**"* **[S-FULL]** Margócsy §2.1.6, verbatim.

Sign rule: $(-1)^{\pi_U+\pi_L-n_\eta}$, where $\pi_U,\pi_L$ are the parities of the upper and
lower index permutations counted separately and $n_\eta$ is the number of $\eta_1$ factors in
the term (the extra $(-1)^{n_\eta}$ corrects for $\eta^p_q=\langle a_q a_p\rangle$ having its
indices in the "wrong" order). **[S-FULL]** Margócsy §2.1.6. This is a concrete sign
convention worth diffing against `MKNO.nb`.

**Self-contractions vanish**: only contractions *connecting* distinct normal-ordered groups
survive. **[S-FULL]** Evangelista Eq. (8); Margócsy §2.1.6; Hanauer §3.2.2.

**In commutators**, single pairwise contractions give a Kronecker delta while *single
multi-leg contractions give null contributions* — a useful simplification, and a distinct
statement from "multi-leg contractions never appear". **[S-FULL]** Li & Evangelista 2016,
text after Eq. (23).

### 3.4 Where multi-leg contractions can land

For a **CAS** reference, cumulants vanish if any index is in a fully-occupied or empty
orbital, so multi-leg contractions connect only *active* indices, and $\gamma_1$ is block
diagonal (identity on core, zero on virtual). **[S-FULL]** Evangelista 2022, Table I and
the paragraph before it; Li & Evangelista 2016; Margócsy §2.1.4.

**This escape hatch does not apply to `MKNO.nb`.** The engine normal-orders with respect to a
correlated $\gamma$ over the *whole* orbital space with no core/active/virtual partition — in
Evangelista's Table I taxonomy every orbital is of type *General*. So every cumulant block is
nontrivial and the missing multi-leg contractions are not confined to a small active space;
they are everywhere. **[INF]** (the taxonomy is [S-FULL]; the application to `MKNO.nb` is mine).

### 3.5 So: is a $\gamma/\eta$-only engine "GNO with $\lambda_2=\lambda_3=0$"?

**Yes — and it has a name in the primary literature.** Kong, Nooijen and Mukherjee's
algebraic proof of the GWT explicitly builds this object:

> *"it follows from the proof that the normal order theory and corresponding contraction
> rules hold equally well if the contractions are only defined up to a certain level. These
> relaxations enable us to extend the original normal order theory. **As a particular
> example, a quasi-normal-order theory is developed, in which only one-body contractions are
> present. These contractions are based on the one-particle reduced density matrix.**"*
> **[P-ABS]** Kong, Nooijen, Mukherjee, JCP **132**, 234107 (2010), doi:10.1063/1.3439395.

That is precisely `MKNO.nb`'s algebra. It is a real, published, self-consistent operator
algebra — **but it is not MK GNO**, and it pays a specific price:

> *"There are alternative or further generalized extensions of normal ordering for MR wave
> functions, such as the **quasi-normal ordering of Kong et al.**, or the generalized
> antisymmetric ordered product of Sinha et al. When applied to a single determinant, the
> traditional normal ordering is recovered by these different orderings, but **they are not
> equivalent in general. For example, the quasi-normal ordering does not satisfy Eq. (2.1.37)
> [i.e. $\langle\{S\}\rangle=0$] for an MR wave function.**"*
> **[S-FULL]** Margócsy §2.1.5, verbatim.

**[NUM]** I quantified exactly how it fails. Define the rank-2 "normal-ordered" operator the
way a $\gamma/\eta$-only engine must (fully-contracted piece $=\gamma\wedge\gamma$ rather than
$\gamma_2$). Then

$$\langle\Psi_0|\{\hat a^{pq}_{rs}\}_{\text{quasi}}|\Psi_0\rangle = \lambda^{pq}_{rs}
\quad\text{exactly}$$

(max deviation $1.6\times10^{-16}$; max residual magnitude $0.136$ on the test reference).
**The defining property of normal ordering fails by exactly $\lambda_2$.**

### 3.6 The sharp answer to the brief's question

The brief asks whether a $\gamma/\eta$-only engine is "setting $\lambda_2=\lambda_3=0$ at the
Wick level, a stronger and different approximation than truncating the cumulant expansion of
the 3-RDM afterwards." **It is, and here is the precise statement:**

1. Setting $\lambda_3=0$ *only in a 3-RDM reconstruction* (Valdemoro/Nakatsuji/Mazziotti
   style) leaves $\gamma_3 \approx \lambda_2\wedge\gamma_1 + \tfrac16\gamma_1^{\wedge3}$ —
   $\lambda_2$ survives, the reference stays correlated. **[S-FULL]** Hanauer Eq. (6.4)/(6.8).
2. Killing $\lambda_2$ *at the Wick level* is strictly stronger: it removes $\lambda_2$ from
   every contraction pattern at every rank, and it changes the *operator algebra*, not just a
   reconstruction formula.
3. And there is an exactness boundary that is not about idempotency:

> **$\lambda_k \equiv 0$ for all $k$ $\iff$ the reference is quasi-free (fermionic Gaussian).**
> A Slater determinant is one such state; so is a *mean-field thermal / grand-canonical*
> state, whose $\gamma$ is emphatically **not** idempotent. **[NUM]**: for
> $\rho \propto e^{-\sum_p \varepsilon_p \hat n_p}$ I get non-idempotent $\gamma$ with
> $\|\lambda_2\| = 3.7\times10^{-16}$, $\|\lambda_3\| = 7.8\times10^{-16}$.

And from the primary abstract, for **pure** states:

> *"Necessary and sufficient for $\lambda_2=0$ is the well-known idempotency condition
> $\gamma^2=\gamma$ for $\gamma\equiv\gamma_1$."* **[P-ABS]** KM, JCP **110**, 2800 (1999).

Reconciling those two (that KM's statement is scoped to pure $n$-fermion states, while the
Gaussian counterexample is a mixed state) is **[INF]** — I have only the abstract, not the
paper's hypotheses.

**Consequence for `MKNO.nb`.** Its contraction algebra is the *exact* Wick algebra of a
quasi-free state, and of nothing else. Non-idempotency of $\gamma$ is *not* the problem;
correlation is. If the reference is intended to be a correlated **pure** state, the engine is
internally inconsistent: its $\gamma$ says "correlated" and its Wick rules say "quasi-free".
If the reference is a mean-field **mixed** state (the rank-1 Fermi–Dirac lineage, and the
target of a Stinespring/ancilla purification), the engine is exact and there is nothing to
fix at the Wick level. **[INF]**, built on [P-ABS] + [NUM].

---

## 4. (c) Cumulant decomposition of the 2- and 3-RDM, with conventions

### 4.1 MK convention as used by Evangelista / Hanauer / Margócsy / Hergert

Normalisation: $\gamma^{pq\cdots}_{rs\cdots} = \langle\Psi_0|\hat a^\dagger_p\hat a^\dagger_q
\cdots \hat a_s\hat a_r|\Psi_0\rangle$ (**no** $1/k!$, **no** Löwdin normalisation), so
$\mathrm{Tr}\,\gamma_1 = N$ and $\sum_{pq}\gamma^{pq}_{pq} = N(N-1)$.

$$\gamma^{pq}_{rs} = \gamma^{p}_{r}\gamma^{q}_{s} - \gamma^{p}_{s}\gamma^{q}_{r} + \lambda^{pq}_{rs}$$

**[S-FULL]** Evangelista Eq. (6) (as the inverse, $\lambda^{pq}_{rs}=\gamma^{pq}_{rs}-\gamma^p_r\gamma^q_s+\gamma^p_s\gamma^q_r$);
Hanauer Eq. (2.25); Margócsy Eq. (2.1.28); Hergert Eq. (26). **All four agree exactly,
including signs.** **[NUM]** verified.

In wedge (Grassmann) notation, with $\wedge$ carrying the "modified prefactor" convention
that excludes redundant permutations within each antisymmetric tensor:

$$\gamma_1=\lambda_1,\qquad
\gamma_2=\lambda_2+\tfrac12\lambda_1\wedge\lambda_1,\qquad
\gamma_3=\lambda_3+\lambda_2\wedge\lambda_1+\tfrac16\lambda_1^{\wedge3}$$
$$\gamma_4=\lambda_4+\lambda_3\wedge\lambda_1+\tfrac12\lambda_2\wedge\lambda_2
+\tfrac12\lambda_2\wedge\lambda_1^{\wedge2}+\tfrac1{24}\lambda_1^{\wedge4}$$

**[S-FULL]** Margócsy Eqs. (2.1.24)–(2.1.26); Hanauer Eqs. (6.2)–(6.5) (identical, with
$\otimes_A$ for $\wedge$). General formula: Margócsy Eq. (2.1.27) / Hanauer Eq. (6.1),
$\gamma_K=\sum_{\nu}\sum_{(K_1\ge\cdots\ge K_\nu)\vdash K} w^{-1}_{K_1\cdots K_\nu}
\bigwedge_\mu \lambda_{K_\mu}$ with $w = n_1!\cdots n_m!$ counting repeats in the partition.

Inverted:
$$\lambda_2=\gamma_2-\tfrac12\gamma_1\wedge\gamma_1,\qquad
\lambda_3=\gamma_3-\gamma_2\wedge\gamma_1+\tfrac13\gamma_1^{\wedge3}$$
$$\lambda_4=\gamma_4-\gamma_3\wedge\gamma_1-\tfrac12\gamma_2\wedge\gamma_2
+\gamma_2\wedge\gamma_1^{\wedge2}-\tfrac14\gamma_1^{\wedge4}$$
**[S-FULL]** Margócsy Eq. (2.1.29); Hanauer Eqs. (6.7)–(6.9). **[NUM]** verified for
$\lambda_2,\lambda_3$ by explicit index-level construction (nine $\lambda_2\!\otimes\!\gamma_1$
placements with sign $(-1)^{i+j}$, plus the $3\times3$ determinant of $\gamma_1$).

### 4.2 Convention warnings — this is where a hand-built engine drifts

- **Normalisation changes the coefficients.** *"as Eq. (2.1.27) is a non-linear relation, the
  coefficients of the expansion depend on the normalization convention used for $\gamma_K$.
  See for example Ref. [17] for the same formulae with Löwdin normalization."*
  **[S-FULL]** Margócsy, note after Eq. (2.1.28). If `MKNO.nb` uses $\tfrac14\sum v\,a^{pq}_{rs}$
  style prefactors, check whether its $\gamma_2$ carries a $1/2$.
- **Wedge prefactor.** Hanauer's $\otimes_A$ is *"a Grassmann (or wedge) product with a
  modified prefactor"* that implicitly excludes redundant permutations; his footnote 13 spells
  out $(\gamma_1\otimes_A\tau_2)^{rst}_{opq}$ term by term (9 terms). **[S-FULL]** Using the
  textbook $\wedge$ with $1/(p!q!)$ instead will introduce exactly the factor-of-2 errors the
  brief worries about.
- **Spin-free vs spin-orbital.** Spin-free cumulants require a *separate definition* —
  *"the definition of the cumulants of the reduced densities has to be generalized … making
  use of irreducible tensor operators with respect to SU(2). **There are two alternative
  definitions of spin-free cumulants**, of which the one in terms of spin-free reduced density
  matrices appears to be preferable."* **[P-ABS]** KM, JCP **116**, 4787 (2002). The
  systematic treatment is Kutzelnigg/Shamasundar/Mukherjee, *Mol. Phys.* **108**, 433 (2010)
  and Shamasundar, JCP **131**, 174109 (2009) — both paywalled, both unread here. **If
  `MKNO.nb` is spin-summed anywhere, this is an unaudited gap.**

### 4.3 A free consistency test for the engine's cumulants

$$\sum_r \lambda^{pr}_{qr} = \delta^{p}_{q}\,(n_q^2 - n_q)\quad\text{(natural-orbital basis)}$$
$$\sum_t \lambda^{prt}_{qst} = \lambda^{pr}_{qs}\,(n_p+n_r+n_q+n_s-2)$$

**[S-FULL]** Margócsy Eqs. (2.1.31)ff, who notes explicitly *"these relations can serve as
tests to check the numerical implementation of cumulants."* **[NUM]** the $\lambda_2$ one is
verified here. Note the RHS $n^2-n \le 0$: **partial traces of $\lambda_2$ are negative
semidefinite**, and their magnitude is exactly the occupation-number variance.

---

## 5. (d) N-representability, positivity, particle–hole symmetry under GNO truncation

### 5.1 What the literature actually says

Honest summary: **I found no paper that states "λ₂ is *the* obstruction to positivity of a
GNO-truncated 2-RDM" as a theorem.** What exists is (i) the exact-conditions literature,
(ii) empirical failure reports, and (iii) the p–h equivalence result.

- **Cumulant truncation demonstrably destabilises methods.** *"Numerical experiments show
  that the rigorous evaluation of contributions containing four-body density cumulants is
  critical to avoid **artificial variational collapse** of the EOM energies, as has been
  reported before in cumulant-truncated ic-MRCI methods."* **[S-FULL]**
  arXiv:2506.13693 (Evangelista group, IP-EOM-DSRG), citing Saitow, Kurashige, Yanai
  (DMRG-ic-MRCI). This is the closest published analogue to your failure mode.
- Hanauer's own numerics: truncating at $n_{\max}=2$ (i.e. keeping only $\lambda_2$ and
  discarding $\lambda_{\ge3}$) degrades ic-MRCCSD energies by up to 3.4 mE_h on N₂ *and*
  causes the solver to fail to converge — *"this problem affected especially calculations with
  a low $n_{\max}$."* **[S-FULL]** Hanauer §6.4, Figs. 36–37. Note this is the *mild* case;
  `MKNO.nb` is effectively $n_{\max}=1$.
- **Cumulants have their own N-representability problem, and it is open.** Margócsy, on
  dressed cumulants: *"Finding reasonable (e.g. N-representable) approximations for
  $\bar\lambda_K$ could…"* **[S-FULL]** — i.e. treated in 2022 as an unsolved design question,
  not a solved one.
- **Particle–hole equivalence in N-representability** is the classical result:
  M. B. Ruskai, *N-Representability Problem: Particle-Hole Equivalence*, J. Math. Phys.
  **11**, 3218 (1970) — cited by Margócsy as ref. [16]. **[S-FULL]** (citation only; I did not
  read Ruskai).
- KM's own N-representability work is the IBC/ICSE hierarchy: *"a hierarchy of $k$-particle
  approximations for the direct calculation of the cumulants $\lambda_k$"*, and Part I notes
  that expressing BC$_k$/CSE$_k$ in cumulants *"is not sufficient to make [them] separable
  (extensive) … however, in a formulation based on the … general normal ordering with respect
  to arbitrary wave functions, the irreducible counterparts IBC$_k$ and ICSE$_k$ can be
  defined."* **[P-ABS]** JCP 114, 2047 (2001) / JCP 116, 4787 (2002). **These abstracts claim
  extensivity, not N-representability or positivity.** Do not conflate the two.

### 5.2 What I can state and have verified: λ₂ *is* the positivity obstruction

Not found verbatim in the literature, so **[INF] + [NUM]**, with a proof:

**Claim.** Let $0 \preceq \gamma \preceq I$. Then the $\lambda_2=0$ surrogate
$D = \tfrac12\gamma\wedge\gamma$ (i.e. $D^{pq}_{rs}=\gamma^p_r\gamma^q_s-\gamma^p_s\gamma^q_r$)
is **automatically positive semidefinite** on the antisymmetric two-particle space, and so is
$Q=\tfrac12\eta\wedge\eta$.

*Proof.* $\gamma\otimes\gamma \succeq 0$ on the full two-particle space. Let $\mathcal A$ be
the orthogonal projector onto the antisymmetric subspace. Then
$\mathcal A(\gamma\otimes\gamma)\mathcal A \succeq 0$, and on that subspace it equals
$\tfrac12\gamma\wedge\gamma$ up to the positive normalisation. Same for $\eta$. ∎

**[NUM]** 300 random $0\preceq\gamma\preceq I$: zero violations. On the correlated test
reference, $\min\mathrm{eig}(D_{\lambda_2=0}) = +8.3\times10^{-3}$ and
$\min\mathrm{eig}(Q_{\lambda_2=0}) = +8.3\times10^{-3}$.

**Conversely**, adding an arbitrary antisymmetric, Hermitian $\lambda_2$ to
$\gamma\wedge\gamma$ readily produces indefinite $D$: worst $\min$ eigenvalue $-1.85$ over 400
random draws. **[NUM]**

**Therefore: in the decomposition $D = \tfrac12\gamma\wedge\gamma + \lambda_2$, with
$\gamma$ constrained to $0\preceq\gamma\preceq I$, $\lambda_2$ is the *only* place where
$D\succeq0$ can fail.** Same statement for $Q$ with $\eta$ (see §6).

**And a stronger, cleaner statement.** $\lambda_k\equiv0\ \forall k$ means the whole RDM
hierarchy is that of a fermionic Gaussian state. A Gaussian state is a genuine density
operator, so its RDMs are **ensemble-N-representable by construction** — every positivity
condition ($D$, $Q$, $G$, $T_1$, $T_2$, all of them) holds automatically. **[INF]**

**The uncomfortable corollary for `MKNO.nb`.** Its Wick algebra can never *produce* a
positivity violation from the reference side — it is structurally incapable of representing a
correlated reference at all. So if the propagated 2-RDM goes non-positive, the mechanism is
**not** "the reference wasn't N-representable"; it is that the EOM is missing the
$\lambda_2$-mediated terms that couple the correlated part of $^2D$ back into its own
evolution. Those terms are exactly the ones that keep $\lambda_2$ inside the admissible set.
Derived here, **[INF]** — but it is the direct reading of §3.5 + §5.2 and it is testable:
re-derive one EOM term with and without the $\lambda_2$ contractions and watch the sign of
$\min\mathrm{eig}(^2D)$.

---

## 6. (e) Particle–hole conjugation of GNO — **there is a clean published statement**

Define the hole RDMs $\eta_K$ by the reversed string
$\eta^{p_1\cdots p_K}_{q_1\cdots q_K}=\langle a_{q_1}\cdots a_{q_K}\, a^\dagger_{p_K}\cdots a^\dagger_{p_1}\rangle$,
so $\eta^p_q=\delta^p_q-\gamma^p_q$. Then:

$$\eta_2 = \lambda_2 + \tfrac12\,\eta_1\wedge\eta_1$$
$$\eta_3 = -\lambda_3 + \lambda_2\wedge\eta_1 + \tfrac16\,\eta_1^{\wedge3}$$
$$\eta_4 = \lambda_4 - \lambda_3\wedge\eta_1 + \tfrac12\lambda_2\wedge\lambda_2
+ \tfrac12\lambda_2\wedge\eta_1^{\wedge2} + \tfrac1{24}\eta_1^{\wedge4}$$

> *"In the general case, cumulant expansion of $\eta_K$ can be formally obtained from that of
> $\gamma_K$ by replacing $\lambda_1=\gamma_1$ with $\eta_1$, and $\lambda_{2K'-1}$ with
> $-\lambda_{2K'-1}$ for $K'\ge2$."* **[S-FULL]** Margócsy §2.1.3, verbatim, with the three
> equations above.

**[NUM]** Both verified exactly: hole 2-cumulant $=$ particle 2-cumulant; hole 3-cumulant
$=-$ particle 3-cumulant.

**So GNO has a clean but *not perfect* $\gamma\leftrightarrow\eta$ covariance, and it is
stated in the literature.** Precisely:

- **Even-rank cumulants are p–h invariant.** In particular $\lambda_2$ is *the same object*
  in the particle and hole channels.
- **Odd-rank cumulants ($k\ge3$) flip sign.**
- $\lambda_1$ maps to $\eta_1 = I-\gamma_1$ — the only "non-covariant" piece, and it is exactly
  the $\delta$ shift.
- Margócsy also proves the operator-level statement behind it: *"operators anticommute with
  each other within the cumulant expectation value $\langle\cdots\rangle_c$ for $K\ge2$"*,
  and *"any product of $a_p,a_r,a_q,a_s$ has the same cumulant part in its expansion"*
  (e.g. $\langle a_s a_q a_p a_r\rangle_c = \langle a_p a_r a_s a_q\rangle_c = \lambda^{pr}_{qs}$).
  **[S-FULL]** Eqs. (2.1.32)–(2.1.33). *This* is the algebraic root of the p–h covariance:
  the cumulant part of a contraction does not care about operator order beyond a sign.
- Margócsy also notes $\gamma_K$ and $\eta_K$ have symmetric roles and near-identical partial
  trace relations, $\sum_r \eta^{p_1\cdots r}_{q_1\cdots r} = (M-N-K+1)\eta^{p_1\cdots}_{q_1\cdots}$.
  **[S-FULL]** Eq. (2.1.23).

### 6.1 Bearing on your particle-hole covariance conjecture

Your conjecture is that p–h covariance of the term list yields simultaneous $D$ and $Q$
preservation. The literature above says something quite specific about that:

$$D = \tfrac12\gamma_1\wedge\gamma_1 + \lambda_2, \qquad
Q = \tfrac12\eta_1\wedge\eta_1 + \lambda_2$$

— **the same $\lambda_2$, with the same sign, in both.** Combined with §5.2
($\gamma\wedge\gamma\succeq0$ and $\eta\wedge\eta\succeq0$ automatically whenever
$0\preceq\gamma\preceq I$), this gives a sharp version of your conjecture:

> $D\succeq0$ and $Q\succeq0$ are two constraints on **one and the same object** $\lambda_2$,
> sitting on the two different positive backgrounds $\gamma\wedge\gamma$ and
> $\eta\wedge\eta$. A term list that is p–h covariant therefore cannot preserve one and break
> the other *by construction* — but covariance alone does not bound $\lambda_2$, so it does
> not by itself guarantee either. **[INF]**, from [S-FULL] Margócsy + [NUM].

This also matches the memory note that the rank-2 pair-hole factor is
$(I\wedge I) - 2D = 2Q$: that identity is exactly the $\eta_2$ relation above at $\lambda_2$
fixed. And it explains why the rank-1 case was easy — at rank 1 the "cumulant" is $\gamma_1$
itself and the $\gamma\leftrightarrow\eta$ map is a pure affine shift with no cumulant to
worry about.

---

## 7. Verdict on `MKNO.nb`, and what to do

1. **The engine is not implementing MK GNO.** It implements the Kong–Nooijen–Mukherjee
   *quasi-normal ordering* (one-body contractions only). That is a published, internally
   consistent algebra — but it is a different algebra, with a different defining property.
2. **What it costs, exactly.** $\langle\{\hat a^{pq}_{rs}\}\rangle = \lambda^{pq}_{rs}$
   instead of $0$. Every derived EOM term is missing the $\lambda_2$ (and $\lambda_3$)
   contraction patterns of the true GWT. In the rank-1×rank-1 product that is one missing
   term out of five; in the rank-2×rank-1 product it is 9 missing terms out of 19.
3. **Non-idempotency was never the issue.** The algebra is exact for any quasi-free state,
   idempotent or not — including the finite-temperature mean-field reference of the rank-1
   Fermi–Dirac work. It is wrong only for *correlated* references. This is likely the precise
   reason the rank-1 result was clean and rank 2 is not.
4. **This is a stronger approximation than 3-RDM cumulant truncation**, and the brief's
   suspicion is confirmed on primary-abstract + multiple-full-text-secondary evidence.
5. **Concrete regression tests to run against `MKNO.nb`** (all in
   `tools/verify_gwt.py`, all passing there):
   - $\langle\{\hat a^{pq}_{rs}\}\rangle \stackrel{?}{=} 0$ on a correlated reference;
   - the rank-1×rank-1 product rule of §3.2, term by term, including the $\lambda^{pr}_{qs}$;
   - Margócsy Eq. (2.1.53), $\{\hat a^{pr}_{qs}\}\{\hat a^{t}_{u}\}$, 19 terms — the strongest
     available derivation-to-code test in the open literature;
   - the cumulant partial-trace identities of §4.3;
   - the p–h conjugation identities of §6.
6. **Two coherent ways forward** (both consistent with everything above):
   - *(A) Restore the cumulant contractions.* Extend the engine to true GWT. This is the
     Evangelista/Hanauer route and it is what every working multireference code does.
   - *(B) Keep the quasi-free algebra and make the reference honestly quasi-free.* Purify
     with an ancilla so the reference is a Gaussian mixed state. Then the engine is *exact*,
     all $\lambda_k$ genuinely vanish, and $D,Q,G$ positivity is automatic — at the price that
     all correlation must live in the ancilla coupling rather than in the reference. This is
     the principled version of the Stinespring idea already on the plan.

   These are not interchangeable: (A) keeps a correlated reference and buys back
   the missing physics; (B) keeps the current algebra and moves the correlation somewhere the
   algebra can see it. **[INF]**

## 8. Open gaps I did not close

- No full text of any KM primary. Everything attributed to KM-1997/1999/2001/2002 above is
  either from the publisher abstract **[P-ABS]** or from a secondary source that cites it
  **[S-FULL]**. If a term-level convention dispute arises, we need the actual PDFs (ILL, an
  institutional subscription, or an author request to the Mukherjee group).
- Mukherjee, CPL **274**, 561 (1997) — nothing obtained, not even an abstract. Its relation to
  KM-1997 (the KM-1997 abstract says the formulation is *"different (but essentially
  equivalent)"* to *"that suggested previously by one of the present authors"* — i.e. this
  paper) is therefore unverified beyond that sentence.
- Shamasundar 2009 and Kutzelnigg/Shamasundar/Mukherjee 2010: the **spin-free** GWT proofs and
  spin-summation rules are entirely unaudited. If `MKNO.nb` is spin-summed, this is a live risk.
- Misiewicz/Turney/Schaefer 2020 (size-consistency combinatorics of GNO): not obtained.
- Hanauer & Köhn, *Chem. Phys.* **401**, 50 (2012): the journal article is paywalled, but
  Hanauer's Mainz thesis (mirrored, full text) is its extended version and is what I cited.
- I did not find any published theorem of the form "GNO truncation preserves / destroys
  N-representability". §5.2 is my own derivation, numerically verified, and is labelled as such.
