# Literature archive

PDFs + `pypdf`-extracted text. Fetched from arXiv 2026-08-02.
Analysis and cross-references: `docs/03_literature_and_positivity.md`.

## N-representability during propagation (the (P1) workstream)

| file | citation |
|---|---|
| `1411.0495v1` | Lackner, Březinová, Sato, Ishikawa, Burgdörfer — *Propagating two-particle reduced density matrices without wavefunctions*, Phys. Rev. A **91**, 023412 (2015) |
| `1611.00888v1` | Lackner et al. — *Accurate high-harmonic spectra from TD-2RDM theory*, Phys. Rev. A **95**, 033414 (2017) |
| `2303.04576v1` | Donsa, Lackner, Burgdörfer, Bonitz, Kloss, Rubio, Březinová — *Non-equilibrium correlation dynamics in the 1D Fermi-Hubbard model*, Phys. Rev. Research **5**, 033022 (2023) |
| `2412.13566v1` | Pescoller, Eder, Březinová — *Projective purification of correlated reduced density matrices* (2024). **The principled replacement for `Posify()`.** |
| `1204.4395v1` | Akbari, Hashemi, Nieminen, van Leeuwen, Rubio — *Challenges in truncating the hierarchy of TD-RDM equations*, Phys. Rev. B **85**, 235121 (2012). Failure-mode catalogue. |

## TDDM / ERPA — ground-state stability (P2) and strength positivity (P3)

| file | citation |
|---|---|
| `2009.00591v2` | Schuck, Delion, Dukelsky, Jemai, Litvinova, Röpke, Tohyama — *EOM method for strongly correlated Fermi systems and extended RPA*, Phys. Rep. **929**, 1 (2021). **Read first.** |
| `1003.0246v1` | Tohyama & Schuck — *Density-matrix formalism with three-body ground-state correlations*, Eur. Phys. J. A **45**, 257 (2010). The **(P2)** stability-matrix result; see §0 of the analysis doc for what it does and does not say. |
| `nucl-th_0305088v1` | Tohyama & Schuck — *Eigenstates of the TDDM theory*, Eur. Phys. J. A **19**, 215 (2004). Non-Hermitian ERPA matrix; **(P3)** strength-function remark. |
| `nucl-th_0305087v1` | Tohyama & Schuck — *Spurious modes in extended RPA theories*, Eur. Phys. J. A **19**, 203 (2004). Zero-mode / Kohn / continuity conditions. Standing regression test. |
| `nucl-th_0311016v1` | Tohyama, Takahara, Schuck — *Extended RPA with ground-state correlations*, Eur. Phys. J. A **21**, 217 (2004). STDDM. |
| `1603.08276v1` | Schuck & Tohyama — *TDDM meets Self-Consistent RPA*, Phys. Rev. B **93**, 165117 (2016). STDDM\*-b ≡ SCRPA in the one-body sector. |
| `1301.2026v1` | Tohyama & Schuck — *Odd-particle-number RPA and extensions*, Phys. Rev. C **87**, 044316 (2013). Parallel particle/hole sectors. |
| `1708.07219v1` | Tohyama — *Role of small-norm components in extended RPA* (2017). **Do not prune the term list by norm.** |

## TDDM truncation schemes (the non-MK closure, for comparison)

| file | citation |
|---|---|
| `1312.3703v1` | Tohyama & Schuck — *Truncation scheme of TDDM approach*, Eur. Phys. J. A (2014). 3-RDM ≈ antisymmetrized product of 2-RDMs. |
| `1708.00569v1` | — *…II* (2017). Normalization patch. |
| `1904.06780v1` | — *…III* (2019). Lipkin, large-N. |

## Not mirrored (no arXiv preprint) — cite from DOI

- Head-Gordon, Rico, Oumi, Lee, *Chem. Phys. Lett.* **219**, 21 (1994) — CIS(D). doi:10.1016/0009-2614(94)00070-0
- Head-Gordon, Oumi, Maurice, *Mol. Phys.* **96**, 593 (1999) — CIS(D₀)/CIS(D∞). doi:10.1080/00268979909482996
- Schirmer, *Phys. Rev. A* **26**, 2395 (1982) — ADC. doi:10.1103/PhysRevA.26.2395
- Trofimov & Schirmer, *J. Phys. B* **28**, 2299 (1995). doi:10.1088/0953-4075/28/12/003
- Hättig, *Adv. Quantum Chem.* **50**, 37 (2005) — CIS(D∞) ≡ P-ADC(2). doi:10.1016/S0065-3276(05)50003-0
- Christiansen, Koch, Jørgensen, *Chem. Phys. Lett.* **243**, 409 (1995) — CC2. doi:10.1016/0009-2614(95)00841-Q
- Dreuw & Wormit, *WIREs CMS* **5**, 82 (2015) — ADC review. doi:10.1002/wcms.1206
- Chatterjee & Pernal, *J. Chem. Phys.* **137**, 204109 (2012) — ERPA from 1-/2-RDM; exact-2-RDM-≠-exact-poles. doi:10.1063/1.4766934
- Coleman, *Rev. Mod. Phys.* **35**, 668 (1963); Garrod & Percus, *J. Math. Phys.* **5**, 1756 (1964) — D/Q conjugacy.
- Mazziotti, *Phys. Rev. E* **65**, 026704 (2002); *Chem. Rev.* **112**, 244 (2012); *PRL* **97**, 143002 (2006).

## Also in `refs/papers/` (JP's own work, fetched earlier)

`1411.5324v2` = Nguyen, Nanguneri, Parkhill, *How electronic dynamics with Pauli
exclusion produces Fermi–Dirac statistics*, J. Chem. Phys. **142**, 134113
(2015) — the rank-1 result this project builds on; LaTeX source in
`src_1411.5324/`. Also `1112.1665v2`, `1206.2567v2`, `1411.1110v1`.

---

## MK generalized normal ordering / generalized Wick theorem (`refs/papers/mkno/`)

Fetched 2026-08-02. PDFs + `pypdf`-extracted text, same convention as above.
Analysis: `docs/07_mkno_primary_sources.md`. Numerical audit: `tools/verify_gwt.py`.

### Full text obtained (open access)

| file | citation | why it's here |
|---|---|---|
| `2205.01178v1` | F. A. Evangelista — *Automatic derivation of many-body theories based on general Fermi vacua* (Wick&d), **J. Chem. Phys. 157, 062502 (2022)**, doi:10.1063/5.0097858 | **Read first.** §II is a complete modern restatement of MK GNO + GWT. Eq. (2) = GNO defining property; Eq. (3) = GWT with multi-leg terms; Eqs. (4)–(6) = γ/η/λ₂ contractions; Table I = orbital-subspace taxonomy. |
| `margocsy_thesis_elte2022` | Á. Margócsy — *Dressing of Interaction Vertices in Multi-Reference Coupled Cluster Theory*, PhD thesis, ELTE (2022), doi:10.15476/elte.2022.042 | **The most complete open derivation.** §2.1.3 cumulants + **particle–hole conjugation** (η₂=λ₂+½η∧η, η₃=−λ₃+…); §2.1.5 GNO; §2.1.6 GWT with worked products Eqs. (2.1.52)/(2.1.53). Also names the KNM *quasi-normal ordering* and its failure of ⟨{S}⟩=0. |
| `hanauer_thesis_mainz2013` | M. Hanauer — *Internally contracted multireference coupled-cluster methods*, PhD thesis, Mainz (2013), doi:10.25358/openscience-2706 | Extended version of Hanauer & Köhn, *Chem. Phys.* **401**, 50 (2012) (paywalled). §3.2.2 GNO/GWT (Eqs. 3.49–3.56); ch. 6 cumulant expansions to γ₅ and the **numerical cost of truncating cumulants** (Figs. 35–37). |
| `1602.05667v1` | Li & Evangelista — *MR-LDSRG(2)*, **J. Chem. Phys. 144, 164114 (2016)** | Third independent statement of the GWT rules, Eqs. (19)–(21); CAS simplifications. |
| `2506.13693v2` | Evangelista group — *Multireference EOM-DSRG: theoretical foundations and ionized states* (2025) | Published report that dropping high-rank cumulants causes **artificial variational collapse** of EOM energies. |
| `1406.0114v2` | Evangelista — *A driven similarity renormalization group approach to quantum many-body problems*, JCP **141**, 054109 (2014) | Context for the DSRG family. |
| `2106.07097v2` | Li & Evangelista — *Spin-free MR-DSRG* (2021) | Spin-free conventions for comparison. |
| `1512.06956` | Hergert, Bogner, Morris, Schwenk, Tsukiyama — *The In-Medium SRG*, **Phys. Rep. 621, 165 (2016)** | §3.2 single-reference normal ordering / Wick's theorem, as the λ=0 baseline. |
| `hergert_frontiers2020_tour` | H. Hergert — *A Guided Tour of ab initio Nuclear Many-Body Theory*, **Front. Phys. 8, 379 (2020)** | Nuclear restatement of MK GNO with irreducible densities λ^(k), Eqs. (25)–(26) — independent sign-convention cross-check. |
| `1408.6555` | Hergert et al. — *Ab initio MR-IMSRG for even Ca and Ni isotopes*, **Phys. Rev. C 90, 041302(R) (2014)** | MR-IMSRG use case citing KM-1997/CPL-1997 directly. |
| `1911.06836` | *Multireference Electron Correlation Methods: Journeys along Potential Energy Surfaces*, **Chem. Rev. 120, 5878 (2020)** | Review context. |

### NOT mirrored — paywalled, no open copy exists (checked Unpaywall/OpenAlex/S2/arXiv 2026-08-02)

**All four requested MK primaries are closed access. None was read in full text.**

- Kutzelnigg & Mukherjee, *Normal order and extended Wick theorem for a multiconfiguration reference wave function*, **J. Chem. Phys. 107, 432–449 (1997)**. doi:10.1063/1.474405 — *abstract obtained*.
- Mukherjee, *Normal ordering and a Wick-like reduction theorem for fermions with respect to a multi-determinantal reference state*, **Chem. Phys. Lett. 274, 561–566 (1997)**. doi:10.1016/S0009-2614(97)00714-8 — *no abstract obtainable*.
- Kutzelnigg & Mukherjee, *Cumulant expansion of the reduced density matrices*, **J. Chem. Phys. 110, 2800–2809 (1999)**. doi:10.1063/1.478189 — *abstract obtained; contains the λ₂=0 ⇔ γ²=γ statement.*
- Mukherjee & Kutzelnigg, *Irreducible Brillouin conditions and contracted Schrödinger equations… **I. The equations satisfied by the density cumulants***, **J. Chem. Phys. 114, 2047–2061 (2001)**. doi:10.1063/1.1337058 — *abstract obtained.* **This is the paper to cite for IBC_k/ICSE_k, not Part II.**
- Kutzelnigg & Mukherjee, *…**II. Spin-free formulation***, **J. Chem. Phys. 116, 4787–4801 (2002)**. doi:10.1063/1.1448827 — *abstract obtained*.
- Kutzelnigg & Mukherjee, *Chem. Phys. Lett.* **317**, 567 (2000) — the earlier IBC companion.
- Kong, Nooijen, Mukherjee, *An algebraic proof of generalized Wick theorem*, **J. Chem. Phys. 132, 234107 (2010)**. doi:10.1063/1.3439395 — *abstract obtained; defines the **quasi-normal-order theory** with one-body contractions only. Directly describes `MKNO.nb`'s algebra.*
- Shamasundar, *Cumulant decomposition of RDMs, multireference normal ordering, and Wick's theorem: a spin-free approach*, **J. Chem. Phys. 131, 174109 (2009)**. doi:10.1063/1.3256237 — *abstract obtained*.
- Kutzelnigg, Shamasundar, Mukherjee, *Spinfree formulation of RDMs, density cumulants and generalised normal ordering*, **Mol. Phys. 108, 433–451 (2010)**. doi:10.1080/00268970903547926.
- Hanauer & Köhn, *Meaning and magnitude of the reduced density matrix cumulants*, **Chem. Phys. 401, 50–61 (2012)**. doi:10.1016/j.chemphys.2011.09.024 — superseded here by Hanauer's Mainz thesis (mirrored).
- Neuscamman, Yanai, Chan, *Quadratic canonical transformation theory and higher order density matrices*, **J. Chem. Phys. 130, 124102 (2009)**. doi:10.1063/1.3086932.
- Yanai & Chan, *Canonical transformation theory for multireference problems*, **J. Chem. Phys. 124, 194106 (2006)**. doi:10.1063/1.2196410.
- Ruskai, *N-Representability Problem: Particle-Hole Equivalence*, **J. Math. Phys. 11, 3218 (1970)** — cited but not read.

### Could not obtain at all

- Misiewicz, Turney, Schaefer, *Reduced Density Matrix Cumulants: The Combinatorics of Size-Consistency and Generalized Normal Ordering*, **JCTC 16, 6150 (2020)**, doi:10.1021/acs.jctc.0c00422. Journal paywalled; the ChemRxiv preprint (doi:10.26434/chemrxiv.12197886) is behind a Cloudflare challenge that rejects scripted clients. Nothing in `docs/07` depends on it.
- Evangelista, *Multireference Theories of Electron Correlation Based on the DSRG*, **Annu. Rev. Phys. Chem. 70, 349 (2019)**, doi:10.1146/annurev-physchem-042018-052416 — Annual Reviews blocks scripted download.
