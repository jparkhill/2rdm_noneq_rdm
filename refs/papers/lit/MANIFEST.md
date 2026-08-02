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
