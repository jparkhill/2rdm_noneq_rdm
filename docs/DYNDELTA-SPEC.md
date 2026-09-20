# The dynamical-δ correction: complete specification and cost

**Purpose.** Define every ingredient of the closure propagated in episodes 8–9 (`tools/f0delta.py`,
`tools/dyndelta_gpu.py`), give the cost of each step as implemented, and show which steps can be reformulated
to bring the total from O(M¹⁰) to O(M⁶). Nothing here is new physics; it is the accounting needed to make the
method usable beyond M ≈ 12.

## 0. Notation

| symbol | meaning | size |
|---|---|---|
| `M` | number of spin-orbitals | |
| `N` | number of electrons (enters only as a parameter) | |
| `p,v,q,w,…` | spin-orbital indices, `0..M−1` | |
| `P = {(p,v): p<v}` | sorted pairs, `nP = M(M−1)/2` | |
| `γ₁[p,q] = ⟨a†ₚa_q⟩` | 1-RDM | `M×M` |
| `γ₂[p,v,q,w] = ⟨a†ₚa†ᵥa_w a_q⟩` | 2-RDM, ket `(p,v)`, bra `(q,w)` | `M⁴` |
| `X = Mat₂` | `X[(p,v),(q,w)] = 2γ₂[p,v,q,w]` on sorted pairs; `X = X†`, `X ≥ 0` | `nP×nP` |
| `λ₂ = γ₂ − γ₁∧γ₁` | 2-cumulant, `(a∧b)[pvqw] = a[pq]b[vw] − a[pw]b[vq]` | `M⁴` |
| `Tr₂` | `(Tr₂ T)[p,q] = Σ_v T[p,v,q,v]`; `γ₁ = Tr₂γ₂/(N−1)` | |
| `h`, `V` | one-body integrals; antisymmetrised two-body integrals, `H = Σ h a†a + ¼ Σ V a†a†aa` | |
| `A⊕A` | pair-space lift of a one-body matrix: `(A⊕A)[(pv),(qw)] = A[pq]δ_vw + A[vw]δ_pq − A[pw]δ_vq − A[vq]δ_pw` | `nP×nP` |
| `K_h = h⊕h`, `K_V[(pv),(qw)] = V[pvqw]` | drift generators | `nP×nP` |
| `N_op = γ₁⊕γ₁` | pair-space occupation | `nP×nP` |
| `W[p,q] = Σ_rs V[p,r,q,s]γ₁[r,s]` | mean field | `M×M` |
| `C[p,q] = Σ_rtu V[r,t,p,u]λ₂[r,t,q,u]` | rank-1 collision matrix | `M×M` |
| `ph(T)[(p,q),(v,w)] = T[p,v,q,w]` | particle-hole (G-matrix) arrangement | `M²×M²` |
| `Sw` | swaps the two entries of a ph pair index | `M²×M²` permutation |
| `D, Q, G` | pp, hh, ph Gram matrices (§3) | `nP, nP, M²` |
| `{B_k}` | real orthonormal basis of Hermitian `nP×nP` matrices, `k = 1..n2`, `n2 = nP²` | |
| `δ = Σ_k c_k B_k` | the correction, Hermitian on pair space; `c ∈ ℝⁿ²` | |

Every propagated quantity is `γ₂` (equivalently `X`) and the vector `c`. `γ₁` is never propagated.

## 1. The equations of motion

```
d/dt X   =  F0(X)  +  δ                                                  (1)
τ d/dt c =  −(c − c_hard(X))  −  κ Zᵀ ∇_c Φ(X, F0 + δ)  −  κ_neg Zᵀ ∇_c Φ_neg(X, F0 + δ)     (2)
```

`F0` is the exact λ₃=0 flow (§2). `c_hard` is the particular solution enforcing the two hard linear constraints
(§4). `Z` projects onto the null space of those constraints so that (2) never violates them. `Φ` is the smooth
tangency potential and `Φ_neg` the cone barrier (§5). Parameters used: `τ = 0.5 a.u.`, `κ = 100`, `κ_neg = 10³–10⁴`,
`σ = 10⁻³`, `ε = 10⁻³`.

Interpretation: δ relaxes on a timescale τ toward the minimiser of `½‖c − c_hard‖² + κΦ + κ_neg Φ_neg`. As τ→0 this
is the instantaneous projection of episode 7 (non-smooth, no spectrum); for finite τ it is a smooth ODE. Both
potentials are functions of the *flow* `F0 + δ` (through the cone rates), so the forces are linear in `c`.

## 2. `F0`: the exact λ₃ = 0 flow (episode 5c, verified 1e-15)

```
X_d   = Mat₂(γ₁∧γ₁),  L = X − X_d = Mat₂(λ₂)
S     = (Wᵀ⊕Wᵀ) − K_V N_op
Q_C   = ½(N_op K_C − K_{γC}) + ½(K_C N_op − K_{Cγ})        K_C = C⊕C, K_{γC} = (γ₁C)⊕(γ₁C), K_{Cγ} = (Cγ₁)⊕(Cγ₁)
P₁ = Sw(γ₁ᵀ⊗I)V_d L_d,  P₄ = Sw(I⊗γ₁)V_d L_d,  P₂ = L_d V_d(Sw(γ₁⊗I))†,  P₃ = L_d V_d(Sw(I⊗γ₁ᵀ))†
        with V_d = ph(V), L_d = ph(λ₂)
Φ_ph  = 2i[ back(P₁+P₂−P₃−P₄) − h.c. ]     back = ph⁻¹, antisymmetrise, → Mat₂/2
F0    = i[K_h + K_V, X] + i(S†X − XS) + i(Q_C − Q_C†) + Φ_ph
```

Cost: pair-space matmuls `O(nP³) = O(M⁶)`; ph-space matmuls `O((M²)³) = O(M⁶)`; `W`, `C` einsums `O(M⁴)`, `O(M⁵)`.
**`F0` is O(M⁶) time, O(M⁴) memory.** No `γ₃` is ever formed.

## 3. The cones and their rates

All three are affine in `(γ₁, γ₂)` and are what must stay PSD:

```
D[(pv),(qw)] = γ₂[p,v,q,w]                                   = X/2
Q[(pv),(qw)] = δ_pq δ_vw − δ_pw δ_vq − δ_vw γ₁[q,p] + δ_vq γ₁[w,p] + δ_pw γ₁[q,v] − δ_pq γ₁[w,v] + γ₂[w,q,v,p]
G[(pq),(rs)] = δ_qs γ₁[p,r] − γ₂[p,s,r,q]                                                        (M²×M²)
```

For a flow `F` (a Hermitian pair matrix, `F = d/dt X`), its induced rates are the same maps applied to
`(γ̇₁, γ̇₂) = (Tr₂ F/2(N−1), F/2)` with the constant terms dropped: `Ḋ = F/2`, `Q̇ = Q_lin(F)`, `Ġ = G_lin(F)`.
These maps are **fixed linear maps**, independent of the state. Cost of applying one: O(M⁴).

## 4. Hard constraints and the projector `Z`

Two properties of the exact flow are imposed exactly at all times:

(b) contraction: `Tr₂(F0 + δ)/(N−1) = γ̇₁^exact(γ₁, γ₂)` — the rank-1 EOM has no γ₃, so it is known exactly
    from γ₂; this is `M²` real equations (Hermitian `M×M`).
(c) energy: `Re tr(∇E · (F0 + δ)) = 0`, `∇E = K_h/(2(N−1)) + K_V/2` — one equation.

Both are linear in `c`: `C c = r(X)` with `C` an `(M²+1)×n2` **constant** matrix (rows = the constraint
functionals applied to each `B_k`), and `r(X) = [hvec(γ̇₁^exact − Tr₂F0/(N−1)); −Re tr(∇E F0)]`.
Precomputed once: SVD of `C`, its pseudo-inverse `C⁺`, and `Z` = orthonormal basis of `ker C` (`n2×(n2−M²−1)`).
Per step: `c_hard = C⁺ r(X)` — O(M² · n2) = O(M⁶); `Zᵀ(·)` and `Z(·)` — O(n2²) = **O(M⁸) as stored densely**.

**Reduction.** `Z Zᵀ = I − C⁺C` is a projector onto the complement of an (M²+1)-dimensional space, so
`Z Zᵀ g = g − C⁺(C g)` costs O(M²·n2) = O(M⁶). Never form `Z`. (Not yet done in code.)

## 5. The potentials

### 5.1 Tangency, `Φ`

For each cone `P ∈ {D, Q, G}` with eigendecomposition `P = U diag(x) U†`, define the spectral weight
`f(x) = σ²/(x² + σ²)` (≈1 for |x| ≪ σ, ≈0 for |x| ≫ σ; analytic) and `W_P = U diag(√f(x))`, so `W_P W_P† = f(P)`
is a smooth matrix function. Then

```
Φ = Σ_P ‖ W_P† Ṗ(F0+δ) W_P ‖_F²
```

is the squared norm of the cone rate restricted (smoothly) to the near-null subspace; the exact flow has
`Φ = O(σ)` because its rate entries vanish on null vectors (XR theorem). `Φ` is quadratic in `c`.

Gradient: with `R_P = W_P† Ṗ(F0+δ) W_P` (an `n_P×n_P` matrix, `n_P` = dim of the cone) and
`R_P^{(k)} = W_P† Ṗ(B_k) W_P`,
```
∂Φ/∂c_k = 2 Σ_P Re ⟨R_P^{(k)}, R_P⟩_F .
```

### 5.2 Barrier, `Φ_neg`

```
b(x) = ε² softplus(−x/ε)²        (C^∞; → x² for x ≪ −ε; → 0 exponentially for x ≫ ε)
Φ_neg = Σ_P Σ_i b(x_i)
d/dt Φ_neg = Σ_P Σ_i b'(x_i) u_i† Ṗ u_i              (first-order perturbation of eigenvalues)
∂Φ_neg/∂c_k := Σ_P Σ_i b'(x_i) u_i† Ṗ(B_k) u_i        (force on δ that raises negative eigenvalues)
```
`b'(x) = −2ε softplus(−x/ε) sigmoid(−x/ε) ≤ 0`, so the force pushes each negative eigenvalue's rate upward in
proportion to its depth. It is not literally a gradient of a potential in `c` (it is the gradient of `Φ_neg` with
respect to the *state*, chained through the flow); it is a smooth restoring force. Only eigenvalues with
`b'(x_i)` non-negligible contribute — typically a handful.

## 6. Where the O(M¹⁰) comes from, and how to remove it

The expensive object is `R_P^{(k)} = W_P† Ṗ(B_k) W_P` **for every basis element k** (n2 = O(M⁴) of them),
each an `n_P×n_P` matrix, giving O(M⁴ · M⁶) = O(M¹⁰) per cone per RHS, and it requires storing `Ṗ(B_k)` for all k:
O(M⁸) memory (51 MB at M=8, 16 GB at M=16). This was done to make the gradient a simple contraction.

**It is unnecessary.** Both forces have the form `Σ_k c-independent ⟨A_P, Ṗ(B_k)⟩_F` for some Hermitian
`A_P` (for Φ: `A_P = W_P R_P W_P†`; for Φ_neg: `A_P = Σ_i b'(x_i) u_i u_i†`). Since `Ṗ(·)` is a fixed linear map
and `{B_k}` is an orthonormal basis, the vector `g_k = ⟨A_P, Ṗ(B_k)⟩` is exactly the coordinate vector of the
**adjoint** `Ṗ†(A_P)` in the basis `{B_k}`:

```
∇_c (Φ-type term)  =  coords_B [ Σ_P Ṗ_P†( A_P ) ]
```

`Ṗ_D† = ½·id`; `Ṗ_Q†` and `Ṗ_G†` are the transposes of the affine maps in §3 — explicit index reshuffles plus a
partial trace, each O(M⁴). So the gradient costs: three eigendecompositions O(M⁶) (or O(k M⁴) iteratively,
since only the |x| ≲ σ band and the negative eigenvalues are needed), forming `A_P` O(M⁶), applying three adjoints
O(M⁴), and `c_hard` O(M⁶). **Total O(M⁶) per RHS, O(M⁴) memory**, with no stored basis images and no `Z`.

Remaining O(M⁶) pieces are the same order as `F0` itself, so the corrected closure costs a small constant times
the uncorrected one. This is the formulation to implement next (`tools/dyndelta_fast.py`); it must reproduce the
current RHS to round-off before use.

## 7. Cost table, per RHS evaluation (4 per RK4 step), N-independent

| ingredient | as implemented | reformulated |
|---|---|---|
| `F0` | O(M⁶) | O(M⁶) |
| exact rank-1 EOM | O(M⁵) | O(M⁵) |
| cone matrices + rates of `F0+δ` | O(M⁴) | O(M⁴) |
| eigendecompositions D, Q, G | O(M⁶) | O(M⁶), or O(k M⁴) iterative |
| basis-image rotation `W†Ṗ(B_k)W` ∀k | **O(M¹⁰), O(M⁸) mem** | eliminated |
| gradient assembly | O(M⁸) | O(M⁴) via adjoints |
| `c_hard`, projector `ZZᵀ` | O(M⁶) + **O(M⁸)** | O(M⁶) |
| **total** | **O(M¹⁰) time, O(M⁸) memory**; 0.10 s at M=8 (RTX 3090) | **O(M⁶) time, O(M⁴) memory** |

## 8. What is and is not structural (for the record)

Exact by construction: Hermiticity, antisymmetry, trace, contraction to the exact rank-1 EOM, free limit,
λ₃=0 exactness. Approached dynamically: energy (to the relaxation lag), cone positivity (bounded by κ_neg, or
inside with softplus). Not held: exact extensivity of δ (1–2% cross-fragment; the min-norm structure couples
fragments — a block-local norm would fix it and is compatible with §6), time-reversal symmetry (δ is dissipative),
Hamiltonian structure.
