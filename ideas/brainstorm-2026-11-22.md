# Brainstorm 2026-11-22: non-probabilistic levers on the exceptional set E

Setting (from I-008): E = {N : orbit never reaches 1} is completely invariant (T⁻¹E = E, T(E) = E),
closed under predecessors, contains the orbit of y* = min E > 2^71, and (Shaik Thm I.2) has positive
lower density ≥ d_z with z = Syr(y*). The question is whether any *structural* property of E can force
E = ∅, instead of the density race f(y) < D(y).

Scripts live in the session scratchpad (`.../scratchpad/ideas/`): `basin_rigidity.py`, `corr.py`,
`fixedpoint_rank.py`, `fixedpoint_lp.py`, plus `tools/density/predcount.c` run along the orbit of 27.

## 0. The unifying kill test (read this first)

**Sparse-modification test (M′).** Take the loneliest z₀ ≈ y, so d_{z₀} ≈ D(y) ≈ y^(−1.35). Define M′ = T
everywhere except on a sparse path S = {z₀, z₁, …} (zⱼ growing geometrically, each chosen lonely),
where M′(zⱼ) := zⱼ₊₁. Then:
- pred_{M′}(z₀) = pred_T(z₀), up to a negligible set, so M′ has a divergent orbit whose basin has positive
  density ≈ D(y). Shaik's lower bound and Tao's (log y)^(−c) bound both hold for M′, because a random
  orbit hits S with probability ≈ Σ D(zⱼ).
- So **any argument that is insensitive to changing T on a set hit with probability ≤ η cannot prove
  f(y) < η.** With η = D(y), this is *exactly* the I-008 threshold. That is why the exponent-1 (really
  1.35) barrier is sharp. It is not an artefact of Tao's method.
- What separates T from M′: M′'s E′ is T-backward-closed, but **not** T-forward-invariant at the
  points of S, which have density 0. So only **pointwise** use of the exact orbit of y* can work. Every
  density, residue or correlation property of E is shared by E′.

**Second structural remark.** Kill test 2 (3x−1 / negatives) only bites on *cycle* arguments. 3x−1 has
no known divergent orbits, so a no-divergence argument may legitimately also cover 3x−1 and every 3x+k.
For divergence the binding controls are 5x+1 and M′.

## Ideas

### 1. Residue-density rigidity (is the transfer fixed point forced?)
- **Statement.** Let ν(B) = density of E in a 3-adic class B. Complete invariance gives the closed
  equation ν(B) = ½ν(B/2) + (3/2)ν((3B+1)/2). It needs no parity data, because E ∩ evens = 2E. Since
  ν ≤ Haar, write ν = h·Haar. Then h is a *bounded harmonic function* of the 3-adic random affine walk
  y ↦ y/2 or (3y+1)/2 with probability ½ each. A Liouville theorem (from Tao-type fine-scale mixing)
  would give h ≡ d(E): E is equidistributed mod 3^k for every k.
- **Experiments.**
  - At finite level the fixed point is *not* forced. For each k there are 3^{k−1} independent equations
    on 3^k unknowns, so the solution space has dimension 2·3^{k−1} (`fixedpoint_rank.py`, k ≤ 6).
  - With all levels below K and positivity, the LP range for ν(0 mod 3) shrinks slowly toward 1/3:
    [0, .50] at K = 2 and [.26, .43] at K = 6 (`fixedpoint_lp.py`). So it is forced only in the
    projective limit, via mixing.
  - The three basins of 3x−1 (densities .327, .325, .348 at 10⁷) are equidistributed mod 3, 9, 27,
    81, 4, 8, 16, 5, 25 and 7, within 1.3% (`basin_rigidity.py`).
- **Needs.** A Liouville property for bounded harmonic h, i.e. Tao Prop 1.14-type mixing.
- **Kill tests.**
  - 5x+1: n/a.
  - 3x−1: ✘. Real nonempty invariant sets satisfy it.
  - 2-adic: ~ (3-adic only).
  - Conway: ~.
  - Heuristic: ✔.
  - M′: ✘.
- **Verdict: dead as a route.** It is probably a true lemma: every invariant set with a density is
  equidistributed in every residue class, and a renewal argument then gives a natural density. But it
  cannot separate d(E) > 0 from E = ∅.

### 2. Pair-correlation rigidity
- **Statement.** c(h) = density{n : n, n+h ∈ E} satisfies the exact relation
  c(2h) = ½c(h) + (3/2)·c_{[2]₃}(3h). Hope: positive-definiteness plus this relation forces c to be
  constant, i.e. E pseudo-random, which might clash with E being tiny.
- **Evidence.** In lattice coordinates (a, b) for h = 2^a 3^b h₀, the bounded solutions include
  non-constant modes z^a(2z−1)^b. Real 3x−1 basins have c(h)/d² = 1.96 at h = 1, falling to 1.19 at
  h = 81 (`corr.py`). Invariant sets are strongly clustered, not random.
- **Kill tests.** 3x−1 ✘ · M′ ✘.
- **Verdict: dead.**

### 3. Ergodicity / zero-one law on Ẑ
- **Statement.** T acts continuously on Ẑ (parity from Z₂), Haar is T-quasi-invariant, and an invariant
  clopen-measurable set has Haar measure 0 or 1. On Z_p for p ≥ 5 the two affine maps generate all
  translations; on Z₃, g₁ forces dependence to drop a level; Z₂ is Bernoulli. So if E were *rationally
  almost periodic* (Besicovitch-approximable by periodic sets), then d(E) ∈ {0, 1}. Since
  d(pred 1) > 0, that gives d(E) = 0, and Shaik then gives E = ∅.
- **Kill tests.** The 3x−1 basins are invariant, have density ≈ 1/3 and are **not** RAP. So "E is RAP"
  is exactly the 3x+1-specific content: a restatement.
- **Verdict: dead** (a clean reformulation: Collatz ⟺ E is RAP).

### 4. Generating function: mixed (2,3) Mahler equation and Pólya–Carlson
- **Statement.** F(x) = Σ_{n∈E} xⁿ satisfies F(y³) = F(y⁶) + y^(−1)·⅓Σ_j ω^(−2j) F(ω^j y²). By
  Pólya–Carlson, F is either rational (E eventually periodic, which rigidity plus d(pred 1) > 0 in
  every class excludes) or has |x| = 1 as a natural boundary. Hope: the equation continues F across an
  arc, which would force rationality.
- **Why it fails.** Every substitution in the equation maps the unit circle to itself (θ ↦ 2θ, 3θ, 6θ),
  so there is no continuation mechanism. The 3x−1 analog (y^(+1), the 3t+1 section) has non-rational
  {0,1}-solutions. Cobham-type theorems for Mahler functions (Schäfke–Singer, Adamczewski–Bell) need
  *separate* k- and l-Mahler equations, and this one is mixed.
- **Verdict: dead.**

### 5. Minimal element and bottom-layer constraints
- **Facts.** y* ≡ 3 (mod 4) and y* ≢ 2 (mod 3), and y* has infinite stopping time. Every element of E
  below 2y* is odd. Elements below (3/2)^j·y* avoid −1 (mod 3^j). For a divergent orbit, z_n < 3^{m_n}
  eventually, so z_n is the least residue of c_n·2^(−n) mod 3^{m_n}.
- **Why it fails.** All of these are tautological consequences of "orbit min ≥ y*", and that is
  exactly f(y*). The quantitative version also fails:
  - Exact 2-adic (Terras) counting bounds #E ∩ [N, 2N] ≤ N^(1−0.05) only for N ≤ y*^1.26 (the Chernoff
    rate is I₀ = 0.0349 per step).
  - Beating the floor X/y* would need an exponent of about 1/β ≈ 0.8: a factor of about 16.
  - Shaik's density is asymptotic, with no effective X₀.
- **Verdict: dead.**

### 6. Basin density is monotone along a divergent orbit
- **Statement.** pred(z_j) ⊂ pred(z_{j+1}), so a divergent orbit has d̄(pred z_j) ≥ d_{z₀} > 0 for all j.
  Hence no divergence ⟸ sup_{z≥Z} d̄_z → 0.
- **Evidence.** Along the orbit of 27 (predcount, X = 10⁸), d_z climbs from 0.056 (z = 82) to 0.39
  (z = 9232, z·d_z = 3583). "Trunk" points carry huge basins.
- **Why it fails.** The hypothesis is *stronger* than (B): it also forbids backward trunks with
  non-vanishing density, and it amounts to polynomial excursion bounds. Upper bounds on d̄_z from
  first-passage estimates are blind to returns after first passage, so M′ ✘.
- **Verdict: dead.**

### 7. Divergent orbits as 2-adic rationals (Lagarias periodicity conjecture)
- **Facts (checked).**
  - A divergent orbit has lim inf m_n/n ≥ log₃2. The proof: log z_n = log z₀ + m log 3 − n log 2 + o(n).
  - The series Σ_i 2^{k_i}·3^(−i−1) converges both in R (to z₀(C−1), where C = Π(1+1/3zᵢ)) and in Q₂
    (to −z₀).
  - So no divergence for 3x+1 on Z_(2), and hence for all 3x+k, is equivalent to the periodicity
    conjecture: Φ(aperiodic) ∉ Q.
- **Tools.** Subspace and Ridout need S-unit denominators, but the periodic-point denominators are
  2^a − 3^b. All that survives is a Liouville-type lemma:

  **Lemma.** If the parity sequence has W^{1+w} at position n, then 2^{w|W|} ≤ |z_n − z_{n+|W|}|, so
  w|W| ≤ log₂ max(z_n, z_{n+|W|}). The reason is that z_n ≡ z_{n+|W|} (mod 2^{w|W|}).

  For balanced (Sturmian) parity of slope p with partial quotients ≤ C, squares s_k s_k appear at
  positions ≤ (C+2)q_k. This excludes divergence when λ := p·log₂3 − 1 < 1/(C+3).
- **Kill tests.**
  - 5x+1: ✔ (partial only).
  - 3x−1: ✔ (aperiodic only).
  - 2-adic: ✔ (uses size).
  - Heuristic: ✔.
  - M′: ✔ (M′ is not affine).
  - Conway: generic, but partial.
- **Verdict: survives only as a partial-result generator** ("no structured divergent orbit").
  Random-looking parity sequences have Diophantine exponent 1, so it cannot reach (B).
- **First test.** Search the literature ("3x+1" Sturmian / low-complexity parity vectors). Then compute
  the exact threshold p₀(C) and check whether repetitions at linear scale are forced in every
  sequence of complexity ≤ Cn (Adamczewski–Cassaigne Diophantine exponent).

### 8. Ranking functions and automata (I-003, Yolcu–Aaronson–Heule)
- **Statement.** A certificate W with W(Tn) < W(n) is pointwise, so it passes every kill test,
  including M′.
- **New obstruction (sketch).**
  - **Block-additive W fails.** Take W = α·len + Σ φ(k-windows). Here ×3 is a 3-state carry
    transducer, and ΔW at odd steps is an additive functional with mean O(1). By the CLT/Livšic
    argument it must be a coboundary, otherwise some n have ΔW ≈ +√L. A coboundary makes
    W = κ·log n + O(1), and the climb 2^L−1 → 3^L−1 then raises W by κL·log(3/2).
  - **Matrix (weighted-automaton) interpretations evaluated on canonical binary strings** face the
    same issue. log W(Tn) − log W(n) has CLT fluctuations of order √L (Furstenberg cocycle) unless
    there is a cohomological miracle.
  - YAH's mixed-base intermediate strings may carry history and evade this.
- **Verdict: the only pointwise route left, with a new heuristic obstacle.**
- **First test.** Take YAH's published interpretations for the variants they *did* prove, and check
  two things: (a) whether W(Tn) < W(n) holds on canonical binary strings; (b) the empirical variance of
  log W(Tn) − log W(n) as L grows. If (b) grows like L for every interpretation they found, the
  obstruction is real and I-003 should target mixed-base certificates only.

### 9. Density plus Szemerédi plus exact affine action on progressions
- **Statement.** T^k(r + 2^k t) = T^k(r) + 3^{m(r)}·t. So 1_E(r + 2^k t) = 1_E(T^k r + 3^m t)
  *exactly*: E along 2-power progressions is affinely identical to E along 3-power progressions.
  Szemerédi then gives long APs in E, which T moves affinely.
- **Why it fails.** It is just the definition of invariance, rewritten. APs inside E move to APs inside
  E with positive increments, and nothing ever drops below y* unless f(y*) says so. 3x−1 basins share
  it ✘, and M′ ✘.
- **Verdict: dead.**

## Bottom line
- Nothing structural about E (residues, correlations, ergodicity, generating functions, minimality,
  APs) can force E = ∅. Each property is shared by real 3x−1 basins or by the M′ basin.
- The M′ test shows the I-008 threshold D(y) is *exactly* the limit of every sparse-robust method.
- What is left is pointwise:
  - Diophantine arguments for cycles (done to m ≤ 133);
  - 2-adic rationality for divergence (partial results only, idea 7);
  - ranking-function certificates (idea 8). Check the fluctuation obstruction before investing.
