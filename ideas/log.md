# Ideas Log

Every idea gets an entry, including the ones that die. A failure recorded with its
reason is progress. The kill tests are defined in `../README.md`.

## Template

```
### I-XXX: <name>
- Status: open | dead | partial result | parked
- Claim: <one precise sentence>
- Kill tests: 5x+1 [ ] · negatives [ ] · 2-adic [ ] · Conway [ ] · heuristic [ ]
- Evidence / experiments: <commands run, results>
- Notes / why it failed:
```

---

### I-001: Prove convergence from residues mod 2ᵏ
- **Status:** dead (known obstruction). Kept as a reference example.
- **Claim:** For a suitable k, every residue class mod 2ᵏ forces a drop below n within k steps.
- **Kill tests:** fails the 2-adic test.
- **Evidence:** `python3 tools/collatz.py bijection 12` shows every parity pattern of length k
  occurs for exactly one residue mod 2ᵏ, for both 3x+1 and 5x+1. In particular the
  all-odd pattern always occurs (the class n ≡ −1 mod 2ᵏ), and it increases for k steps.
- **Notes:** Residue-class arguments can say "density 1" (Terras) but never "all".
  Any further step has to use the archimedean size of n.

### I-002: Raise the no-m-cycles bound past Hercher's m ≤ 91 (attack line 1)
- **Status:** open, with a first round of results (2026-11-20, see "Results" below).
  The direct route does not work: 2⁷¹ alone does **not** raise m.
- **Claim:** No nontrivial cycle has m ≤ M for some M > 91.
- **Kill tests:** negatives ✓ by construction. The argument uses the inequality 2ᴷ > 3ᵐ,
  which fails for the negative cycles, so they're correctly excluded. 5x+1 ✓, since the
  same Diophantine setup with 5 has small solutions (the cycles found by `controls`).
- **Background:** A cycle with parity vector v of length K and m odd steps has
  n = c(v) / (2ᴷ − 3ᵐ) (`rational_fixed_point` in the toolkit). A positive cycle needs
  2ᴷ > 3ᵐ, so K/m > log₂3, and a small n means 2ᴷ − 3ᵐ is small relative to c(v).
  That is controlled by the convergents of log₂3 and by Baker / linear forms in logs.
  The m-cycle results combine this with the verification bound (all n < 2⁷¹), so the
  2⁷¹ result (Barina 2025) may already imply a bound above 91. **Check whether Hercher's
  method plus 2⁷¹ gives a new m bound. That could be the first concrete contribution.**
- **Tool:** `tools/cycle_bounds.py` re-implements Hercher's Theorem 21, Theorem 27, the
  Lemma 22 continued-fraction step, and the Theorem 23 iteration, with conservative
  decimal arithmetic at 150 digits.
  - `selftest` checks the continued fraction against brute force and against the
    published numbers.
  - `mcycles --x0 <units of 2^60>` runs the m-cycle argument.
  - `allm` computes the bound that holds for every m.
- **Validation** (everything reproduces exactly):
  - The Sinisalo-style bound at 2075·2⁶⁰ gives K ≥ 72,057,431,991 (matches the Breeden repo).
  - Theorem 27 jumps to 137,528,045,312 at exactly 2836·2⁶⁰ (matches Hercher's Remark 28).
  - Theorem 23 step 1 gives m₂ = 47, ε = 6.8·10⁻³², K > 5.27·10¹⁵ (paper: 47, 6.9·10⁻³², 5.2·10¹⁵).
  - At X0 = 704·2⁶⁰ the result is **exactly m ≤ 91**. The m = 91 chain is
    2.74e19 → 7.77e19 → 2.06e20 → 7.94e21, matching the paper, and m = 92 stalls.
- **Results:**
  1. **2⁷¹ (Barina 2025) leaves the m-bound at 91.** For m = 92 the K lower bound now
     climbs to 2.06·10²⁰ (it stalled at 7.77·10¹⁹ before). Using m₂ = m needs about
     (δ⁹²−1)/(δ−1)·log₂(1.67·X0) ≈ 3.1·10²⁰, so the chain stops short.
  2. **Verification needed for each m** (from a binary search; both the interval and
     the direct continued-fraction methods agree):

     | exclude all m ≤ | need X0 ≥ | vs. Barina's 2075·2⁶⁰ |
     |---|---|---|
     | 92 | ≈ 14953·2⁶⁰ ≈ 2^73.87 | 7.2× |
     | 93 | ≈ 21899·2⁶⁰ ≈ 2^74.42 | 10.6× |
     | 94 | ≈ 27735·2⁶⁰ ≈ 2^74.76 | 13.4× |
     | 95 | ≈ 33588·2⁶⁰ ≈ 2^75.04 | 16.2× |
     | 96 | ≈ 39425·2⁶⁰ ≈ 2^75.27 | 19.0× |
     | 101 | 2⁷⁷ | 64× |

     Angeltveit (arXiv:2602.10466) estimates 2⁷⁷ is feasible with Barina-scale resources.
     **If that verification happens, Hercher's method gives m ≤ 101 with no new theory.**
  3. **All-m bound:** Theorem 27 alone still gives K ≥ 72,057,431,991 at 2075·2⁶⁰.
     Hercher's Corollary 29 (backed by a 5-week computation) needs only X0 ≥ 1536·2⁶⁰,
     which Barina has now passed. So **every nontrivial cycle has K > 1.375·10¹¹ odd
     elements**, conditional on Barina's run and Hercher's computation. The Breeden repo
     reports only 7.2·10¹⁰, so it seems to miss Cor. 29.
  4. The direct continued-fraction method (N − Kδ < η instead of the Lemma 22 interval)
     converges in fewer steps but reaches the same fixed point, so it gives no extra m.
- **Round 2 (2026-11-20): primary-source check.** Text extracted with pypdf from the
  arXiv v3 PDF, the JIS corrigendum, and Simons–de Weger v1.44 (2010).
  1. **The Theorem 21 transcription is exact**, including the special cases m₂ = m−1 and m₂ = m.
  2. **The corrigendum fixes a real gap.** The original proof bounded T for the later
     window minima using Lemma 20, but Lemma 20 gives *upper* bounds on later minima,
     which can't bound T from above. Xinjun Wang caught this. The corrigendum's fix:
     - Set xᵢ = log₂(nᵢ+1) and pick the window that maximizes Σxᵢ, which gives the
       wrap-around constraint x_first < δ·x_last.
     - S = Σ 3/(2^{xᵢ}−1) is convex, so its maximum is at a vertex.
     - That vertex is the geometric chain v, δv, ….

     **I checked the fix and it holds.** One step it leaves implicit: the feasible region
     is unbounded. That's fine, because recession directions are nonnegative (from xᵢ ≥ 1)
     and S is decreasing, so the maximum is still at a vertex.
  3. **The premise is stronger than needed.** Validity only needs (δ^{m₂}−1)/(δ−1) < (m₂/m)K.
     The log₂(162/97·X0) factor only makes the bound beat Theorem 14. Relaxing the
     premise in the code changes nothing (tested at 704 and 2075).
  4. **The bottleneck idea from round 1 is dead.** Corollary 13 is a *disjunction*: one
     of three partial-sum bounds holds. It is not a per-term bound of 97/54·1/X0. The
     single-term bound 3/n (Remark 7) is essentially sharp: take n = a·2ᵏ − 1 just above
     X0 with k ≥ 60, which gives T(n) ≈ 3/X0.
  5. **The Simons–de Weger bound is confirmed.** Theorem 3(d) gives K < 1.4784·m·δᵐ for
     91 ≤ m ≤ 515619 (and 1.3993 for 78 ≤ m ≤ 90). It comes from their Lemma 7
     (Λ < m·cₘ·2^{−(δ−1)K/(δᵐ−1)}, the same vertex idea) together with a lower bound on Λ
     from the largest relevant partial quotient A(m) = 55. They say its dependence on X0
     is negligible.
  6. **Exhaustive survivor enumeration** (`cycle_bounds.py survivors`, checked against
     brute force) finds every K in [iterated lower bound, S–dW upper bound] that passes
     all the conditions: Theorem 16/21 for every m₂ plus continued fractions.
     - At X0 = 2075·2⁶⁰, **for m = 92 exactly one K survives:**
       **K\* = 205,632,218,873,398,596,256** (a convergent denominator of log₂3), with
       **L\* = 120,287,136,981,023,372,109** and gap ⌈K\*δ⌉ − K\*δ = 1.284·10⁻²².
     - For m = 93…98 the survivors are exactly K\*, 2K\*, …, ⌊U/K\*⌋·K\*. Multiples
       appear at m = 99 and above.
     - **So "no 92-cycles" is equivalent, within this framework, to ruling out a
       92-cycle with exactly (K\*, L\*) = (205632218873398596256, 120287136981023372109).**
  7. **Consequence:** sharpening the S–dW upper bound can't help with m = 92, since K\*
     sits below it and satisfies every known inequality. The only routes are:
     - (a) raise X0 to about 2^73.9 (the table above);
     - (b) a better single-minimum bound in the m₂ = m−1 case (shown sharp for one term,
       so it would need the window's structure);
     - (c) a new argument aimed at the single pair (K\*, L\*, m = 92). A cycle there
       would need 91 consecutive near-maximal climbs (each xᵢ₊₁ ≈ δxᵢ) followed by one
       crash back to about X0. That is very rigid, which suggests a targeted residue or
       structure argument could work.

  **Next step: try (c).** Model the vertex-extremal shape of a (K\*, L\*, 92)-cycle and
  see whether the congruence conditions nᵢ ≡ −1 (mod 2^{kᵢ}) plus the exact totals give
  a contradiction. Pass the 2-adic kill test first: it can't rest on residues alone.
- **Round 3 (2026-11-20): attacking K\* directly. Result: no 92-cycles (claimed; needs
  independent review).** Tool: `tools/kstar_sieve.py`. Driver: `tools/run_kstar_all.py`.
  Toy validation: `tools/test_kstar_sieve.py`.
  - **Setup.** Let n₁ be the cycle's smallest element. Write each local minimum as
    nᵢ + 1 = aᵢ·2^{kᵢ} with aᵢ odd, and set yᵢ = log₂(nᵢ+1), αᵢ = log₂aᵢ = yᵢ − kᵢ ≥ 0.
    Then yᵢ₊₁ = yᵢ + (δ−1)kᵢ − lᵢ + εᵢ, with 0 < εᵢ < 1/(X0 ln 2), and
    Eᵢ := δyᵢ − yᵢ₊₁ = (δ−1)αᵢ + lᵢ − εᵢ ≥ 1 − ε.
  - **Fact 1 (budget).** Σᵢ<ₘ G(m−i)Eᵢ + Σᵢ αᵢ = G(m)·y₁ − K exactly, with
    G(n) = (δⁿ−1)/(δ−1). Every term on the left is ≥ 0, so any partial sum is at most
    G(m)·y₁ − K, which is linear in y₁ once the first runs (kᵢ, lᵢ) are fixed.
    For m = 92 and K = K\*, κ = K/G(92) = 47.69, while y₁ ≈ 73. So the first minima must be
    "almost all ones" in binary: each unit of αᵢ costs about 0.37·y₁ of a budget of only y₁ − κ ≈ 25.
  - **Fact 2 (window).** Σ over odd x of log(1 + 1/(3x)) = g·ln 2 exactly, where
    g = N − Kδ = 1.284·10⁻²². With T(n₁) ∈ [3c/(n₁+1), 3c/n₁] and c = 1 − (2/3)^{k₁},
    and with Fact 1 forcing y₂ ≥ 77.3, y₃ ≥ 121.5, …, this pins
    n₁ ∈ [3c/A_hi − 1, 3c/(A − R)] where A = 3g ln 2 and R ≤ 0.06A. So y₁ ∈ [71.02, 73.36]
    and each k₁-window is about 2^69.3 wide. Slicing in y₁ is needed: above y₁ ≈ 75.6
    the budget allows n₉₂ to be small, and those slices are inconsistent.
  - **Sieve.** A DFS over n₁'s exact run lengths (k₁, l₁, k₂, …), pruned by Fact 1 (plus
    αᵢ ≥ 0 and nᵢ ≥ n₁). It stops once the forced parity prefix is long enough to leave at
    most about one candidate in the window. Every candidate is then checked exactly
    (big-integer orbit, Fact 1 with exact yᵢ, and Fact 2's T-sum).
  - **Proof draft for reviewers:** https://claude.ai/artifact/RWYbwNv2UHcxLnyGYSwK11
    (LaTeX version: `paper/no-92-cycles.tex` → `paper/no-92-cycles.pdf`, 6 pages;
    build with `tectonic no-92-cycles.tex`)
  - **Result for (m, K, L) = (92, K\*, L\*):** 1,268 DFS nodes and 1,257 leaves, with forced prefixes of
    **87–188 bits** against a 2^69.3 window. The expected number of candidates is
    3·10⁻⁴, and the actual number is **0**. The same holds at X0 = 2⁷¹ exactly.
    Combined with round 2 (K\* is the only K left for m = 92), **there is no Collatz
    92-cycle**, conditional on the inputs below.
  - **m = 93…98:** only the largest surviving multiple jK\* closes the window and is ruled
    out (m=93: 2K\*; 94: 3K\*; 95: 4K\*; 96: 6K\*, 8 candidates all killed at minimum 3).
    The smaller multiples have too much budget slack (κ much less than y₁), so R ≥ A and
    the window never closes. **Open: K\* for m = 93, and so on.** A sharper bound on the
    later minima (an LP instead of the single-ρ bound in `_ylin`) might close more cases.
  - **Validation:**
    - Toy brute force: 6 cases with m = 3–6 and n₁ up to 2²¹. The sieve's survivors exactly
      equal the brute-force set.
    - Stress test: the real window with a budget weakened to 0.92·K\* produces 4
      candidates, and with 0.90·K\* produces 9. All are killed by the exact check, so the
      emit path works at full scale.
    - Instrumented prefix-length histogram (above).
  - **Kill tests.** 2-adic: passes, because Fact 2 uses n₁'s archimedean size (the
    T-window) and the pigeonhole compares prefix length to window size. Conway, 5x+1,
    negatives: the argument is a computation for one specific (3, 2, m, K, L), with no
    generic step.
  - **Before telling anyone:** (1) have the derivation checked by someone else; (2) port
    the sieve to exact rational arithmetic, or audit every Decimal rounding direction;
    (3) check whether the ResearchGate "m ≤ 95" paper already covers this.
- **Round 4 (2026-11-21): prior work found. The result is not new; hardening stopped.**
  - Wang (m ≤ 95 preprint, June 2026), Alick (m ≤ 101, Lean, kernel-checked), Palmqvist
    (m ≤ 104, claimed), Marshall (m ≤ 111 Lean; m ≤ 133 computed). See bibliography 9.3–9.6.
  - Alick's "two-rise bound" plus "chain bound" is the same mechanism as our Corollary
    A1 search plus Lemma A. His worst case K = 77,692,117,359,936,589,403 is the K\*
    our tools found at X0 = 704·2⁶⁰. So the independent rediscovery is consistent.
  - Possible value left: an *independent* cross-check of the computed (not Lean) ranges,
    m = 112–133 in Marshall's repo, using our separately written sieve. Their README says
    those ranges have no human referee and only a 0.37% margin at m = 127.
- **Caveats:** everything is conditional on Barina's computation (X0 = 2075·2⁶⁰, or 2⁷¹
  from the paper), Hercher's Theorems 14, 21 and 27 (with the corrigendum) and his Cor. 19
  program, and the Simons–de Weger upper bound K < 1.4784·m·δᵐ.
- **Unvetted competition:** a ResearchGate item titled "Non-existence of Collatz m-cycles
  for m ≤ 95" (authors unknown; the page returns 403). Check it before claiming anything.
  If it's correct, it must get past the K\* obstruction, which would be informative.

### I-005: Independent cross-check of "no m-cycles for m ≤ 133" (Marshall et al.)
- **Status:** done for now (2026-11-21).
  - Chain and shape stages: independently confirmed for m = 92–133.
  - Tables: no unsafe entry found.
  - **All 3 entries where +1 alone breaks m = 133 are confirmed exact**, including
    B₄(317) = 337, by the C meet-in-the-middle solver (see below).
- **Code:** `tools/xcheck/`.
  - `rigorous.py`: exact rational interval arithmetic with series and tail bounds; no floats.
  - `rises.py`: own B_r(N) solver (closed-form 2-adic ν plus Terras least representative).
  - `mchain.py`: own chain theorem, own n₀ bound, own exhaustive search for dangerous shapes.
  - `test_mchain.py`: controls.
  - Tables from `marshall_data/` (Apache-2.0, unmodified) are the only shared input.
- **What is independent:**
  - The n₀ bound is re-derived from Hercher's Thm 21 with the corrigendum. Its window is
    the one maximising Σx, with m₂ = 0…m.
  - The chain lemma is re-derived from 2^(k+1)(n'+1) ≤ 3^k(n+1) + 2^k. Our rule 1 is
    δY − 1 + ε, slightly stronger than Marshall's δY.
  - All arithmetic is rigorous intervals at 320 bits.
  - Shapes: instead of enumerating Marshall's 84M survivors, each dyadic K-range below
    the Rhin ceiling gets a gap threshold η. F(Y) is monotone, so every shape with
    gap ≥ η is provably closed. The shapes below η are listed exhaustively with a
    unimodular convergent basis and tested exactly.
- **Results:**
  - m = 92…133: **all closed** in both modes. For m ≤ 132 no shape falls below η. At
    m = 133, 2 shapes do, and both close; the worst is K = 27,444,133,206,411,171,953 with
    ratio **0.9236463** (own) / **0.9236464** (compat).
  - Compat ratios equal Marshall's refchain to 7 digits on his three published hard shapes
    (0.9236464, 0.5171293, 0.7649742).
  - The full sweep takes about 2.5 min, against about 4 h for 84M shapes.
- **Controls:**
  - lattice() vs brute force: 40 random ranges, no misses.
  - Negative control: m = 134 and 135 report open shapes, with the hard shape at 1.464
    and 2.320, i.e. growth ×1.585 per m as predicted.
  - η inflated ×10⁶ at m = 120: 3,909 extra candidates, all closed.
  - All 5,174 (m, K) pairs from one-sided best approximations of δ are closed by the
    exact test.
- **Tables** (own solver, `rises.py`, validated against brute force for r = 1–4):
  - B₂ at all 88 grid points with N ≤ 10,000: 74 exact, and 14 where the table is 1–2 too
    *high* (safe). These are exactly the 14 that Marshall discloses for `grid_own.json`.
  - B₃ at all 39 grid points with N ≤ 400: all exact.
  - B₄ at all 7 grid points with N ≤ 80: all exact.
  - **0 unsafe entries.** About 6.5 CPU-h in total.
- **Sensitivity** (`sensitivity.py`, `sensitivity_points.py`). The hard shape at m = 133
  uses 41 distinct table entries. Raised by +1 one at a time, **only three break m = 133**:
  - B₄(73) = 87: independently checked, exact.
  - B₄(119) = 136: independently checked, exact (3,570 s; witness runs (58,1,22,3,23,2)).
  - B₄(317) = 337: **independently checked, exact.** It was out of reach for the Python
    solver (about 10 CPU-days), so `b4mitm.c` handles it. That solver splits
    ν = A(k₀,l₀,k₁,l₁) + 2^(l₀+l₁)·R(K₂,l₂), joins on a w-bit window key (carry 0/1) and
    tests every match exactly mod 2^512.
    - Bounds: K₂ ≤ B₃(325) = 340 (already verified), and l₂ < N + (δ−1)K₂.
    - Validated against rises.py at N = 14…40 (even), 72, 80 and 119: max = B₄, nothing
      at B₄ + 1, witnesses re-simulated.
    - N = 317: T = 338 gives no hits (418M left halves, 1.69B exact checks, 117 s).
      T = 337 gives exactly one hit, runs (103, 1, 117, 1, 102, 2). Re-simulated as an
      actual 317-bit odd n with rises 103, 117, 102, 15 = 337.
  - All B₂ and B₃ entries survive +1 individually. B₂ tolerates at least +63 on every
    unchecked entry at once. A uniform +1 on all unchecked B₃ (or B₄) entries breaks m = 133,
    with the hard shape going 0.92 → 1.46.
- **Not checked independently:** the other large table entries (B₂ above N = 10,000, B₃
  above 400, B₄ above 119 except 317) and Rhin's bound. None of them breaks m = 133 when
  raised by +1 alone. `b4mitm.c` could check the remaining B₄ points (all N ≤ 852) in
  roughly hours of CPU. B₂ and B₃ at large N would need their own fast solvers.

### I-003: SAT and rewriting-system termination (attack line 2)
- **Status:** open. Reproduction done in iteration 7 (see I-015); a long run on the new system H is under way.
- **Claim:** Find a matrix or arctic interpretation proving termination of a stronger
  Collatz variant than they reached.
- **Kill tests:** Conway is the main risk. The rewriting system encodes the specific map,
  so a termination proof is specific to it, but check that the same encoding for 5x+1 is
  *not* provably terminating.
- **Notes:** The best fit for CS skills. Needs a SAT solver (kissat / cadical) and a
  termination tool (e.g. matchbox, AProVE, TTT2).

### I-004: Push the Krasikov–Lagarias exponent above 0.84 (attack line 3)
- **Status:** dead (superseded), 2026-11-21. The 2026 positive-density claims (bibliography
  9.8–9.10) give ≥ cX of the numbers below X reaching 1, which beats x^c for every c < 1.
  Reopen only if those claims fail review.
- **Claim:** At least x^c of the numbers ≤ x reach 1, for some c > 0.84.
- **Notes:** The method is a system of difference inequalities over residue classes mod 3ᵏ,
  solved as a linear program. More compute and a larger modulus might give a small gain.
  Check whether the exponent saturates as k grows (their paper may say).

### I-006: Tao's method: what is left (reading 2026-11-21)
- **Status:** open (study done; no attack started).
- **Claim:** Prove a *good* explicit density of numbers reaching 1 (for example ≥ 1/2), using
  polynomial error terms plus verification to 2^71. The current claimed constant is
  2^(−2^(2^163)).
- **Kill tests:** 5x+1 [✔ uses negative drift] · negatives [✘ Tao's estimates hold verbatim
  for 3x−1; only the finite verification is specific to 3x+1] · 2-adic [✔ uses size and
  3-adics] · Conway [~ the arithmetic enters only via offset injectivity and the triangle
  lemma] · heuristic [✔ rigorous].
- **Notes:** Full reading in `sources/notes/tao2019.md`. Key facts:
  1. Collatz ⟺ the set reaching 1 has density **exactly** 1 (via predecessor density, 9.10).
     Density 1 − ε proves nothing, so this line yields quantitative theorems, not a proof.
  2. The cheapest first step is to compute max_ξ |E e(ξ·Syrac(Z/3ⁿ)/3ⁿ)| for n ≤ ~15 with an
     FFT transfer operator, to see the true decay rate (Tao asks whether it is exp(−cn)).
  3. Active AI-driven competition (ProofAtlas, AxiomProver, Shaik).

### I-007: True decay rate of the Syracuse characteristic function (cycle 1, 2026-11-22)
- **Status:** done. The result is a confirmation, not new: it matches a 2026 conjecture by WakauraH.
- **Claim tested:** Tao asks (after Prop 1.17) whether max_{3∤ξ} |E e(ξ·Syrac(Z/3ⁿ)/3ⁿ)| decays like
  exp(−cn).
- **Method:**
  - `tools/fourier/syrac_decay.py`: exact law in discrete-log coordinates. Every unit ξ is 2^k, so a
    single cross-correlation FFT gives every unit frequency. Validated against Tao's Z/9 table and
    brute force. Exact for n ≤ 17.
  - `tools/fourier/band_recursion.py`: Φₙ(k) = Σ_a 2^(−a)·e(2^(k−a)/3ⁿ)·Φₙ₋₁(k−a). Paths only move
    down in k, so a band k ≥ −J0 suffices. It matches the exact values for n ≤ 17 to all printed
    digits and does not change with J0 ∈ {30, 60, 120}. Run to n = 1200.
  - `tools/fourier/mc_decompose.py`: Monte Carlo split of the extremal coefficient by trajectory.
- **Result:**
  - max|Φₙ| ≈ C·e^(−I₀n)·n^(−c), with I₀ = I(log₂3) = (x−1)·ln(2−2/x) + ln(2/x) at x = log₂3,
    which is 0.05498 nats (0.0793 bits), and c ≈ 1–1.3.
  - The extremal frequency is ξ = 2^(k*) with k* ≈ ⌊n·log₂3⌋ − 7, so 2^(k*) ≈ 3ⁿ: the cycle-equation
    balance.
  - I₀ is the Cramér rate for a total valuation s_n ≈ n·log₂3, i.e. a trajectory that neither grows
    nor shrinks. The mass of the extremal coefficient comes from walks that stay inside Tao's main
    black triangle until the end.
- **Retraction:** for n ≤ 17 the data look like exp(−1.06√n). That is a transient; −log/√n climbs
  from 1.04 to 1.65 by n = 600.
- **Prior work:** github.com/WakauraH/collatz-barrier-repro conjectures δ_m = Θ(2^(−I₀m)·m^(−c)),
  with I₀ = (1−H(θ))/θ in bits. That equals ours for θ = log₃2. They report certified numerics to
  m = 1200 and a conditional proof. Our computation is independent and agrees.
  - Possibly new: the explicit extremal ξ, the offset of 7, and the large-deviation interpretation.
    Their README defines none of these, so this is unverified.
- **Kill tests:** this is a statement about a random variable, not a proof attempt. 3x−1 has the
  conjugate characteristic function, so the same rate holds there. No route to Collatz.
- **Value:** if the rate e^(−0.055n) were proved, it would feed I-006 item 1 (good constants). But
  0.055 is small, so explicit constants stay weak.

### I-008: The density race, and the exponent-1 barrier (cycle 2, 2026-11-22)
- **Status:** open. It is a reformulation that locates the difficulty precisely. It is not a proof.
- **Lemma (elementary, given Shaik's Theorem I.2):**
  - Definitions. E = {N : the orbit of N never reaches 1}. f(y) = upper density of
    {N : Col_min(N) ≥ y}. For y ≥ 2 this set lies inside E. d_z = lower density of the
    predecessor basin of z.
  - Claim. Suppose E ≠ ∅, and let y* = min E. Then y* is odd and y* > 2^71. Also z = Syr(y*)
    satisfies 3 ∤ z and z ≤ (3y*+1)/2, and every N in pred(z) lies in E. So Col_min(N) ≥ y* and
    f(y*) ≥ d_z > 0.
  - **Consequence.** Collatz follows if f(y) < D(1.5y+1) for every odd y > 2^71, where
    D(Z) = min over z ≤ Z with 3 ∤ z of d_z. Collatz implies f ≡ 0, so this is an equivalence.
  - This is a "zero-or-gap" principle: the exceptional density is either 0 or at least the basin
    floor.
- **Scale of the floor:**
  - Numerically (`tools/density/predcount.c`, X = 10⁸ and 10⁹), z·d_z is O(1) for typical z: median
    about 7 for z ~ 10⁴, 1% quantile about 0.6, and stable in X.
  - Heuristically, z·d_z is the limit of a branching-tree martingale, so D(Z) ≈ Z^(−1−o(1)).
  - Shaik's proven floor is much weaker. It has W_z ≍ 1/z, but the density is ≈ W·2^(−C·W^(−0.011)),
    which is stretched-exponentially small in z.
- **Barrier:** Tao-type bounds give f(y) ≪ (log y)^(−c). Proving f(y) = o(y^(−1−ε)) is *equivalent*
  to Collatz (modulo a finite check). So the quantitative Tao program (I-006) can never pass rate
  exponent 1 without proving the conjecture. Any rate y^(−A) with A < 1 stays merely quantitative.
- **Kill tests:**
  - 3x−1: ✔. Its cycles have min 1, 5 and 17, so f ≡ 0 above 17 is consistent.
  - 5x+1: ✔. Tao fails there because the drift is positive.
  - 2-adic: ✔. The argument uses size.
  - Conway: ~. The predecessor theorem uses 3x+1 structure.
  - Heuristic: ✔. Both sides are rigorous statements.
- **Floor exponent (cycle 3; details in `sources/notes/predecessor_floor.md`):**
  - Heuristic truth: the left tail of W = z·d_z is P(W < w) ≈ w^2.9 (6,000 samples, X = 10⁸). So
    min over z ≤ Z of d_z ≈ Z^(−1.35).
  - Shaik as written: d_z ≈ exp(−c·z^0.011). The cause is polynomial mixing, from both the Fourier
    input and his typicality windows.
  - Conditional on exponential decay at rate I (I-007 gives I ≈ 0.055), and with the windows
    replaced by a linear envelope (a subagent's derivation, unrefereed): d_z ≥ c₀·z^(−B) with
    B ≈ 3 + 4.15/I ≈ 76. The best possible for that route is B ≈ 14.5.
  - So a proven race needs f(y) = o(y^(−76)), while the truth is about y^(−1.35). Both are far
    beyond Tao's (log y)^(−c).
- **Perturbation barrier (cycle 3 brainstorm; `ideas/brainstorm-2026-11-22.md`):**
  - Construct M′: equal to T except on a sparse path that starts at the z₀ with the smallest
    basin density, and is forced to diverge.
  - Every estimate that is insensitive to changes on a set hit with probability η still holds for
    M′. So no such method can prove f(y) < η.
  - With η = D(y), this is exactly the race threshold, so the barrier is sharp. Closing the race
    needs pointwise use of the actual orbit of y*.
  - The brainstorm also killed: residue-class rigidity, pair correlations, profinite zero-one
    laws, generating functions, minimal-element constraints, and Szemerédi-type invariance.
    Real 3x−1 basins are the counterexamples.
- **Paper:** `paper/density-race.tex` (built to PDF).
- **Next:**
  1. Improve the floor to D(Z) ≥ Z^(−B). This is a concrete target, plausibly attainable given
     exponential decay.
  2. The other side, A > B, is the conjecture itself.

### I-009: Count m-cycles with exponential sums mod q = 2^K − 3ⁿ (cycle 3, 2026-11-22)
- **Status:** dead.
- **Claim:**
  - A cycle with n odd steps and K halvings solves N·q = F(a) with q = 2^K − 3ⁿ, where
    F(a) = Σ 3^(n−i)·2^(s_{i−1}).
  - So the number of patterns with q | F is #pat/q + (#pat/q)·Σ_{t≠0} E e(tF/q).
  - Since #pat ≈ C(K−1, n−1) ≈ 3ⁿ·e^(−I₀n) and q ≥ 3ⁿ/n^C (Rhin), the main term is below 1. So
    enough cancellation would rule out every long cycle.
- **Kill tests:**
  - 3x−1 and 5x+1: ✔. Their known cycles are short, and the heuristic only claims n ≥ n₀.
  - **Square-root barrier: ✘.** Getting a count below 1 needs Σ_{t≠0} |E e(tF/q)| < e^(I₀n). But even
    square-root cancellation (|E e| ≈ #pat^(−1/2) ≈ 3^(−n/2)) summed over q ≈ 3ⁿ frequencies gives
    about 3^(n/2). That misses by an exponential factor. Fourier equidistribution cannot show that a
    set whose expected size is below 1 is empty.
- **Notes:** This is the standard heuristic for why there are no long cycles. Turning it into a
  proof would need an algebraic, non-averaged argument, which is what Baker-type bounds plus
  computation provide.

### I-010: Low-complexity parity sequences cannot diverge (cycles 3–4, 2026-11-22)
- **Status:** done, and largely known (literature check in cycle 5).
  - [arXiv 2101.12747](https://arxiv.org/abs/2101.12747) (the 3x+1 periodicity conjecture in ℝ)
    shows that no rational odd 2-adic number has a Sturmian parity vector of slope ≠ ln2/ln3.
  - At slope exactly ln2/ln3, a positive integer's orbit is bounded, hence eventually periodic,
    which a Sturmian word is not.
  - So for integers the Sturmian case is covered. Our general complexity criterion
    (p(L) ≤ L/ε) is at most a small generalization. Related: Fernández–Ibáñez,
    [arXiv 2607.24844](https://arxiv.org/abs/2607.24844) (Christoffel words as extremal
    structures; about cycles).
- **Lemma (repetition):** if the parity words of z_n and z_{n+P} (accelerated map) agree for L
  steps, then z_n ≡ z_{n+P} mod 2^L. A divergent orbit never repeats, so
  2^L ≤ max(z_n, z_{n+P}). We also have log₂ z_j ≤ log₂(N+1) + o_j·log₂3 − j + o(1), where o_j
  is the number of odd steps among the first j.
- **Corollary (complexity criterion):** suppose N diverges. Let G_m = max_{j≤m}(o_j·log₂3 − j).
  Then every parity factor of length L > log₂(N+1) + G_m + c that starts at a position ≤ m − L
  occurs only once. By pigeonhole, a parity word with
  liminf over L of (L − ε·p(L)) = ∞, where ε = α·log₂3 − 1 ≤ 0.585 and p(L) is the factor
  complexity, cannot be the parity sequence of a divergent orbit.
- **Sturmian case:** every Sturmian word (p(L) = L+1) is excluded. Rational slopes give periodic
  words, which are already excluded (2-adic rationals). For irrational slopes we checked
  numerically (`tools/parity/sturmian_repetition.py`): over 9 slopes with bounded and unbounded
  partial quotients, the excess L(n,P) − growth grows linearly in M. The mechanism: at P = q_k
  (a convergent), the agreement runs have length about q_{k+1} ≥ q_k.
- **Kill tests:**
  - Perturbation barrier: ✔. The argument is pointwise.
  - 2-adic: ✔. It combines the 2-adic size with the real growth.
  - 5x+1: ✔. The same exclusion holds there and is consistent, since its divergent orbits are
    expected to have random-like parity.
- **Why it can't finish the job:** random-like words have p_m(L) = m − L + 1 once L > 2·log₂m,
  so they satisfy the criterion trivially. A real divergent orbit, if one exists, has
  random-looking parity.

### I-011: No sliding-window ranking function on binary digits (cycle 5, 2026-11-22)
- **Status:** done. A method-class no-go, backed by numerics. Not a proof idea, and probably
  folklore.
- **Test (`tools/ranking/window_lp.py`):** an LP for weights w on b-bit windows (LSB-first binary,
  zero-padded above), with W(n) = Σ windows and W(T(n)) ≤ W(n) − 1 for all 3 ≤ n < 2^Lmax.
  - Control: the terminating map n ↦ n/2 | (n+1)/2 is feasible for every b.
  - Collatz is infeasible for every b ≤ 8 already at Lmax = 12.
  - The smallest feasible window b* is 8, 10, 11 and 13 for Lmax = 8, 10, 12 and 14. So
    b* ≈ Lmax − 1: a window potential works only by memorizing whole numbers.
- **Reading:** a ranking function with fixed finite memory over plain binary does not exist,
  consistent with the subagent's √L-fluctuation argument and with YAH's need for a mixed
  binary/ternary encoding. This closes the plain-binary variant of I-003. Matrix interpretations
  over mixed bases remain untested.

### I-012: Pointwise sieve on the minimal counterexample (cycle 6, 2026-11-22)
- **Status:** dead. Known material.
- **Claim:** y* = min E must avoid forward stopping classes mod 2^k (Terras) and backward
  classes, since an odd predecessor (2w−1)/3 < w exists when w ≡ 2 mod 3.
- **Result (`tools/sieve/min_counterexample_sieve.py`):**
  - The forward survivor density among odd classes is 0.148, 0.110, 0.065 and 0.052 at
    k = 8, 12, 16 and 20.
  - The backward test removes exactly the y ≡ 2 (mod 3) third. Deeper chains add nothing,
    because the first odd predecessor is already smaller.
  - The density tends to 0, but no level is empty.
- **Kill test:** 2-adic/3-adic ✘. Every finite modulus admits surviving classes; standard.

### Loop summary (cycles 1–6, 2026-11-22)
- The stop rule fired: cycles 4–6 produced nothing new toward a proof. I-010 is largely known,
  I-011 is a method no-go, and I-012 is known.
- What the session established:
  1. **Fourier decay.** Tao's coefficient decays like e^(−0.055n), with extremal ξ ≈ 3ⁿ in
     2-power form (I-007; matches the WakauraH conjecture).
  2. **Density race.** Collatz ⟺ the exceptional density is below the basin floor, which is
     about y^(−1.35). Rate exponent 1 is the barrier, and the perturbation argument makes it
     sharp (I-008, `paper/density-race.pdf`).
  3. **Small exclusions.** Low-complexity parity sequences cannot diverge; finite-window binary
     ranking functions do not exist.
- **Honest assessment:**
  - Every density-insensitive method is blocked by the perturbation barrier.
  - Pointwise methods reduce to known partial results: cycle bounds (I-002, I-005), parity
    complexity (I-010), and sieves (I-012).
  - No candidate for a novel proof survived.

### I-013: Two-sided parity coding of orbit elements (iteration 7, 2026-11-22)
- **Status:** done. A small refinement of I-010, probably implicit in the Bernstein–Lagarias
  2-adic/3-adic conjugacy.
- **Identity (checked, 66,000 cases, 0 mismatches):** z_n mod 3^k depends only on the last k
  Syracuse valuations, since 2^s·z_n = 3^k·z_{n−k} + F(a). And z_n mod 2^F is fixed by the next
  F parities.
- **Corollary:** if a divergent orbit has the same parity window [n−B, n+F) at two positions,
  then 2^F·3^(o_B) ≤ max(z_n, z_{n'}), where o_B is the number of odd steps in the past part.
  That strengthens I-010 by a factor α·log₂3 ∈ (1, 1.585] per past symbol.
- **Why it can't finish the job:** the same reason as I-010. Random-like parity words never
  repeat windows of logarithmic length.

### I-014: Track B, toward a proof of exponential Fourier decay (iteration 7, 2026-11-22)
- **Status:** open. Reduced to one conjecture, not proved. Notes and proofs in
  `tools/fourier/proof/NOTES.md`; scripts in the same folder.
- **Proved:**
  - Lemma 1 (exact pairing majorant F_m ≥ |Φ_m|).
  - Theorem A: |Φ_n(k)| ≤ μ^t + P_k(V < t), with μ = 0.95225, where V is the number of white
    points (‖θ‖ ≥ 1/5) the paired walk visits.
  - Triangle structure lemmas, also checked exactly for m ≤ 16 and |j| ≤ 300: unique
    decomposition, white rows below each triangle, slanted exits land white, separation, and
    height ≤ m·log₂3.
  - The exit lemma (white with probability ≥ 1/4 after a bottom exit), except at the left tip.
- **Reduction:**
  - **(★)** P_k(V < εn/2) ≤ C·e^(−cn) uniformly in k. This implies max|Φ_n| ≤ ρⁿ; with ε = 0.05,
    ρ = 0.99878.
  - Numerically, about 60% of visited points are white. The window m ≤ 600 contains 96,000
    triangles with heights following a 2^(−H) law (tallest 15.7). The majorant F decays like
    e^(−0.035n).
- **Obstacle:**
  - Triangles of height ≍ m that chain together would give only O(log n) exits.
  - Ruling them out needs 2^d·r₂ ≢ r₁ (mod 3^(m₂)) for |r_i| up to 2^(ηm): a 3-adic Diophantine
    statement beyond Baker/Yu. This looks like WakauraH's "×2×3 rigidity hypothesis"
    (bibliography 9.11).
  - A uniform sup-norm block contraction (F1) is false.
  - "Every L-block after the first white point has a white point" (F2) is false globally, since
    there are tall triangles at every large level. Only an averaged statement can hold.
- **Kill tests:** this is not a Collatz proof route (3x−1 has the same Φ). Its value is
  quantitative, feeding I-006 and I-008.

### I-015: Track A, reproducing Yolcu–Aaronson–Heule (iteration 7, 2026-11-22)
- **Status:** reproduction done. Nothing open proved. Code and notes in `tools/rewriting/`
  (NOTES.md).
- **Built:** all systems from arXiv 2105.14697v3, with arithmetic semantics checked by
  simulation (n ≤ 3000, and the 5x+1 loops from 13 and 17). Direct CNF encoding (CaDiCaL/Kissat)
  agrees with z3 on 144/144 queries. Every proof is re-checked with exact integers
  (`check.py`). Tampering control: 15 of 20 random edits rejected; the other 5 left the proof
  valid.
- **Reproduced:** all 11 one-rule-removed subsystems (Table 1), W′, S3, and 2 of the 4
  Thm 4.9 cases. Not reproduced within budget: Farkas F, 2 Thm 4.9 cases, and the mod-8 cases.
- **Controls:**
  - Full T: UNSAT for removing any rule at natural d ≤ 3 and arctic d ≤ 3.
  - The 5x+1 analogues are unprovable, including the false 5x+1 version of Conj 4.14.
- **Push:** S1, S2, the Conj 4.14 subsystem, and the new accelerated system **H**
  (ff⊿→0⊿, ttt⊿→22⊿, plus X) are UNSAT up to d ≈ 3 and time out at d = 4.
- **Limits:** UNSAT covers only the shapes searched, and timeouts mean nothing. The search used
  about 1 CPU-hour, against YAH's thousands, and the prover is about 10× slower than theirs. The
  arctic top-termination strictness argument is unformalized.
- **Next:** a long background run on H (`logs/H_long.txt`).

### I-016: Track C, no m-cycles for m ≤ 135 (computed; iteration 7, 2026-11-22)
- **Status:** computed, not proved. m ≤ 134 is cross-checked by two solvers; the m = 135 cross-check is running.
  Code in `tools/xcheck/b5/` (NOTES.md).
- **Method:** a five-rise bound B₅(N), computed exactly by meet-in-the-middle (`b5mitm.c`, exact
  mod 2^512), added as rule 5 to our independent Marshall-style chain (`mchain.py`).
  - Sensitivity analysis first: m = 134 needs only B₅(73) ≤ 92 and B₅(121) ≤ 146. m = 135 also
    needs B₅(199) ≤ 230 plus cheap points. m = 136 needs points near N ≈ 320, at about 4–9 h each.
- **Values:**
  - Needed points: B₅(73) = 91, B₅(121) = 143, B₅(199) = 225, with margins +1, +3 and +5.
  - Also: 74 → 91, 79 → 101, 114 → 137, 116 → 139, 119 → 141, 120 → 143, 122 → 145,
    124 → 146, 128 → 152, 130 → 152, 132 → 158.
- **Chain result:** every shape for m = 92–135 is closed (worst ratio 0.92365), in both own and
  compat modes. **m = 136 and 137 stay open** (negative controls).
- **Validation:**
  - The solver matches brute force and `rises.py` for N = 6–36.
  - A thread-count check gives identical results.
  - The slow, differently structured `b5mitm_v0.c` gives identical hit lists at N = 44, 50 and 73.
  - Every maximizer is re-simulated as an integer.
  - **N = 121 cross-checked:** v0 on 8 workers gives the same 4 maximizers (total 143) as
    `b5mitm.c` (`logs/v0_hits_N121_T143.txt`). So m = 134 rests on two independent solvers.
  - N = 199, needed only for m = 135: a v0 check of 'no total ≥ 226' is running, about 13 h on
    6 workers.
- **Relies on:** Barina's 2^71 verification, Rhin's bound (unchecked), Marshall's B₂–B₄ tables
  (spot-checked; see `tools/xcheck/README.md`), and the rule-5 chain lemma (rules 2–4 composed).
- **Significance:** this goes two steps past the best published computed bound (m ≤ 133,
  Marshall) and past the Lean-checked m ≤ 111. It is still finite: the cycle side can never be
  closed this way.
