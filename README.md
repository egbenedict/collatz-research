# Collatz Project

Goal: attack the Collatz (3x+1) conjecture seriously. That means understanding the
literature, building tools to test ideas, and pursuing attack lines that can produce
real partial results even if a full proof stays out of reach.

**The conjecture.** Let `C(n) = n/2` if n is even and `3n+1` if n is odd. For every
positive integer n, iterating `C` eventually reaches 1.

It splits into two independent open claims:

- **(A) No nontrivial cycles.** The only cycle on positive integers is 1 → 4 → 2 → 1.
- **(B) No divergent orbits.** No orbit grows forever.

Most work here uses the accelerated map `T(n) = n/2` (even), `(3n+1)/2` (odd), because
its parity vectors are much cleaner.

## Layout

| Path | What's there |
|---|---|
| `sources/bibliography.md` | Annotated bibliography by theme, each entry flagged verified or unverified |
| `tools/collatz.py` | Toolkit: orbits, stopping times, parity vectors, cycle search, control maps |
| `tools/cycle_bounds.py` | Hercher's cycle-bound machinery, re-runnable with any verification bound (I-002) |
| `tools/kstar_sieve.py` | Rules out an m-cycle with a specific K (budget identity + T-window + prefix sieve) |
| `tools/run_kstar_all.py` | Runs the sieve on every surviving (m, K) |
| `tools/test_kstar_sieve.py` | Brute-force validation of the sieve on toy cases |
| `tools/test_collatz.py` | Known facts the toolkit must reproduce (run before trusting any experiment) |
| `tools/fourier/` | Exact Syracuse law and characteristic function (`syrac_decay.py`), banded transfer recursion to n ≈ 1000 (`band_recursion.py`), Monte Carlo decomposition (I-007) |
| `tools/density/predcount.c` | Predecessor-basin counts #(P(z) ∩ [1, X]) by inverse-tree DFS (I-008) |
| `paper/density-race.tex` | The zero-or-gap lemma and the exponent-one barrier (I-008) |
| `ideas/log.md` | Every idea, its kill-test results, and why it failed or survived |

## Kill tests (run every idea through these first)

1. **The 5x+1 test.** If the argument also proves 5x+1 converges, it's wrong. 5x+1 has
   cycles {1,3,8,4,2}, {13,33,83,…}, and {17,43,…} (see `python3 tools/collatz.py controls`).
2. **The negative-integers test.** If the argument rules out cycles for 3x+1 on negative
   integers, it's wrong. Those cycles exist: {−1}, {−5,−7,−10}, and {−17,…} with 11 elements.
3. **The 2-adic test.** If the argument only uses residues mod 2ᵏ, it's dead. Every parity
   pattern of length k occurs for exactly one residue mod 2ᵏ
   (`python3 tools/collatz.py bijection 12` checks this). A proof has to use the
   ordinary size of n together with its 2-adic behavior.
4. **The Conway test.** If the argument would work for any map of the same shape (affine
   on residue classes), it's wrong. That class is undecidable (Conway 1972; Kurtz–Simon
   2007), so a proof has to use something specific about 3 and 2.
5. **The heuristic test.** If the argument assumes parity steps behave like independent
   coin flips, that's the heuristic, not a proof.

## Phases

1. **Literature** (`sources/`): read and annotate. Priority: Lagarias's 1985 survey,
   Tao 2019, Simons–de Weger 2005 and Hercher 2023, Yolcu–Aaronson–Heule 2021.
2. **Tooling** (`tools/`): every claim gets a numerical check, always with the 5x+1 and
   negative-integer controls.
3. **Attack lines**, ranked by how tractable they are:
   1. **Cycle bounds.** Push "no m-cycles" past m ≤ 91 (Hercher), using Baker-type bounds
      on |2ᴷ − 3ᵐ| plus computation.
   2. **Automated termination.** Extend the Yolcu–Aaronson–Heule rewriting-system + SAT
      approach with better encodings, new interpretation classes, and larger searches.
   3. **Counting.** Raise the Krasikov–Lagarias exponent: at least x^0.84 of the numbers
      up to x are known to reach 1.
   4. **Extending Tao.** Go from logarithmic density to natural density, or make the
      result quantitative.
   5. **Moonshot.** Find a new invariant that ties a number's ordinary size to its 2-adic
      behavior and has to decrease.
4. **Proof hygiene.** Formalize key lemmas in Lean, run every lemma through the kill
   tests, and log every failure in `ideas/log.md`.

## Current status

- Verified computationally for all n < 2⁷¹ (Barina 2025).
- No nontrivial cycle has m ≤ 91 local minima (Hercher 2023). Our re-run shows 2⁷¹
  alone doesn't raise this; m ≤ 92 needs about 2^73.9 and m ≤ 101 needs 2⁷⁷ (see I-002).
- Every nontrivial cycle has more than 1.375·10¹¹ odd elements (Hercher Cor. 29 + Barina 2⁷¹).
- **The cycle frontier has moved past us.** Kernel-checked Lean proofs (unrefereed) now
  give no m-cycles for m ≤ 101 (Alick) and m ≤ 111 (Marshall); m ≤ 133 is computed.
  Our independent m ≤ 92 result (I-002 round 3, `paper/`) uses essentially the same idea.
  It is a valid cross-check, but not new.
- Almost all orbits (logarithmic density) attain almost bounded values (Tao 2019/2022;
  notes in `sources/notes/tao2019.md`). In 2026 there are unrefereed, AI-generated claims of:
  - a natural-density version;
  - a positive density of integers reaching 1 (Lean-checked, constant 2^(−2^(2^163)));
  - a positive predecessor density for every y with 3 ∤ y. This would make Collatz
    equivalent to "density exactly 1". See bibliography 9.7–9.10 and I-006.
- **Our computed extension: no m-cycles for m ≤ 137** (I-016, I-017). It adds exact five-rise
  bounds B₅ to an independent Marshall-style chain. For m = 137 it also uses a new fall-weighted
  block rule, C_r(N) = max[(δ−1)Σk − Σl]. That rule was needed, because B₅(532) = 563 is too large
  for the standard rule. The values up to N = 199 are confirmed by two solvers. B₅(326), B₅(327),
  B₅(328) and the exhaustive N = 532 list (all five-rise patterns ≥ 557, 256 shards, about 600
  runner-hours) come from one solver plus positive controls. It still relies on 2⁷¹, Rhin's bound,
  and Marshall's B₂–B₄ tables. m = 138 would need B₅ near N ≈ 850, about 16× the N = 532 cost,
  which is out of reach.
- Fourier decay of Tao's Syracuse variable: computed exactly to n = 17 and to n = 1200 with a
  validated banded recursion. It decays exponentially, but at the tiny rate I₀ = I(log₂3) ≈ 0.055,
  with extremal ξ = 2^k and 2^k ≈ 3ⁿ (I-007). This agrees with a 2026 conjecture by WakauraH.
- The probabilistic line, fully mapped (I-008): Collatz ⟺ the exceptional density
  f(y) = d{Col_min ≥ y} falls below the predecessor-basin floor, which is about y^(−1−o(1)).
  Tao-type methods give (log y)^(−c). Rate exponent 1 is exactly the barrier.
- Both (A) and (B) are open.
