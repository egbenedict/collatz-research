# Track C: a five-rise bound B₅ for the m-cycle chain (2026-11-22)

Everything here lives in `tools/xcheck/b5/`. No existing file was edited. `mchain.py`, `rises.py`,
`rigorous.py`, `b4mitm` are imported/run unchanged; `b5mitm.c` started as a copy of `b4mitm.c`.

## Definitions (same as the existing code)

- Rises: runs of o-steps of T(n) = n/2, (3n+1)/2. B₅(N) = max over odd n < 2^N of k₀+…+k₄
  (as in `rises.py`, `Solver.B(5)`).
- Chain rule r = 5: Y_t ≤ Y_{t−5} + (δ−1)·b₅(⌈Y_{t−5}⌉) − 5 + 5ε, where b₅(N) = min of the table over
  grid points ≥ N. Same derivation as rules 2–4 (block lemma composed five times; the five rises after
  minimum n_{t−5} are its first five rises). It plugs in by adding `M.Table(b5)` to `Checker.tab`
  (`b5lib.checker`); `chain()` already loops over every rule in `tab`.

## Task 1: which B₅ values are needed (before computing)

Hard-shape chain at m = 134 (`trace134.py`): a saw-tooth. In each block, rules 2–4 from the first
value of the block keep the next three values flat; then rule 1 forces a jump ×δ. After Y ≈ 9·10⁶
(t ≈ 73) only rule 1 is left, so every step saved early divides the final sum by about δ = 1.585.
The m = 134 ratio is 1.464, so **one saved step closes m = 134**. Each m needs one more step.

The B₄-derived bound on B₅ is very weak: the chain from Y₀ = N implies B₅(73) ≤ 206, B₅(121) ≤ 335,
B₅(199) ≤ 535. The true values are about B₄ + 1…8. So a B₅ rule makes each block one value longer.

Hypothetical table B₅ := B₄ + c on B₄'s grid N ≤ Nmax (`sens5.py`, `sens5_nmax.py`; c = 1 is the most
optimistic possible, since B₅ ≥ B₄ + 1):

| m | smallest Nmax that closes the hard shape, c = 1 / 3 / 5 / 8 |
|---|---|
| 134 | 73 / 120 / 121 / 201 |
| 135 | 119 / 198 / 201 / 328 |
| 136 | 194 / 321 / 328 / 852 |

With B₅ = B₄ + 5 on all B₄ grid points ≤ 852, the hard-shape ratio stays below 1 up to m = 138.

## Task 2: solver and cost

`b5mitm.c`: meet in the middle, ν = A(k₀,l₀,k₁,l₁,k₂,l₂) + 2^{L₃}·R(K₃,l₃). The joining key is the
w = T+4−N bit window (carry 0/1), and every match is tested exactly mod 2^512. Bounds: K₃ ≤ KMAX ≥ B₄(N);
l₃ < N + (δ−1)K₃; L₃ < N + 0.585·KMAX + 4. That last one is a hard assertion: if a feasible left goes
past it, the run aborts with a non-zero exit. Right-side tables are precomputed for every shift L₃.
For fixed (k₀…l₁, k₂), V = A₁ + 2^{L₁}u(K₂) fixes the bits that decide feasibility, so the l₂ loop
bound is exact. The work is shared by 8 pthreads.

`b5mitm_v0.c` is the first version: L₃ in the outer loop, as in b4mitm. It is about 40× slower and
kept only as a second implementation for the cross-checks.

The cost is set by the left side: about 1.05·C(N,6) prefixes. At 8 threads on this machine
(4P + 4E cores) the wall time per prefix is 7 ns at N = 60, 13 ns at N = 73 and about 20 ns at N = 121–132.

| N | left prefixes | wall time (8 threads) |
|---|---|---|
| 73 | 1.9·10⁸ | 2.5 s |
| 121 | 4.1·10⁹ | 80 s |
| 132 | 6.9·10⁹ | 145 s |
| 199 | 8.2·10¹⁰ | 885 s (10.7 ns/prefix) |
| 220 (est.) | ≈1.6·10¹¹ | ≈ 30–55 min |
| 320 (est.) | ≈1.5·10¹² | ≈ 4–9 h per point |

So N ≤ ~210 fits the budget and N ≈ 320 (needed for m = 136) does not. Total compute this session:
about 1 h (B₅ runs ≈ 27 min, the killed v0 cross-check ≈ 13 min, chain/sensitivity runs ≈ 20 min, single core).

## Validation

- `test_b5mitm.py` (`logs/test_b5mitm.log`): for N = 6…36 (even), the solver's maximum at T = B₅ equals
  `rises.py` B₅ (brute force for N ≤ 18), and nothing is found at T = B₅+1. There are no UNRESOLVED
  matches, the 4-thread hit count equals the 1-thread count, and every witness is re-simulated as an
  actual odd n < 2^N (`rises.realise`).
- v0 and v1 give identical hit lists at (N,T) = (44,54), (50,62) and (73,91) (`logs/v*_hits_N73_T91.txt`).
- Every computed value: KMAX = b₄(N) from Marshall's table, first confirmed by `../b4mitm` (no n < 2^N with
  a four-rise total ≥ b₄(N)+1, using KMAX = b₃(N); B₃ is independently exact for N ≤ 400). The maximiser
  is re-simulated as an actual integer.
- Not done: a v0 run at N = 121 (killed after 11 min, too slow). There is no independent solver
  above N = 73.

## Computed B₅ values (exact; `b5_values.json`, `logs/compute.log`)

| N | 73 | 74 | 79 | 114 | 116 | 119 | 120 | 121 | 122 | 124 | 128 | 130 | 132 | 199 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| B₅ | 91 | 91 | 101 | 137 | 139 | 141 | 143 | 143 | 145 | 146 | 152 | 152 | 158 | 225 |
| b₄(N) | 87 | 89 | 93 | 135 | 135 | 136 | 138 | 139 | 144 | 144 | 147 | 147 | 156 | 221 |

## Which points matter (computed-value sensitivity)

- m = 134 needs only B₅(73) and B₅(121). It closes when B₅(73) ≤ 92 and B₅(121) ≤ 146, so the exact values
  91 and 143 have margins of +1 and +3 (`tol134.py`). Dropping either opens it. The other ≤ 132 points
  are not needed for m = 134 (`sens_b5_points.py`).
- m = 135 also needs one point near 199 (`minset135.py`, `single135.py`; B₄+c is hypothetical at the
  points not yet computed). N ≤ 197 does not work even with B₅ = B₄+1. N = 199 or 201 works when
  B₅ ≤ 230. **Computed B₅(199) = 225, a margin of +5.** The cheap points 74, 116 and 122–132 are what
  close the other two hard shapes, K = 7.77·10¹⁹ and 2.06·10²⁰.

## Task 3: chain results with the computed table (`final_runs.sh`, logs/)

The table is {73: 91, 74: 91, 79: 101, 114: 137, 116: 139, 119: 141, 120: 143, 121: 143, 122: 145, 124: 146,
128: 152, 130: 152, 132: 158, 199: 225}.

- **Own mode, m = 92…135: all closed** (`logs/sweep_own_92_136.log`). At m = 135 there are 3 dangerous
  candidates, all closed. The worst ratio is **0.9236498**, at K = 27,444,133,206,411,171,953, the same
  hard shape as at m = 133, now one step later.
- Compat mode (rule 1 = δY): m = 134 and 135 closed, worst ratio 0.9236500 (`logs/chain_compat_134_136.log`).
- Negative controls: m = 136 is open in both modes (2 shapes; hard shape 1.4640, K = 2.06·10²⁰ at 1.2243).
  m = 137 is open in own mode (4 shapes; hard shape 3.678).
- The ratio still grows about ×1.585 per m, so m = 136 needs points near N ≈ 320 (`sens5_nmax.py`).

**Status: computed only, not proved.** "No m-cycles for m ≤ 135" depends on:
- E1: Barina's verification to 2^71;
- E3: Rhin's bound (Simons–de Weger Lemma 12);
- Marshall's B₂, B₃, B₄ tables, which are independently checked only in the ranges listed in `../README.md`;
- b₄(N) at the 14 B₅ points (re-confirmed here with b4mitm);
- the new B₅ values (one implementation, validated only up to N = 36 against rises.py, and to N = 73
  against v0);
- the rule-5 extension of the chain lemma.

Adding a rule only lowers the chain, so the existing m ≤ 133 results are unaffected.

## Files

`b5mitm.c` (solver), `b5mitm_v0.c` (slow second implementation), `test_b5mitm.py`, `compute_b5.py` (driver
→ `b5_values.json`), `b5lib.py` (Checker + rule 5), `run_chain.py`, `final_runs.sh`, and the sensitivity
scripts `sens5*.py`, `plan135.py`, `minset135.py`, `single135.py`, `tol134.py`, `sens_b5_points.py`,
`trace134.py`, `binding5.py`, `small_b5.py`.

Reproduce: `cc -O3 -o b5mitm b5mitm.c && python3 test_b5mitm.py 36 && python3 compute_b5.py 73 121 199:218 && ./final_runs.sh`.
