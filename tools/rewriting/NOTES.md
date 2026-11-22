# Track A: rewriting-system termination (YAH 2021/2023 reproduction + push)

Source: Yolcu, Aaronson, Heule, *An Automated Approach to the Collatz Conjecture*, arXiv 2105.14697v3
(JAR 2023). All systems below were transcribed from the v3 PDF (Sections 3.2, 4.1-4.5).

## Files

| File | What |
|---|---|
| `srs.py` | All rewriting systems (ASCII: `<`=◁, `>`=⊿, `*`=•, `f,t`=binary 0/1, `0,1,2`=ternary) |
| `simulate.py` | Control: rules have the right arithmetic meaning; random-strategy rewriting reproduces T, F, S, C, H trajectories for n ≤ N; 5x+1 systems really loop (13-cycle, 17-cycle) |
| `matrix.py` | z3 bit-vector encoding (reference) + **exact verifier** (`verify`) |
| `cnf.py` | Direct CNF encoding (binary circuits for natural, order encoding for arctic) + SAT portfolio (CaDiCaL 1.9.5, Kissat, Glucose, MiniSat; pysat) |
| `crosscheck.py` | z3 vs CNF encodings agree on sat/unsat (144/144 small queries) |
| `prover.py` | Stepwise rule-removal prover; writes `proofs/NAME.json` |
| `check.py` | Independent proof checker: bookkeeping + exact matrix check + direct numeric evaluation in random contexts |
| `run_batch.py` | Parallel driver; logs in `logs/` |

Setup: `python3 -m venv .venv && .venv/bin/pip install z3-solver python-sat pypdf`.

## Proof scheme (what a "PROVED" means)

For a system R = D ∪ X' (D = dynamic rules ending in `>`/`*`, Bf = rules starting with `<`):

1. ob1: SN(R \ Bf), ob2: SN(R \ D), each by stepwise rule removal with extended monotone
   algebras (natural: M₁₁ ≥ 1; arctic: v = −∞, M₁₁ finite). By type introduction
   (YAH Lemma 3.16, generalised: the typing argument needs exactly ob1 and ob2), only canonical
   strings `<w>` can be nonterminating.
2. ob3: SN(Bf_top / R) (forward) or SN(rev(D)_top / rev(R)) (reverse) by stepwise removal with
   weakly monotone algebras (YAH Lemma 3.18 / Thm 2.15). On canonical strings these rules only
   fire at the end, so this gives SN(Bf/R) (resp. SN(D/R)); Thm 2.6 + ob1 (resp. ob2) gives SN(R).

Orders. Natural: `[l] > [r]` iff M_l ≥ M_r entrywise, (v_l)₁ > (v_r)₁, (v_l)ᵢ ≥ (v_r)ᵢ (YAH (2)).
Arctic ext mode: KW09 entrywise ≫. Arctic top mode: I use a slightly relaxed strictness —
row 1 of M and v₁ must be ≫, all other entries only ≥. Soundness: on the domain N × A^{d−1}
(closed because every symbol has M₁₁ or v₁ finite) the order "x₁ > y₁ and xᵢ ≥ yᵢ" is
well-founded, arctic affine maps are monotone for ≥, and if row 1 is ≫ then the first
component of [l](x) strictly exceeds that of [r](x) (the right side's first component is a
finite max attained by some finite term, and the matching left term is strictly larger).
Weak steps anywhere are fine by ≥-monotonicity. Natural entries are in [0, 2^w); arctic entries
in {−∞, 0, …, w} (so arctic `w` = max entry; the paper's V counts distinct values).

Every SAT model is re-verified with exact integer arithmetic before it is used, and every saved
proof is re-checked by `check.py`, so a CNF bug can cost completeness but not soundness.
Tamper test: 15/20 random single-entry corruptions of a saved proof are rejected (the other 5
remain valid interpretations).

## Results (2026-11-22, one laptop, 8 cores, ~1.5 h total)

### Controls
- `simulate.py 3000`: every X rule preserves value, every dynamic rule applies its map; random
  rewriting reproduces T, F, S, C(mod-8 systems), H trajectories for all n ≤ 3000. PASSED.
- 5x+1: `T5` (binary/quinary analogue of T) loops from 13; `ODD8_5-tft>-ttt>` (5x+1 analogue of
  Conj 4.14, false because the 5x+1 cycle 17→43→108→54→27→68→34 avoids 5, 7 mod 8) loops from 17.
  Prover on both: **no proof**; all queries with d ≤ 3 (nat w ≤ 1, arctic max ≤ 2) are genuine
  UNSAT in both directions, larger ones time out (60 s). (Any proof here would expose a bug.)
- Full Collatz T: no first-step interpretation exists (UNSAT, both directions) for nat d ≤ 3, w = 1;
  nat d ≤ 2, w = 2; arctic d ≤ 3 max ≤ 2 (d = 2 max 3). nat d=3 w=2, arctic d=3 max 3, arctic
  d=4 max 2: timeout at 180 s each, both directions. (`logs/smallT.log`, `logs/pushT*_T.log`)

### Reproduced (proof JSONs in `proofs/`, all pass `check.py`)
| System | YAH | ours (ob3 step sizes) |
|---|---|---|
| T \ {each of the 11 rules} (Table 1) | all 11 proved | **all 11 proved** (largest: T\{t1→2f} arctic d=5 max 2, 22 min; T\{f0→0f} arctic d=5 max 2, 12 min; others < 2 min) |
| W′ (Thm 4.1) | proved | proved (YAH: nat d=2; ours: fwd arctic d=3) |
| S3 = S \ {t•→2•} (Prop 4.6) | proved | proved, nat d ≤ 2, 0.5 s |
| C^Even/4 \ {tf⊿→t⊿}, T^Odd/4 \ {tt⊿→2t⊿} (Thm 4.9, 2 of 4) | proved | proved |
| C^Even/4 \ {t⊿→2f⊿}, T^Odd/4 \ {f⊿→⊿} | trivial | proved |

Not reproduced within budget (killed or timed out, not refuted): C^Even/4\{ff⊿}, T^Odd/4\{ft⊿},
mod-8 cases 2, 3 (EVEN8-ftf>, ODD8-ftt>), {1,5}/{1,7} variants (never started), and the Farkas
system F (arctic d=5 max 4 reverse: timeout at 300 s with a 2-solver portfolio; YAH needed
negative branching + 8 parallel shuffles to get 5-90 s).

### Push on open systems (nothing proved)
| System | meaning | outcome |
|---|---|---|
| T | Collatz | see controls: UNSAT at small sizes, timeouts at nat 3/w2, arctic 3/3, arctic 4/2 |
| S1, S2 (Prop 4.5) | S-trajectories must hit 1 / 5 mod 8 | fwd: UNSAT up to nat d=3 w=1, arctic d=3 max 2; timeouts up to arctic d=4 max 3 (90 s); rev: UNSAT through arctic d=3 max 2 |
| ODD8 \ {tft⊿, ttt⊿} (Conj 4.14) | must hit 5 or 7 mod 8 | UNSAT ≤ nat d=3 w=1 / arctic d=2 max 3; timeouts to arctic d=4 max 3 |
| **H** (new encoding) | Conj 4.14 via YAH's partial map H(n)=3n/4 (n≡0 mod 4), (9n+1)/8 (n≡7 mod 8): rules ff⊿→0⊿, ttt⊿→22⊿ plus X | UNSAT ≤ nat d=3 w=1 / arctic d=3 max 2; timeouts nat d=3 w=2 to arctic d=4 max 3 |
| G, single residues 0,1,5,7 mod 8 | — | defined in `srs.py`, not run (budget) |

H is, as far as the paper says, a new SRS formulation of Conj 4.14 (the paper states H as a
function only). It has 2 dynamic rules instead of 4 and is accelerated (each H step = several
Collatz steps), so its derivational complexity is lower, which is the paper's own heuristic for
provability. It deserves a long run (hours, 8-way portfolio) — the cheapest real shot here.

## Honest limits
- Everything new is negative and bounded: UNSAT results only exclude the specific interpretation
  shapes (dimension, entry range, my order variants, one-rule-at-least first step); timeouts mean
  nothing. YAH spent thousands of CPU-hours on these open cases; we spent ~1 CPU-hour.
- My CNF + CaDiCaL is ~10x slower than YAH's prover on the hardest steps (they report seconds for
  T\{t1→2f}; we needed 22 min of sweep), so the push reached only d ≤ 4.
- A proof here is only as good as the reduction lemmas (type introduction, YAH Lemma 3.18, rule
  removal), which I re-derived for general D ∪ X systems but did not formalize.

## Rerun
```
cd tools/rewriting
.venv/bin/python simulate.py 3000
.venv/bin/python crosscheck.py
CFG=nat:1:1,arc:1:1,nat:1:2,arc:1:2,nat:2:1,arc:2:1,nat:2:2,arc:2:2,arc:2:3,nat:3:1,arc:3:1,arc:3:2,nat:3:2,arc:3:3,nat:4:1,arc:4:1,arc:4:2,nat:4:2,nat:5:1,arc:5:2
PAR=6 .venv/bin/python run_batch.py table1 "T-f> T-t> T-f0 T-f1 T-f2 T-t0 T-t1 T-t2 T-<0 T-<1 T-<2 Wp" -- --configs $CFG --timeout 60
.venv/bin/python prover.py S3 --dmax 2 --wmax 2
.venv/bin/python check.py proofs/*.json
# open problems / kill tests
.venv/bin/python prover.py T --first-step-only --configs nat:3:2,arc:3:3,arc:4:2 --timeout 180
.venv/bin/python prover.py H --configs nat:3:2,arc:3:3,arc:4:2,arc:4:3 --timeout 600 --jobs 4
.venv/bin/python prover.py T5 --configs nat:2:2,arc:3:2 --timeout 60
```
