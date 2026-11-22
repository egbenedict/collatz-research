# Independent cross-check of "no Collatz m-cycles for m ≤ 133"

This checks the computed (not Lean-verified) result in
[LeeMarshall1113/collatz-m-cycles](https://github.com/LeeMarshall1113/collatz-m-cycles)
with code written separately. The only shared input is the published rise-bound tables in
`marshall_data/`, copied unmodified under Apache-2.0 (see `marshall_data/SOURCE.md`).
They are spot-checked here. The derivations are in the docstring of `mchain.py`.

## Results (2026-11-21)

| Part | Result |
|---|---|
| Chain + all shapes, m = 92–133 | **All closed.** Worst ratio 0.9236463 at K = 27,444,133,206,411,171,953, m = 133 |
| Agreement with Marshall's referee chain | Equal to 7 digits on his three published hard shapes |
| B₂ table, 88 points N ≤ 10,000 | 74 exact; 14 are 1–2 too high (safe; the ones Marshall discloses) |
| B₃ table, 39 points N ≤ 400 | All exact |
| B₄ table, 7 points N ≤ 80, plus N = 119 and 317 | All exact |
| Entries where +1 alone breaks m = 133 | **All three checked exact**: B₄(73) = 87, B₄(119) = 136, B₄(317) = 337 |

Not checked: the other table entries above the ranges listed (each tolerates +1 on its own),
and Rhin's bound (an external input).

## Files

| File | What it does |
|---|---|
| `rigorous.py` | Exact rational interval arithmetic: ln 2, ln 3, log₂3, ln, exp, 2^x with tail bounds. `python3 rigorous.py` self-tests |
| `b4mitm.c` | Meet-in-the-middle B₄ solver in C (exact mod 2^512). `cc -O2 -o b4mitm b4mitm.c; ./b4mitm N T KMAX` lists every odd n < 2^N whose first four rises total ≥ T |
| `test_b4mitm.py` | Validates `b4mitm` against `rises.py` (N = 14–40, 72, 80) with witness re-simulation |
| `rises.py` | Own B_r(N) solver. `selftest` checks it against brute force; `compare --r R --nmax N` checks it against the tables |
| `mchain.py` | Own n₀ bound, chain and exhaustive search for dangerous shapes. `--mmin/--mmax`, `--mode own/compat`, `--hard` |
| `test_mchain.py` | Controls: lattice vs brute force, negative control m = 134/135, widened thresholds, semiconvergents |
| `sensitivity.py`, `sensitivity_points.py` | How large an error in the unchecked table entries would break m = 133 |
| `logs/` | Outputs of the runs above |

## Reproduce

```bash
python3 rigorous.py
python3 rises.py selftest
python3 test_mchain.py
python3 mchain.py --mmin 92 --mmax 133
python3 sensitivity_points.py
```
