"""Adaptive-grid idea: would exact B_r at the chain's own (off-grid) queries help? HYPOTHETICAL dense
tables: for rules in RULES, fill every integer N between Marshall grid points by linear interpolation of
the table (rounded down = optimistic, or up = pessimistic), keep the real B5 table. An estimate of what
a dense computation might give, not a bound. Usage: adaptive.py m down|up RULES [N:v extra B5 ...]"""
import io, contextlib, json, sys
from math import floor, ceil
from b5lib import *
import rises
TB = rises.load_marshall(ROOT)
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
for a in sys.argv[4:]:
    k, v = map(int, a.split(":")); real[k] = v   # hypothetical extra B5 values
m, mode, rules = int(sys.argv[1]), sys.argv[2], [int(c) for c in sys.argv[3]]
import os
NCUT = int(os.environ.get("NCUT", "10000000"))   # densify only N <= NCUT
ch = checker(real)
for r in rules:
    g = TB[r]; Ns = sorted(g); d = dict(g)
    for a, b in zip(Ns, Ns[1:]):
        if a >= NCUT: break
        for N in range(a + 1, min(b, NCUT + 1)):
            x = g[a] + (g[b] - g[a]) * (N - a) / (b - a)
            d[N] = max(g[a], floor(x) if mode == "down" else ceil(x))
    ch.tab[r] = M.Table(d)
with contextlib.redirect_stdout(io.StringIO()):
    fails, worst = M.check_m(ch, m)
print(f"m={m} dense-interpolated rules {rules} ({mode}) B5 extra {sys.argv[4:]} NCUT={NCUT}: {'closed' if not fails else 'open'} (worst {worst[0]:.4f})")
