"""Chain with fall-weighted constants where computed (logs/weighted_r3.log etc.): rule-r table value at N
replaced by b_eff = (C_r(N) + r - 1)/Dh1 (C_r an exact rational upper bound); other points unchanged.
Usage: wchain.py m [EXTRA5 N:v ...]   (reads C values recomputed from weighted.py output files)"""
import io, contextlib, json, sys, re
from fractions import Fraction as F
from b5lib import *
import rises
TB = rises.load_marshall(ROOT)
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
for a in sys.argv[2:]:
    k, v = map(int, a.split(":")); real[k] = v
ch = checker(real); Dh1 = ch.Dh1
def load(path, r):
    out = {}
    for line in open(path):
        mm = re.match(rf"r={r} N=(\d+) .*C{r}<?=([\d.]+)", line)
        if mm: out[int(mm.group(1))] = mm.group(2)
    return out
import weighted_vals as WV   # exact rationals recomputed
for r, vals in WV.load().items():
    d = dict(TB[r]) if r < 5 else dict(real)
    for N, C in vals.items(): d[N] = min(d.get(N, 10**9), (C + r - 1) / Dh1)
    ch.tab[r] = M.Table(d)
m = int(sys.argv[1])
with contextlib.redirect_stdout(io.StringIO()):
    f, w = M.check_m(ch, m)
print(f"m={m} weighted tables {sorted((r, len(v)) for r, v in WV.load().items())} extra5 {sys.argv[2:]}: {'closed' if not f else 'open'} (worst {w[0]:.4f})")
