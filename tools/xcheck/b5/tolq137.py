"""m=137 with rule-5 points placed at the chain's own queries (not Marshall's grid). HYPOTHETICAL values.
For B5(N1) = v1 (N1 near 315), report the next query ceil(Y) in the 500s and the largest v2 at that N2
that closes m (scan v2 down from b4-ish). Usage: tolq137.py m N1 N2a,N2b,... v1 [v1 ...]"""
import io, contextlib, json, sys
from math import ceil
from b5lib import *
import rises
B4 = rises.load_marshall(ROOT)[4]
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
import os
for a in os.environ.get("EXTRA", "").split():
    k, v = map(int, a.split(":")); real[k] = v   # e.g. EXTRA='532:557' (a proven lower bound used as value: optimistic)
m, N1 = int(sys.argv[1]), int(sys.argv[2])
N2S = [int(x) for x in sys.argv[3].split(',')]
def closes(t):
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t), m)
    return not fails, worst[0]
for v1 in map(int, sys.argv[4:]):
    t = dict(real); t[N1] = v1
    ok, w = closes(t)
    _, Y, _, _ = hard_trace(checker(t), m)
    qs = [ceil(y) for y in Y if 400 < float(y) < 700]
    if ok:
        print(f"m={m} B5({N1})={v1}: closes ALONE (worst {w:.4f})"); continue
    for N2 in N2S:
        best = None
        lo, hi = N2, N2 + 70          # binary search: closing is monotone (smaller v2 is a stronger rule)
        t2 = dict(t); t2[N2] = lo
        if not closes(t2)[0]:
            print(f"m={m} B5({N1})={v1}: N2={N2} never closes", flush=True); continue
        while lo < hi:
            mid = (lo + hi + 1) // 2; t2 = dict(t); t2[N2] = mid
            if closes(t2)[0]: lo = mid
            else: hi = mid - 1
        print(f"m={m} B5({N1})={v1}: open alone; largest closing B5({N2}) = {lo}  (b4 table {B4.get(N2)})", flush=True)
