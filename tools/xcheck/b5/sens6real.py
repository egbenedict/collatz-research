"""Rule-6 sensitivity on the computed-B5 points only. HYPOTHETICAL B6(N) := B5(N) + c (B6 >= B5 + 1, so
c = 1 is the most optimistic) at points of b5_values.json. Real B5 table for rule 5. Reports whether m
closes with all those points and a greedy minimal set (drops largest N first: B6 cost ~ C(N,8)).
Usage: sens6real.py m c [c ...]"""
import io, contextlib, json, sys
from b5lib import *
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
m = int(sys.argv[1])
def closes(b6):
    ch = checker(real); ch.tab[6] = M.Table(b6)
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(ch, m)
    return not fails, worst[0]
for c in map(int, sys.argv[2:]):
    S = sorted(real)
    ok, w = closes({N: real[N] + c for N in S})
    if not ok:
        print(f"m={m} c={c}: open even with B6 at all {len(S)} points (worst {w:.4f})", flush=True); continue
    for N in sorted(S, reverse=True):
        T = [x for x in S if x != N]
        if T and closes({x: real[x] + c for x in T})[0]: S = T
    print(f"m={m} c={c}: closes; greedy minimal set {S} -> B6 {[real[x]+c for x in S]}", flush=True)
