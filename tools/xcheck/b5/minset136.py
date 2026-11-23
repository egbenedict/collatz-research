"""m=136: which B5 points beyond the computed table are needed? Hypothetical value b4(N)+c at
uncomputed B4 grid points N in (132, Nmax]; computed values where known. Prints the smallest Nmax that
closes m=136, then a greedy minimal set (drops the most expensive points first)."""
import io, contextlib, json, sys
from b5lib import *
import rises
B4 = rises.load_marshall(ROOT)[4]
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
c = int(sys.argv[1]); m = int(sys.argv[2]) if len(sys.argv) > 2 else 136
grid = sorted(N for N in B4 if N > 132 and N not in real)
def run(pts):
    t = {N: B4[N] + c for N in pts}; t.update(real)
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t), m)
    return not fails, worst[0]
Nmax = None
for hi in [x for x in grid if x >= 240]:
    ok, w = run([x for x in grid if x <= hi])
    print(f"c={c} m={m} all grid points <= {hi}: {'closed' if ok else 'open'} (worst {w:.4f})", flush=True)
    if ok: Nmax = hi; break
if Nmax is None: sys.exit()
S = [x for x in grid if x <= Nmax]
for N in sorted(S, reverse=True):
    trial = [x for x in S if x != N]
    ok, w = run(trial)
    if ok: S = trial
print(f"c={c} m={m} greedy minimal set: {S} -> hypothetical B5 {[B4[N]+c for N in S]}", flush=True)
