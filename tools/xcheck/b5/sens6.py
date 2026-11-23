"""Would a six-rise rule help? Real B5 table (b5_values.json, rule 5) plus a HYPOTHETICAL rule-6 table
B6(N) := B5base(N) + c at grid points N <= Nmax, where B5base is the real B5 if computed, else b4(N)+1
(the most optimistic lower bound, since B6 >= B5 + 1 >= b4 + 2). Prints the smallest Nmax closing m.
Hypothetical: no B6 value is computed here. Usage: sens6.py c m"""
import io, contextlib, json, sys
from b5lib import *
import rises
B4 = rises.load_marshall(ROOT)[4]
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
c = int(sys.argv[1]); m = int(sys.argv[2])
grid = sorted(set(B4) | set(real))
def run(hi):
    b6 = {N: (real[N] if N in real else B4[N] + 1) + c for N in grid if N <= hi}
    ch = checker(real); ch.tab[6] = M.Table(b6)
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(ch, m)
    return not fails, worst[0]
for hi in [x for x in grid if x >= 40]:
    ok, w = run(hi)
    if ok or hi in (60, 80, 100, 121, 150, 201, 250, 328, 400, 532):
        print(f"c={c} m={m} B6 on grid <= {hi}: {'closed' if ok else 'open'} (worst {w:.4f})", flush=True)
    if ok: break
