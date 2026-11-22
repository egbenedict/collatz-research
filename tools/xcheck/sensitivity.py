"""How large an error in the *unchecked* table entries would break m <= 133?
Adds +c to every B_r(N) above the independently checked range and reruns the full m = 133 check."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mchain as M
import rises

LIM = {2: 10000, 3: 400, 4: 80}
orig = rises.load_marshall

def run(bump):
    def patched(root):
        t = orig(root)
        return {r: {n: v + (bump.get(r, 0) if n > LIM[r] else 0) for n, v in g.items()} for r, g in t.items()}
    M.load_marshall = patched
    ch = M.Checker(Path(__file__).resolve().parent / "marshall_data")
    fails, worst = M.check_m(ch, 133)
    M.load_marshall = orig
    return worst[0], len(fails)

for r in (2, 3, 4):
    lo, hi = 0, 64
    while hi - lo > 1:                     # largest c that still closes m = 133
        mid = (lo + hi) // 2
        w, nf = run({r: mid})
        lo, hi = (mid, hi) if nf == 0 else (lo, mid)
    print(f"B{r}: m=133 still closes with every unchecked B{r} entry raised by +{lo}" + ("" if hi == 64 and lo == 63 else f"; +{hi} opens it") + (" (search capped at 63)" if lo == 63 else ""), flush=True)
for c in (1, 2, 3):
    w, nf = run({2: c, 3: c, 4: c})
    print(f"all three +{c}: worst ratio {w:.4f}, open shapes {nf}", flush=True)
