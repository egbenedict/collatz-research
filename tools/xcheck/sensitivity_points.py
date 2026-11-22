"""Which single table entries is m = 133 sensitive to? Raise one binding entry at a time by +1."""
import bisect
import sys
from math import ceil
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mchain as M
import rigorous as R
import rises

root = Path(__file__).resolve().parent / "marshall_data"
LIM = {2: 10000, 3: 400, 4: 80}
base = rises.load_marshall(root)
ch = M.Checker(root)
K = 27444133206411171953
p = ceil(K * M.D.hi)
Lam = ch.Lambda(K, p)
Ub, m2 = min(((133 - m2) / (Lam.lo - th), m2) for m2 in range(0, 133)
             if (th := ch.tail_hi(133, m2, K)) is not None and Lam.lo > th)
Y0 = R.log2(ceil(Ub)).hi
s, Y, use = ch.chain(133, Y0)
binding = set()
for t in range(1, 133):
    best_c, how, pt = ch.Dh * Y[t - 1] - 1 + M.EPS, 1, None
    for r, tb in ch.tab.items():
        if t >= r:
            i = bisect.bisect_left(tb.N, ceil(Y[t - r]))
            if i < len(tb.N):
                c = Y[t - r] + ch.Dh1 * tb.suf[i] - r + r * M.EPS
                if c < best_c:
                    q = next(q for q in range(i, len(tb.N)) if base[r][tb.N[q]] == tb.suf[i])
                    best_c, how, pt = c, r, tb.N[q]
    if how > 1:
        binding.add((how, pt))
print(f"{len(binding)} distinct binding table entries for the hard shape at m=133")
orig = rises.load_marshall
fragile = []
for r, N in sorted(binding):
    def patched(rt, r=r, N=N):
        t = orig(rt)
        t[r] = dict(t[r]); t[r][N] += 1
        return t
    M.load_marshall = patched
    c2 = M.Checker(root)
    fails, worst = M.check_m(c2, 133)
    M.load_marshall = orig
    checked = N <= LIM[r]
    if fails:
        fragile.append((r, N, base[r][N], checked))
    print(f"  B{r}({N}) = {base[r][N]} {'[independently checked]' if checked else ''}: +1 -> worst {worst[0]:.4f} {'OPENS m=133' if fails else 'still closed'}", flush=True)
print(f"\n{len(fragile)} single entries where +1 alone breaks m=133:")
for r, N, v, checked in fragile:
    print(f"  B{r}({N}) = {v}  {'independently checked' if checked else 'NOT independently checked'}")
