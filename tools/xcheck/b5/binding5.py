"""List the rule-5 table entries that bind in the hard-shape chain for a hypothetical B5 = B4 + c."""
import sys
from math import ceil
from sens5 import *
def binding(ch, Y):
    out = []
    for t in range(1, len(Y)):
        best, how, pt = ch.Dh * Y[t-1] - 1 + M.EPS, 1, None
        for r, tb in ch.tab.items():
            if t >= r:
                b = tb.b(ceil(Y[t-r]))
                if b is not None:
                    c = Y[t-r] + ch.Dh1 * b - r + r * M.EPS
                    if c < best:
                        best, how, pt = c, r, (ceil(Y[t-r]), b)
        out.append((t, float(Y[t]), how, pt))
    return out
m, c, Nmax = map(int, sys.argv[1:4])
ch = checker(hyp(c, Nmax))
ratio, Y, use, m2 = hard_trace(ch, m)
print(f"m={m} c={c} Nmax={Nmax} ratio={ratio:.4f} use={use}")
for t, y, how, pt in binding(ch, Y)[:40]:
    print(t, f"{y:.2f}", how, pt if how == 5 else "")
