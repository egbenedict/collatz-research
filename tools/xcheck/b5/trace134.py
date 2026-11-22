import bisect, sys
from math import ceil
from b5lib import *
ch = checker()
for m in (133, 134, 135):
    ratio, Y, use, m2 = hard_trace(ch, m)
    print(f"m={m} ratio={ratio:.5f} m2={m2} use={use}")
m = 134
ratio, Y, use, m2 = hard_trace(ch, m)
for t in range(m):
    # which rule bound Y_t
    how = "-"
    if t:
        best, how = ch.Dh * Y[t-1] - 1 + M.EPS, 1
        for r, tb in ch.tab.items():
            if t >= r:
                b = tb.b(ceil(Y[t-r]))
                if b is not None:
                    c = Y[t-r] + ch.Dh1 * b - r + r * M.EPS
                    if c < best: best, how = c, r
    print(t, f"{float(Y[t]):.3f}", how)
