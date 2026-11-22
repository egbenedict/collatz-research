"""Which B5 points bind for the open m-shapes? Real computed values + hypothetical B4(N)+c elsewhere."""
import json, sys
from math import ceil
from b5lib import *
from binding5 import binding
import rises
B4 = rises.load_marshall(ROOT)[4]
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
m, c, Nmax = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
extra = [int(x) for x in sys.argv[4:]]       # extra hypothetical points (value b4(N)+c)
T4 = M.Table(B4)
hyp = {N: v + c for N, v in B4.items() if N <= Nmax and N not in real}
hyp.update({N: T4.b(N) + c for N in extra if N not in real})
tab = dict(hyp); tab.update(real)
ch = checker(tab)
for K in HARD:
    p = ceil(K * M.D.hi)
    ok, why, ratio = ch.closes(m, K, p)
    Lam = ch.Lambda(K, p)
    Ub, m2 = min(((m - m2) / (Lam.lo - th), m2) for m2 in range(0, m) if (th := ch.tail_hi(m, m2, K)) is not None and Lam.lo > th)
    s, Y, use = ch.chain(m, R.log2(ceil(Ub)).hi)
    pts = sorted({pt for t, y, how, pt in binding(ch, Y) if how == 5})
    print(f"K={K:.3e}: ratio {ratio:.4f} closed={ok}; rule-5 points used: " +
          ", ".join(f"{N}:{v}{'' if N in real else '*'}" for N, v in pts))
print("(* = hypothetical B4+c)")
