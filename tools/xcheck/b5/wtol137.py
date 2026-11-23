"""Fall-weighted rule (I-017, iteration 10). Block lemma with the actual fall lengths:
x_t <= x_{t-r} + C_r(N) - 1 + r eps, C_r(N) = max over odd n < 2^N of (delta-1)*sum(k_0..k_{r-1}) - sum(l_0..l_{r-2})
(rigorous with the upper rational bound Dh1 of delta-1). Plugged into the chain as b_eff = (C + r - 1)/Dh1.
This script: largest C_5(532) (as a float threshold, bisection on rationals) that closes m with the real table.
Usage: wtol137.py m"""
import io, contextlib, json, sys
from fractions import Fraction as F
from b5lib import *
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
m = int(sys.argv[1]) if len(sys.argv) > 1 else 137
ch0 = checker(real); Dh1 = ch0.Dh1
def closes(C):
    t = dict(real); t[532] = (C + 4) / Dh1
    with contextlib.redirect_stdout(io.StringIO()):
        f, w = M.check_m(checker(t), m)
    return not f
std = lambda b: Dh1 * b - 4
print("Dh1 =", float(Dh1), " standard C for b=560:", float(std(560)), " b=561:", float(std(561)))
lo, hi = std(557), std(600)
assert closes(lo)
for _ in range(30):
    mid = (lo + hi) / 2
    mid = F(round(mid * 2**20), 2**20)
    if closes(mid): lo = mid
    else: hi = mid
print(f"m={m}: closes iff C_5(532) <= ~{float(lo):.5f} (equivalent b = {float((lo + 4) / Dh1):.3f})")
