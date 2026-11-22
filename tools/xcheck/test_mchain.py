"""Controls for mchain.py. Run: python3 test_mchain.py"""
import math
import random
import sys
from fractions import Fraction as F
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mchain as M  # noqa: E402

# 1. lattice enumeration == brute force on small ranges
random.seed(3)
d = math.log2(3)
for _ in range(40):
    hi = random.randint(50, 60000)
    lo = random.randint(1, hi)
    eta = F(10 ** random.uniform(-5.5, -1)).limit_denominator(10 ** 9)
    got = {K for K, N in M.lattice(lo, hi, eta)}
    want = {K for K in range(lo, hi + 1) if 0 < math.ceil(K * d) - K * d < float(eta)}
    # float brute force can disagree only within 1e-9 of the threshold; check there are no misses
    miss = {K for K in want - got if abs((math.ceil(K * d) - K * d) - float(eta)) > 1e-9}
    assert not miss, (lo, hi, float(eta), sorted(miss)[:5])
print("ok  lattice() finds every K the brute force finds (40 random ranges)")

ch = M.Checker(Path(__file__).resolve().parent / "marshall_data")

# 2. negative control: beyond m = 133 the tables are too short, so something must stay open
for m in (134, 135):
    fails, worst = M.check_m(ch, m)
    assert fails, f"m={m} unexpectedly closed"
print("ok  negative control: m = 134, 135 report open shapes")

# 3. widened thresholds: inflate eta, every extra candidate must still be closed by the exact test
orig = M.lattice
count = [0]
def wide(lo, hi, eta):
    c = orig(lo, hi, min(F(1), eta * 10 ** 6))
    count[0] += len(c)
    return c
M.lattice = wide
fails, worst = M.check_m(ch, 120)
M.lattice = orig
assert not fails and count[0] > 1000, count[0]
print(f"ok  widened eta x1e6 at m=120: {count[0]} candidates, all closed by the exact per-shape test")

# 4. consistency: every one-sided best approximation K (smallest gaps) is closed by the exact test
conv = M.convergents()
semis = set()
for idx in range(1, len(conv) - 1):
    (p0, q0), (p1, q1), (p2, q2) = conv[idx - 1], conv[idx], conv[idx + 1]
    a = (q2 - q0) // q1
    for j in range(1, a + 1):
        semis.add((q0 + j * q1, p0 + j * p1))
semis = sorted((K, N) for K, N in semis if F(N, K) > M.D.hi)     # above delta: gap > 0
tested = 0
worst = (0, None, None)
for m in range(92, 134):
    Kc = ch.Kceil(m)
    for K, N in semis:
        if m <= K <= Kc:
            ok, why, ratio = ch.closes(m, K, N)
            tested += 1
            assert ok, (m, K, why, ratio)
            if ratio and ratio > worst[0]:
                worst = (ratio, m, K)
print(f"ok  {tested} (m, K) pairs from one-sided best approximations all closed; worst ratio {worst[0]:.6f} (m={worst[1]}, K={worst[2]})")
