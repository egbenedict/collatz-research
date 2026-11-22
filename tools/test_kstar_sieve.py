"""Toy validation of the Fact-1 sieve against brute force.

Toy statement (the same algebra as the cycle case, minus the cycle closure):
for odd n_1 in [lo, hi] whose next m-1 local minima are all >= n_1, list those with
sum_{i<=m} k_i >= K_t. The sieve must find every one (its output is a superset).

Run: python3 test_kstar_sieve.py
"""
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kstar_sieve import Problem, sieve, _trailing_ones  # noqa: E402


def minima_runs(n, m):
    out = []
    for _ in range(m):
        k = _trailing_ones(n)
        peak = ((n + 1) >> k) * 3 ** k - 1
        l = (peak & -peak).bit_length() - 1
        out.append((n, k, l))
        n = peak >> l
    return out


def brute(lo, hi, m, Kt):
    hits = set()
    for n in range(lo | 1, hi + 1, 2):
        r = minima_runs(n, m)
        if all(x[0] >= n for x in r[1:]) and sum(x[1] for x in r) >= Kt:
            hits.add(n)
    return hits


def run_case(bits, m, quantile, seed):
    random.seed(seed)
    lo, hi = (1 << bits) + 1, (1 << (bits + 1)) - 1
    # choose K_t near the top of the distribution so the sieve has to work
    sample = []
    for _ in range(20000):
        n = random.randrange(lo | 1, hi, 2)
        r = minima_runs(n, m)
        if all(x[0] >= n for x in r[1:]):
            sample.append(sum(x[1] for x in r))
    sample.sort()
    Kt = sample[int(quantile * (len(sample) - 1))]
    t = time.time()
    B = brute(lo, hi, m, Kt)
    tb = time.time() - t
    P = Problem(m, Kt, 0, lo - 1, toy_window=(lo, hi))
    t = time.time()
    st = sieve(P)
    ts = time.time() - t
    S = {n for n, _ in st["survivors"]}
    missing = B - S
    print(f"bits={bits} m={m} K_t={Kt}: brute {len(B)} hits ({tb:.1f}s); sieve checked "
          f"{st['candidates']} candidates, {len(S)} survivors ({ts:.1f}s); missing {len(missing)}")
    assert not missing, sorted(missing)[:10]
    return len(B), st["candidates"]


if __name__ == "__main__":
    for bits, m, q, seed in [(16, 3, 0.999, 1), (16, 4, 0.999, 2), (18, 4, 0.9995, 3),
                             (18, 5, 0.999, 4), (20, 4, 0.9999, 5), (14, 6, 0.99, 6)]:
        run_case(bits, m, q, seed)
    print("ok: sieve never missed a brute-force hit")
