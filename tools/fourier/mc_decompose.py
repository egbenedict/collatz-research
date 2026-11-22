"""Monte Carlo decomposition of phi_n(2^k) = E e(2^k Syrac(Z/3^n)/3^n) by trajectory features.

X = sum_i 3^{i-1} 2^{-s_i} mod 3^n, s_i = a_1+...+a_i (Tao 1.29).
We group the contribution E[e(.) ; feature = f] by: total valuation s_n, and the first pair-index where
s_i > k (exit of the main black triangle). Usage: python3 mc_decompose.py n k samples
"""
import sys
import numpy as np

n, k, S = int(sys.argv[1]), int(sys.argv[2]), int(float(sys.argv[3]))
mod = 3 ** n
inv2 = pow(2, -1, mod)
rng = np.random.default_rng(1)
tot = 0j
by_sn = {}
by_exit = {}
chunk = 2_000_000
inv2pows = np.array([pow(inv2, t, mod) for t in range(0, 400)], dtype=object)
for c0 in range(0, S, chunk):
    m = min(chunk, S - c0)
    a = rng.geometric(0.5, size=(m, n))
    s = np.cumsum(a, axis=1)
    X = np.zeros(m, dtype=np.int64)
    # 2^{-s} mod 3^n via table of 2^{-t} for t < 400 (s_n < 400 w.h.p.)
    tab = np.array([pow(inv2, t, mod) for t in range(400)], dtype=np.int64)
    s = np.minimum(s, 399)
    p3 = 1
    for i in range(n):
        X = (X + (p3 * tab[s[:, i]]) % mod) % mod
        p3 *= 3
    ph = np.exp(2j * np.pi * ((pow(2, k, mod) * X.astype(object)) % mod).astype(np.float64) / mod)
    tot += ph.sum()
    sn = s[:, -1]
    ex = np.argmax(s > k, axis=1)            # first i (0-based) with s_i > k
    ex[~(s > k).any(axis=1)] = n
    for key, arr, d in (("sn", sn, by_sn), ("ex", ex, by_exit)):
        u, inv = np.unique(arr, return_inverse=True)
        sums = np.zeros(len(u), complex)
        np.add.at(sums, inv, ph)
        cnt = np.bincount(inv)
        for uu, ss, cc in zip(u, sums, cnt):
            v = d.get(int(uu), [0j, 0])
            v[0] += ss; v[1] += cc
            d[int(uu)] = v
print(f"n={n} k={k} samples={S}: phi ~ {tot/S:.5f}  |phi| ~ {abs(tot/S):.5f}  (noise ~ {1/np.sqrt(S):.1e})")
for name, d in (("total valuation s_n", by_sn), ("first i with s_i > k", by_exit)):
    print(f"\ncontribution by {name}:  value  P(value)  E[e;value]  |E[e|value]|")
    for key in sorted(d):
        ssum, cnt = d[key]
        if cnt < S * 1e-4:
            continue
        print(f"  {key:4d}  {cnt/S:.4f}  {ssum/S:+.5f}  {abs(ssum/cnt):.4f}")
