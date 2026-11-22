"""Worst-case block gray persistence q_L = max P(Q-walk from (m,J) is gray at all of its next L positions).

Exact DP over the band J in [-J0, KMAX] (positions below -J0 are treated as gray = worst case, so the
reported q_L is an upper bound for starts in the band, up to float rounding).  Q-steps (-2, -b), P(b)=(b-1)2^-b.
Reported separately for the main triangle region (0 <= J < m log2 3) and the 2-adic region (J < 0).
Usage: python3 block_gray.py NMAX J0 L1,L2,...
"""
import sys
import math
import numpy as np

BMAX = 60


def run(nmax, J0, Ls):
    KMAX = int(1.6 * nmax) + 4
    L = J0 + KMAX + 1
    gray = {}
    for m in range(1, nmax + 1):
        mod = 3 ** m
        r = pow(2, -J0, mod)
        g = np.empty(L, bool)
        for i in range(L):
            x = r / mod
            g[i] = min(x, 1 - x) < 0.2
            r = 2 * r % mod
        gray[m] = g
    pb = np.array([(b - 1) * 2.0 ** -b for b in range(BMAX + 1)])
    Js = np.arange(-J0, KMAX + 1)
    for Lb in Ls:
        # q[m] = prob all of next Lb positions (including current) gray, starting at level m
        worst_adic, worst_main, where = 0.0, 0.0, None
        q = {}
        for m in range(1, nmax + 1):
            # base: q^{(1)}_m = gray[m]
            pass
        # iterate t = 1..Lb : Q^{(t)}_m(J) = gray_m(J) * sum_b pb Q^{(t-1)}_{m-2}(J-b)
        prev = {m: gray[m].astype(float) for m in range(1, nmax + 1)}
        for t in range(2, Lb + 1):
            cur = {}
            for m in range(1, nmax + 1):
                if m - 2 < 1:
                    cur[m] = np.zeros(L)  # walk ended: count as not-all-gray? treat as 1 (worst)
                    cur[m][:] = gray[m]
                    continue
                s = np.zeros(L)
                pv = prev[m - 2]
                for b in range(2, BMAX + 1):
                    sh = np.ones(L)          # below band -> gray (worst case)
                    sh[b:] = pv[:-b]
                    s += pb[b] * sh
                cur[m] = gray[m] * s
            prev = cur
        for m in range(2 * Lb + 1, nmax + 1):
            v = prev[m]
            mt = (Js >= 0) & (Js < m * math.log2(3))
            ad = (Js < 0) & (Js >= -J0 + 2 * Lb * 4)  # stay away from the artificial floor
            wa = v[ad].max()
            if wa > worst_adic:
                worst_adic, where = wa, (m, int(Js[ad][np.argmax(v[ad])]))
            worst_main = max(worst_main, v[mt].max())
        print(f"L={Lb:3d}: max over main-triangle starts = {worst_main:.4f};  max over 2-adic starts "
              f"(m<={nmax}, {-J0 + 8*Lb}<=J<0) = {worst_adic:.3e} at (m,J)={where}", flush=True)


if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), [int(x) for x in sys.argv[3].split(",")])
