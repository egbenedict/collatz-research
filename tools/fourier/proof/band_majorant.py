"""Pair majorant F_m(J) on a band J in [-J0, KMAX], for large m.

F_m(J) = sum_{b>=2} 2^-b |sum_{a=1}^{b-1} e(th(m,J-a))| F_{m-2}(J-b),  F_0 = 1, F_1 = |Phi_1| (<=1).
Below the band (J < -J0) we use F = floor_val:
   floor_val = 1  -> rigorous majorant (up to float rounding) for all J in the band,
   floor_val = 0  -> heuristic (drops deep 2-adic paths, like band_recursion.py).
Usage: python3 band_majorant.py NMAX J0 [floor_val]
"""
import sys
import math
import numpy as np

BMAX = 60


def theta_band(m, J0, KMAX):
    mod = 3 ** m
    r = pow(2, -J0, mod)
    out = np.empty(J0 + KMAX + 1)
    for i in range(J0 + KMAX + 1):
        out[i] = r / mod
        r = (2 * r) % mod
    return out


def run(nmax, J0, floor_val, report=None):
    KMAX = int(1.6 * nmax) + 4
    L = J0 + KMAX + 1
    w = 2.0 ** -np.arange(0, BMAX + 1)
    pad = BMAX + 2
    F = {0: np.ones(L)}
    # Phi_1(J) = sum_a 2^-a e(2^(J-a)/3): exact periodic values, |Phi_1| = 1/sqrt3 for all J
    F[1] = np.full(L, 1 / math.sqrt(3))
    res = {}
    for m in range(2, nmax + 1):
        th = theta_band(m, J0, KMAX)
        E = np.exp(2j * np.pi * th)
        # extend below band: E for J in [-J0-pad, -J0) ; compute exactly
        mod = 3 ** m
        r = pow(2, -J0 - pad, mod)
        Elow = np.empty(pad, complex)
        for i in range(pad):
            Elow[i] = np.exp(2j * np.pi * (r / mod))
            r = 2 * r % mod
        Eext = np.concatenate([Elow, E])            # index i <-> J = -J0-pad+i
        Fext = np.concatenate([np.full(pad, floor_val), F[m - 2]])
        # also everything with J-b < -J0-pad uses floor_val (tail weight tiny for BMAX=60)
        cur = np.zeros(L)
        Sb = np.zeros(L, complex)
        base = pad  # Eext index of J=-J0
        idx = np.arange(L)
        for b in range(2, BMAX + 1):
            Sb += Eext[base + idx - (b - 1)] if b - 1 <= pad else 0
            j = base + idx - b
            fv = np.where(j >= 0, Fext[np.maximum(j, 0)], floor_val)
            cur += w[b] * np.abs(Sb) * fv
        F[m] = cur
        sel = slice(J0, None)  # J >= 0
        kb = int(np.argmax(cur[sel]))
        res[m] = (kb, cur[J0 + kb], cur[:J0].max())
        if report and m % report == 0:
            print(f"m={m:4d}  max_(J>=0) F = {cur[J0+kb]:.4e} at k={kb} (k - m log2 3 = {kb - m*math.log2(3):.1f})"
                  f"   max_(J<0) F = {cur[:J0].max():.4e}   -log/m = {-math.log(cur[J0+kb])/m:.4f}", flush=True)
        F.pop(m - 3, None)
    return res


if __name__ == "__main__":
    nmax = int(sys.argv[1]); J0 = int(sys.argv[2]); fl = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    run(nmax, J0, fl, report=int(sys.argv[4]) if len(sys.argv) > 4 else 20)
