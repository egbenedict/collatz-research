"""Banded transfer recursion for Phi_n(k) := E e(2^k Syrac(Z/3^n) / 3^n), every unit frequency xi = 2^k.

From Syrac(Z/3^n) = 2^{-a} (3 Syrac(Z/3^{n-1}) + 1):

    Phi_n(k) = sum_{a >= 1} 2^{-a} e(2^{k-a} / 3^n) Phi_{n-1}(k - a),       Phi_0 = 1,

where 2^{k-a} is reduced mod 3^n (negative exponents are inverses). The paths only move down in k.
We keep exact values (full FFT tables from syrac_decay.run) for levels m <= m_exact. Above that we keep
only k in the band [-J0, KMAX] and treat Phi_m(k) = 0 for k < -J0. That drops the paths that go deep
into the white region, where the phases are essentially random. We validate against exact n <= 17.

Usage: python3 band_recursion.py NMAX [J0] [m_exact]
"""
import sys
import math
import numpy as np
from syrac_decay import run

TWO_PI = 2 * math.pi


def phase(j, m):
    """e(2^j / 3^m) with 2^j reduced mod 3^m (j may be negative)."""
    mod = 3 ** m
    r = pow(2, j, mod)
    return complex(math.cos(TWO_PI * r / mod), math.sin(TWO_PI * r / mod)) if m < 60 else \
        _phase_big(r, mod)


def _phase_big(r, mod):
    # r/mod as a float loses nothing important: we need it to ~1e-16 relative, fine.
    from fractions import Fraction
    x = float(Fraction(r, mod))
    return complex(math.cos(TWO_PI * x), math.sin(TWO_PI * x))


def band(nmax, J0=60, m_exact=11, kmax=None, amax=200):
    kmax = kmax or int(2.1 * nmax) + 10
    _, laws = run(m_exact, verbose=False, keep=True)
    exact = {}
    for m in range(1, m_exact + 1):
        vals, p, _ = laws[m]
        # complex Phi_m(k) for all k mod M: cross-correlation with e(2^e/3^m)
        mod = 3 ** m
        g = np.exp(2j * np.pi * vals.astype(float) / mod)
        exact[m] = np.fft.ifft(np.conj(np.fft.fft(p)) * np.fft.fft(g))
    ks = np.arange(-J0, kmax + 1)
    off = J0
    prev = np.array([exact[m_exact][k % len(exact[m_exact])] for k in ks])
    out = {}
    w = 2.0 ** -np.arange(1, amax + 1)
    for m in range(m_exact + 1, nmax + 1):
        psi = np.array([phase(int(k), m) for k in ks])
        # Psi(j) = e(2^j/3^m) Phi_{m-1}(j); then Phi_m(k) = sum_a 2^{-a} Psi(k - a)
        Psi = psi * prev
        if m - 1 == m_exact:
            # below the band we still know Phi_{m-1} exactly: add the tail sum_{k-a < -J0} once
            tail_vals = None
        cur = np.zeros(len(ks), complex)
        for a in range(1, amax + 1):
            cur[a:] += w[a - 1] * Psi[:-a]
        if m - 1 == m_exact:
            Mx = len(exact[m_exact])
            for idx, k in enumerate(ks):
                # contributions with k - a < -J0
                amin = int(k) + J0 + 1
                s = 0j
                for a in range(max(1, amin), amin + amax):
                    j = int(k) - a
                    s += 2.0 ** -a * phase(j, m) * exact[m_exact][j % Mx]
                cur[idx] += s
        prev = cur
        mag = np.abs(cur)
        sel = ks >= 0
        kb = int(ks[sel][np.argmax(mag[sel])])
        out[m] = (kb, float(mag[off + kb]), mag.copy())
    return ks, out


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    J0 = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    me = int(sys.argv[3]) if len(sys.argv) > 3 else 11
    ks, out = band(nmax, J0, me)
    L = math.log2(3)
    x = L
    I0 = (x - 1) * math.log(2 - 2 / x) + math.log(2 / x)
    print(f"I0 = I(log2 3) = {I0:.5f}")
    print(" n    k*   k*/n   max|Phi|      -log/n   -log/sqrt(n)  (-log - I0 n)/n^(2/3)")
    for m, (kb, v, _) in out.items():
        lg = -math.log(v)
        print(f"{m:3d} {kb:5d}  {kb/m:.3f}  {v:.4e}   {lg/m:.4f}   {lg/math.sqrt(m):.4f}    {(lg - I0*m)/m**(2/3):.4f}", flush=True)
