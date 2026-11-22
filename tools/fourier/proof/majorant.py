"""Exact (float) comparison of |Phi_n(k)| with the pairing majorants, all k mod 2*3^(n-1).

Phi_m(J) = sum_{a>=1} 2^-a e(th(m,J-a)) Phi_{m-1}(J-a),  th(m,j) = (2^j mod 3^m)/3^m.

Pair majorant (Lemma 1 in NOTES.md): grouping the steps at levels m, m-1 with b = a+a',
    |Phi_m(J)| <= sum_{b>=2} 2^-b |sum_{a=1}^{b-1} e(th(m,J-a))| |Phi_{m-2}(J-b)|.
F: fixed pairing (m,m-1),(m-2,m-3),...;  F_0 = 1, F_1 = |Phi_1|.
U: adaptive DP  U_m = min(single-step bound from U_{m-1}, pair bound from U_{m-2}).
Both are rigorous majorants of |Phi_m| (up to float rounding and the b<=BMAX truncation, tail < 1e-15).

Usage: python3 majorant.py NMAX
"""
import sys
import numpy as np

BMAX = 64


def theta_row(m):
    P = 2 * 3 ** (m - 1)
    mod = 3 ** m
    t = np.empty(P, dtype=np.float64)
    r = 1
    # vectorised doubling in chunks of int64-safe size (mod < 2^62 for m <= 39)
    vals = np.empty(P, dtype=np.int64)
    vals[0] = 1
    k = 1
    while k < P:
        s = min(k, P - k)
        vals[k:k + s] = (vals[:s].astype(object) * pow(2, k, mod) % mod).astype(np.int64) if mod > 2**31 else (vals[:s] * pow(2, k, mod)) % mod
        k += s
    return vals.astype(np.float64) / mod


def lift(f, P):
    """periodic f of period len(f) dividing P, extended to length P"""
    return np.tile(f, P // len(f))


def shift(f, s):
    """g(J) = f(J - s) for periodic f"""
    return np.roll(f, s)


def run(nmax):
    rows = {}
    Phi = {0: np.ones(1, complex)}
    F = {0: np.ones(1)}
    U = {0: np.ones(1)}
    w = 2.0 ** -np.arange(0, BMAX + 1)
    out = []
    for m in range(1, nmax + 1):
        P = 2 * 3 ** (m - 1)
        th = theta_row(m)
        E = np.exp(2j * np.pi * th)
        # Phi_m(J) = sum_a 2^-a E(J-a) Phi_{m-1}(J-a)
        Pm1 = lift(Phi[m - 1], P)
        G = E * Pm1
        cur = np.zeros(P, complex)
        Gs = G.copy()
        for a in range(1, BMAX + 1):
            Gs = np.roll(Gs, 1)          # Gs(J) = G(J-a)
            cur += w[a] * Gs
        Phi[m] = cur
        # single-step bound from U_{m-1}
        Um1 = lift(U[m - 1], P)
        single = np.zeros(P)
        s = Um1.copy()
        for a in range(1, BMAX + 1):
            s = np.roll(s, 1)
            single += w[a] * s
        if m == 1:
            F[1] = np.abs(cur)
            U[1] = np.minimum(single, np.abs(cur) * 0 + single)
        else:
            Fm2 = lift(F[m - 2], P)
            Um2 = lift(U[m - 2], P)
            Sb = np.zeros(P, complex)      # Sb(J) = sum_{a=1}^{b-1} E(J-a)
            Eb = E.copy()
            pairF = np.zeros(P)
            pairU = np.zeros(P)
            Fs = np.roll(Fm2, 1)
            Us = np.roll(Um2, 1)
            for b in range(2, BMAX + 1):
                Eb = np.roll(Eb, 1)        # E(J-(b-1))
                Sb += Eb
                Fs = np.roll(Fs, 1)        # F(J-b)
                Us = np.roll(Us, 1)
                ab = np.abs(Sb)
                pairF += w[b] * ab * Fs
                pairU += w[b] * ab * Us
            F[m] = pairF
            U[m] = np.minimum(single, pairU)
        assert np.all(np.abs(cur) <= U[m] + 1e-12) and np.all(np.abs(cur) <= F[m] + 1e-12)
        res = (m, np.abs(cur).max(), int(np.abs(cur).argmax()), F[m].max(), int(F[m].argmax()), U[m].max(), int(U[m].argmax()))
        out.append(res)
        print("m=%2d  max|Phi|=%.4e (k=%d)  maxF=%.4e (k=%d)  maxU=%.4e (k=%d)" % res, flush=True)
        for d in (Phi, F, U):
            d.pop(m - 3, None)
    return out


if __name__ == "__main__":
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 12)
