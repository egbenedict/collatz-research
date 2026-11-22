"""Exact law and characteristic function of Tao's Syracuse random variable Syrac(Z/3^n).

Definition (Tao, Lemma 1.12):  Syrac(Z/3^{n+1}) = (3*Syrac(Z/3^n) + 1) * 2^{-a},  a ~ Geom(2) iid,
with Syrac(Z/3^0) = 0.

Every value is a unit mod 3^n, and 2 is a primitive root mod 3^n. So we store the law in
discrete-log coordinates e mod M = 2*3^{n-1} (value 2^e). Multiplying by 2^{-a} is a cyclic
convolution with the folded geometric kernel w[-r mod M] += 2^{-r}.

Every frequency xi with 3 not dividing xi is also a power of 2: xi = 2^k. So
    phi_n(2^k) = sum_e p[e] * e(2^{e+k} / 3^n)
is a cyclic cross-correlation of p with g[e] = e(2^e/3^n). One FFT of size M gives every unit xi.

Usage:  python3 syrac_decay.py NMAX         (NMAX <= 18 fits in 16 GB)
        python3 syrac_decay.py --selftest
"""
import sys
import time
import numpy as np
from fractions import Fraction


def pow2_table(M, mod, dtype=np.int64):
    """t[e] = 2^e mod `mod` for e in [0, M)."""
    t = np.empty(M, dtype=np.int64)
    t[0] = 1
    k = 1
    while k < M:
        step = min(k, M - k)
        t[k:k + step] = (t[:step] * pow(2, k, mod)) % mod
        k += step
    return t.astype(dtype)


def geom_kernel_rfft(M):
    w = np.zeros(M)
    rmax = min(M, 1100)          # 2^{-1100} underflows anyway
    r = np.arange(1, rmax + 1)
    np.add.at(w, (-r) % M, np.exp2(-r.astype(float)))
    if M <= 1100:
        w /= (1.0 - 2.0 ** (-M))
    return np.fft.rfft(w)


def unit_char(p, vals, mod):
    """Return |phi(2^k)| for k in [0, M): cross-correlation of p with e(2^e/mod)."""
    g = np.exp(2j * np.pi * (vals.astype(np.float64) / mod))
    c = np.fft.ifft(np.conj(np.fft.fft(p)) * np.fft.fft(g))
    del g
    return np.abs(c)


def run(nmax, verbose=True, keep=False):
    results, laws = [], {}
    prev_vals = prev_p = None
    for n in range(1, nmax + 1):
        t0 = time.time()
        mod = 3 ** n
        M = 2 * 3 ** (n - 1)
        pw = pow2_table(M, mod)
        log = np.full(mod, -1, dtype=np.int32)
        log[pw] = np.arange(M, dtype=np.int32)
        q = np.zeros(M)
        if n == 1:
            q[log[1]] = 1.0
        else:
            y = (3 * prev_vals.astype(np.int64) + 1) % mod
            np.add.at(q, log[y], prev_p)
            del y
        del log
        p = np.fft.irfft(np.fft.rfft(q) * geom_kernel_rfft(M), n=M)
        del q
        p = np.maximum(p, 0.0)
        p /= p.sum()
        mag = unit_char(p, pw, mod)
        k = int(np.argmax(mag))
        top = np.argsort(mag)[-4:][::-1]
        xi = pow(2, k, mod)
        xi_signed = xi if xi <= mod // 2 else xi - mod
        r = dict(n=n, max=float(mag[k]), k=k, xi=xi_signed,
                 top=[(int(t), float(mag[t])) for t in top],
                 mean_sq=float(np.mean(mag ** 2)))
        results.append(r)
        if verbose:
            print(f"n={n:2d}  max|phi|={r['max']:.6e}  k={k:6d}  xi={xi_signed:>10d}  "
                  f"rms={np.sqrt(r['mean_sq']):.3e}  ({time.time()-t0:.1f}s)", flush=True)
        if keep:
            laws[n] = (pw.copy(), p.copy(), mag.copy())
        del mag
        prev_vals, prev_p = pw, p
    return results, laws


def exact_small(n, amax):
    mod = 3 ** n
    inv2 = pow(2, -1, mod)
    dist = {0: Fraction(1)}
    for _ in range(n):
        new = {}
        for x, pr in dist.items():
            base = (3 * x + 1) % mod
            for a in range(1, amax + 1):
                y = base * pow(inv2, a, mod) % mod
                new[y] = new.get(y, 0) + pr * Fraction(1, 2 ** a)
        dist = new
    return dist


def selftest():
    res, laws = run(5, verbose=False, keep=True)
    tao = [0, 8, 16, 0, 11, 4, 0, 2, 22]            # Tao, after Lemma 1.12: 63 * P(Syrac(Z/9) = x)
    vals, p, _ = laws[2]
    law = np.zeros(9)
    law[vals] = p
    assert np.allclose(law * 63, tao, atol=1e-9), law * 63
    print("Z/9 law x63 matches Tao:", tao)
    for n, amax in ((3, 40), (4, 25)):
        d = exact_small(n, amax)
        vals, p, mag = laws[n]
        law = np.zeros(3 ** n)
        law[vals] = p
        bf = np.zeros(3 ** n)
        for x, pr in d.items():
            bf[x] = float(pr)
        err = np.abs(law - bf).max()
        assert err < 2 ** -(amax - 6), err
        # characteristic function against the additive FFT
        full = np.abs(np.fft.fft(law))
        M = len(vals)
        for k in range(M):
            xi = pow(2, k, 3 ** n)
            assert abs(full[(-xi) % 3 ** n] - mag[k]) < 1e-12 or abs(full[xi] - mag[k]) < 1e-12
        print(f"n={n}: law vs brute force {err:.1e}; unit char fn matches additive FFT")
    print("selftest OK")


if __name__ == "__main__":
    if sys.argv[1:] == ["--selftest"]:
        selftest()
    else:
        nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
        res, _ = run(nmax)
        print("\n n   max|phi|      ratio    k     k/n    log max / n   log max / sqrt n")
        for i, r in enumerate(res):
            rat = r['max'] / res[i - 1]['max'] if i else float('nan')
            lm = np.log(r['max'])
            print(f"{r['n']:2d}  {r['max']:.6e}  {rat:.4f}  {r['k']:5d}  {r['k']/r['n']:.3f}  "
                  f"{lm/r['n']:.4f}   {lm/np.sqrt(r['n']):.4f}")
