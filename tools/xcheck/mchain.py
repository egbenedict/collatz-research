#!/usr/bin/env python3
"""Independent cross-check of "no Collatz m-cycle for m <= 133" (Marshall et al.).

Only the rise-bound tables B_2, B_3, B_4 are taken from the published data
(marshall_data/, spot-checked by rises.py). Everything else is re-derived here and
implemented from scratch in rigorous rational interval arithmetic (rigorous.py).

Setting: T(n) = n/2 or (3n+1)/2. A cycle has K odd steps, period p, m local minima
n_0 (the least element), n_1, ..., n_{m-1} in cyclic order, x_t = log2(n_t + 1),
k_t = rise after n_t (k_t <= x_t, so sum x_t >= K), Lambda = p ln2 - K ln3 > 0.
Inputs: E1 every n <= X0 = 2^71 reaches 1; E3 (Rhin, via Simons-de Weger Lemma 12)
Lambda > exp(-13.3 (0.46057 + ln K)); the tables B_r(N) >= max first-r-rise total of odd n < 2^N.

(a) Upper bound on n_0. From prod(3 + 1/x) = 2^p, Lambda = sum log(1 + 1/(3x)) < (1/3) sum_t T(n_t),
    with T(n) < 3/n. Take a window of m2 consecutive minima maximising sum x_t (so its sum is
    >= (m2/m) K). By Lemma 20 (x_{t+1} <= delta x_t) and the corrigendum's vertex argument, when
    G(m2) = (delta^m2 - 1)/(delta - 1) < (m2/m) K the window's (1/3) sum T is
    < tail(m2) = 1/(2^v - 1) + (m2 - 1)/(2^(delta v) - 1), v = (m2/m) K / G(m2).
    The other m - m2 minima give (1/3) T(n_t) < 1/n_t <= 1/n_0. Hence
        n_0 < Ub := (m - m2)/(Lambda - tail(m2)),
    with m2 = 0 allowed (tail 0), and m2 = m giving the contradiction Lambda < tail(m).
    If Ub <= X0 + 1 the shape is impossible (E1). Universal: Lambda < K/(3 X0).
(b) Chain. Block (k, l) from odd n to n': 2^(k+1)(n'+1) <= 3^k (n+1) + 2^k (since l >= 1), so
    x' <= x + (delta-1) k - 1 + eps with eps <= 1/((n+1) ln 2) < 2^-70. With k <= x this is
    x' <= delta x - 1 + eps. Over r blocks, x_t <= x_{t-r} + (delta-1) s - r + r eps, where
    s <= B_r(N) whenever n_{t-r} < 2^N, and B_r is nondecreasing in N, so any grid point
    N >= ceil(Y_{t-r}) may be used. Define Y_0 >= log2 ceil(Ub) and
        Y_t = min(delta Y_{t-1} - 1 + eps, Y_{t-r} + (delta-1) b_r(Y_{t-r}) - r + r eps, r = 2, 3, 4).
    If sum_t Y_t < K the shape is impossible. ("compat" mode uses delta Y_{t-1}, as Marshall does.)
(c) Shapes. p = ceil(K delta) and 0 < p - K delta < Lambda/ln2. Above Kceil, E3 contradicts
    Lambda < tail(m). Below it, split K into dyadic ranges [lo, hi]. The chain sum F(Y) is
    nondecreasing in Y, so with Y* = sup{Y : F(Y) < lo}, every K in the range whose Ub is
    <= U* = max(X0 + 1, 2^Y* - 1) is closed. That holds whenever Lambda >= tail(m2, lo) + (m - m2)/U*,
    so only K with gap < eta (the minimum over m2) can survive. Those are enumerated
    exhaustively with a unimodular convergent basis and each is tested exactly.
"""
from __future__ import annotations

import argparse
import bisect
import sys
import time
from fractions import Fraction as F
from math import ceil, floor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rigorous as R                     # noqa: E402
from rises import load_marshall          # noqa: E402

X0 = 2 ** 71
D = R.DELTA
EPS = F(1, 2 ** 70)
Q64 = 2 ** 64


def ru64(x: F) -> F:
    return F(ceil(x * Q64), Q64)


def rd64(x: F) -> F:
    return F(floor(x * Q64), Q64)


class Table:
    def __init__(self, g: dict[int, int]):
        self.N = sorted(g)
        vals = [g[n] for n in self.N]
        self.suf = vals[:]
        for i in range(len(vals) - 2, -1, -1):
            self.suf[i] = min(self.suf[i], self.suf[i + 1])

    def b(self, Nmin: int) -> int | None:
        """min table value over grid points N >= Nmin (None if none)."""
        i = bisect.bisect_left(self.N, Nmin)
        return self.suf[i] if i < len(self.N) else None


class Checker:
    def __init__(self, data_root: Path, mode: str = "own", rules=(2, 3, 4)):
        tabs = load_marshall(data_root)
        self.tab = {r: Table(tabs[r]) for r in rules}
        self.mode = mode
        self.Dh = D.hi
        self.Dh1 = D.hi - 1
        self._G = {}

    # ----- chain -----
    def chain(self, m: int, Y0: F, stop_above: F | None = None, trace=False):
        """Upper bound F(Y0) on sum_t x_t; returns (sum, Ys, usage)."""
        Y = [ru64(Y0)]
        s = Y[0]
        use = {}
        for t in range(1, m):
            if self.mode == "own":
                best, how = self.Dh * Y[t - 1] - 1 + EPS, 1
            else:
                best, how = self.Dh * Y[t - 1], 1
            for r, tb in self.tab.items():
                if t >= r:
                    b = tb.b(ceil(Y[t - r]))
                    if b is not None:
                        c = Y[t - r] + self.Dh1 * b - r + r * EPS
                        if c < best:
                            best, how = c, r
            y = ru64(best)
            Y.append(y)
            s += y
            use[how] = use.get(how, 0) + 1
            if stop_above is not None and s >= stop_above:
                return s, Y, use
        return s, Y, use

    def Ystar(self, m: int, Klo: int, ylo: F, yhi: F = F(4000)) -> F | None:
        """Largest Y (on a 2^-40 grid) with chain(m, Y) < Klo, or None if even ylo fails."""
        if self.chain(m, ylo, stop_above=F(Klo))[0] >= Klo:
            return None
        if self.chain(m, yhi, stop_above=F(Klo))[0] < Klo:
            return yhi
        a, b = ylo, yhi
        while b - a > F(1, 2 ** 40):
            mid = rd64((a + b) / 2) if (a + b) / 2 - rd64((a + b) / 2) < F(1, 2 ** 41) else (a + b) / 2
            mid = F(floor(mid * 2 ** 40), 2 ** 40)
            if mid <= a:
                mid = a + F(1, 2 ** 40)
            if self.chain(m, mid, stop_above=F(Klo))[0] < Klo:
                a = mid
            else:
                b = mid
        return a

    # ----- Hercher/corrigendum bound on n_0 -----
    def G(self, n: int) -> R.I:
        if n not in self._G:
            self._G[n] = ((D ** n) - 1) / (D - 1)
        return self._G[n]

    def tail_hi(self, m: int, m2: int, K: int) -> F | None:
        """Upper bound on tail(m2) for this K, or None if m2 is not admissible."""
        if m2 == 0:
            return F(0)
        Gm2 = self.G(m2)
        if not (Gm2.hi < F(m2, m) * K):
            return None
        v_lo = F(m2, m) * K / Gm2.hi
        if v_lo > 4000:
            return F(m2) / F(2) ** 3990
        a = R.pow2(v_lo).lo - 1
        t = 1 / a
        if m2 > 1:
            b = R.pow2(D.lo * v_lo).lo - 1
            t += (m2 - 1) / b
        return R.rup(t)

    def Lambda(self, K: int, p: int) -> R.I:
        gap = R.I(p - K * D.hi, p - K * D.lo)
        return R.LN2 * gap

    def closes(self, m: int, K: int, p: int):
        """Exact test for one shape. Returns (closed, reason, best ratio or None)."""
        Lam = self.Lambda(K, p)
        if Lam.hi <= 0:
            return True, "p < K delta", None
        if Lam.lo >= F(K, 3 * X0):
            return True, "Lambda >= K/(3 X0)", None
        best = None
        for m2 in range(0, m + 1):
            th = self.tail_hi(m, m2, K)
            if th is None:
                continue
            if m2 == m:
                if Lam.lo >= th:
                    return True, "Lambda >= tail(m)", None
                continue
            if Lam.lo <= th:
                continue
            Ub = (m - m2) / (Lam.lo - th)
            if Ub <= X0 + 1:
                return True, f"n_0 < Ub <= X0+1 (m2={m2})", None
            if best is None or Ub < best[0]:
                best = (Ub, m2)
        if best is None:
            return False, "no usable m2", None
        Ub, m2 = best
        Y0 = R.log2(ceil(Ub)).hi
        s, _, use = self.chain(m, Y0)
        return s < K, f"chain (m2={m2}, Y0={float(Y0):.4f}, use={use})", float(s / K)

    # ----- Rhin ceiling -----
    def rhin_lo(self, K: int) -> F:
        z = F(133, 10) * (F(46057, 100000) + R.ln(K).hi)
        return 1 / R.exp(R.rup(z)).hi

    def Kceil(self, m: int) -> int:
        Gm = self.G(m).hi
        K = 1 << max(1, ceil(F(133, 10) * Gm / R.LN2.lo).bit_length())
        while True:
            v = F(K) / Gm
            bound = 2 * m / R.pow2(min(v, F(4000))).lo if v < 4000 else F(2 * m) / F(2) ** 3990
            if self.rhin_lo(K) > bound:
                return K
            K <<= 1


# ----- exhaustive lattice enumeration of K in [lo, hi] with 0 < ceil(K d) - K d < eta -----
def convergents(depth: int = 200):
    """Convergents of delta from the common continued-fraction prefix of D.lo and D.hi."""
    def cf(x: F, n):
        out = []
        for _ in range(n):
            a = floor(x)
            out.append(a)
            if x == a:
                break
            x = 1 / (x - a)
        return out
    A, Bq = cf(D.lo, depth), cf(D.hi, depth)
    common = []
    for a, b in zip(A, Bq):
        if a != b:
            break
        common.append(a)
    common = common[:-1]                     # drop the last agreed term to be safe
    out = []
    p0, q0, p1, q1 = 1, 0, common[0], 1
    out.append((p1, q1))
    for a in common[1:]:
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
        out.append((p1, q1))
    return out


CONV = None


def lattice(lo: int, hi: int, eta: F):
    """All (K, N) with lo <= K <= hi and 0 < N - K delta < eta (a superset, by interval delta)."""
    global CONV
    if CONV is None:
        CONV = convergents()
    idx = max(i for i, (_, q) in enumerate(CONV) if q <= hi)
    (pa, qa), (pb, qb) = CONV[idx], CONV[idx + 1]
    assert abs(pa * qb - pb * qa) == 1
    ea = R.I(pa - qa * D.hi, pa - qa * D.lo)
    eb = R.I(pb - qb * D.hi, pb - qb * D.lo)
    J = floor(hi * max(abs(ea.lo), abs(ea.hi)) + eta * qa) + 2
    out = []
    for j in range(-J, J + 1):
        i_lo_k = ceil(F(lo - j * qb, qa))
        i_hi_k = floor(F(hi - j * qb, qa))
        if i_lo_k > i_hi_k:
            continue
        # e = i ea + j eb in (0, eta): bound i using all endpoint combinations (conservative)
        cands = []
        for a in (ea.lo, ea.hi):
            for b in (eb.lo, eb.hi):
                if a == 0:
                    continue
                cands += [(-j * b) / a, (eta - j * b) / a]
        if cands:
            i_lo = max(i_lo_k, floor(min(cands)) - 1)
            i_hi = min(i_hi_k, ceil(max(cands)) + 1)
        else:
            i_lo, i_hi = i_lo_k, i_hi_k
        for i in range(i_lo, i_hi + 1):
            K, N = i * qa + j * qb, i * pa + j * pb
            if lo <= K <= hi and N - K * D.lo > 0 and N - K * D.hi < eta:
                out.append((K, N))
    return out


def check_m(ch: Checker, m: int, verbose=False):
    t0 = time.time()
    Kc = ch.Kceil(m)
    ylo = R.log2(X0 + 2).lo                       # n_0 >= X0 + 1, so x_0 >= log2(X0 + 2)
    lo = m
    n_ranges = n_cands = 0
    fails, worst = [], (0.0, None)
    while lo <= Kc:
        hi = min(2 * lo - 1, Kc)
        n_ranges += 1
        eta_univ = F(hi, 3 * X0) / R.LN2.lo
        Ys = ch.Ystar(m, lo, ylo)
        Ustar = F(X0 + 1)
        if Ys is not None:
            Ustar = max(Ustar, R.pow2(Ys).lo - 1)
        eta_chain = None
        for m2 in range(0, m + 1):
            th = ch.tail_hi(m, m2, lo)
            if th is None:
                continue
            lam = th if m2 == m else th + (m - m2) / Ustar
            e = lam / R.LN2.lo
            if eta_chain is None or e < eta_chain:
                eta_chain = e
        eta = min(eta_univ, eta_chain) if eta_chain is not None else eta_univ
        if eta >= 1:
            eta = F(1)
        cands = lattice(lo, hi, eta)
        n_cands += len(cands)
        for K, N in cands:
            ok, why, ratio = ch.closes(m, K, N)
            if ratio is not None and ratio > worst[0]:
                worst = (ratio, K)
            if not ok:
                fails.append((K, N, why, ratio))
            if verbose:
                print(f"   m={m} K={K} p={N}: {'closed' if ok else 'OPEN'} by {why}, ratio {ratio}")
        lo = hi + 1
    print(f"m={m:3d}: Kceil={float(Kc):.4e}, {n_ranges} K-ranges, {n_cands} dangerous candidates tested, "
          f"{len(fails)} open, worst chain ratio {worst[0]:.7f} at K={worst[1]}  [{time.time()-t0:.0f}s]", flush=True)
    for f_ in fails:
        print("   OPEN:", f_, flush=True)
    return fails, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent / "marshall_data"))
    ap.add_argument("--mode", choices=["own", "compat"], default="own")
    ap.add_argument("--mmin", type=int, default=92)
    ap.add_argument("--mmax", type=int, default=133)
    ap.add_argument("--hard", action="store_true", help="only print ratios of known hard shapes")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args()
    ch = Checker(Path(a.root), a.mode)
    if a.hard:
        for m in range(a.mmin, a.mmax + 1):
            for K in (27444133206411171953, 77692117359936589403, 205632218873398596256):
                p = ceil(K * D.hi)
                assert ceil(K * D.lo) == p
                ok, why, ratio = ch.closes(m, K, p)
                print(f"m={m} K={K}: {'closed' if ok else 'OPEN'} ratio={ratio} ({why})", flush=True)
        return
    allopen = []
    for m in range(a.mmin, a.mmax + 1):
        fails, _ = check_m(ch, m, a.verbose)
        allopen += [(m,) + f for f in fails]
    print("ALL CLOSED" if not allopen else f"{len(allopen)} OPEN SHAPES")


if __name__ == "__main__":
    main()
