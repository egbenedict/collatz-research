#!/usr/bin/env python3
"""Rule out a Collatz m-cycle with a *specific* number K of odd terms (idea I-002, round 3).

Setting (Hercher's accelerated map C(n) = n/2 or (3n+1)/2). A hypothetical cycle has
N = K + L elements, K odd, and m local minima n_1 < (all others), WLOG n_1 is the cycle's
smallest element. Write each minimum as n_i + 1 = a_i * 2^{k_i} (a_i odd), so the cycle
takes k_i o-steps after n_i, then l_i >= 1 e-steps down to n_{i+1}.

    y_i     := log2(n_i + 1),   alpha_i := log2(a_i) = y_i - k_i >= 0
    y_{i+1}  = y_i + (delta-1) k_i - l_i + eps_i,   0 < eps_i < eps_max (tiny)
    E_i     := delta*y_i - y_{i+1} = (delta-1) alpha_i + l_i - eps_i

Fact 1 (budget). Summing y over the m minima,
    sum_{i<m} G(m-i) E_i + sum_{i<=m} alpha_i = G(m) y_1 - K,      G(n) = (delta^n - 1)/(delta - 1).
Every term on the left is nonnegative (E_i >= 1 - eps_max), so any partial sum of the left
side is <= G(m) y_1 - K. That is linear in y_1 once (k_i, l_i) for i <= j are fixed.

Fact 2 (window). Taking logs of the cycle equation, sum over odd x of log(1 + 1/(3x)) = g ln 2
with g = N - K delta. Hence A <= sum 1/x <= A_hi, where A = 3 g ln 2. sum 1/x is the sum
of T(n_i) over runs, and (3/(n+1)) c_k <= T(n) <= (3/n) c_k with c_k = 1 - (2/3)^k.
Fact 1 forces the minima after n_1 to be huge, so their T values are tiny (R below). This
pins n_1 into [3c/A_hi - 1, 3c/(A - R)].

Sieve. Do a DFS over the exact run lengths (k_1, l_1, k_2, l_2, ...) of n_1's orbit, pruned
by Fact 1. As soon as the forced parity prefix of n_1 has more than log2(window) bits, the
residue class mod 2^P has at most one or two members in the window. Each candidate is then
checked exactly against Facts 1 and 2.

All real arithmetic is Decimal at 120 digits, and every bound is rounded so that it can only
*keep* a candidate (never wrongly discard one).
"""
from __future__ import annotations

import argparse
import sys
import time
from decimal import Decimal as D, getcontext
from itertools import count
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from collatz import COLLATZ, residue_for_parity  # noqa: E402

getcontext().prec = 120
LN2 = D(2).ln()
DELTA = D(3).ln() / LN2
TWO60 = 2 ** 60


def log2D(n: int) -> D:
    return D(n).ln() / LN2


def three_over(y: D) -> D:
    """Upper bound on 3/(2^y - 1) = 3/n for a minimum with y = log2(n+1); safe for huge y."""
    if y > 4000:
        return D(10) ** -1000            # 3/(2^y - 1) < 2^-3990 < 10^-1000
    return 3 / (D(2) ** y - 1)


class Problem:
    """A hypothetical m-cycle with exactly K odd and L even elements; all n <= X0 verified."""

    def __init__(self, m: int, K: int, L: int, X0: int, toy_window: tuple[int, int] | None = None):
        self.m, self.K, self.L, self.X0 = m, K, L, X0
        self.G = [(DELTA ** n - 1) / (DELTA - 1) for n in range(m + 1)]
        # eps_i < log2(1 + 2^{l_i}/(a_i 3^{k_i})) <= log2(1 + 1/n_{i+1}) < 1/(n_{i+1} ln 2) <= 1/(X0 ln 2)
        self.eps = 2 / (D(X0) * LN2)
        self.tail = [D(0)] * (m + 2)          # tail[j] = sum_{i=j}^{m-1} G(m-i) * (1 - eps)
        for j in range(m - 1, 0, -1):
            self.tail[j] = self.tail[j + 1] + self.G[m - j] * (1 - self.eps)
        self.toy = toy_window
        if toy_window is None:
            self.g = D(K + L) - D(K) * DELTA
            assert self.g > 0, "need 2^N > 3^K"
            self.A = 3 * self.g * LN2
            u = 1 / (3 * D(X0))
            self.A_hi = self.A / (1 - u / 2)
            self._init_window()

    # ---------------- Fact 1 machinery ----------------
    def _ylin(self, j: int, y1: D) -> D:
        """delta^{j-1} y1 - rho_j * B_j(y1): a lower bound on y_j (j >= 2), linear in y1.

        B_j(y1) = G(m) y1 - K - sum_{i=j}^{m-1} G(m-i)(1-eps) >= sum_{i<j} G(m-i) E_i, and
        sum_{i<j} delta^{j-1-i} E_i <= rho_j * B_j with rho_j = max_{i<j} delta^{j-1-i}/G(m-i),
        which is valid because every E_i >= 0. (If B_j < 0 that y1 is impossible anyway.)
        """
        m = self.m
        B = self.G[m] * y1 - self.K - self.tail[j]
        rho = max(DELTA ** (j - 1 - i) / self.G[m - i] for i in range(1, j))
        return DELTA ** (j - 1) * y1 - rho * B

    def ylow_from_budget(self, j: int, y1: D) -> D:
        """Lower bound on y_j for a cycle with exactly this y_1 (also y_j >= y_1)."""
        return max(y1, self._ylin(j, y1))

    def ylow_range(self, j: int, y_lo: D, y_hi: D) -> D:
        """Lower bound on y_j valid for every y_1 in [y_lo, y_hi] (the linear part is minimized at an endpoint)."""
        return max(y_lo, min(self._ylin(j, y_lo), self._ylin(j, y_hi)))

    def feasible_interval(self, runs: list[tuple[int, int]], lo: D, hi: D) -> tuple[D, D] | None:
        """Sub-interval of y_1 in [lo, hi] consistent with Fact 1 for completed runs (k_i, l_i)."""
        m, eps = self.m, self.eps
        j = len(runs)
        if j > m:
            return None
        coef, const = D(0), D(0)       # LHS(y1) = coef*y1 + const  (lower bound, eps_i -> 0)
        off = D(0)
        a, b = lo, hi
        for i, (k, l) in enumerate(runs, start=1):
            # y_i = y1 + off (+ up to (i-1) eps)
            if i >= 2:
                # n_i >= n_1  =>  y_i >= y_1  =>  off + (i-1) eps >= 0
                if off + (i - 1) * eps < 0:
                    return None
            # alpha_i >= 0: y1 + off + (i-1) eps >= k
            a = max(a, k - off - (i - 1) * eps)
            # alpha_i term
            coef += 1
            const += off - k
            if i <= m - 1:
                w = self.G[m - i]
                coef += w * (DELTA - 1)
                const += w * ((DELTA - 1) * (off - k) + l - eps)
            off += (DELTA - 1) * k - l
        if j + 1 <= m and j >= 1 and off + j * eps < 0:
            return None                              # n_{j+1} < n_1
        const += self.tail[j + 1] if j + 1 <= m - 1 else 0
        # need coef*y1 + const <= G(m) y1 - K   <=>   (G - coef) y1 >= const + K
        d = self.G[m] - coef
        rhs = const + self.K
        if d > 0:
            a = max(a, rhs / d)
        elif d < 0:
            b = min(b, rhs / d)
        elif rhs > 0:
            return None
        return (a, b) if a <= b else None

    # ---------------- Fact 2 machinery ----------------
    def R_bound(self, y_lo: D, y_hi: D) -> D:
        """Upper bound on sum_{j>=2} T(n_j) for any cycle with y_1 in [y_lo, y_hi]."""
        R = D(0)
        for j in range(2, self.m + 1):
            R += three_over(self.ylow_range(j, y_lo, y_hi))
        return R

    def _slice_feasible(self, a: D, b: D) -> bool:
        """Can n_1 have y_1 in [a, b]?  Needs n_1 <= 3/(A - R(a, b)) (since T(n_1) < 3/n_1)."""
        R = self.R_bound(a, b)
        if R >= self.A:
            return True
        return D(2) ** a - 1 <= 3 / (self.A - R)

    def _init_window(self):
        y_lo = log2D(self.X0 + 2)
        y_top = log2D(int(3 * self.m / self.A) + 2)           # sum T < 3m/n_1 <= ... crude
        # Slice [y_lo, y_top] and keep slices that are self-consistent. For large y_1 the
        # budget lets the last minima be small again, so the slicing must be fine enough.
        step = D("0.02")
        y_hi = None
        a = y_lo
        while a < y_top:
            b = min(a + step, y_top)
            if self._slice_feasible(a, b):
                y_hi = b
            a = b
        assert y_hi is not None
        self.y1_lo, self.y1_hi = y_lo, y_hi
        self.R = self.R_bound(y_lo, y_hi)
        if self.R >= self.A:
            raise RuntimeError("window did not close: R >= A on the final range")

    def n1_window(self, k1: int) -> tuple[int, int]:
        if self.toy:
            return self.toy
        c = 1 - (D(2) / 3) ** k1
        lo = max(self.X0 + 1, int(3 * c / self.A_hi) - 2)
        if self.R >= self.A:
            hi = int(D(2) ** self.y1_hi)
        else:
            hi = int(3 * c / (self.A - self.R)) + 2
        return lo, hi

    # ---------------- exact check of one candidate ----------------
    def check_candidate(self, n1: int, max_minima: int = 8) -> tuple[bool, str]:
        lo, hi = self.n1_window(_trailing_ones(n1))
        if not (lo <= n1 <= hi) or n1 % 2 == 0:
            return False, "outside window"
        runs, ns = [], []
        n = n1
        y1 = log2D(n1 + 1)
        for _ in range(min(max_minima, self.m)):
            k = _trailing_ones(n)
            peak = ((n + 1) >> k) * 3 ** k - 1
            l = (peak & -peak).bit_length() - 1
            ns.append((n, k))
            runs.append((k, l))
            n = peak >> l
            if n < n1 and len(runs) < self.m:
                return False, f"minimum {len(runs)+1} below n_1"
            if self.feasible_interval(runs, y1, y1) is None:
                return False, f"budget fails at minimum {len(runs)}"
        if self.toy:
            return True, "survived budget"
        # Fact 2: A - R_{>j} <= sum_{i<=j} T(n_i) <= A_hi
        S = D(0)
        for (nn, k) in ns:
            base = D(nn + 1)
            S += sum(1 / (base * (D(3) / 2) ** t - 1) for t in range(k))
        j = len(ns)
        Rrest = D(0)
        for jj in range(j + 1, self.m + 1):
            Rrest += three_over(self.ylow_from_budget(jj, y1))
        if S > self.A_hi:
            return False, f"T-sum too large ({float(S/self.A):.6f} A)"
        if S < self.A - Rrest:
            return False, f"T-sum too small ({float(S/self.A):.6f} A)"
        return True, "SURVIVOR"


def _trailing_ones(n: int) -> int:
    return ((n + 1) & -(n + 1)).bit_length() - 1


def sieve(P: Problem, verbose: bool = False, k1_max: int | None = None):
    """Enumerate every n_1 consistent with Fact 1 (and Fact 2's window). Returns stats + survivors."""
    stats = {"nodes": 0, "candidates": 0, "eliminated": {}, "survivors": []}
    lo_all = P.y1_lo if not P.toy else log2D(P.toy[0] + 1)
    hi_all = P.y1_hi if not P.toy else log2D(P.toy[1] + 1)
    k1_top = int(hi_all) + 1 if k1_max is None else k1_max

    def emit(prefix_bits: list[int], k1: int):
        lo, hi = P.n1_window(k1)
        Pbits = len(prefix_bits)
        r = residue_for_parity(tuple(prefix_bits), COLLATZ)
        mod = 1 << Pbits
        start = r + ((lo - r + mod - 1) // mod) * mod if r < lo else r
        n = start
        while n <= hi:
            stats["candidates"] += 1
            ok, why = P.check_candidate(n)
            if ok:
                stats["survivors"].append((n, why))
            else:
                stats["eliminated"][why.split(" (")[0]] = stats["eliminated"].get(why.split(" (")[0], 0) + 1
            n += mod

    def dfs(runs: list[tuple[int, int]], bits: list[int], k1: int, ylo: D, yhi: D):
        stats["nodes"] += 1
        lo, hi = P.n1_window(k1)
        width = max(hi - lo + 1, 1)
        # smallest next run length that keeps Fact 1 feasible (with l = 1)
        j = len(runs)
        kmin = None
        ymax_next = yhi + sum((DELTA - 1) * k - l for k, l in runs) + j * P.eps
        for k in range(1, int(ymax_next) + 2):
            if P.feasible_interval(runs + [(k, 1)], ylo, yhi) is not None:
                kmin = k
                break
        if kmin is None:
            return
        forced = bits + [1] * kmin
        if (1 << len(forced)) >= 4 * width or j + 1 >= P.m:
            emit(forced, k1)
            return
        for k in range(kmin, int(ymax_next) + 2):
            for l in count(1):
                iv = P.feasible_interval(runs + [(k, l)], ylo, yhi)
                if iv is None:
                    break
                dfs(runs + [(k, l)], bits + [1] * k + [0] * l, k1, iv[0], iv[1])

    t0 = time.time()
    for k1 in range(1, k1_top + 1):
        lo, hi = P.n1_window(k1)
        if lo > hi:
            continue
        ylo, yhi = max(lo_all, log2D(lo + 1)), min(hi_all, log2D(hi + 1))
        for l1 in count(1):
            iv = P.feasible_interval([(k1, l1)], ylo, yhi)
            if iv is None:
                break
            dfs([(k1, l1)], [1] * k1 + [0] * l1, k1, iv[0], iv[1])
        if verbose:
            print(f"  k1={k1:3d}: nodes {stats['nodes']}, candidates {stats['candidates']}, "
                  f"survivors {len(stats['survivors'])}  [{time.time()-t0:.1f}s]", flush=True)
    return stats


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--m", type=int, default=92)
    p.add_argument("--K", type=int, default=205632218873398596256)
    p.add_argument("--L", type=int, default=120287136981023372109)
    p.add_argument("--x0", type=int, default=2075, help="X0 in units of 2^60")
    p.add_argument("-v", "--verbose", action="store_true")
    a = p.parse_args(argv)
    P = Problem(a.m, a.K, a.L, a.x0 * TWO60)
    print(f"m={a.m}, K={a.K:,}, L={a.L:,}, X0={a.x0}*2^60")
    print(f"gap g = {float(P.g):.6e},  A = 3 g ln2 = {float(P.A):.6e}")
    print(f"y_1 = log2(n_1+1) in [{float(P.y1_lo):.4f}, {float(P.y1_hi):.4f}],  R = {float(P.R):.3e} ({float(P.R/P.A):.3f} A)")
    print(f"kappa = K/G(m) = {float(D(a.K)/P.G[a.m]):.4f}")
    st = sieve(P, verbose=a.verbose)
    print(f"\nDFS nodes {st['nodes']}, candidates checked {st['candidates']}")
    for why, c in sorted(st["eliminated"].items()):
        print(f"  eliminated: {c:6d}  {why}")
    print(f"SURVIVORS: {len(st['survivors'])}")
    for n, why in st["survivors"][:20]:
        print(f"  n_1 = {n}  ({why})")


if __name__ == "__main__":
    main()
