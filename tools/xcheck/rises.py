#!/usr/bin/env python3
"""Independent computation of the r-rise bounds B_r(N), for cross-checking published tables.

B_r(N) = max over odd n < 2^N of k_0 + ... + k_{r-1}, the total length of the first r
rises (runs of o-steps) of n under T(n) = n/2, (3n+1)/2.

Method (own derivation, see notes):
* For run lengths sigma = (k_0, l_0, ..., k_{r-2}, l_{r-2}), let nu(sigma) be the 2-adic
  integer whose parity vector is 1^k0 0^l0 ... 1^k_{r-2} 0^l_{r-2} 1^inf. Undoing the steps from
  the fixed point -1 (x -> 2x for an e-step, x -> (2x-1)/3 for an o-step, and
  g^k(x) + 1 = (2/3)^k (x+1)) gives
      nu = -1 + sum_i 2^(l_0+...+l_{i-1}) * u(k_0+...+k_i) * (1 - 2^(l_i)),   u(K) = (2/3)^K in Z_2.
* Terras: the first P parities are in bijection with residues mod 2^P, so n with parity
  prefix pi lies in the class nu mod 2^P. That class has an odd member below 2^N iff
  bits N..P-1 of nu are zero, i.e. iff P <= z_N(nu) := lowest set bit of nu at position >= N.
* So the last rise can be extended to j* = z_N(nu) - P, and
      B_r(N) = max over feasible sigma of (k_0 + ... + k_{r-2}) + z_N(nu) - P.
  The DFS stops each k- and l-loop as soon as the bare prefix becomes infeasible.

Validation: brute force for N <= 18 (all r), and literal re-simulation of every reported
maximiser (if the bound is attained, an actual integer n < 2^N attains it).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path


def lowest_set_bit_at_or_above(x: int, N: int, W: int) -> int | None:
    y = (x % (1 << W)) >> N
    if y == 0:
        return None
    return N + ((y & -y).bit_length() - 1)


class Solver:
    def __init__(self, N: int, extra: int = 256):
        self.N = N
        self.W = N + extra
        M = 1 << self.W
        inv3 = pow(3, -1, M)
        # u(K) = 2^K * 3^-K mod 2^W, for K = 0..W (u(K) = 0 mod 2^W once K >= W)
        self.u = [0] * (self.W + 1)
        p = 1
        for K in range(self.W + 1):
            self.u[K] = ((p << K) % M) if K < self.W else 0
            p = p * inv3 % M
        self.M = M

    def uK(self, K: int) -> int:
        return self.u[K] if K <= self.W else 0

    def z(self, nu: int) -> int:
        zz = lowest_set_bit_at_or_above(nu, self.N, self.W)
        if zz is None:
            raise OverflowError("increase W")
        return zz

    def B(self, r: int) -> tuple[int, tuple]:
        """(B_r(N), maximising sigma). r >= 1."""
        N, M = self.N, self.M
        best = [-1, None]

        def rec(nu: int, P: int, Ksum: int, Lsum: int, depth: int, sigma: tuple):
            zz = self.z(nu)
            if depth == r - 1:
                if zz > P:
                    tot = Ksum + zz - P
                    if tot > best[0]:
                        best[0], best[1] = tot, sigma
                return
            # next pair (k, l): bare prefix ...1^k feasible iff P + k <= zz
            for k in range(1, zz - P + 1):
                K2 = Ksum + k
                base = self.uK(K2)
                l = 1
                while True:
                    P2 = P + k + l
                    term = ((base * (1 - (1 << l))) << Lsum) % M
                    nu2 = (nu + term) % M
                    z2 = self.z(nu2)
                    if P2 > z2:            # bare prefix ...1^k 0^l infeasible; so is every longer l
                        break
                    rec(nu2, P2, K2, Lsum + l, depth + 1, sigma + (k, l))
                    l += 1

        rec(M - 1, 0, 0, 0, 0, ())
        return best[0], best[1]


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def first_rises(n: int, r: int) -> list[int]:
    out = []
    for _ in range(r):
        k = 0
        while n % 2 == 1:
            n = T(n)
            k += 1
        out.append(k)
        while n % 2 == 0:
            n = T(n)
    return out


def brute(N: int, r: int) -> int:
    return max(sum(first_rises(n, r)) for n in range(1, 1 << N, 2))


def realise(N: int, r: int, sigma: tuple, total: int) -> int:
    """Return an odd n < 2^N with parity prefix sigma + 1^j*, and check its first r rises sum to >= total."""
    S = Solver(N)
    nu = S.M - 1
    P = Ksum = Lsum = 0
    for i in range(0, len(sigma), 2):
        k, l = sigma[i], sigma[i + 1]
        Ksum += k
        nu = (nu + ((S.uK(Ksum) * (1 - (1 << l))) << Lsum)) % S.M
        Lsum += l
        P += k + l
    zz = S.z(nu)
    n = nu % (1 << zz)                       # least representative of the class mod 2^zz
    assert n % 2 == 1 and n < (1 << N), (n, N)
    got = sum(first_rises(n, r))
    assert got >= total, (got, total, sigma)
    return n


def load_marshall(root: Path):
    import glob
    def load(fs, key):
        g = {}
        for f in fs:
            for k, v in json.load(open(f)).items():
                g[int(k)] = v[key] if isinstance(v, dict) else v
        return g
    c = root / "computed"
    g2 = load([c / "target1a/chain/grid_own.json", c / "target7_threerise/b2_extra.json"]
              + sorted(glob.glob(str(c / "target7b_push/b2_new*.json"))) + [c / "target13_m128/b2_new_g13.json"], "B2")
    g3 = load(sorted(glob.glob(str(c / "target7_threerise/grid3*.json")))
              + sorted(glob.glob(str(c / "target7b_push/b3_new*.json"))) + [c / "target13_m128/b3_new_g13.json"], "B3")
    g4 = load([c / "target11_fourrise/b4_grid_snapshot_le852.json"], "B4")
    return {2: g2, 3: g3, 4: g4}


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["selftest", "compare"])
    ap.add_argument("--root", help="Marshall repo checkout (for compare)")
    ap.add_argument("--r", type=int, default=2)
    ap.add_argument("--nmax", type=int, default=400)
    ap.add_argument("--points", type=int, nargs="*", help="check only these grid points")
    a = ap.parse_args()
    if a.cmd == "selftest":
        for r in (1, 2, 3, 4):
            for N in range(2, 19 if r <= 3 else 17):
                b, sig = Solver(N).B(r)
                bb = brute(N, r)
                assert b == bb, (r, N, b, bb)
                if sig is not None:
                    realise(N, r, sig, b)
            print(f"ok  r={r}: DFS == brute force for N = 2..{18 if r <= 3 else 16}, maximisers re-simulated")
        return
    tabs = load_marshall(Path(a.root))
    g = tabs[a.r]
    pts = sorted(a.points) if a.points else sorted(n for n in g if n <= a.nmax)
    bad = 0
    t0 = time.time()
    for N in pts:
        t = time.time()
        b, sig = Solver(N).B(a.r)
        n = realise(N, a.r, sig, b)
        flag = "OK" if b <= g[N] else "TABLE TOO SMALL"
        if b != g[N]:
            flag += f" (table {g[N]}, exact {b})"
        bad += b > g[N]
        print(f"B{a.r}({N}) = {b}  table {g[N]}  {flag}  witness n has {n.bit_length()} bits, "
              f"runs {sig}  [{time.time()-t:.1f}s]", flush=True)
    print(f"{len(pts)} points, {bad} unsafe, {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
