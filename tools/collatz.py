#!/usr/bin/env python3
"""Collatz toolkit: orbits, parity vectors, cycle search, and control maps.

Everything is parameterized by a generalized accelerated map

    T_{a,b}(n) = n / 2          if n is even
               = (a*n + b) / 2  if n is odd        (a, b odd)

so T_{3,1} is the usual accelerated Collatz map, and T_{5,1} and T_{3,-1}
are the standard controls. 3x+1 on negative integers is T_{3,1} with n < 0
(equivalent to 3x-1 on positives under n -> -n).

Usage:
    python3 collatz.py orbit 27 [--a 3 --b 1] [--standard]
    python3 collatz.py stats 1000000
    python3 collatz.py cycles --kmax 18 [--a 5 --b 1]
    python3 collatz.py bijection 12
    python3 collatz.py controls
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from fractions import Fraction
from itertools import product


# --------------------------------------------------------------------------
# Maps
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Map:
    a: int = 3
    b: int = 1

    def __post_init__(self):
        if self.a % 2 == 0 or self.b % 2 == 0:
            raise ValueError("a and b must be odd so odd steps land on integers")

    def __call__(self, n: int) -> int:
        return n // 2 if n % 2 == 0 else (self.a * n + self.b) // 2

    def __str__(self) -> str:
        sign = "+" if self.b >= 0 else "-"
        return f"T({self.a}x{sign}{abs(self.b)})"


COLLATZ = Map(3, 1)


def standard_step(n: int) -> int:
    """The unaccelerated map C(n) = n/2 or 3n+1."""
    return n // 2 if n % 2 == 0 else 3 * n + 1


# --------------------------------------------------------------------------
# Orbits
# --------------------------------------------------------------------------

@dataclass
class Orbit:
    start: int
    path: list[int]          # path[0] = start; ends at first repeated value or on bailout
    cycle: list[int] | None  # the cycle reached, if any
    bailed: bool             # True if max_steps or max_abs was hit

    @property
    def steps_to_cycle(self) -> int | None:
        if self.cycle is None:
            return None
        return self.path.index(self.cycle[0])

    @property
    def peak(self) -> int:
        return max(self.path, key=abs)


def orbit(n: int, f=COLLATZ, max_steps: int = 10**6, max_abs: int | None = None) -> Orbit:
    """Iterate f from n until a value repeats (cycle found) or a bound is hit."""
    seen: dict[int, int] = {}
    path: list[int] = []
    x = n
    for _ in range(max_steps):
        if x in seen:
            i = seen[x]
            return Orbit(n, path, path[i:], False)
        seen[x] = len(path)
        path.append(x)
        if max_abs is not None and abs(x) > max_abs:
            return Orbit(n, path, None, True)
        x = f(x)
    return Orbit(n, path, None, True)


def canonical_cycle(cyc: list[int]) -> tuple[int, ...]:
    """Rotate a cycle so it starts at its smallest absolute value (ties: smallest value)."""
    i = min(range(len(cyc)), key=lambda j: (abs(cyc[j]), cyc[j]))
    return tuple(cyc[i:] + cyc[:i])


def stopping_time(n: int, f=COLLATZ, max_steps: int = 10**5) -> int | None:
    """Smallest k >= 1 with T^k(n) < n (Terras). None if not reached."""
    x = n
    for k in range(1, max_steps + 1):
        x = f(x)
        if x < n:
            return k
    return None


def total_stopping_time(n: int, f=COLLATZ, target: int = 1, max_steps: int = 10**6) -> int | None:
    """Number of steps to reach `target`. None if not reached within max_steps."""
    x, k = n, 0
    while x != target:
        if k >= max_steps:
            return None
        x = f(x)
        k += 1
    return k


# --------------------------------------------------------------------------
# Parity vectors and the mod 2^k bijection
# --------------------------------------------------------------------------

def parity_vector(n: int, k: int, f=COLLATZ) -> tuple[int, ...]:
    """The first k parities n, T(n), ..., T^{k-1}(n) (1 = odd)."""
    v = []
    x = n
    for _ in range(k):
        v.append(x & 1)
        x = f(x)
    return tuple(v)


def affine_after(v: tuple[int, ...], f: Map = COLLATZ) -> tuple[int, int, int]:
    """Return (A, C, D) with T^k(n) = (A*n + C) / D for every n with parity vector v.

    A = a^m, D = 2^k, where m = number of 1s in v.
    """
    A, C, D = 1, 0, 1
    for bit in v:
        if bit:
            A *= f.a
            C = f.a * C + f.b * D
        D *= 2
    return A, C, D


def residue_for_parity(v: tuple[int, ...], f: Map = COLLATZ) -> int:
    """The unique residue r mod 2^k whose first k parities are v.

    Built one bit at a time: the parity of T^i(n) depends on n mod 2^{i+1},
    and flipping bit i of n flips that parity (because a is odd).
    """
    r = 0
    for i, bit in enumerate(v):
        if parity_vector(r, i + 1, f)[i] != bit:
            r += 1 << i
    return r


def check_bijection(k: int, f: Map = COLLATZ) -> bool:
    """Check that n mod 2^k -> (first k parities) is a bijection onto {0,1}^k."""
    seen = {parity_vector(r, k, f) for r in range(1 << k)}
    if len(seen) != 1 << k:
        return False
    return all(residue_for_parity(v, f) == r
               for r in range(1 << k)
               for v in [parity_vector(r, k, f)])


# --------------------------------------------------------------------------
# Cycles (Böhm–Sontacchi)
# --------------------------------------------------------------------------

def rational_fixed_point(v: tuple[int, ...], f: Map = COLLATZ) -> Fraction | None:
    """The rational n with T_v(n) = n, where T_v is the affine map for parity vector v.

    n = C / (D - A). Integer solutions with the right parity vector are integer cycles.
    """
    A, C, D = affine_after(v, f)
    if D == A:
        return None
    return Fraction(C, D - A)


def find_cycles(kmax: int, f: Map = COLLATZ) -> dict[tuple[int, ...], int]:
    """Every integer cycle of f whose length (in T-steps) is at most kmax.

    Exhaustive: each cycle of length k has a periodic parity vector of length k,
    and its elements are integer fixed points of the matching affine map.
    Returns {canonical cycle: length}.
    """
    found: dict[tuple[int, ...], int] = {}
    for k in range(1, kmax + 1):
        for v in product((0, 1), repeat=k):
            q = rational_fixed_point(v, f)
            if q is None or q.denominator != 1:
                continue
            n = q.numerator
            # Check it really is a cycle with this parity vector.
            if parity_vector(n, k, f) != v:
                continue
            o = orbit(n, f, max_steps=k + 1)
            if o.cycle is None:
                continue
            cyc = canonical_cycle(o.cycle)
            found.setdefault(cyc, len(cyc))
    return found


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------

def log_drift_stats(N: int, k: int = 40, f: Map = COLLATZ) -> dict[str, float]:
    """Average of log(T^k(n)/n)/k over N/2 <= n < N, against the coin-flip prediction.

    Under the heuristic, each T-step multiplies by 1/2 or ~a/2 with probability 1/2,
    so the drift per step is (1/2) log(a/4).
    """
    total, count = 0.0, 0
    for n in range(N // 2, N):
        x = n
        for _ in range(k):
            x = f(x)
            if x == 0:
                break
        if x > 0:
            total += math.log(x / n) / k
            count += 1
    return {
        "empirical_drift": total / count,
        "heuristic_drift": 0.5 * math.log(f.a / 4),
        "samples": count,
    }


def stopping_time_stats(N: int) -> dict:
    """Records for the standard map C on 1 <= n < N."""
    longest, longest_n = 0, 1
    highest, highest_n = 1, 1
    for n in range(1, N):
        x, steps, peak = n, 0, n
        while x != 1:
            x = standard_step(x)
            steps += 1
            if x > peak:
                peak = x
        if steps > longest:
            longest, longest_n = steps, n
        if peak > highest:
            highest, highest_n = peak, n
    return {"N": N, "longest_total_stopping_time": (longest_n, longest),
            "highest_peak": (highest_n, highest)}


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _cmd_orbit(args):
    if args.standard:
        x, path = args.n, [args.n]
        while x != 1 and len(path) < 10**6:
            x = standard_step(x)
            path.append(x)
        print(" → ".join(map(str, path)))
        print(f"\nsteps: {len(path) - 1}, peak: {max(path)}")
        return
    f = Map(args.a, args.b)
    o = orbit(args.n, f, max_abs=args.max_abs)
    print(" → ".join(map(str, o.path[:200])) + (" …" if len(o.path) > 200 else ""))
    if o.cycle is not None:
        print(f"\n{f}: reached cycle {canonical_cycle(o.cycle)} after {o.steps_to_cycle} steps; peak {o.peak}")
    else:
        print(f"\n{f}: no cycle within bounds (len {len(o.path)}, last |x| ≈ 2^{abs(o.path[-1]).bit_length()})")


def _cmd_cycles(args):
    f = Map(args.a, args.b)
    cycles = find_cycles(args.kmax, f)
    print(f"{f}: integer cycles with length ≤ {args.kmax} T-steps")
    for cyc, length in sorted(cycles.items(), key=lambda kv: (kv[1], abs(kv[0][0]))):
        side = "positive" if cyc[0] > 0 else ("zero" if cyc[0] == 0 else "negative")
        print(f"  length {length:3d}  [{side}]  {cyc}")


def _cmd_bijection(args):
    for k in range(1, args.k + 1):
        ok3 = check_bijection(k, COLLATZ)
        ok5 = check_bijection(k, Map(5, 1))
        print(f"k={k:2d}: 3x+1 bijection {'OK' if ok3 else 'FAIL'}, 5x+1 bijection {'OK' if ok5 else 'FAIL'}")
    print("\nMod 2^k data cannot tell 3x+1 apart from 5x+1: both give every parity pattern exactly once.")


def _cmd_stats(args):
    s = stopping_time_stats(args.N)
    print(f"standard map, n < {s['N']}:")
    print(f"  longest total stopping time: n={s['longest_total_stopping_time'][0]} "
          f"({s['longest_total_stopping_time'][1]} steps)")
    print(f"  highest peak: n={s['highest_peak'][0]} reaches {s['highest_peak'][1]}")
    for f in (COLLATZ, Map(5, 1)):
        d = log_drift_stats(min(args.N, 200_000), k=30, f=f)
        print(f"{f}: empirical log-drift/step {d['empirical_drift']:+.4f}, "
              f"coin-flip heuristic {d['heuristic_drift']:+.4f}  ({d['samples']} samples)")


def _cmd_controls(args):
    k = args.kmax
    print(f"Control systems: every integer cycle with ≤ {k} T-steps (exhaustive search)\n")
    for f, label in [(COLLATZ, "3x+1 (positives and negatives)"),
                     (Map(5, 1), "5x+1"),
                     (Map(3, -1), "3x-1")]:
        cycles = find_cycles(k, f)
        print(f"{label}:")
        for cyc, length in sorted(cycles.items(), key=lambda kv: (kv[1], abs(kv[0][0]))):
            if cyc == (0,):
                continue
            shown = cyc if len(cyc) <= 12 else cyc[:12] + ("…",)
            print(f"  length {length:3d}: {shown}")
        print()
    print("Any argument that rules out the 5x+1 cycles, or the negative 3x+1 cycles, is wrong.")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    o = sub.add_parser("orbit", help="print an orbit")
    o.add_argument("n", type=int)
    o.add_argument("--a", type=int, default=3)
    o.add_argument("--b", type=int, default=1)
    o.add_argument("--standard", action="store_true", help="use C(n)=n/2, 3n+1 instead of T")
    o.add_argument("--max-abs", type=int, default=10**60)
    o.set_defaults(func=_cmd_orbit)

    c = sub.add_parser("cycles", help="exhaustive cycle search by parity vector")
    c.add_argument("--kmax", type=int, default=16)
    c.add_argument("--a", type=int, default=3)
    c.add_argument("--b", type=int, default=1)
    c.set_defaults(func=_cmd_cycles)

    b = sub.add_parser("bijection", help="check the mod 2^k parity-vector bijection")
    b.add_argument("k", type=int, nargs="?", default=12)
    b.set_defaults(func=_cmd_bijection)

    s = sub.add_parser("stats", help="stopping-time records and drift")
    s.add_argument("N", type=int, nargs="?", default=10**6)
    s.set_defaults(func=_cmd_stats)

    k = sub.add_parser("controls", help="cycles of the control systems")
    k.add_argument("--kmax", type=int, default=16)
    k.set_defaults(func=_cmd_controls)

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
