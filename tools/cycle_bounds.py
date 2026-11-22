#!/usr/bin/env python3
"""Lower bounds on Collatz cycles from a verification bound X0 (idea I-002).

Re-implements the continued-fraction machinery of
  C. Hercher, "There are no Collatz m-Cycles with m <= 91", JIS 26 (2023), arXiv:2201.00406
so we can rerun it with newer verification bounds (Barina 2025: all n < 2^71;
reported exactly as 2075 * 2^60).

Notation (Hercher): a nontrivial cycle has K odd and L even elements and m local
minima. delta = log2(3). For every such cycle

    delta < (K+L)/K < delta + eps,

where eps comes from bounds on sum(1/x) over the cycle's odd elements. The smallest
denominator of a rational in (delta, delta+eps) is a lower bound on K. The m-cycle
argument raises this lower bound until it passes the Simons-de Weger upper bound
K < 1.4784 * m * delta^m.

The arithmetic is decimal at PREC digits. Every rounding goes in the direction that
makes the K lower bound *smaller* (conservative): intervals are widened, eps is
rounded up, and delta is replaced by its rational enclosure.

Usage:
    python3 cycle_bounds.py selftest
    python3 cycle_bounds.py allm  --x0 2075        # m-independent bound (X0 = 2075*2^60)
    python3 cycle_bounds.py mcycles --x0 704 --mmax 100 [--method hercher|direct]
"""
from __future__ import annotations

import argparse
from decimal import Decimal as D, getcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction as F
from math import floor

PREC = 150
_ctx = getcontext()
_ctx.prec = PREC
_ctx.Emax = 10 ** 9
_ctx.Emin = -(10 ** 9)
V_CAP = D(10) ** 6                          # 2^-1e6 is far below anything that matters

LN2 = D(2).ln()
DELTA = D(3).ln() / LN2                     # log2(3), PREC digits
_ULP = D(10) ** (-(PREC - 10))              # generous error bar on DELTA
DELTA_LO = F(DELTA - _ULP)
DELTA_HI = F(DELTA + _ULP)
TWO60 = 2 ** 60
SAFETY = D(1) + D(10) ** -40                # inflate eps a little: conservative


# --------------------------------------------------------------------------
# Continued fractions
# --------------------------------------------------------------------------

def cf_terms(x: F, n: int) -> list[int]:
    out = []
    for _ in range(n):
        a = floor(x)
        out.append(a)
        if x == a:
            break
        x = 1 / (x - a)
    return out


def simplest_in_open(a: F, b: F) -> F:
    """The rational with smallest denominator in the open interval (a, b), 0 < a < b."""
    assert 0 < a < b
    fl = floor(a)
    if fl + 1 < b:
        return F(fl + 1)
    if a == fl:
        n = floor(1 / (b - fl)) + 1          # fl + 1/n < b
        return fl + F(1, n)
    r = simplest_in_open(1 / (b - fl), 1 / (a - fl))
    return fl + 1 / r


def min_K_interval(eps: D) -> int:
    """Lower bound on K from delta < (K+L)/K < delta + eps (Hercher's Lemma 22 step).

    Searches the widened interval (DELTA_LO, DELTA_HI + eps), so the answer can only be
    at or below the true smallest denominator.
    """
    upper = DELTA_HI + F(eps * SAFETY)
    return simplest_in_open(DELTA_LO, upper).denominator


def _semiconvergents_above(max_q: int):
    """Rationals p/q > delta that are one-sided best approximations from above,
    in increasing q (a superset: all semiconvergents above delta)."""
    terms = cf_terms(F(DELTA), 120)
    p_prev, q_prev, p, q = 1, 0, terms[0], 1
    cands = []
    for a in terms[1:]:
        for j in range(1, a + 1):
            pj, qj = p_prev + j * p, q_prev + j * q
            if qj > max_q:
                break
            cands.append((qj, pj))
        p_prev, q_prev, p, q = p, q, a * p + p_prev, a * q + q_prev
        if q > max_q:
            break
    cands = sorted(set(cands))
    return [(qq, pp) for qq, pp in cands if F(pp, qq) > DELTA_HI]


def min_K_direct(eta: D, max_q: int = 10 ** 60) -> int:
    """Smallest K with 0 < ceil(K delta) - K delta < eta (a stronger, K-aware version).

    Uses N - K*delta < eta directly instead of an interval fixed at the previous
    lower bound. The first K that works is a one-sided best approximation, so it
    appears among the semiconvergents above delta.
    """
    eta_up = F(eta * SAFETY)
    for q, p in _semiconvergents_above(max_q):
        gap_lo = p - q * DELTA_HI           # <= true gap
        if gap_lo < eta_up:
            return q
    raise ValueError("max_q too small")


def min_K_direct_bruteforce(eta: float, limit: int) -> int | None:
    import math
    d = math.log2(3)
    for K in range(1, limit):
        g = math.ceil(K * d) - K * d
        if 0 < g < eta:
            return K
    return None


# --------------------------------------------------------------------------
# Hercher's bounds
# --------------------------------------------------------------------------

def eps_all_m(X0: D) -> D:
    """Theorem 27: (K+L)/K - delta < (1/(3 ln 2)) * (3/4) / X0, for every cycle."""
    return (D(3) / 4) / (3 * LN2 * X0)


def eps_breeden(X0: D) -> D:
    """The simpler Sinisalo-style bound log2(3 + 1/X0) - log2(3)."""
    return (D(3) + 1 / X0).ln() / LN2 - DELTA


def bracket_theorem21(m: int, m2: int, K: D, X0: D) -> D | None:
    """Theorem 21: returns S with (K+L)/K - delta < S / (3 K ln 2), or None if m2 is not allowed.

    Condition: (delta^m2 - 1)/(delta - 1) * log2(162/97 * X0) <= (m2/m) * K
    v = (m2/m) * K * (delta - 1)/(delta^m2 - 1)
    S = (97(m-m2)+73)/(54 X0) + 3/(2^v - 1) + 3(m2-1)/(2^v - 1)^delta
    Special cases: m2 = m-1: first term 3/X0 and last term 3(m-2)/(2^v-1)^delta.
                   m2 = m:   first term 0.
    """
    geo = (DELTA ** m2 - 1) / (DELTA - 1)
    if geo * ((D(162) / 97 * X0).ln() / LN2) > D(m2) / m * K:
        return None
    v = D(m2) / m * K / geo
    # The terms decrease in v, so capping v only weakens (never invalidates) the bound.
    v = min(v, V_CAP)
    t = D(2) ** v - 1
    if m2 == m:
        first, last = D(0), 3 * D(m - 1) / t ** DELTA
    elif m2 == m - 1:
        first, last = 3 / X0, 3 * D(m - 2) / t ** DELTA
    else:
        first, last = D(97 * (m - m2) + 73) / (54 * X0), 3 * D(m2 - 1) / t ** DELTA
    return first + 3 / t + last


def bracket_cor17(m: int, X0: D) -> D:
    """Corollary 17: S = (97/54) * m / X0."""
    return D(97) / 54 * m / X0


def best_bracket(m: int, K: int, X0: D) -> tuple[D, int]:
    best, best_m2 = bracket_cor17(m, X0), 0
    for m2 in range(1, m + 1):
        s = bracket_theorem21(m, m2, D(K), X0)
        if s is not None and s < best:
            best, best_m2 = s, m2
    return best, best_m2


def sdw_upper(m: int) -> D:
    """Simons-de Weger upper bound on K for an m-cycle, as quoted by Hercher."""
    return D("1.4784") * m * DELTA ** m


def mcycle_chain(m: int, X0: D, K0: int, method: str = "hercher", max_iter: int = 200):
    """Iterate Theorem 21 + continued fractions. Returns (excluded?, chain).

    chain = list of (K_lower, m2, eps). Excluded when K_lower exceeds sdw_upper(m).
    """
    K = K0
    chain = []
    upper = sdw_upper(m)
    for _ in range(max_iter):
        if D(K) > upper:
            return True, chain
        S, m2 = best_bracket(m, K, X0)
        eps = S / (3 * LN2 * K)
        if method == "hercher":
            K_new = min_K_interval(eps)
        else:
            K_new = min_K_direct(S / (3 * LN2))
        chain.append((K_new, m2, eps))
        if K_new <= K:
            return D(K) > upper, chain
        K = K_new
    return D(K) > upper, chain


def convergents(n: int = 80) -> list[tuple[int, int]]:
    terms = cf_terms(F(DELTA), n)
    out = [(terms[0], 1)]
    p0, q0, p1, q1 = 1, 0, terms[0], 1
    for a in terms[1:]:
        p0, q0, p1, q1 = p1, q1, a * p1 + p0, a * q1 + q0
        out.append((p1, q1))
    return out


def surviving_K(m: int, X0: D, K0: int = 7 * 10 ** 11) -> tuple[int, int, list[tuple[int, int, float, int]]]:
    """Every K in [iterated lower bound, Simons-de Weger upper bound] that passes
    gap(K) = ceil(K delta) - K delta < eta(K) = min_m2 S(m2, K) / (3 ln 2).

    Exhaustive. (K, N) is written in the unimodular basis of two consecutive convergents,
    K = i q_a + j q_b, so the conditions cut out a thin parallelogram with few lattice points.
    Bounds: |j| <= K |e_a| + e q_a, and i lies in a window of width eta/|e_a| for each j.
    Returns (K_lower, K_upper, [(K, L, gap, best m2)]).
    """
    U = int(sdw_upper(m))
    ok, chain = mcycle_chain(m, X0, K0)
    if ok:
        return U + 1, U, []
    Klo = chain[-1][0]
    conv = convergents()
    k = max(i for i, (_, q) in enumerate(conv) if q <= U)
    (pa, qa), (pb, qb) = conv[k], conv[k + 1]
    d = F(DELTA)
    ea, eb = pa - qa * d, pb - qb * d
    S0, _ = best_bracket(m, Klo, X0)
    eta_max = F(S0 / (3 * LN2)) * 2                 # eta(K) is nonincreasing in K
    jmax = int(U * abs(ea) + eta_max * qa) + 3
    out = set()
    for j in range(-jmax, jmax + 1):
        lo, hi = sorted(((-j * eb) / ea, (eta_max - j * eb) / ea))
        for i in range(floor(lo) - 1, floor(hi) + 3):
            K, N = i * qa + j * qb, i * pa + j * pb
            if not (Klo <= K <= U) or N - K * DELTA_LO <= 0:
                continue
            S, m2 = best_bracket(m, K, X0)
            if N - K * DELTA_HI < F(S * SAFETY / (3 * LN2)):
                out.add((K, N - K, float(N - K * d), m2))
    return Klo, U, sorted(out)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def x0_from(arg: str) -> D:
    """'704' means 704 * 2^60; '2^71' means 2^71."""
    if arg.startswith("2^"):
        return D(2) ** int(arg[2:])
    return D(arg) * TWO60


def cmd_selftest(_):
    # 1. Continued fraction of log2(3) (OEIS A028507).
    known = [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 15, 1, 9, 2, 5, 7, 1, 1, 4, 8]
    got = cf_terms(F(DELTA), len(known))
    assert got == known, got
    print("ok  continued fraction of log2(3) matches A028507 (30 terms)")

    # 2. Interval search against brute force on small cases.
    import random
    random.seed(1)
    for _ in range(300):
        a = F(random.randint(1, 10 ** 6), random.randint(1, 10 ** 5))
        b = a + F(1, random.randint(2, 10 ** 7))
        r = simplest_in_open(a, b)
        assert a < r < b
        for q in range(1, r.denominator):
            assert not any(a < F(p, q) < b for p in range(floor(a * q), floor(b * q) + 2)), (a, b, q)
    print("ok  simplest_in_open agrees with brute force (300 random intervals)")

    # 3. Direct K-search against brute force.
    for eta in [0.3, 0.1, 0.03, 0.01, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]:
        assert min_K_direct(D(eta)) == min_K_direct_bruteforce(eta, 10 ** 6), eta
    print("ok  min_K_direct agrees with brute force (eta down to 1e-5)")

    # 4. Published numbers.
    b = min_K_interval(eps_breeden(x0_from("2075")))
    print(f"    Sinisalo-style bound at 2075*2^60: K >= {b:,}  (Breeden repo: 72,057,431,991)")
    assert b == 72_057_431_991
    lo = min_K_interval(eps_all_m(x0_from("2835")))
    hi = min_K_interval(eps_all_m(x0_from("2836")))
    print(f"    Theorem 27 at 2835*2^60: K >= {lo:,};  at 2836*2^60: K >= {hi:,}")
    print("    (Hercher Remark 28: Theorem 27 reaches the next level once X0 >= 2836*2^60)")
    assert lo < 137_528_045_312 <= hi
    # Hercher step 1 of Theorem 23: K > 7e11, X0 = 704*2^60, m = 91 -> m2 = 47, eps ~ 6.9e-32.
    S, m2 = best_bracket(91, 7 * 10 ** 11, x0_from("704"))
    eps = S / (3 * LN2 * 7 * 10 ** 11)
    print(f"    Theorem 23 step 1: m2 = {m2}, eps = {float(eps):.3g}  (paper: m2 = 47, eps = 6.9e-32)")
    assert m2 == 47 and 6.5e-32 < eps < 7.0e-32
    K1 = min_K_interval(eps)
    print(f"    Theorem 23 step 1: K > {float(K1):.3g}  (paper: 5.2e15)")
    print("ok  published numbers reproduced")


def cmd_allm(args):
    X0 = x0_from(args.x0)
    for name, eps in [("Sinisalo-style log2(3+1/X0)", eps_breeden(X0)),
                      ("Hercher Thm 27 (3/4 factor)", eps_all_m(X0))]:
        print(f"{name:32s}: every nontrivial cycle has K >= {min_K_interval(eps):,} odd elements")


def cmd_mcycles(args):
    X0 = x0_from(args.x0)
    K0 = args.k0 if args.k0 else min_K_interval(eps_all_m(X0))
    print(f"X0 = {args.x0} ({float(X0):.4g}), method = {args.method}, starting K >= {K0:,}\n")
    last_ok = 0
    for m in range(1, args.mmax + 1):
        ok, chain = mcycle_chain(m, X0, K0, args.method)
        if args.verbose or (m >= args.mmax - 15):
            steps = " → ".join(f"{float(k):.3g}(m2={m2})" for k, m2, _ in chain[-4:])
            print(f"m={m:3d}  {'EXCLUDED' if ok else 'not excluded':12s} upper {float(sdw_upper(m)):.3g}  chain …{steps}")
        if ok and last_ok == m - 1:
            last_ok = m
    print(f"\nNo m-cycles for all m <= {last_ok}  (given X0 and the quoted Simons-de Weger bound)")


def cmd_survivors(args):
    X0 = x0_from(args.x0)
    conv_q = {q for _, q in convergents()}
    base = None
    for m in range(args.mmin, args.mmax + 1):
        Klo, U, surv = surviving_K(m, X0)
        if not surv:
            print(f"m={m:3d}: excluded (no K survives)")
            continue
        base = base or surv[0][0]
        desc = []
        for K, L, gap, m2 in surv:
            tag = f"{K // base}·K*" if K % base == 0 else f"{K:.4g}"
            desc.append(tag + (" [convergent]" if K in conv_q else ""))
        print(f"m={m:3d}: K in [{Klo:.4g}, {U:.4g}] → {len(surv)} survivor(s): " + ", ".join(desc))
    if base:
        K, L, gap, m2 = next(s for s in surviving_K(args.mmin, X0)[2] if s[0] == base)
        print(f"\nK* = {K:,}   L* = {L:,}   gap(K*) = {gap:.4e}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(func=cmd_selftest)
    a = sub.add_parser("allm")
    a.add_argument("--x0", default="2075")
    a.set_defaults(func=cmd_allm)
    mc = sub.add_parser("mcycles")
    mc.add_argument("--x0", default="704")
    mc.add_argument("--mmax", type=int, default=100)
    mc.add_argument("--k0", type=int, default=0)
    mc.add_argument("--method", choices=["hercher", "direct"], default="hercher")
    mc.add_argument("-v", "--verbose", action="store_true")
    mc.set_defaults(func=cmd_mcycles)
    sv = sub.add_parser("survivors", help="every K still allowed for an m-cycle")
    sv.add_argument("--x0", default="2075")
    sv.add_argument("--mmin", type=int, default=92)
    sv.add_argument("--mmax", type=int, default=100)
    sv.set_defaults(func=cmd_survivors)
    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
