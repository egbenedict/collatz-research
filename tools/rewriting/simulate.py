"""Control: the rewriting systems faithfully simulate their maps.

1. Symbolic check: every auxiliary rule preserves value (Gamma_l == Gamma_r) and every
   dynamic rule maps Gamma_l(x) to MAP(Gamma_l(x)), for x = 0..200.
2. Dynamic check: from the canonical (binary) string of every n in [2, N], rewrite with a
   random strategy until a normal form; the values seen at dynamic steps must be exactly the
   map's trajectory, and the normal form must be a sink value.
3. 5x+1 control: from 13, the T5 rewriting never terminates (value returns to 13 repeatedly).

Usage: .venv/bin/python simulate.py [N]
"""
import random
import sys

from srs import (SYSTEMS, X, X5, D_T, D_F, D_S, D_EVEN8, D_ODD8, D_EVEN4, D_ODD4, D_T5, D_W,
                 gamma, value)


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def C(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


def T5(n):
    return n // 2 if n % 2 == 0 else (5 * n + 1) // 2


def F(n):
    if n % 3 == 1:
        return (n - 1) // 3
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def S(n):  # odd n
    if n % 8 == 1:
        return (3 * n + 1) // 4
    if n % 8 == 5:
        return (n - 1) // 4
    return (3 * n + 1) // 2


def W(n):  # W': f> -> 0> : 2x -> 3x on evens, odd fixed
    return 3 * n // 2 if n % 2 == 0 else n


def check_rules(aux, dyn, fmap, base=3, end=">"):
    for l, r in aux:
        for x in range(201):
            assert gamma(l, x, base) == gamma(r, x, base), (l, r, x)
    for l, r in dyn:
        for x in range(201):
            gl, gr = gamma(l, x, base), gamma(r, x, base)
            if end == "*":  # S works on odd numbers; '*' already included in gamma
                pass
            assert fmap(gl) == gr, (l, r, x, gl, gr)


def canonical(n, end=">"):
    if end == "*":
        m = (n - 1) // 2            # n = 2m + 1, m >= 1
        bits = bin(m)[3:]
    else:
        bits = bin(n)[3:]
    return "<" + bits.replace("0", "f").replace("1", "t") + end


def ternary_start(n):
    """F starts in ternary (YAH Sec 4.2): leading 1 -> '<', leading 2 -> '<f'."""
    d = []
    while n:
        d.append(str(n % 3))
        n //= 3
    d = d[::-1]
    return ("<" if d[0] == "1" else "<f") + "".join(d[1:]) + ">"


def run(s, rules, dyn_lhs, base, rng, max_steps=10**6, rightmost=False):
    """Random-strategy rewriting. Returns (list of values after each dynamic step, normal form)."""
    vals = []
    for _ in range(max_steps):
        redexes = []
        for l, r in rules:
            i = s.find(l)
            while i != -1:
                redexes.append((i, l, r))
                i = s.find(l, i + 1)
        if not redexes:
            return vals, s
        i, l, r = max(redexes) if rightmost else rng.choice(redexes)
        s = s[:i] + r + s[i + len(l):]
        if l in dyn_lhs:
            vals.append(value(s, base))
    return vals, None


def check_dynamic(name, dyn, fmap, sinks, N, end=">", base=3, seed=1, start=None,
                  rightmost=False):
    rules = SYSTEMS[name]
    dyn_lhs = {l for l, _ in dyn}
    rng = random.Random(seed)
    starts = range(3, N + 1, 2) if end == "*" else range(2, N + 1)
    seen = set()
    for n in starts:
        s = start(n) if start else canonical(n, end)
        assert value(s, base) == n
        vals, nf = run(s, rules, dyn_lhs, base, rng, rightmost=rightmost)
        assert nf is not None, (name, n)
        x = n
        for v in vals:
            x = fmap(x)
            assert v == x, (name, n, vals)
        assert sinks is None or value(nf, base) in sinks, (name, n, nf)
        seen.add(value(nf, base))
    print(f"  {name}: simulates its map faithfully on n <= {N} ({'rightmost' if rightmost else 'random'} strategy); "
          f"normal-form values {sorted(seen) if len(seen) < 12 else str(len(seen)) + ' distinct (H undefined there)'}")


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    check_rules(X, D_T, T)
    check_rules(X, D_W, W)
    check_rules(X, D_F, F)
    check_rules(X, D_S, S, end="*")
    check_rules(X, D_EVEN8 + D_EVEN4, C)
    check_rules(X, D_ODD8 + D_ODD4, T)
    check_rules(X5, D_T5, T5, base=5)
    print("rule semantics: all auxiliary rules value-preserving, all dynamic rules apply the map")
    check_dynamic("T", D_T, T, {1}, N)
    # F: pure-binary strings are normal forms, so start in ternary and rewrite rightmost-first
    # (pure-binary strings ARE normal forms of F, so the sink set is not checked, only that
    # every dynamic step applies F)
    check_dynamic("F", D_F, F, None, N, start=ternary_start, rightmost=True)
    check_dynamic("S", D_S, S, {3, 5}, N, end="*")
    check_dynamic("EVEN8", D_EVEN8, C, {1, 2, 4, 6}, N)
    check_dynamic("ODD8", D_ODD8, T, {1, 3, 5, 7}, N)
    # 5x+1: T5 converges for some n, cycles for others; check value faithfulness on a convergent n
    rng = random.Random(2)
    vals, nf = run(canonical(13), SYSTEMS["T5"], {"f>", "t>"}, 5, rng, max_steps=20000)
    assert nf is None, "T5 terminated from 13?!"
    x = 13
    for v in vals:
        x = T5(x)
        assert v == x
    print(f"  T5 (5x+1): from 13, {len(vals)} dynamic steps in 20000 rewrites, value hit 13 "
          f"{sum(v == 13 for v in vals)} times -> nonterminating, as it must be")
    # H (Conj 4.14 as a partial map): dynamic steps apply H; normal forms are where H is undefined
    from srs import D_H, D_ODD8_5
    def H(n):
        return 3 * n // 4 if n % 4 == 0 else (9 * n + 1) // 8
    check_rules(X, D_H, H)
    check_dynamic("H", D_H, H, None, N)
    # 5x+1 analogue of Conj 4.14 is FALSE: the T5-cycle 17 -> 43 -> 108 -> 54 -> 27 -> 68 -> 34 -> 17
    # avoids 5 and 7 mod 8, so ODD8_5 - {tft>, ttt>} must not terminate from 17
    check_rules(X5, D_ODD8_5, T5, base=5)
    vals, nf = run(canonical(17), SYSTEMS["ODD8_5-tft>-ttt>"], {l for l, _ in D_ODD8_5}, 5, rng,
                   max_steps=20000)
    assert nf is None and vals.count(17) > 10
    print(f"  ODD8_5-tft>-ttt> (5x+1, avoid 5,7 mod 8): from 17 nonterminating "
          f"(value 17 revisited {vals.count(17)} times)")
    print("SIMULATION CONTROL PASSED")


if __name__ == "__main__":
    main()
