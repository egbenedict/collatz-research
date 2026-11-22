"""Matrix-interpretation search (natural and arctic) via z3 bit-vectors, plus an exact verifier.

One search step: given rule lists `strict` (candidates to remove) and `weak` (everything that
must weakly decrease), find interpretations such that every rule in strict+weak weakly
decreases and at least one rule of `strict` strictly decreases.

mode = "ext": extended monotone algebra -> proves SN(removed / all) (relative termination).
              natural: (M_s)_11 >= 1 for every symbol; arctic: v_s = -inf, (M_s)_11 finite.
mode = "top": weakly monotone algebra -> proves SN(removed_top / all) (relative top termination).
              natural: no extra condition; arctic: (M_s)_11 or (v_s)_1 finite.

Orders (x, y vectors):
  natural: x > y  iff x1 > y1 and xi >= yi (i >= 2);   [l] > [r] iff M_l >= M_r and v_l > v_r.
  arctic ext: KW09 Lemma 6.5, [l] >> [r] iff M_l >> M_r entrywise (a >> b: a > b or both -inf).
  arctic top: first row strict (M_l[1,:] >> M_r[1,:], v_l1 >> v_r1), all other entries >=.
             Sound for weakly monotone algebras on N x A^{d-1} (argument in NOTES.md).
Entries: natural in [0, 2^w); arctic in {-inf} u [0, 2^w).
"""
import itertools
import time

import z3

NEG = None  # arctic -infinity


# ============================================================ exact (python) semantics
def nat_compose(interp, s, d):
    M = [[int(i == j) for j in range(d)] for i in range(d)]
    v = [0] * d
    for c in reversed(s):            # [s] = [s1] o ... o [sn]
        Mc, vc = interp[c]
        v = [sum(Mc[i][k] * v[k] for k in range(d)) + vc[i] for i in range(d)]
        M = [[sum(Mc[i][k] * M[k][j] for k in range(d)) for j in range(d)] for i in range(d)]
    return M, v


def a_add(a, b):
    return NEG if a is NEG or b is NEG else a + b


def a_max(xs):
    xs = [x for x in xs if x is not NEG]
    return max(xs) if xs else NEG


def arc_compose(interp, s, d):
    M = [[0 if i == j else NEG for j in range(d)] for i in range(d)]
    v = [NEG] * d
    for c in reversed(s):
        Mc, vc = interp[c]
        v = [a_max([a_add(Mc[i][k], v[k]) for k in range(d)] + [vc[i]]) for i in range(d)]
        M = [[a_max([a_add(Mc[i][k], M[k][j]) for k in range(d)]) for j in range(d)]
             for i in range(d)]
    return M, v


def a_gg(a, b):   # a >> b
    return (a is NEG and b is NEG) or (a is not NEG and (b is NEG or a > b))


def a_ge(a, b):
    return b is NEG or (a is not NEG and a >= b)


def verify(interp, kind, mode, d, rules_strict_claimed, rules_weak):
    """Exact check. Returns the list of rules that are strictly decreasing (must include all of
    rules_strict_claimed) or raises AssertionError."""
    syms = sorted({c for l, r in rules_strict_claimed + rules_weak for c in l + r})
    for c in syms:
        M, v = interp[c]
        if kind == "nat":
            assert all(x >= 0 for row in M for x in row) and all(x >= 0 for x in v)
            if mode == "ext":
                assert M[0][0] >= 1, c
        else:
            if mode == "ext":
                assert all(x is NEG for x in v) and M[0][0] is not NEG, c
            else:
                assert M[0][0] is not NEG or v[0] is not NEG, c
    strict = []
    for l, r in rules_strict_claimed + rules_weak:
        if kind == "nat":
            Ml, vl = nat_compose(interp, l, d)
            Mr, vr = nat_compose(interp, r, d)
            weak = all(Ml[i][j] >= Mr[i][j] for i in range(d) for j in range(d)) and \
                all(vl[i] >= vr[i] for i in range(d))
            st = weak and vl[0] > vr[0]
        else:
            Ml, vl = arc_compose(interp, l, d)
            Mr, vr = arc_compose(interp, r, d)
            weak = all(a_ge(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d)) and \
                all(a_ge(vl[i], vr[i]) for i in range(d))
            if mode == "ext":
                st = weak and all(a_gg(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d))
            else:
                st = weak and all(a_gg(Ml[0][j], Mr[0][j]) for j in range(d)) and a_gg(vl[0], vr[0])
        assert weak, ("not weakly decreasing", l, r)
        if st:
            strict.append((l, r))
    for rule in rules_strict_claimed:
        pass
    return strict


# ============================================================ z3 encoding
class Encoder:
    def __init__(self, kind, mode, d, w, rules):
        self.kind, self.mode, self.d, self.w = kind, mode, d, w
        self.syms = sorted({c for l, r in rules for c in l + r})
        K = max(max(len(l), len(r)) for l, r in rules)
        self.K = K
        self.s = z3.Solver()
        self.memo = {}
        dd = range(d)
        if kind == "nat":
            self.W = w * K + (K * d ** K).bit_length() + 2
            self.var = {}
            for c in self.syms:
                M = [[z3.BitVec(f"M_{c}_{i}_{j}", w) for j in dd] for i in dd]
                v = [z3.BitVec(f"v_{c}_{i}", w) for i in dd]
                self.var[c] = (M, v)
                if mode == "ext":
                    self.s.add(M[0][0] != 0)
            self.base = {c: ([[self.ext(x) for x in row] for row in M], [self.ext(x) for x in v])
                         for c, (M, v) in self.var.items()}
        else:
            self.W = w + (K + 1).bit_length() + 1
            self.var = {}
            for c in self.syms:
                M = [[(z3.Bool(f"F_{c}_{i}_{j}"), z3.BitVec(f"M_{c}_{i}_{j}", w)) for j in dd]
                     for i in dd]
                v = [(z3.Bool(f"G_{c}_{i}"), z3.BitVec(f"v_{c}_{i}", w)) for i in dd]
                self.var[c] = (M, v)
                if mode == "ext":
                    self.s.add(M[0][0][0])
                    for f, _ in v:
                        self.s.add(z3.Not(f))
                else:
                    self.s.add(z3.Or(M[0][0][0], v[0][0]))
            self.base = {c: ([[(f, self.ext(x)) for f, x in row] for row in M],
                             [(f, self.ext(x)) for f, x in v]) for c, (M, v) in self.var.items()}

    def ext(self, x):
        return z3.ZeroExt(self.W - self.w, x)

    # ---------------- arctic helpers
    @staticmethod
    def amax(xs):
        xs = [x for x in xs if x is not None]
        if not xs:
            return None
        f, v = xs[0]
        for g, u in xs[1:]:
            take = z3.And(g, z3.Or(z3.Not(f), z3.UGT(u, v)))
            v = z3.If(take, u, v)
            f = z3.Or(f, g)
        return (z3.simplify(f), v)

    @staticmethod
    def aplus(a, b):
        if a is None or b is None:
            return None
        return (z3.And(a[0], b[0]), a[1] + b[1])

    def compose(self, s):
        if s in self.memo:
            return self.memo[s]
        d = self.d
        if len(s) == 1:
            res = self.base[s]
        else:
            Mc, vc = self.base[s[0]]
            M, v = self.compose(s[1:])
            if self.kind == "nat":
                v2 = [sum((Mc[i][k] * v[k] for k in range(d)), vc[i]) for i in range(d)]
                M2 = [[sum(Mc[i][k] * M[k][j] for k in range(1, d)) + Mc[i][0] * M[0][j]
                       for j in range(d)] for i in range(d)]
            else:
                v2 = [self.amax([self.aplus(Mc[i][k], v[k]) for k in range(d)] + [vc[i]])
                      for i in range(d)]
                M2 = [[self.amax([self.aplus(Mc[i][k], M[k][j]) for k in range(d)])
                       for j in range(d)] for i in range(d)]
            res = (M2, v2)
        self.memo[s] = res
        return res

    # arctic comparisons
    @staticmethod
    def a_ge(a, b):
        return z3.Or(z3.Not(b[0]), z3.And(a[0], z3.UGE(a[1], b[1])))

    @staticmethod
    def a_gg(a, b):
        return z3.Or(z3.And(z3.Not(a[0]), z3.Not(b[0])),
                     z3.And(a[0], z3.Or(z3.Not(b[0]), z3.UGT(a[1], b[1]))))

    def weak_and_strict(self, l, r):
        Ml, vl = self.compose(l)
        Mr, vr = self.compose(r)
        d = self.d
        if self.kind == "nat":
            weak = z3.And([z3.UGE(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d)] +
                          [z3.UGE(vl[i], vr[i]) for i in range(d)])
            strict = z3.UGT(vl[0], vr[0])
        else:
            weak = z3.And([self.a_ge(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d)] +
                          [self.a_ge(vl[i], vr[i]) for i in range(d)])
            if self.mode == "ext":
                strict = z3.And([self.a_gg(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d)])
            else:
                strict = z3.And([self.a_gg(Ml[0][j], Mr[0][j]) for j in range(d)] +
                                [self.a_gg(vl[0], vr[0])])
        return weak, strict

    def extract(self, m):
        out = {}
        for c, (M, v) in self.var.items():
            if self.kind == "nat":
                out[c] = ([[m.eval(x, model_completion=True).as_long() for x in row] for row in M],
                          [m.eval(x, model_completion=True).as_long() for x in v])
            else:
                def g(p):
                    f, x = p
                    if not z3.is_true(m.eval(f, model_completion=True)):
                        return NEG
                    return m.eval(x, model_completion=True).as_long()
                out[c] = ([[g(p) for p in row] for row in M], [g(p) for p in v])
        return out


def search(strict_cands, weak_rules, kind, mode, d, w, timeout_s=60, seed=0):
    """Return (status, interp, removed_rules, seconds). status in sat/unsat/unknown."""
    rules = strict_cands + weak_rules
    t0 = time.time()
    enc = Encoder(kind, mode, d, w, rules)
    sel = []
    for l, r in rules:
        weak, strict = enc.weak_and_strict(l, r)
        enc.s.add(weak)
        if (l, r) in strict_cands:
            b = z3.Bool(f"sel_{l}_{r}")
            enc.s.add(z3.Implies(b, strict))
            sel.append(b)
    enc.s.add(z3.Or(sel))
    enc.s.set("timeout", int(timeout_s * 1000))
    enc.s.set("random_seed", seed)
    res = enc.s.check()
    dt = time.time() - t0
    if res == z3.sat:
        interp = enc.extract(enc.s.model())
        removed = verify(interp, kind, mode, d, [], rules)
        removed = [r for r in removed if r in strict_cands]
        assert removed, "solver model removes nothing?!"
        return "sat", interp, removed, dt
    return ("unsat" if res == z3.unsat else "unknown"), None, [], dt


def fmt_interp(interp, kind):
    def e(x):
        return "-" if x is NEG else str(x)
    lines = []
    for c, (M, v) in sorted(interp.items()):
        lines.append(f"  [{c}]: M={[[e(x) for x in row] for row in M]} v={[e(x) for x in v]}")
    return "\n".join(lines)
