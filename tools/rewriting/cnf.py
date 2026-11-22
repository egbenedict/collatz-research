"""Direct CNF encoding of the matrix-interpretation search (faster than z3 bit-vectors).

  natural: unsigned binary bit-vectors with exact widths (no overflow possible), ripple adders,
           shift-and-add multipliers, constant propagation.
  arctic : order encoding: X represented by bits [X >= 0], [X >= 1], ..., [X >= U];
           -inf = all false. max = bitwise OR; X + Y >= k = OR_{i+j=k} ([X>=i] and [Y>=j]).
Parameter w: natural entries in [0, 2^w); arctic entries in {-inf, 0, 1, ..., w}.
Same orders/modes as matrix.py (see its docstring). Models are decoded and then re-verified
with matrix.verify (exact integer arithmetic), so encoding bugs cannot produce false proofs.

Solving: a portfolio of SAT processes (CaDiCaL 1.9.5 / Kissat / Glucose / MiniSat, optional
negative branching, clause shuffles), first answer wins, hard timeout per query.
"""
import multiprocessing as mp
import random
import time

from pysat.solvers import Solver

from matrix import verify, NEG

TRUE = "T"
FALSE = "F"


class CNF:
    def __init__(self):
        self.nv = 0
        self.clauses = []
        self.cache = {}

    def new(self):
        self.nv += 1
        return self.nv

    # literals: int (var / -var) or TRUE / FALSE constants
    @staticmethod
    def neg(a):
        if a == TRUE:
            return FALSE
        if a == FALSE:
            return TRUE
        return -a

    def add(self, cl):
        out = []
        for x in cl:
            if x == TRUE:
                return
            if x == FALSE:
                continue
            out.append(x)
        if not out:
            self.clauses.append([])  # UNSAT marker
            return
        self.clauses.append(out)

    def AND(self, xs):
        xs = [x for x in xs if x != TRUE]
        if any(x == FALSE for x in xs):
            return FALSE
        xs = sorted(set(xs))
        if not xs:
            return TRUE
        if len(xs) == 1:
            return xs[0]
        if any(-x in xs for x in xs):
            return FALSE
        key = ("and", tuple(xs))
        if key in self.cache:
            return self.cache[key]
        z = self.new()
        for x in xs:
            self.add([-z, x])
        self.add([z] + [-x for x in xs])
        self.cache[key] = z
        return z

    def OR(self, xs):
        return self.neg(self.AND([self.neg(x) for x in xs]))

    def XOR(self, a, b):
        if a == FALSE:
            return b
        if b == FALSE:
            return a
        if a == TRUE:
            return self.neg(b)
        if b == TRUE:
            return self.neg(a)
        if a == b:
            return FALSE
        if a == -b:
            return TRUE
        key = ("xor", min(a, b), max(a, b))
        if key in self.cache:
            return self.cache[key]
        z = self.new()
        self.add([-z, a, b]); self.add([-z, -a, -b]); self.add([z, -a, b]); self.add([z, a, -b])
        self.cache[key] = z
        return z

    def MAJ(self, a, b, c):
        return self.OR([self.AND([a, b]), self.AND([a, c]), self.AND([b, c])])

    # ---------------- binary unsigned numbers: list of literals, LSB first
    def badd(self, x, y):
        n = max(len(x), len(y))
        x = x + [FALSE] * (n - len(x))
        y = y + [FALSE] * (n - len(y))
        out, c = [], FALSE
        for a, b in zip(x, y):
            out.append(self.XOR(self.XOR(a, b), c))
            c = self.MAJ(a, b, c)
        out.append(c)
        return self.trim(out)

    @staticmethod
    def trim(x):
        while x and x[-1] == FALSE:
            x = x[:-1]
        return x

    def bmul(self, x, y):
        acc = []
        for i, b in enumerate(y):
            if b == FALSE:
                continue
            pp = [FALSE] * i + [self.AND([a, b]) for a in x]
            acc = self.badd(acc, self.trim(pp)) if acc else self.trim(pp)
        return acc

    def bsum(self, xs):
        xs = [x for x in xs if x]
        while len(xs) > 1:          # balanced tree
            nxt = []
            for i in range(0, len(xs) - 1, 2):
                nxt.append(self.badd(xs[i], xs[i + 1]))
            if len(xs) % 2:
                nxt.append(xs[-1])
            xs = nxt
        return xs[0] if xs else []

    def bge(self, x, y, strict=False):
        """literal for x >= y (or x > y)."""
        n = max(len(x), len(y))
        x = x + [FALSE] * (n - len(x))
        y = y + [FALSE] * (n - len(y))
        res = FALSE if strict else TRUE   # comparison of empty suffix
        for a, b in zip(x, y):           # LSB to MSB: res_i = (a>b) or (a==b and res_{i-1})
            gt = self.AND([a, self.neg(b)])
            eq = self.neg(self.XOR(a, b))
            res = self.OR([gt, self.AND([eq, res])])
        return res

    # ---------------- arctic order-encoded: list u where u[k] = [X >= k]; [] = -inf
    def amax(self, xs):
        n = max((len(x) for x in xs), default=0)
        out = []
        for k in range(n):
            out.append(self.OR([x[k] for x in xs if len(x) > k]))
        return self.atrim(out)

    def aplus(self, x, y):
        if not x or not y:
            return []
        out = []
        for k in range(len(x) + len(y) - 1):
            terms = [self.AND([x[i], y[k - i]]) for i in range(len(x)) if 0 <= k - i < len(y)]
            out.append(self.OR(terms))
        return self.atrim(out)

    @staticmethod
    def atrim(x):
        while x and x[-1] == FALSE:
            x = x[:-1]
        return x

    def age(self, x, y):
        """x >= y (arctic)"""
        cl = []
        for k in range(len(y)):
            xk = x[k] if k < len(x) else FALSE
            cl.append(self.OR([self.neg(y[k]), xk]))
        return self.AND(cl)

    def agg(self, x, y):
        """x >> y: x > y or both -inf  <=>  for all k >= 0: y >= k -> x >= k+1"""
        cl = []
        for k in range(len(y)):
            xk1 = x[k + 1] if k + 1 < len(x) else FALSE
            cl.append(self.OR([self.neg(y[k]), xk1]))
        return self.AND(cl)


class Enc:
    def __init__(self, kind, mode, d, w, rules):
        self.kind, self.mode, self.d, self.w = kind, mode, d, w
        self.c = CNF()
        self.syms = sorted({ch for l, r in rules for ch in l + r})
        self.var = {}
        self.memo = {}
        c, dd = self.c, range(d)
        U = 2 ** w - 1 if kind == "nat" else w   # arctic: parameter w IS the max finite entry
        for s in self.syms:
            if kind == "nat":
                M = [[[c.new() for _ in range(w)] for j in dd] for i in dd]
                v = [[c.new() for _ in range(w)] for i in dd]
                if mode == "ext":
                    c.add(M[0][0])          # M11 >= 1
            else:
                def avar():
                    bits = [c.new() for _ in range(U + 1)]
                    for k in range(U):
                        c.add([-bits[k + 1], bits[k]])
                    return bits
                M = [[avar() for j in dd] for i in dd]
                v = [avar() for i in dd]
                if mode == "ext":
                    for x in v:
                        for b in x:
                            c.add([-b])
                    c.add([M[0][0][0]])
                else:
                    c.add([M[0][0][0], v[0][0]])
            self.var[s] = (M, v)

    def compose(self, s):
        if s in self.memo:
            return self.memo[s]
        if len(s) == 1:
            res = self.var[s]
        else:
            c, d = self.c, self.d
            Mc, vc = self.var[s[0]]
            M, v = self.compose(s[1:])
            if self.kind == "nat":
                v2 = [c.bsum([c.bmul(Mc[i][k], v[k]) for k in range(d)] + [vc[i]]) for i in range(d)]
                M2 = [[c.bsum([c.bmul(Mc[i][k], M[k][j]) for k in range(d)]) for j in range(d)]
                      for i in range(d)]
            else:
                v2 = [c.amax([c.aplus(Mc[i][k], v[k]) for k in range(d)] + [vc[i]]) for i in range(d)]
                M2 = [[c.amax([c.aplus(Mc[i][k], M[k][j]) for k in range(d)]) for j in range(d)]
                      for i in range(d)]
            res = (M2, v2)
        self.memo[s] = res
        return res

    def weak_strict(self, l, r):
        c, d = self.c, self.d
        Ml, vl = self.compose(l)
        Mr, vr = self.compose(r)
        if self.kind == "nat":
            weak = [c.bge(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d)] + \
                   [c.bge(vl[i], vr[i]) for i in range(d)]
            strict = c.bge(vl[0], vr[0], strict=True)
        else:
            weak = [c.age(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d)] + \
                   [c.age(vl[i], vr[i]) for i in range(d)]
            if self.mode == "ext":
                strict = c.AND([c.agg(Ml[i][j], Mr[i][j]) for i in range(d) for j in range(d)])
            else:
                strict = c.AND([c.agg(Ml[0][j], Mr[0][j]) for j in range(d)] + [c.agg(vl[0], vr[0])])
        return weak, strict

    def decode(self, model):
        pos = set(x for x in model if x > 0)

        def val(lit):
            return lit == TRUE or (lit != FALSE and lit in pos)
        out = {}
        for s, (M, v) in self.var.items():
            if self.kind == "nat":
                f = lambda bits: sum(1 << i for i, b in enumerate(bits) if val(b))
            else:
                def f(bits):
                    if not val(bits[0]):
                        return NEG
                    return max(k for k in range(len(bits)) if val(bits[k]))
            out[s] = ([[f(x) for x in row] for row in M], [f(x) for x in v])
        return out


def build(strict_cands, weak_rules, kind, mode, d, w):
    rules = strict_cands + weak_rules
    e = Enc(kind, mode, d, w, rules)
    sel = []
    for l, r in rules:
        weak, strict = e.weak_strict(l, r)
        for x in weak:
            e.c.add([x])
        if (l, r) in strict_cands:
            sel.append(strict)
    e.c.add(sel)
    return e


SOLVER_CYCLE = [("cadical195", False), ("cadical195", True), ("kissat404", False),
                ("glucose4", True), ("minisat22", True), ("cadical195", False),
                ("cadical195", True), ("glucose4", False)]


def _worker(args):
    clauses, nv, solver, negbranch, seed, q = args
    cl = list(clauses)
    if seed:
        random.Random(seed).shuffle(cl)
    with Solver(name=solver, bootstrap_with=cl) as s:
        if negbranch:
            try:
                s.set_phases([-v for v in range(1, nv + 1)])
            except Exception:
                pass
        ok = s.solve()
        q.put((ok, s.get_model() if ok else None, solver, negbranch, seed))


def solve_portfolio(clauses, nv, timeout_s, jobs=1):
    if any(len(c) == 0 for c in clauses):
        return False, None, "trivial"
    ctx = mp.get_context("fork")
    q = ctx.Queue()
    procs = []
    for i in range(jobs):
        solver, neg = SOLVER_CYCLE[i % len(SOLVER_CYCLE)]
        p = ctx.Process(target=_worker, args=((clauses, nv, solver, neg, i, q),))
        p.start()
        procs.append(p)
    try:
        ok, model, solver, neg, seed = q.get(timeout=timeout_s)
        tag = f"{solver}{'-neg' if neg else ''}#{seed}"
    except Exception:
        ok, model, tag = None, None, "timeout"
    for p in procs:
        if p.is_alive():
            p.kill()
        p.join()
    return ok, model, tag


def search(strict_cands, weak_rules, kind, mode, d, w, timeout_s=60, jobs=1):
    t0 = time.time()
    e = build(strict_cands, weak_rules, kind, mode, d, w)
    tb = time.time() - t0
    ok, model, tag = solve_portfolio(e.c.clauses, e.c.nv, timeout_s, jobs)
    dt = time.time() - t0
    info = {"vars": e.c.nv, "clauses": len(e.c.clauses), "build_s": round(tb, 2), "solver": tag}
    if ok:
        interp = e.decode(model)
        allr = strict_cands + weak_rules
        removed = [r for r in verify(interp, kind, mode, d, [], allr) if r in strict_cands]
        assert removed, "model removes nothing?!"
        return "sat", interp, removed, dt, info
    return ("unsat" if ok is False else "unknown"), None, [], dt, info
