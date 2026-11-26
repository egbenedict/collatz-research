"""Semiclassical (hbar) expansion of the PURE (p,q) string equation [P,(L^q)_+] = kappa at a branch over t = infinity.

Variables: t = z^{-N}, N = p+q-1, d/dt z^e = (-e/N) z^{e+N}.  An operator of weight w has, at slot d^m, a coefficient
sum_h c_h z^{(m-w) + h(N+1)};  h is the hbar-order (each derivative costs one unit of h).  With no lower times and no
integration constants the branch solution only has relative orders that are multiples of N+1 (see notes), so a
coefficient is a short vector over h.  We keep h <= HMAX = 3, enough to get the hbar^2 term of u = u_{p-2}, which is
exactly [t^-2]u.

Usage:  python pure_branch.py p q          (all branches, numerically, 40 digits)
Prints, for every leading solution (kappa normalised to 1): det J, prod of critical values of Lambda, and [t^-2]u,
compared with p*gamma1 = p(p-1)(q-1)/(24(p+q-1))."""
import sys, itertools
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
HMAX = 3

class Ops:
    def __init__(self, p, q, mpmode):
        self.p, self.q, self.N = p, q, p + q - 1
        self.mp = mpmode
        self.zero = mp.mpc(0) if mpmode else sp.Integer(0)
    def c(self, x):
        return mp.mpc(x) if self.mp else sp.nsimplify(x) if not isinstance(x, sp.Basic) else x
    def binom(self, m, r):
        num = 1
        for s in range(r): num *= (m - s)
        den = 1
        for s in range(1, r + 1): den *= s
        return (mp.mpf(num) / den) if self.mp else sp.Rational(num, den)
    def dfac(self, e, r):
        f = mp.mpf(1) if self.mp else sp.Integer(1)
        for s in range(r):
            f = f * ((mp.mpf(-(e + s * self.N)) / self.N) if self.mp else sp.Rational(-(e + s * self.N), self.N))
        return f
    def mul(self, A, wA, B, wB, minslot):
        """A*B for operators A (weight wA), B (weight wB); dict slot -> [coef_h]."""
        N = self.N; out = {}
        for m, ca in A.items():
            for n, cb in B.items():
                for r in range(0, HMAX + 1):
                    s = m + n - r
                    if s < minslot: break
                    bn = self.binom(m, r)
                    if bn == 0: continue
                    for ha, xa in enumerate(ca):
                        if xa == 0: continue
                        for hb, xb in enumerate(cb):
                            h = ha + hb + r
                            if h > HMAX: break
                            if xb == 0: continue
                            eb = (n - wB) + hb * (N + 1)
                            val = bn * xa * xb * self.dfac(eb, r)
                            row = out.setdefault(s, [self.zero] * (HMAX + 1))
                            row[h] = row[h] + val
        return out
    def sub(self, A, B):
        out = {k: v[:] for k, v in A.items()}
        for k, v in B.items():
            row = out.setdefault(k, [self.zero] * (HMAX + 1))
            for h in range(HMAX + 1): row[h] = row[h] - v[h]
        return out

def string_residual(O, U):
    """U: dict i -> [u_i at h=0..HMAX] for i=0..p-2.  Returns commutator slots 0..p-2 (lists over h)."""
    p, q = O.p, O.q
    one = [O.c(1)] + [O.zero] * HMAX
    zero = [O.zero] * (HMAX + 1)
    P = {p: one[:], p - 1: zero[:]}
    for i in range(p - 1): P[i] = U[i][:]
    D = q + 1
    # L = P^{1/p} = d + sum_{r>=1} l_r d^{-r}
    L = {1: one[:]}
    for r in range(1, D + 1):
        X = L
        for _ in range(p - 1): X = O.mul(X, 1, L, 1, -D - p)
        diff = O.sub(P, {k: v for k, v in X.items()})
        tgt = diff.get(p - 1 - r, zero)
        L[-r] = [x / p for x in tgt]
        L[1 - 1] = L.get(0, zero[:])
    # Q = (L^q)_+
    X = L
    for _ in range(q - 1): X = O.mul(X, 1, L, 1, -D)
    Q = {k: v for k, v in X.items() if k >= 0}
    PQ = O.mul(P, p, Q, q, 0); QP = O.mul(Q, q, P, p, 0)
    C = O.sub(PQ, QP)
    return {j: C.get(j, zero[:]) for j in range(p - 1)}

def top_system(p, q):
    O = Ops(p, q, False)
    Cs = sp.symbols('C0:%d' % (p - 1))
    U = {i: [Cs[i]] + [sp.Integer(0)] * HMAX for i in range(p - 1)}
    R = string_residual(O, U)
    eqs = [sp.expand(R[j][1] - (1 if j == 0 else 0)) for j in range(p - 1)]   # kappa = 1
    for j in range(p - 1):
        for h in (0,): assert sp.expand(R[j][h]) == 0
    return Cs, eqs

def solve_top(Cs, eqs, tol=mp.mpf(10) ** -25):
    """All solutions of the top system, numerically (60 digits), via a lex Groebner basis and triangular back-substitution."""
    mp.mp.dps = 60
    gb = list(sp.groebner(eqs, *Cs, order='lex'))
    n = len(Cs)
    partial = [dict()]
    for k in reversed(range(n)):
        var = Cs[k]; later = set(Cs[k + 1:])
        polys = [g for g in gb if g.free_symbols <= later | {var} and var in g.free_symbols]
        new = []
        for sol in partial:
            subs = {c: sp.Float(sp.re(v), 70) + sp.I * sp.Float(sp.im(v), 70) for c, v in sol.items()} if sol else {}
            uni = []
            for g in polys:
                gg = sp.Poly(sp.expand(g.subs(subs)), var)
                cs = [complex(sp.N(c, 60)) for c in gg.all_coeffs()]
                if gg.degree() >= 1 and max(abs(c) for c in cs) > 1e-30: uni.append(gg)
            if not uni: raise RuntimeError('positive-dimensional fibre')
            g0 = min(uni, key=lambda P: P.degree())
            co = [mp.mpc(mp.mpf(str(sp.re(sp.N(c, 60)))), mp.mpf(str(sp.im(sp.N(c, 60))))) for c in g0.all_coeffs()]
            roots = mp.polyroots(co, maxsteps=500, extraprec=400) if len(co) > 2 else [-co[1] / co[0]]
            for r in roots:
                ok = True
                for P in uni:
                    cP = [mp.mpc(mp.mpf(str(sp.re(sp.N(c, 60)))), mp.mpf(str(sp.im(sp.N(c, 60))))) for c in P.all_coeffs()]
                    if abs(mp.polyval(cP, r)) > mp.mpf(10) ** -20 * (1 + max(abs(c) for c in cP)): ok = False; break
                if ok:
                    d = dict(sol); d[var] = r; new.append(d)
        partial = new
    out = []
    for sol in partial:
        row = [sol[c] for c in Cs]
        if not any(max(abs(x - y) for x, y in zip(row, o)) < mp.mpf(10) ** -15 for o in out): out.append(row)
    mp.mp.dps = 40
    return out

def branch_value(p, q, Cvals):
    O = Ops(p, q, True)
    m = p - 1
    def U_of(a, b):
        return {i: [Cvals[i], a[i], b[i], mp.mpc(0)] for i in range(m)}
    z = [mp.mpc(0)] * m
    # check top
    R0 = string_residual(O, U_of(z, z))
    top_err = max(abs(R0[j][1] - (1 if j == 0 else 0)) for j in range(m))
    # J for h=1 unknowns: residual at h=2 is affine in a
    base = [R0[j][2] for j in range(m)]
    J = mp.matrix(m, m)
    for i in range(m):
        e = [mp.mpc(0)] * m; e[i] = mp.mpc(1)
        Ri = string_residual(O, U_of(e, z))
        for j in range(m): J[j, i] = Ri[j][2] - base[j]
    a = mp.lu_solve(J, mp.matrix([-x for x in base]))
    a = [a[i] for i in range(m)]
    R1 = string_residual(O, U_of(a, z))
    base2 = [R1[j][3] for j in range(m)]
    J2 = mp.matrix(m, m)
    for i in range(m):
        e = [mp.mpc(0)] * m; e[i] = mp.mpc(1)
        Ri = string_residual(O, U_of(a, e))
        for j in range(m): J2[j, i] = Ri[j][3] - base2[j]
    b = mp.lu_solve(J2, mp.matrix([-x for x in base2]))
    b = [b[i] for i in range(m)]
    Rf = string_residual(O, U_of(a, b))
    res_err = max(abs(Rf[j][h]) for j in range(m) for h in (2, 3))
    # critical values of Lambda(eta) = eta^p + sum C_i eta^i
    eta = sp.Symbol('eta')
    coeffs = [1, 0] + [complex(Cvals[i]) for i in reversed(range(m))]   # eta^p, eta^{p-1}=0, then C_{p-2}..C_0
    dcoeffs = [c * (p - k) for k, c in enumerate(coeffs[:-1])]
    Lam = lambda x: sum(mp.mpc(c) * x ** (p - k) for k, c in enumerate(coeffs))
    try:
        crit = mp.polyroots([mp.mpc(c) for c in dcoeffs], maxsteps=400, extraprec=400)
        prodcv = mp.fprod([Lam(x) for x in crit])
    except mp.libmp.libhyper.NoConvergence:   # diagnostic only (repeated critical points)
        prodcv = mp.nan
    return dict(a=a, b=b, detJ=mp.det(J), detJ2=mp.det(J2), prodcv=prodcv, val=b[m - 1], a_top=a[m - 1], top_err=top_err, res_err=res_err)

if __name__ == '__main__':
    p, q = int(sys.argv[1]), int(sys.argv[2])
    N = p + q - 1
    pred = sp.Rational(p * (p - 1) * (q - 1), 24 * N)
    Cs, eqs = top_system(p, q)
    sols = solve_top(Cs, eqs)
    sols = [s for s in sols if any(abs(x) > 1e-20 for x in s) or True]
    print('(p,q)=(%d,%d)  N=%d  #leading solutions (kappa=1): %d   predicted [t^-2]u = p*gamma1 = %s' % (p, q, N, len(sols), pred))
    worst = 0
    for s in sols:
        r = branch_value(p, q, s)
        dev = abs(r['val'] - mp.mpf(pred.p) / pred.q)
        worst = max(worst, dev)
        print('  C=%s  |  [t^-2]u=%s  dev=%.2e  detJ=%s  prod crit vals=%s  ratio=%s  top_err=%.1e res_err=%.1e' % (
            [mp.nstr(x, 6) for x in s], mp.nstr(r['val'], 12), float(dev), mp.nstr(r['detJ'], 8), mp.nstr(r['prodcv'], 8),
            mp.nstr(r['detJ'] / r['prodcv'], 8) if abs(r['prodcv']) > 1e-30 else 'inf', float(r['top_err']), float(r['res_err'])))
    print('max deviation from p*gamma1: %.3e' % float(worst))
