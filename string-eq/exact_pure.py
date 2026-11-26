"""EXACT proof of Hypothesis P (one loop at the pure point) for a given coprime (p,q), on ALL leading symbols at once.

Method.
 1. Top system (kappa = 1) for the leading coefficients C_0..C_{p-2}; adjoin w = sum r_i C_i (random integers) and take a
    lex Groebner basis with w smallest.  We require shape position: GB = {C_i - g_i(w), f(w)} with f squarefree.
    Then the leading symbols are in bijection with the roots of f, and C_i = g_i(w).
 2. Work in the ring A = Q[w]/(f(w)) (a product of number fields, one factor per Galois orbit of branches).
    Run the hbar-recursion of pure_branch.py with coefficients in A: solve the hbar^1 and hbar^2 linear systems by
    adjugate / determinant, checking that det J1, det J2 are UNITS of A (nonsingular on every branch).
 3. Check exactly: residuals at hbar^0..hbar^3 vanish in A, and a_{p-2,2} - p*gamma1 == 0 in A.
 4. Also report on how many branches every res L^k (k<p) vanishes through hbar^2 (where (F) fails at the pure point).
Usage: python exact_pure.py p q [seed]"""
import sys, random, itertools, time
from fractions import Fraction as Fr
import sympy as sp
from flint import fmpq_poly, fmpq
import pure_branch as pb
from pure_flows_check import res_powers

class A:
    f = None
    __slots__ = ('v',)
    def __init__(self, v):
        if not isinstance(v, fmpq_poly): v = fmpq_poly([v]) if not isinstance(v, Fr) else fmpq_poly([fmpq(v.numerator, v.denominator)])
        self.v = v % A.f if v.degree() >= A.f.degree() else v
    @staticmethod
    def lift(x):
        if isinstance(x, A): return x.v
        if isinstance(x, Fr): return fmpq_poly([fmpq(x.numerator, x.denominator)])
        return fmpq_poly([x])
    def __add__(s, o): return A(s.v + A.lift(o))
    __radd__ = __add__
    def __sub__(s, o): return A(s.v - A.lift(o))
    def __rsub__(s, o): return A(A.lift(o) - s.v)
    def __neg__(s): return A(-s.v)
    def __mul__(s, o): return A(s.v * A.lift(o))
    __rmul__ = __mul__
    def __truediv__(s, o):
        if isinstance(o, A): return s * o.inv()
        o = Fr(o); return A(s.v * fmpq_poly([fmpq(o.denominator, o.numerator)]))
    def __eq__(s, o): return (s.v - A.lift(o)) == 0
    def __ne__(s, o): return not s.__eq__(o)
    __hash__ = None
    def is_unit(s):
        g = s.v.gcd(A.f); return g.degree() == 0
    def inv(s):
        g, a, b = s.v.xgcd(A.f)
        if g.degree() != 0: raise ZeroDivisionError('not a unit in A')
        return A(a / g)

class OpsA(pb.Ops):
    def __init__(self, p, q):
        super().__init__(p, q, True)
        self.zero = A(0)
    def c(self, x): return A(Fr(x))
    def binom(self, m, r):
        num = 1
        for s in range(r): num *= (m - s)
        den = 1
        for s in range(1, r + 1): den *= s
        return Fr(num, den)
    def dfac(self, e, r):
        f = Fr(1)
        for s in range(r): f *= Fr(-(e + s * self.N), self.N)
        return f

def det(M):
    n = len(M)
    if n == 1: return M[0][0]
    tot = A(0)
    for perm in itertools.permutations(range(n)):
        sgn = 1
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]: sgn = -sgn
        t = A(sgn)
        for i in range(n): t = t * M[i][perm[i]]
        tot = tot + t
    return tot

def solve(M, b):
    n = len(M); d = det(M)
    if not d.is_unit(): raise ZeroDivisionError('det J is not a unit: singular on some branch')
    di = d.inv(); x = []
    for k in range(n):
        Mk = [[(b[i] if j == k else M[i][j]) for j in range(n)] for i in range(n)]
        x.append(det(Mk) * di)
    return x, d

def lex_basis(polys, gens):
    """Lex Groebner basis (gens[-1] smallest).  Uses Singular (std in dp, then fglm to lp) when available,
    otherwise sympy grevlex + FGLM."""
    import shutil, subprocess, tempfile, os
    if shutil.which('Singular'):
        names = [str(g) for g in gens]
        src = ['ring r = 0, (%s), dp;' % ','.join(names), 'option(redSB);',
               'ideal I = ' + ', '.join(str(sp.expand(f)).replace('**', '^') for f in polys) + ';',
               'ideal J = std(I);', 'if (dim(J) != 0) { print("NOTZERODIM"); quit; }',
               'ring s = 0, (%s), lp;' % ','.join(names), 'ideal K = fglm(r, J);', 'short = 0;',
               'int i; for (i = 1; i <= size(K); i++) { print("GEN " + string(K[i])); }', 'quit;']
        with tempfile.NamedTemporaryFile('w', suffix='.sing', delete=False) as fh:
            fh.write('\n'.join(src) + '\n'); path = fh.name
        out = subprocess.run(['Singular', '-q', path], capture_output=True, text=True).stdout
        os.unlink(path)
        if 'NOTZERODIM' in out: raise RuntimeError('ideal is not zero-dimensional')
        loc = {n: g for n, g in zip(names, gens)}
        res = [sp.sympify(l[4:].replace('^', '**'), locals=loc) for l in out.splitlines() if l.startswith('GEN ')]
        if res: return res
        print('  Singular returned no generators; falling back to sympy. Singular output head:', out[:400], flush=True)
    G = sp.groebner(polys, *gens, order='grevlex')
    return list(G.fglm('lex'))

def shape_basis(p, q, seed):
    Cs, eqs = pb.top_system(p, q)
    rnd = random.Random(seed)
    for attempt in range(20):
        r = [rnd.randint(-7, 7) for _ in Cs]
        w = sp.Symbol('w')
        tg = time.time()
        gb = lex_basis(eqs + [w - sum(ri * c for ri, c in zip(r, Cs))], list(Cs) + [w])
        print('  Groebner (shape attempt %d): %.1fs, %d generators' % (attempt, time.time() - tg, len(gb)), flush=True)
        fpol = [g for g in gb if g.free_symbols <= {w}]
        lin = {}
        ok = len(fpol) == 1
        for g in gb:
            if g in fpol: continue
            P = sp.Poly(g, *Cs)
            lead = [c for c in Cs if c in g.free_symbols]
            # expect g = c * C_i - h(w) with c rational
            if len(lead) == 1 and sp.Poly(g, lead[0]).degree() == 1:
                cf = sp.Poly(g, lead[0]).all_coeffs()
                if cf[0].free_symbols: ok = False; break
                lin[lead[0]] = sp.expand(-cf[1] / cf[0])
            else:
                ok = False; break
        if ok and len(lin) == len(Cs):
            f = sp.Poly(fpol[0], w)
            if sp.degree(sp.gcd(f, f.diff(w))) == 0:
                return Cs, w, f, lin, r
    raise RuntimeError('no shape position found')

def to_fmpq_poly(expr, w):
    P = sp.Poly(expr, w)
    cs = [sp.Rational(c) for c in reversed(P.all_coeffs())]
    return fmpq_poly([fmpq(int(c.p), int(c.q)) for c in cs])

def main(p, q, seed=1):
    t0 = time.time()
    N = p + q - 1; m = p - 1
    pred = Fr(p * (p - 1) * (q - 1), 24 * N)
    Cs, w, f, lin, r = shape_basis(p, q, seed)
    A.f = to_fmpq_poly(f.as_expr(), w)
    deg = A.f.degree()
    Cv = [A(to_fmpq_poly(lin[c], w)) for c in Cs]
    O = OpsA(p, q)
    Z = A(0)
    def U_of(a, b): return {i: [Cv[i], a[i], b[i], A(0)] for i in range(m)}
    zero = [A(0)] * m
    R0 = pb.string_residual(O, U_of(zero, zero))
    assert all(R0[j][0] == 0 for j in range(m)), 'hbar^0 residual'
    assert all(R0[j][1] == (1 if j == 0 else 0) for j in range(m)), 'top system not satisfied'
    base = [R0[j][2] for j in range(m)]
    J1 = [[None] * m for _ in range(m)]
    for i in range(m):
        e = [A(0)] * m; e[i] = A(1)
        Ri = pb.string_residual(O, U_of(e, zero))
        for j in range(m): J1[j][i] = Ri[j][2] - base[j]
    a, d1 = solve(J1, [-x for x in base])
    R1 = pb.string_residual(O, U_of(a, zero))
    base2 = [R1[j][3] for j in range(m)]
    J2 = [[None] * m for _ in range(m)]
    for i in range(m):
        e = [A(0)] * m; e[i] = A(1)
        Ri = pb.string_residual(O, U_of(a, e))
        for j in range(m): J2[j][i] = Ri[j][3] - base2[j]
    b, d2 = solve(J2, [-x for x in base2])
    Rf = pb.string_residual(O, U_of(a, b))
    assert all(Rf[j][h] == 0 for j in range(m) for h in (2, 3)), 'residual check failed'
    val = b[m - 1]
    ok = (val == A(pred))
    # (F) at the pure point: branches on which res L^k vanishes through hbar^2 for some k < p
    U = U_of(a, b)
    R = res_powers(O, U, p - 1)
    badF = fmpq_poly([0])   # product over k of gcd(f, components) -> count branches with some vanishing res L^k
    nbad = 0; bad_roots = fmpq_poly([1])
    for k in range(1, p):
        g = A.f
        for h in range(3): g = g.gcd(R[k][h].v) if R[k][h].v != 0 else g
        if g.degree() > 0:
            bad_roots = (bad_roots * g) // bad_roots.gcd(g)   # lcm
    nbad = bad_roots.degree()
    print('(p,q)=(%d,%d): %d leading symbols (deg f); det J1, det J2 units: yes; a_{p-2,2} == p*gamma1 = %s EXACTLY on all: %s;'
          ' branches where (F) fails at pure point: %d;  [%.1fs]' % (p, q, deg, pred, ok, nbad, time.time() - t0))
    return ok

if __name__ == '__main__':
    p, q = int(sys.argv[1]), int(sys.argv[2])
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    ok = main(p, q, seed)
    sys.exit(0 if ok else 1)
