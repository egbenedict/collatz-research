"""Side conditions for the flow argument (reduction to the pure point), checked at every pure branch of (p,q):
  (N) nondegeneracy: the linear systems for the hbar^1 and hbar^2 corrections are nonsingular (det J1, det J2 != 0);
  (F) for 1 <= k <= p-1, res L^k is not identically zero on the pure branch (so the flows t_1..t_{p-1} are independent).
Usage: python pure_flows_check.py p q"""
import sys
import mpmath as mp
import pure_branch as pb

def res_powers(O, U, kmax):
    p, q, H = O.p, O.q, pb.HMAX
    one = [O.c(1)] + [O.zero] * H; zero = [O.zero] * (H + 1)
    P = {p: one[:], p - 1: zero[:]}
    for i in range(p - 1): P[i] = U[i][:]
    D = max(q, kmax) + 2
    L = {1: one[:]}
    for r in range(1, D + 1):
        X = L
        for _ in range(p - 1): X = O.mul(X, 1, L, 1, -D - p)
        diff = O.sub(P, X)
        L[-r] = [x / p for x in diff.get(p - 1 - r, zero)]
    out = {}; X = L; w = 1
    for k in range(1, kmax + 1):
        if k > 1: X = O.mul(X, k - 1, L, 1, -D);
        out[k] = X.get(-1, zero)
    return out

if __name__ == '__main__':
    p, q = int(sys.argv[1]), int(sys.argv[2])
    Cs, eqs = pb.top_system(p, q)
    sols = pb.solve_top(Cs, eqs)
    ok = True; mindet = mp.inf; minres = mp.inf
    for s in sols:
        r = pb.branch_value(p, q, s)
        O = pb.Ops(p, q, True)
        # recompute a, b quickly through branch_value internals is not exposed; flows need only h<=2 of res L^k at the
        # solution, and the hbar-0 part already decides non-vanishing except in degenerate cases, so use full solution:
        m = p - 1
        U = {i: [s[i], r['a'][i], r['b'][i], mp.mpc(0)] for i in range(m)}
        R = res_powers(O, U, p - 1)
        lead = [max(abs(x) for x in R[k][:3]) for k in range(1, p)]
        d = min(abs(r['detJ']), abs(r['detJ2']))
        mindet = min(mindet, d); minres = min([minres] + lead)
        print('  C=%s  |detJ1|,|detJ2|=%s,%s  max_h|res L^k|, k=1..%d: %s' % (
            [mp.nstr(x, 5) for x in s], mp.nstr(abs(r['detJ']), 5), mp.nstr(abs(r['detJ2']), 5), p - 1,
            [mp.nstr(x, 5) for x in lead]))
    print('(p,q)=(%d,%d): min |det J| = %s ; min max_{h<=2} |res L^k| = %s' % (p, q, mp.nstr(mindet, 5), mp.nstr(minres, 5)))
