"""I-011: is there a sliding-window ("digit-additive") ranking function for the accelerated map?

W(n) = sum_{i=0}^{len(n)-1} w(s_i s_{i+1} ... s_{i+b-1}), where s = binary digits of n (LSB first) and
the string is padded with zeros above the top bit. We ask for an LP solution of
    W(T(n)) <= W(n) - 1   for all n in [3, 2^Lmax),       T(n) = n/2 or (3n+1)/2.
If no solution exists, no such ranking function exists for that window size b (for all n, a fortiori).
Control: the same LP for the map n -> n/2 (even), (n+1)/2 (odd), which terminates; it should be feasible.
"""
import sys
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix

def window_counts(n, b, nb):
    """vector c with W(n) = c . w"""
    bits = [(n >> i) & 1 for i in range(n.bit_length())] + [0] * (b - 1)
    c = np.zeros(nb)
    for i in range(n.bit_length()):
        idx = 0
        for t in range(b):
            idx |= bits[i + t] << t
        c[idx] += 1
    return c

def feasible(b, Lmax, step):
    nb = 1 << b
    ns = range(3, 1 << Lmax)
    rows = []
    for n in ns:
        m = step(n)
        if m <= 2 and m != 0:
            # reaching the trivial cycle region is fine: only require W(m) <= W(n) - 1 if m >= 3
            continue
        rows.append(window_counts(m, b, nb) - window_counts(n, b, nb))
    A = np.array(rows)
    bvec = -np.ones(len(rows))
    res = linprog(np.zeros(nb), A_ub=A, b_ub=bvec, bounds=[(-1e4, 1e4)] * nb, method="highs")
    return res.status == 0, len(rows)

collatz = lambda n: n // 2 if n % 2 == 0 else (3 * n + 1) // 2
control = lambda n: n // 2 if n % 2 == 0 else (n + 1) // 2

if __name__ == "__main__":
    Lmax = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    for b in range(1, 9):
        fc, k = feasible(b, Lmax, collatz)
        fk, _ = feasible(b, Lmax, control)
        print(f"b={b}: Collatz feasible={fc}   control (n+1)/2 feasible={fk}   ({k} constraints, n < 2^{Lmax})", flush=True)
