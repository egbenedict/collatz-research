"""I-010: can a Sturmian parity word be the parity sequence of a divergent Collatz orbit?

Lemma: if the parity words of z_n and z_{n+P} agree for L steps, then z_n = z_{n+P} mod 2^L. A
divergent orbit never repeats, so 2^L <= max(z_n, z_{n+P}) <= (N+1) * prod over steps of growth.
For odd density alpha, log2 z_m <= log2(N+1) + o_m*log2(3) - m + o(1), where o_m = #odd steps
among the first m.

So the excess  g(n,P) = L(n,P) - max_{m in {n, n+P}} (o_m*log2 3 - m)  must stay <= log2(N+1)+O(1).
If sup g = +infinity, the word is not the parity sequence of any divergent orbit.
We compute max g over n + P <= M for growing M, for mechanical words s_n = floor((n+1)a+r) - floor(na+r).
"""
import math, sys
import numpy as np

L3 = math.log2(3)

def word(a, r, M):
    n = np.arange(M + 1)
    return (np.floor((n + 1) * a + r) - np.floor(n * a + r)).astype(np.int8)

def max_excess(a, r, M):
    s = word(a, r, 2 * M + 10)
    o = np.concatenate([[0], np.cumsum(s)])
    growth = o * L3 - np.arange(len(o))          # log2 z_m - log2 N, roughly
    best = -1e9
    for P in range(1, M):
        # LCP of s[n:] and s[n+P:] for all n, via run lengths of equality
        eq = (s[:-P] == s[P:])
        # run length of trues starting at each index, computed backward
        run = np.zeros(len(eq) + 1, dtype=np.int64)
        for i in range(len(eq) - 1, -1, -1):
            run[i] = run[i + 1] + 1 if eq[i] else 0
        nmax = M - P
        if nmax <= 0:
            break
        n = np.arange(nmax)
        Ls = run[:nmax]
        g = Ls - np.maximum(growth[n], growth[n + P])
        best = max(best, float(g.max()))
    return best

if __name__ == "__main__":
    for a in (0.632, 0.64, 0.66, 0.70, 0.75, 0.80, 0.90, 0.95, (math.sqrt(5) - 1) / 2 + 0.01):
        eps = a * L3 - 1
        row = [max_excess(a, r, M) for M in (250, 500, 1000, 2000) for r in (0.0,)]
        print(f"alpha={a:.4f} eps={eps:.4f}  max excess at M=250,500,1000,2000: " + "  ".join(f"{x:7.1f}" for x in row), flush=True)
