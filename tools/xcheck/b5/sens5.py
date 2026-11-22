"""Task 1: which B5 values are needed? Hypothetical table B5(N) := B4(N) + c on B4's grid points N <= Nmax.
(True B5(N) >= B4(N) + 1 always, so c = 1 is the most optimistic possible table on that grid.)"""
import sys, io, contextlib
from b5lib import *
import rises
B4 = rises.load_marshall(ROOT)[4]

def closes(m, b5):
    ch = checker(b5)
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(ch, m)
    return not fails, worst[0]

def hyp(c, Nmax, Nmin=0):
    return {N: v + c for N, v in B4.items() if Nmin <= N <= Nmax}

if __name__ == "__main__":
    for m in map(int, sys.argv[1:]):
        r = [hard_trace(checker(hyp(c, 852)), m)[0] for c in (1, 3, 5, 8, 12, 20)]
        print(f"m={m}: hard-shape ratio with B5=B4+c on all B4 grid N<=852, c=1,3,5,8,12,20:",
              " ".join(f"{x:.3f}" for x in r), flush=True)
