import sys, time
sys.path.insert(0, '..')
from rises import Solver, brute
for N in range(4, int(sys.argv[1]) + 1, 2):
    t = time.time(); S = Solver(N)
    b = [S.B(r)[0] for r in (2, 3, 4, 5)]
    print(N, "B2..B5 =", b, "B5-B4 =", b[3] - b[2], f"[{time.time()-t:.1f}s]", flush=True)
