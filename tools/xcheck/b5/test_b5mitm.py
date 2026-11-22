"""Validate b5mitm against rises.py (which is validated against brute force) and direct brute force.
For each N: max total found at T = B5 equals B5, nothing found at T = B5 + 1, no UNRESOLVED,
witness re-simulated as an actual odd n < 2^N (rises.realise). Also runs with 4 threads and checks
the hit count equals the single-thread run."""
import re, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from rises import Solver, realise, brute

BIN = str(HERE / "b5mitm")

def run(N, T, KMAX, nw=1):
    best, sig, unres, hits = -1, None, False, 0
    for i in range(1):
        out = subprocess.run([BIN, str(N), str(T), str(KMAX), str(nw)],
                             capture_output=True, text=True, check=True).stdout
        unres |= "UNRESOLVED" in out
        s = re.search(r"hits=(\d+) bad=0 max_total=(-?\d+)(?: at \((.*)\))?", out)
        hits += int(s.group(1))
        if int(s.group(2)) > best:
            best, sig = int(s.group(2)), tuple(map(int, s.group(3).split(",")))
    return best, sig, unres, hits

Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 30
for N in range(6, Nmax + 1, 2):
    S = Solver(N)
    b4, b5 = S.B(4)[0], S.B(5)[0]
    if N <= 18:
        assert brute(N, 5) == b5, N
    if b5 + 4 - N < 1:
        print(f"skip N={N} (window < 1)"); continue
    got, sig, unres, h1 = run(N, b5, b4)
    assert not unres and got == b5, (N, got, b5)
    n = realise(N, 5, sig, b5)
    none, _, unres2, _ = run(N, b5 + 1, b4)
    assert not unres2 and none == -1, (N, none)
    got4, _, _, h4 = run(N, b5, b4, nw=4)
    assert got4 == b5 and h4 == h1, (N, got4, h4, h1)
    print(f"ok  N={N}: B5 = {b5} (rises.py{', brute force' if N <= 18 else ''}); C max {got}, none >= {b5+1}; "
          f"4 threads same ({h4} hits); witness n = {n} runs {sig}", flush=True)
