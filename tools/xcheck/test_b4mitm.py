"""Validate the C meet-in-the-middle solver against rises.py (which matches brute force)."""
import re
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rises import Solver, realise

BIN = str(Path(__file__).resolve().parent / "b4mitm")

def run(N, T, KMAX):
    out = subprocess.run([BIN, str(N), str(T), str(KMAX)], capture_output=True, text=True, check=True).stdout
    s = re.search(r"max_total=(-?\d+)(?: at \((.*)\))?", out)
    unresolved = "UNRESOLVED" in out
    sig = tuple(map(int, s.group(2).split(","))) if s.group(2) else None
    return int(s.group(1)), sig, unresolved

Ns = list(range(14, 41, 2)) + [72, 80]
known = {72: 86, 80: 95}          # exact values from rises.py (logs/b4_compare.log)
for N in Ns:
    b3 = Solver(N).B(3)[0]
    b4 = known.get(N) or Solver(N).B(4)[0]
    got, sig, unres = run(N, b4, b3)
    assert not unres, N
    assert got == b4, (N, got, b4)
    n = realise(N, 4, sig, b4)          # an actual odd n < 2^N whose first 4 rises total >= b4
    none, _, unres2 = run(N, b4 + 1, b3)
    assert not unres2 and none == -1, (N, none)
    print(f"ok  N={N}: C finds max {got} = B4 (python), and nothing >= {b4 + 1}; witness runs {sig}", flush=True)
