"""Fall-weighted rise constants C_r(N) (see wtol137.py). Rigorous upper bounds:
 r=2: exact enumeration of all prefixes 1^k0 0^l0 (feasible), C2 = max Dh1*(z - l0) - l0 (k1 = z - k0 - l0 maximal).
 r=3: from b3mitm hits at threshold T: C3 <= max(max over hits of Dh1*sum k - (l0+l1), Dh1*(T-1) - 2).
Usage: weighted.py 2 N...   |  weighted.py 3 T3OFF N...  (T = B3(N) - T3OFF, B3 from Marshall, KMAX = b2(N))"""
import sys, re, subprocess
from fractions import Fraction as F
from b5lib import *
import rises
TB = rises.load_marshall(ROOT)
Dh1 = checker().Dh1
def C2(N):
    S = rises.Solver(N); M_ = S.M; best = None; arg = None
    for k0 in range(1, N + 1):
        uk = S.uK(k0)
        for l0 in range(1, 4 * N):
            nu = (M_ - 1 + uk * (1 - (1 << l0))) % M_
            z = S.z(nu)
            if k0 + l0 > z: break
            if z <= k0 + l0: continue
            v = Dh1 * (z - l0) - l0
            if best is None or v > best: best, arg = v, (k0, l0, z - k0 - l0)
    return best, arg
def C3(N, off):
    b3 = TB[3][N]; T = b3 - off; kmax = M.Table(TB[2]).b(N)
    LIM = 8 if T + 2 < 510 else 36
    exe = "../b3mitm" if LIM == 8 else "../b3mitm36"
    out = subprocess.run([exe, str(N), str(T), str(kmax)], capture_output=True, text=True).stdout
    assert "unresolved=0" in out, out[-300:]
    best = Dh1 * (T - 1) - 2; arg = "floor"
    for h in re.finditer(r"HIT total=(\d+) k0=(\d+) l0=(\d+) k1=(\d+) l1=(\d+)", out):
        t, k0, l0, k1, l1 = map(int, h.groups())
        v = Dh1 * t - l0 - l1
        if v > best: best, arg = v, (t, l0, l1)
    return best, arg, b3
if sys.argv[1] == "2":
    for N in map(int, sys.argv[2:]):
        b2 = TB[2].get(N); c, a = C2(N); std = Dh1 * b2 - 1 if b2 else None
        print(f"r=2 N={N} B2table={b2} C2={float(c):.3f} at {a}; standard {float(std):.3f}; gain {float(std - c):.3f}", flush=True)
else:
    off = int(sys.argv[2])
    for N in map(int, sys.argv[3:]):
        c, a, b3 = C3(N, off); std = Dh1 * b3 - 2
        print(f"r=3 N={N} B3={b3} C3<={float(c):.3f} ({a}); standard {float(std):.3f}; gain {float(std - c):.3f}", flush=True)
