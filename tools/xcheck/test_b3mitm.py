"""b3mitm vs rises.py for small N: max at T = B3 must equal B3, nothing at T = B3+1, witnesses realised."""
import subprocess, sys, re
import rises
ok = True
for N in range(int(sys.argv[1]), int(sys.argv[2]) + 1, int(sys.argv[3]) if len(sys.argv) > 3 else 1):
    S = rises.Solver(N); b3, _ = S.B(3); b2, _ = S.B(2)
    out = subprocess.run(["./b3mitm", str(N), str(b3), str(b2)], capture_output=True, text=True).stdout
    mx = int(re.search(r"max_total=(-?\d+)", out).group(1))
    out2 = subprocess.run(["./b3mitm", str(N), str(b3 + 1), str(b2)], capture_output=True, text=True).stdout
    mx2 = int(re.search(r"max_total=(-?\d+)", out2).group(1))
    for h in re.finditer(r"HIT total=(\d+) k0=(\d+) l0=(\d+) k1=(\d+) l1=(\d+)", out):
        t, k0, l0, k1, l1 = map(int, h.groups())
        rises.realise(N, 3, (k0, l0, k1, l1), t)
    good = mx == b3 and mx2 == -1 and "UNRESOLVED" not in out + out2
    ok &= good
    print(f"N={N} B3={b3} b3mitm max={mx} at T+1: {mx2} {'OK' if good else 'FAIL'}", flush=True)
print("ALL OK" if ok else "FAILURES")
