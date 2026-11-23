"""b6mitm vs rises.py brute force: max at T = B6 equals B6, nothing at T = B6+1, every witness realised,
1-thread and 4-thread hit counts agree. Usage: test_b6mitm.py Nlo Nhi [step]"""
import subprocess, sys, re
sys.path.insert(0, "..")
import rises
ok = True
for N in range(int(sys.argv[1]), int(sys.argv[2]) + 1, int(sys.argv[3]) if len(sys.argv) > 3 else 1):
    S = rises.Solver(N); b6, _ = S.B(6); b5, _ = S.B(5)
    run = lambda T, th: subprocess.run(["./b6mitm", str(N), str(T), str(b5), str(th)], capture_output=True, text=True).stdout
    o1, o4, o2 = run(b6, 1), run(b6, 4), run(b6 + 1, 1)
    mx = int(re.search(r"max_total=(-?\d+)", o1).group(1)); mx2 = int(re.search(r"max_total=(-?\d+)", o2).group(1))
    h1, h4 = o1.count("HIT"), o4.count("HIT")
    for h in re.finditer(r"HIT total=(\d+) k0=(\d+) l0=(\d+) k1=(\d+) l1=(\d+) k2=(\d+) l2=(\d+) k3=(\d+) l3=(\d+) k4=(\d+) l4=(\d+)", o1):
        v = list(map(int, h.groups())); rises.realise(N, 6, tuple(v[1:]), v[0])
    good = mx == b6 and mx2 == -1 and h1 == h4 and "UNRESOLVED" not in o1 + o2 and " bad=0 " in o1 and " bad=0 " in o2
    ok &= good
    print(f"N={N} B5={b5} B6={b6} b6mitm max={mx} (+1: {mx2}) hits 1t/4t {h1}/{h4} {'OK' if good else 'FAIL'}", flush=True)
print("ALL OK" if ok else "FAILURES")
