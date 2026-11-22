"""Compute exact B5(N) at given N with b5mitm. KMAX = B4 upper bound b4(N) from Marshall's table
(min over grid >= N), first confirmed by ../b4mitm (no n < 2^N with 4-rise total >= b4(N)+1, using
KMAX = b3(N), the B3 table being independently checked exact for N <= 400). Every maximiser is
re-simulated as an actual odd n < 2^N (rises.realise). Results appended to b5_values.json."""
import json, re, subprocess, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mchain as M
from rises import load_marshall, realise
tabs = load_marshall(HERE.parent / "marshall_data")
T3, T4 = M.Table(tabs[3]), M.Table(tabs[4])
OUT = HERE / "b5_values.json"
vals = json.loads(OUT.read_text()) if OUT.exists() else {}
NT = 8

def b4check(N, b4, b3):
    out = subprocess.run([str(HERE.parent / "b4mitm"), str(N), str(b4 + 1), str(b3)], capture_output=True, text=True, check=True).stdout
    assert "UNRESOLVED" not in out and "max_total=-1" in out, out[-300:]

for arg in sys.argv[1:]:
    N, Tforce = (int(arg.split(':')[0]), int(arg.split(':')[1])) if ':' in arg else (int(arg), None)
    if str(N) in vals:
        print(f"N={N}: cached {vals[str(N)]}"); continue
    b3, b4 = T3.b(N), T4.b(N)
    assert N <= 400
    t0 = time.time(); b4check(N, b4, b3); tb4 = time.time() - t0
    T = Tforce or b4 + 1       # N:T forces a lower start (must be <= B5(N) to get hits in one run)
    while True:
        t0 = time.time()
        p = subprocess.run([str(HERE / "b5mitm"), str(N), str(T), str(b4), str(NT)], capture_output=True, text=True)
        assert p.returncode == 0, p.stdout[-500:] + p.stderr
        out = p.stdout
        assert "UNRESOLVED" not in out
        s = re.search(r"left=(\d+) matches=(\d+) hits=(\d+) bad=0 max_total=(-?\d+)(?: at \((.*)\))?", out)
        if int(s.group(4)) >= 0:
            break
        T -= 1
    best = int(s.group(4)); sig = tuple(map(int, s.group(5).split(",")))
    n = realise(N, 5, sig, best)
    rec = dict(B5=best, B4tab=b4, B4_confirmed_by_b4mitm=True, T=T, KMAX=b4, left=int(s.group(1)), hits=int(s.group(3)),
               runs=sig, witness=str(n), secs=round(time.time() - t0, 1), b4check_secs=round(tb4, 1))
    vals[str(N)] = rec
    OUT.write_text(json.dumps(vals, indent=1, sort_keys=True))
    print(f"N={N}: B5 = {best} (B4 <= {b4}, +{best-b4}); runs {sig}; left {rec['left']:.3e}; {rec['secs']}s (+{tb4:.0f}s b4 check)", flush=True)
