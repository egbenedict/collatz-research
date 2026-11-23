"""m=137 after the B5(532) >= 557 witness: for a hypothetical set of extra exact values (args N:v), the
largest B5(532) that still closes m (binary search on 557..600; closing is monotone in the table).
Usage: tol137b.py N:v [N:v ...]"""
import io, contextlib, json, sys
from b5lib import *
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
base = dict(map(int, a.split(":")) for a in sys.argv[1:])
def closes(v):
    t = dict(real); t.update(base); t[532] = v
    with contextlib.redirect_stdout(io.StringIO()):
        fails, _ = M.check_m(checker(t), 137)
    return not fails
lo, hi = 557, 600
if not closes(lo): print(f"{base}: open even at B5(532)=557"); sys.exit()
while lo < hi:
    mid = (lo + hi + 1) // 2
    if closes(mid): lo = mid
    else: hi = mid - 1
print(f"{base}: m=137 closes iff B5(532) <= {lo}", flush=True)
