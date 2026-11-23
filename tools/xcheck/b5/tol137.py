"""m=137: for each value V328 of B5(328) (>= B5(326) = 355, monotone in N), the largest B5(532) that closes."""
import io, contextlib, json
from b5lib import *
import rises
T4 = M.Table(rises.load_marshall(ROOT)[4])
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
def run(v328, v532):
    t = dict(real); t[328] = v328; t[532] = v532
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t), 137)
    return not fails
b = T4.b(532)
for v in range(355, 362):
    lo, hi = b, b + 24
    if not run(v, lo + 1): print(f"B5(328)={v}: no B5(532) closes", flush=True); continue
    a, c = lo + 1, hi
    if run(v, c): print(f"B5(328)={v}: closes even with B5(532)={c}", flush=True); continue
    while c - a > 1:
        mid = (a + c) // 2; a, c = (mid, c) if run(v, mid) else (a, mid)
    print(f"B5(328)={v}: closes iff B5(532) <= {a} (= b4 + {a - b}, b4(532)={b})", flush=True)
