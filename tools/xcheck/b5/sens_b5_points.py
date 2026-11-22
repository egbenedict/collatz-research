"""m=134 robustness: raise each computed B5 entry by +1 (one at a time) and rerun the full m=134 check."""
import io, contextlib, json
from b5lib import *
base = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items() if int(k) <= 121}
for N in sorted(base):
    t = dict(base); t[N] += 1
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t), 134)
    print(f"B5({N}) = {base[N]} -> {base[N]+1}: worst {worst[0]:.4f}, {'OPENS m=134' if fails else 'still closed'}", flush=True)
t = {N: v for N, v in base.items()}
for drop in sorted(base):
    t2 = {N: v for N, v in base.items() if N != drop}
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t2), 134)
    print(f"without B5({drop}): worst {worst[0]:.4f}, {'OPENS' if fails else 'still closed'}", flush=True)
