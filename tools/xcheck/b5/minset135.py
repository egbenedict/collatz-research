"""m=135: greedy minimal set of expensive (N > 132) B5 points, hypothetical value b4(N)+c, full check_m."""
import io, contextlib, json, sys
from b5lib import *
import rises
T4 = M.Table(rises.load_marshall(ROOT)[4])
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
c = int(sys.argv[1])
cheap = [122, 128]                      # cheap (N <= 132) points still uncomputed: hypothetical too
S = [189, 192, 193, 195, 197, 198, 199, 203, 205, 210, 213, 220]
def run(pts):
    t = {N: T4.b(N) + c for N in pts + cheap}
    t.update(real)
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t), 135)
    return not fails, worst[0]
ok, w = run(S); print(f"c={c} full set {S}: {'closed' if ok else 'OPEN'} worst {w:.4f}", flush=True)
if ok:
    for N in list(S):
        trial = [x for x in S if x != N]
        ok2, w2 = run(trial)
        if ok2:
            S = trial
        print(f"  drop {N}: {'closed' if ok2 else 'needed'} (worst {w2:.4f})", flush=True)
    print(f"c={c} minimal greedy set: {S} -> hypothetical values {[T4.b(N)+c for N in S]}")
