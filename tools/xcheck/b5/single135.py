"""m=135 with ONE expensive B5 point at N (value V) + computed cheap points (+ hypothetical b4+c at
uncomputed cheap points 122, 128): largest V that still closes m=135, for each candidate N."""
import io, contextlib, json, sys
from b5lib import *
import rises
T4 = M.Table(rises.load_marshall(ROOT)[4])
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
c = int(sys.argv[1]); m = int(sys.argv[2]); Ns = [int(x) for x in sys.argv[3:]]
cheap = {N: T4.b(N) + c for N in (122, 128) if N not in real}
def run(N, V):
    t = dict(cheap); t[N] = V; t.update(real)
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t), m)
    return not fails
for N in Ns:
    lo, hi = T4.b(N), T4.b(N) + 16          # V = lo surely <= true B5 - 1; search largest closing V
    if not run(N, lo + 1):
        print(f"m={m} N={N}: even B5 = b4+1 = {lo+1} does not close", flush=True); continue
    a, b = lo + 1, hi
    if run(N, b):
        print(f"m={m} N={N}: closes even with B5 = {b}", flush=True); continue
    while b - a > 1:
        mid = (a + b) // 2
        a, b = (mid, b) if run(N, mid) else (a, mid)
    print(f"m={m} N={N}: b4(N) = {lo}; closes iff B5({N}) <= {a} (= b4 + {a-lo})", flush=True)
