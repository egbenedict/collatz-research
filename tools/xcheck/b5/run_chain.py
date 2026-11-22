"""Rerun the full chain + shape search (mchain.check_m) with the computed B5 table added as rule r=5."""
import json, sys, time
from b5lib import *
vals = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
mode = "own"
args = sys.argv[1:]
if args and args[0] in ("own", "compat"):
    mode, args = args[0], args[1:]
print(f"B5 table ({mode}):", dict(sorted(vals.items())))
ch = checker(vals, mode)
for m in map(int, args):
    r, Y, use, m2 = hard_trace(ch, m)
    print(f"m={m}: hard shape K=2.744e19 ratio {r:.7f}, rule use {use}", flush=True)
    M.check_m(ch, m)
