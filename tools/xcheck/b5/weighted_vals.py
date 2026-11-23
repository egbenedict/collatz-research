"""Exact rational C_r upper bounds, recomputed from the saved weighted.py logs' arguments (hits re-derived
on demand is too slow, so the logs store (total, l-sum) or 'floor'; C = Dh1*total - lsum or Dh1*(T-1)-(r-1))."""
import re
from b5lib import *
import rises
def load():
    Dh1 = checker().Dh1; TB = rises.load_marshall(ROOT); out = {}
    for line in open(HERE / "logs/weighted_r3.log"):
        mm = re.match(r"r=3 N=(\d+) B3=(\d+) C3<=\S+ \(\(?(.*?)\)?\);", line)
        if not mm: continue
        N, b3, a = int(mm.group(1)), int(mm.group(2)), mm.group(3)
        if 'floor' in a: C = Dh1 * (b3 - 3 - 1) - 2
        else:
            t, l0, l1 = map(int, a.split(", ")); C = Dh1 * t - l0 - l1
        out.setdefault(3, {})[N] = C
    return out
