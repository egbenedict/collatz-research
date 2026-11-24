"""Final verdict for m = 137 from the T = 557 B5(532) search (possibly split over several Actions runs).
Usage: verdict532.py DIR [DIR ...]   (each DIR contains shard-*/shard_*.txt artifacts)
Checks: every shard index 0..NS-1 present exactly once, bad = 0, no UNRESOLVED, per-shard hit counts = HIT lines.
Then: exact B5(532) (if >= 557), rigorous C5(532) = max(max over hits of Dh1*sum k - sum l0..l3, Dh1*556 - 4),
and the m = 137 chain check with b_eff = (C5 + 4)/Dh1 at N = 532 (everything else = real table)."""
import sys, re, glob, io, contextlib, json
from collections import Counter
from b5lib import *
files = [f for d in sys.argv[1:] for f in glob.glob(f"{d}/**/shard_*.txt", recursive=True)]
seen = Counter(); NS = None; bad = 0; hits = []; want_hits = 0; unres = 0
for f in files:
    txt = open(f).read()
    unres += txt.count("UNRESOLVED")
    for s in re.finditer(r"SUMMARY N=532 T=557 KMAX=554 .*shard=(\d+)/(\d+) .* hits=(\d+) bad=(\d+)", txt):
        seen[int(s.group(1))] += 1; NS = int(s.group(2)); want_hits += int(s.group(3)); bad += int(s.group(4))
    for line in txt.splitlines():
        if line.startswith("HIT"): hits.append(dict((a, int(b)) for a, b in re.findall(r"(\w+)=(\d+)", line)))
missing = [i for i in range(NS or 0) if seen[i] == 0]; dup = [i for i, c in seen.items() if c > 1]
print(f"files {len(files)}, NS {NS}, distinct shards {len(seen)}, missing {missing[:10]}{'...' if len(missing) > 10 else ''}, duplicates {dup}")
print(f"bad {bad}, UNRESOLVED {unres}, HIT lines {len(hits)}, summary hit total {want_hits}")
ok = NS == 256 and not missing and not dup and bad == 0 and unres == 0 and len(hits) == want_hits
ch0 = checker(); Dh1 = ch0.Dh1
B5 = max(h["total"] for h in hits) if hits else None
vals = sorted(((Dh1 * h["total"] - sum(h[f"l{i}"] for i in range(4)), h) for h in hits), key=lambda x: x[0], reverse=True)
floor = Dh1 * 556 - 4
C5 = max(vals[0][0], floor) if vals else floor
print(f"B5(532) = {B5} (exact if >= 557); C5(532) <= {float(C5):.4f} (floor {float(floor):.4f}); top pattern {vals[0][1] if vals else None}")
real = {int(k): v["B5"] for k, v in json.loads((HERE / "b5_values.json").read_text()).items()}
real[532] = (C5 + 4) / Dh1
with contextlib.redirect_stdout(io.StringIO()):
    fails, worst = M.check_m(checker(real), 137)
print(f"m = 137 with weighted C5(532): {'CLOSED' if not fails else 'OPEN'} (worst {float(worst[0]):.4f}); data checks {'PASS' if ok else 'FAIL'}")
