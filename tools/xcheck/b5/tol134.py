"""m=134 needs only B5(73) and B5(121): how much slack does each have (alone and jointly)?"""
import io, contextlib
from b5lib import *
def run(t):
    with contextlib.redirect_stdout(io.StringIO()):
        fails, worst = M.check_m(checker(t), 134)
    return not fails, worst[0]
for a in range(0, 4):
    for b in range(0, 4):
        ok, w = run({73: 91 + a, 121: 143 + b})
        print(f"B5(73)={91+a} B5(121)={143+b}: {'closed' if ok else 'OPEN'} worst {w:.4f}", flush=True)
