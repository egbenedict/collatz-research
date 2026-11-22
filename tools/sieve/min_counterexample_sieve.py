"""Cycle 6: residue sieve on the minimal counterexample y* = min E (accelerated map T).
Forward: if the first k parity steps (fixed by y mod 2^k) force T^j(y) < y for some j <= k, then y != y*.
Backward: if y = T^j(x) with x < y (decidable from y mod 3^b via odd predecessors (2w-1)/3), then y != y*.
We count surviving classes mod 2^k * 3^b. The 2-adic/3-adic kill test predicts the surviving
density tends to 0 but no finite level is empty. This is a sanity check, not a new method.
"""
import math

def forward_survivors(k):
    """residues r mod 2^k (odd) not forced below themselves within k steps (Terras coefficient test)."""
    surv = []
    for r in range(1, 1 << k, 2):
        # track y = a*y0 + c symbolically with y0 = r + 2^k t; the parity of step j uses residue mod 2^{k-j}
        a_num, a_den, x = 1, 1, r
        ok = True
        for j in range(k):
            if x % 2: x = (3 * x + 1) // 2; a_num *= 3; a_den *= 2
            else: x //= 2; a_den *= 2
            if a_num < a_den:          # coefficient < 1: T^{j+1}(y) < y for all large y in class
                ok = False; break
        if ok: surv.append(r)
    return surv

def backward_excluded_mod3(y, depth=6):
    """True if some odd predecessor chain of length <= depth gives x < y (x = (2w-1)/3 type, w = 2 mod 3)."""
    w = y
    for _ in range(depth):
        if w % 3 != 2: return False
        x = (2 * w - 1) // 3
        if x < y: return True
        w = x
    return False

for k in (8, 12, 16, 20):
    fs = forward_survivors(k)
    # combine with mod 3 backward test on representatives in [2^k, 2^k * 3^6)
    tot = 0; keep = 0
    for r in fs[:2000]:
        for t in range(3 ** 4):
            y = r + (1 << k) * t + (1 << k) * 3 ** 4 * 10   # large representatives
            tot += 1
            if not backward_excluded_mod3(y): keep += 1
    print(f"k={k:2d}: forward survivors mod 2^k: {len(fs)} (density among odd {len(fs)/2**(k-1):.4f});"
          f" after backward mod-3 test: {keep/tot:.4f} of those kept  -> combined {len(fs)/2**(k-1)*keep/tot:.4f}")
