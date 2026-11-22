"""Helpers for the B5 study: a Checker with an extra rule r=5 (hypothetical or computed B5 table).
Imports the existing tools unchanged; rule r=5 plugs into mchain.Checker.chain because chain()
iterates over self.tab (rule r: Y_t <= Y_{t-r} + (delta-1) b_r(Y_{t-r}) - r + r eps, valid for any r
by the same block lemma)."""
import sys
from math import ceil
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mchain as M      # noqa: E402
import rigorous as R    # noqa: E402

ROOT = HERE.parent / "marshall_data"
HARD = (27444133206411171953, 77692117359936589403, 205632218873398596256)


def checker(b5: dict | None = None, mode="own"):
    ch = M.Checker(ROOT, mode)
    if b5:
        ch.tab[5] = M.Table(b5)
    return ch


def hard_trace(ch, m, K=HARD[0]):
    p = ceil(K * M.D.hi)
    Lam = ch.Lambda(K, p)
    Ub, m2 = min(((m - m2) / (Lam.lo - th), m2) for m2 in range(0, m)
                 if (th := ch.tail_hi(m, m2, K)) is not None and Lam.lo > th)
    Y0 = R.log2(ceil(Ub)).hi
    s, Y, use = ch.chain(m, Y0)
    return float(s / K), Y, use, m2
