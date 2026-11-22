#!/usr/bin/env python3
"""Run the (K, L, m) sieve for every surviving K from `cycle_bounds.py survivors`.

For each m in [mmin, mmax], list the K values that pass every necessary condition
(Hercher Thm 21 + continued fractions + the Simons-de Weger upper bound), then try to rule
out each one with kstar_sieve. Prints one line per (m, K).

Usage: python3 run_kstar_all.py --x0 2075 --mmin 92 --mmax 98
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cycle_bounds as cb  # noqa: E402
import kstar_sieve as ks   # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--x0", default="2075")
    p.add_argument("--mmin", type=int, default=92)
    p.add_argument("--mmax", type=int, default=98)
    a = p.parse_args()
    X0_dec = cb.x0_from(a.x0)
    X0_int = int(a.x0) * 2 ** 60
    all_ok = True
    for m in range(a.mmin, a.mmax + 1):
        Klo, U, surv = cb.surviving_K(m, X0_dec)
        if not surv:
            print(f"m={m}: no K survives the continued-fraction stage", flush=True)
            continue
        for K, L, gap, _ in surv:
            t = time.time()
            try:
                P = ks.Problem(m, K, L, X0_int)
            except RuntimeError as e:
                print(f"m={m} K={K:,}: window did not close ({e})", flush=True)
                all_ok = False
                continue
            st = ks.sieve(P)
            ns = len(st["survivors"])
            all_ok &= ns == 0
            print(f"m={m} K={K:,} L={L:,}: y1 in [{float(P.y1_lo):.3f}, {float(P.y1_hi):.3f}], "
                  f"{st['candidates']} candidates, {ns} survivors  [{time.time()-t:.0f}s]  "
                  f"{dict(st['eliminated'])}", flush=True)
            for n, why in st["survivors"][:5]:
                print(f"    survivor n_1 = {n} ({why})", flush=True)
    print("ALL RULED OUT" if all_ok else "SOME CASES REMAIN", flush=True)


if __name__ == "__main__":
    main()
