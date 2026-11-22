"""Cross-check the CNF encoding against the independent z3 bit-vector encoding on small queries:
both must agree on sat/unsat (SAT models are re-verified exactly by both).
z3 arctic uses w bits (entries < 2^w); cnf arctic uses max entry w, so compare arc with
cnf w = 2^w - 1.  Usage: .venv/bin/python crosscheck.py"""
import itertools
import cnf
import matrix
from srs import SYSTEMS, head_markers, reverse

def main():
  agree = 0
  cases = []
  for name in ["T-f1", "T-t>", "T-<0", "Wp", "T", "T5", "S3", "EVEN4-ff>"]:
      R = SYSTEMS[name]
      Bf, De = head_markers(R)
      cases.append((name, "top-fwd", Bf, [r for r in R if r not in Bf], "top"))
      H, RR = reverse(De), reverse(R)
      cases.append((name, "top-rev", H, [r for r in RR if r not in H], "top"))
      cases.append((name, "ext", R, [], "ext"))
  for (name, lab, H, W, mode), (kind, d, w) in itertools.product(
          cases, [("nat", 1, 1), ("nat", 1, 2), ("nat", 2, 1), ("arc", 1, 1), ("arc", 2, 1), ("arc", 1, 2)]):
      a = matrix.search(H, W, kind, mode, d, w, timeout_s=60)[0]
      b = cnf.search(H, W, kind, mode, d, w if kind == "nat" else 2 ** w - 1, timeout_s=60)[0]
      flag = "" if a == b else "   <-- DISAGREE"
      if a == b:
          agree += 1
      else:
          print(f"{name:10s} {lab:8s} {kind} d={d} w={w}: z3={a} cnf={b}{flag}")
  print(f"agree on {agree} / {len(cases) * 6} queries")


if __name__ == "__main__":
    main()
