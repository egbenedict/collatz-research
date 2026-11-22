"""Minimal stepwise termination prover (rule removal with matrix interpretations).

Proof scheme for a system R over {<, f, t, digits, > or *} (YAH Lemmas 3.15, 3.16, 3.18 + Thm 2.6):
  ob1: SN(R \\ Bf)   where Bf = rules '<s -> <t'           (full termination, ext mode)
  ob2: SN(R \\ De)   where De = rules 's> -> t>' / 's* -> t*' (full termination, ext mode)
       ob1 + ob2 => only canonical strings <w> can be nonterminating (type introduction).
  ob3 (forward): SN(Bf_top / R)                         (top mode, rules stay as weak)
      or (reverse): SN(rev(De)_top / rev(R))
  ob3 => SN(Bf / R) (resp. SN(De/R)) by Lemma 3.18; with ob1 (resp. ob2) and rule removal,
  SN(R).
A step may also remove rules in ext mode (relative termination, any position).

Usage:
  .venv/bin/python prover.py SYSTEM [--dir fwd|rev|both] [--dmax D] [--wmax W]
         [--kinds nat,arc] [--timeout SEC] [--first-step-only] [--configs nat:3:2,arc:4:3]
"""
import argparse
import json
import os
import sys
import time

import cnf
import matrix

BACKEND = {"name": "cnf", "jobs": 1}


def search(strict, weak, kind, mode, d, w, timeout_s):
    if BACKEND["name"] == "z3":
        st, interp, removed, dt = matrix.search(strict, weak, kind, mode, d, w, timeout_s)
        return st, interp, removed, dt, {}
    return cnf.search(strict, weak, kind, mode, d, w, timeout_s, jobs=BACKEND["jobs"])
from srs import SYSTEMS, YAH_STATUS, head_markers, reverse

HERE = os.path.dirname(os.path.abspath(__file__))


EXT_IN_TOP = [False]


def default_configs(kinds, dmax, wmax):
    cfg = []
    for d in range(1, dmax + 1):
        for w in range(1, wmax + 1):
            for k in kinds:
                cfg.append((k, d, w))
    return cfg


def log(msg):
    print(msg, flush=True)


def stepwise(strict, weak, mode, configs, timeout, label, allow_ext_in_top=False,
             first_step_only=False):
    """Remove all rules of `strict`. mode 'ext' with weak == []: full termination of strict.
    mode 'top': SN(strict_top / strict u weak). Returns list of steps or None."""
    strict = list(strict)
    weak = list(weak)
    steps = []
    full = (mode == "ext" and not weak)
    while strict:
        found = False
        tries = [(k, d, w, mode) for k, d, w in configs]
        if mode == "top" and allow_ext_in_top:
            tries += [(k, d, w, "ext") for k, d, w in configs]
        for kind, d, w, md in tries:
            st, interp, removed, dt, info = search(strict, weak if not full else [], kind, md, d, w,
                                                   timeout_s=timeout)
            log(f"    [{label}] {md:3s} {kind} d={d} w={w}: {st} ({dt:.1f}s"
                + (f", {info['vars']}v/{info['clauses']}c, {info['solver']}" if info else "") + ")"
                + (f" removes {[l + '->' + r for l, r in removed]}" if st == "sat" else ""))
            if st == "sat":
                steps.append({"kind": kind, "mode": md, "d": d, "w": w,
                              "strict": strict[:], "weak": weak[:] if not full else [],
                              "removed": removed,
                              "interp": {c: [M, v] for c, (M, v) in interp.items()}})
                strict = [r for r in strict if r not in removed]
                if not full:
                    weak = weak + removed
                found = True
                break
        if not found:
            return None, steps
        if first_step_only:
            return "partial", steps
    return "done", steps


def prove(name, direction, configs, timeout, first_step_only=False):
    R = SYSTEMS[name]
    Bf, De = head_markers(R)
    log(f"System {name}: {len(R)} rules; YAH status: {YAH_STATUS.get(name, '?')}")
    proof = {"system": name, "rules": R, "obligations": {}}
    if not first_step_only:
        for ob, sub in (("ob1", [r for r in R if r not in Bf]), ("ob2", [r for r in R if r not in De])):
            log(f"  {ob}: SN of {len(sub)} rules (R minus {'Bf' if ob == 'ob1' else 'De'})")
            status, steps = stepwise(sub, [], "ext", configs, timeout, ob)
            if status != "done":
                log(f"  {ob} FAILED")
                return False, proof
            proof["obligations"][ob] = {"rules": sub, "steps": steps}
    dirs = ["fwd", "rev"] if direction == "both" else [direction]
    for dr in dirs:
        if dr == "fwd":
            H, all_rules = Bf, R
        else:
            H, all_rules = reverse(De), reverse(R)
        weak = [r for r in all_rules if r not in H]
        log(f"  ob3-{dr}: SN(H_top / R) with H = {[l + '->' + r for l, r in H]}")
        status, steps = stepwise(H, weak, "top", configs, timeout, "ob3-" + dr,
                                 allow_ext_in_top=EXT_IN_TOP[0], first_step_only=first_step_only)
        if status == "done":
            proof["obligations"]["ob3"] = {"dir": dr, "H": H, "R": all_rules, "steps": steps}
            log(f"  ==> SN({name}) PROVED (direction {dr})")
            return True, proof
        if status == "partial":
            proof["obligations"]["ob3-partial"] = {"dir": dr, "H": H, "R": all_rules,
                                                   "steps": steps}
            log(f"  ==> first step found ({dr}): a rule was removed")
            return "partial", proof
        log(f"  ob3-{dr}: no proof within configs")
    return False, proof


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("system")
    ap.add_argument("--dir", default="both")
    ap.add_argument("--dmax", type=int, default=3)
    ap.add_argument("--wmax", type=int, default=2)
    ap.add_argument("--kinds", default="nat,arc")
    ap.add_argument("--configs", default=None, help="explicit list kind:d:w,...")
    ap.add_argument("--timeout", type=float, default=60)
    ap.add_argument("--first-step-only", action="store_true")
    ap.add_argument("--backend", default="cnf", choices=["cnf", "z3"])
    ap.add_argument("--jobs", type=int, default=1, help="parallel SAT portfolio size")
    ap.add_argument("--ext-in-top", action="store_true",
                    help="also try ext-mode (relative) steps while removing ob3 rules")
    a = ap.parse_args()
    BACKEND["name"], BACKEND["jobs"] = a.backend, a.jobs
    EXT_IN_TOP[0] = a.ext_in_top
    if a.configs:
        configs = [(c.split(":")[0], int(c.split(":")[1]), int(c.split(":")[2]))
                   for c in a.configs.split(",")]
    else:
        configs = default_configs(a.kinds.split(","), a.dmax, a.wmax)
    t0 = time.time()
    ok, proof = prove(a.system, a.dir, configs, a.timeout, a.first_step_only)
    proof["result"] = ok
    proof["seconds"] = time.time() - t0
    os.makedirs(os.path.join(HERE, "proofs"), exist_ok=True)
    fn = os.path.join(HERE, "proofs", a.system.replace(">", "E").replace("*", "S").replace("<", "B")
                      + ".json")
    if ok is True:
        with open(fn, "w") as fh:
            json.dump(proof, fh)
        log(f"  proof written to {fn} ({proof['seconds']:.1f}s total)")
    else:
        log(f"  result: {ok} ({proof['seconds']:.1f}s total)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
