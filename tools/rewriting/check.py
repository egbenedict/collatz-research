"""Independent checker for proofs/*.json produced by prover.py.

For each proof it checks
  * the bookkeeping: ob1 removes every rule of R \\ Bf, ob2 every rule of R \\ De, ob3 every
    rule of H (forward: H = Bf over R; reverse: H = rev(De) over rev(R)), each step's
    strict/weak sets being exactly what remains;
  * every step's interpretation with exact integer matrix arithmetic (matrix.verify);
  * every step again by *direct evaluation*: random vectors x and random context strings p, q,
    compare [p l q](x) with [p r q](x) computed by applying the symbol functions one by one
    (no matrix products), against the order of the mode (ext: in any context; top: p empty).

Usage: .venv/bin/python check.py proofs/NAME.json [...]
"""
import json
import random
import sys

from matrix import verify, NEG
from srs import SYSTEMS, head_markers, reverse


def T(rs):
    return [tuple(r) for r in rs]


def apply_fn(interp, kind, c, x):
    M, v = interp[c]
    d = len(v)
    if kind == "nat":
        return [sum(M[i][k] * x[k] for k in range(d)) + v[i] for i in range(d)]
    out = []
    for i in range(d):
        terms = [M[i][k] + x[k] for k in range(d) if M[i][k] is not NEG and x[k] is not NEG]
        if v[i] is not NEG:
            terms.append(v[i])
        out.append(max(terms) if terms else NEG)
    return out


def ev(interp, kind, s, x):
    for c in reversed(s):
        x = apply_fn(interp, kind, c, x)
    return x


def gt(kind, mode, a, b):
    if kind == "nat":
        return a[0] > b[0] and all(x >= y for x, y in zip(a, b))
    ge = all(y is NEG or (x is not NEG and x >= y) for x, y in zip(a, b))
    if mode == "ext":
        return ge and all((x is NEG and y is NEG) or (x is not NEG and (y is NEG or x > y))
                          for x, y in zip(a, b))
    return ge and a[0] is not NEG and (b[0] is NEG or a[0] > b[0])


def ge(kind, a, b):
    if kind == "nat":
        return all(x >= y for x, y in zip(a, b))
    return all(y is NEG or (x is not NEG and x >= y) for x, y in zip(a, b))


def rand_vec(kind, mode, d, rng):
    if kind == "nat":
        return [rng.randrange(0, 50) for _ in range(d)]
    x = [rng.choice([NEG, rng.randrange(0, 50)]) for _ in range(d)]
    x[0] = rng.randrange(0, 50)          # domain N x A^{d-1}
    return x


def numeric_check(step, alphabet, rng, trials=300):
    interp = {c: (M, v) for c, (M, v) in step["interp"].items()}
    kind, mode, d = step["kind"], step["mode"], step["d"]
    removed = T(step["removed"])
    alphabet = sorted(interp)  # contexts over every interpreted symbol (markers included)
    allr = T(step["strict"]) + T(step["weak"])
    for l, r in allr:
        for _ in range(trials):
            x = rand_vec(kind, mode, d, rng)
            q = "".join(rng.choice(alphabet) for _ in range(rng.randrange(0, 4)))
            top_strict = mode == "top" and (l, r) in removed
            p = "" if top_strict else "".join(rng.choice(alphabet) for _ in range(rng.randrange(0, 4)))
            a = ev(interp, kind, p + l + q, x)
            b = ev(interp, kind, p + r + q, x)
            if (l, r) in removed:
                assert gt(kind, mode, a, b), ("strict fails", l, r, p, q, x, a, b)
            else:
                assert ge(kind, a, b), ("weak fails", l, r, p, q, x, a, b)


def check_chain(rules, steps, mode, alphabet, rng, full):
    strict = list(rules) if full else None
    for st in steps:
        if full:
            assert sorted(T(st["strict"])) == sorted(strict), "strict set mismatch"
        interp = {c: (M, v) for c, (M, v) in st["interp"].items()}
        got = verify(interp, st["kind"], st["mode"], st["d"], [], T(st["strict"]) + T(st["weak"]))
        for r in T(st["removed"]):
            assert r in got and r in T(st["strict"]), ("claimed removal not strict", r)
        numeric_check(st, alphabet, rng)
        if full:
            strict = [r for r in strict if r not in T(st["removed"])]
    return strict


def check(fn):
    P = json.load(open(fn))
    name = P["system"]
    R = SYSTEMS[name]
    assert T(P["rules"]) == R, "rule list differs from srs.py"
    Bf, De = head_markers(R)
    rng = random.Random(0)
    alphabet = sorted({c for l, r in R for c in l + r} - set("<>*"))
    ob = P["obligations"]
    for key, sub in (("ob1", [r for r in R if r not in Bf]), ("ob2", [r for r in R if r not in De])):
        assert sorted(T(ob[key]["rules"])) == sorted(sub)
        left = check_chain(sub, ob[key]["steps"], "ext", alphabet, rng, full=True)
        assert not left, (key, "rules left", left)
    o3 = ob["ob3"]
    if o3["dir"] == "fwd":
        H, allr = Bf, R
    else:
        H, allr = reverse(De), reverse(R)
    assert sorted(T(o3["H"])) == sorted(H) and sorted(T(o3["R"])) == sorted(allr)
    strict = list(H)
    for st in o3["steps"]:
        assert sorted(T(st["strict"])) == sorted(strict)
        assert sorted(T(st["strict"]) + T(st["weak"])) == sorted(allr), "weak set incomplete"
        assert st["mode"] == "top" or st["mode"] == "ext"
        check_chain(None, [st], st["mode"], alphabet, rng, full=False)
        strict = [r for r in strict if r not in T(st["removed"])]
    assert not strict, ("ob3 rules left", strict)
    nsteps = sum(len(ob[k]["steps"]) for k in ("ob1", "ob2", "ob3"))
    big = max((s["kind"], s["d"], s["w"]) for s in o3["steps"])
    print(f"OK  {name}: SN verified ({nsteps} steps; ob3 {o3['dir']}, "
          f"ob3 steps {[(s['kind'], s['d'], s['w']) for s in o3['steps']]})")


if __name__ == "__main__":
    for fn in sys.argv[1:]:
        check(fn)
