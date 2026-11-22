"""Brute-force check (exact integers) of the triangle lemmas of NOTES.md, threshold eps = 1/5.

Plane: s(m,j) = signed residue of 2^j mod 3^m, divided by 3^m, in (-1/2,1/2].  gray <=> |s| < 1/5.
Vertex: gray point with s(m,j-1) != s/2 and s(m+1,j) != s/3.
Checks, for 1 <= m <= M-1 (so that the right neighbour exists), j in a window [-W, W]:
  L3a  every gray point reaches a vertex by the parent rule and has the cone value 3^a 2^b delta exactly
  L3b  rows j_v-1, j_v-2 under the bottom edge (columns with 3^a|delta| < 1/5) are not gray; row j_v-1 = 1/2 + 3^a delta/2
  L3c  hypotenuse: a step (-2,-b), b>=2, from a gray point to a point still in the vertex's cone rows (b<=height)
       but outside T(v) lands on a non-gray point
  L4   separation: vertices v1 != v2 with m1<=m2, j1<=j2: 2^(j2-j1)|d1| >= 1/2 or 3^(m2-m1)|d2| >= 1/2
"""
import sys
from fractions import Fraction as Fr

M = int(sys.argv[1]) if len(sys.argv) > 1 else 9
W = int(sys.argv[2]) if len(sys.argv) > 2 else 40
EPS = Fr(1, 5)


def s(m, j):
    mod = 3 ** m
    r = pow(2, j, mod)
    if 2 * r > mod:
        r -= mod
    return Fr(r, mod)


S = {(m, j): s(m, j) for m in range(1, M + 2) for j in range(-W - 3, W + 40)}
gray = lambda p: abs(S[p]) < EPS
verts = {}
for m in range(1, M):
    for j in range(-W, W + 1):
        v = S[(m, j)]
        if abs(v) < EPS and S[(m, j - 1)] != v / 2 and S[(m + 1, j)] != v / 3:
            verts[(m, j)] = v

bad = 0
# L3a
for m in range(1, M - 1):
    for j in range(-W + 2, W - 2):
        if not gray((m, j)):
            continue
        p = (m, j)
        a = b = 0
        while p not in verts:
            v = S[p]
            if S[(p[0], p[1] - 1)] == v / 2:
                p = (p[0], p[1] - 1); b += 1
            elif p[0] + 1 <= M and S[(p[0] + 1, p[1])] == v / 3:
                p = (p[0] + 1, p[1]); a += 1
            else:
                bad += 1; print("L3a no parent", (m, j)); break
            if p[0] >= M or p[1] < -W:
                break
        else:
            if S[(m, j)] != 3 ** a * 2 ** b * verts[p]:
                bad += 1; print("L3a value", (m, j), p)
print(f"vertices: {len(verts)}   L3a done, bad={bad}")
# L3b
for (mv, jv), d in verts.items():
    a = 0
    while mv - a >= 1 and 3 ** a * abs(d) < EPS:
        x1, x2 = S[(mv - a, jv - 1)], S[(mv - a, jv - 2)]
        if gray((mv - a, jv - 1)) or gray((mv - a, jv - 2)) or (x1 - (Fr(1, 2) + 3 ** a * d / 2)) % 1 != 0:
            bad += 1; print("L3b", (mv, jv), a, x1, x2)
        a += 1
print(f"L3b done, bad={bad}")
# L3c
cnt = 0
for (mv, jv), d in verts.items():
    for a in range(0, mv):
        for h in range(0, 60):
            if 3 ** a * 2 ** h * abs(d) >= EPS:
                break
            for bb in range(2, h + 1):
                a2, h2 = a + 2, h - bb
                if mv - a2 < 1:
                    continue
                if 3 ** a2 * 2 ** h2 * abs(d) >= EPS:   # left T(v) but stayed in cone rows
                    cnt += 1
                    if gray((mv - a2, jv + h2)):
                        bad += 1; print("L3c", (mv, jv), a, h, bb)
print(f"L3c checked {cnt} hypotenuse exits, bad={bad}")
# L4
vl = list(verts.items())
pairs = 0
for (p1, d1) in vl:
    for (p2, d2) in vl:
        if p1 != p2 and p1[0] <= p2[0] and p1[1] <= p2[1]:
            pairs += 1
            if not (2 ** (p2[1] - p1[1]) * abs(d1) >= Fr(1, 2) or 3 ** (p2[0] - p1[0]) * abs(d2) >= Fr(1, 2)):
                bad += 1; print("L4", p1, p2, d1, d2)
print(f"L4 checked {pairs} ordered vertex pairs, bad={bad}")
print("ALL OK" if bad == 0 else f"FAILURES: {bad}")
