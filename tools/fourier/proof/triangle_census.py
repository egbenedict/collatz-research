"""Census of gray triangles (eps=1/5) in the 2-adic window m <= NMAX, -J0 <= J < 0, below the main triangle.

Vertex (m,J): |s| < 1/5, s(m,J-1) = s/2 + 1/2 (not s/2), s(m+1,J) != s/3.  Height H = log2(1/(5|s|)) (rows).
Prints the number of vertices with H >= h, against the 'generic' prediction (area) * (2/5) * P(H >= h | vertex)
-- for a random real s the vertex density with H >= h is ~ (1/5)*(1/2)*(2/3)*2^-h*2 per point (see NOTES).
Usage: python3 triangle_census.py NMAX J0
"""
import sys
import math
from collections import Counter


def sgn(r, mod):
    return (r - mod) / mod if 2 * r > mod else r / mod


def main():
    nmax, J0 = int(sys.argv[1]), int(sys.argv[2])
    rows = {}
    for m in range(1, nmax + 2):
        mod = 3 ** m
        r = pow(2, -J0 - 1, mod)
        row = []
        for _ in range(J0 + 1):
            row.append(sgn(r, mod))
            r = 2 * r % mod
        rows[m] = row          # index i <-> J = -J0-1+i
    H = Counter()
    best = []
    area = 0
    for m in range(1, nmax + 1):
        for i in range(1, J0 + 1):
            area += 1
            s = rows[m][i]
            if abs(s) >= 0.2:
                continue
            if abs(rows[m][i - 1] - s / 2) < 1e-9:
                continue
            if abs(rows[m + 1][i] - s / 3) < 1e-9:
                continue
            h = math.log2(1 / (5 * abs(s)))
            H[int(h)] += 1
            best.append((h, m, i - J0 - 1))
    best.sort(reverse=True)
    tot = sum(H.values())
    print(f"window: m<={nmax}, {-J0}<=J<0, area={area}, vertices={tot} (density {tot/area:.4f})")
    cum = 0
    for h in sorted(H, reverse=True):
        cum += H[h]
        print(f"  H>={h:2d}: {cum:7d}   area*dens*2^-h = {tot*2.0**-h:9.1f}")
    print("largest triangles (H, m, J_v):", [(round(a, 1), b, c) for a, b, c in best[:8]])


if __name__ == "__main__":
    main()
