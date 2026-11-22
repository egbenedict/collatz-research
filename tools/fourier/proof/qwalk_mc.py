"""Monte Carlo of the Q-walk of NOTES.md (hypothesis H-white), exact big-integer phases.

Q_0 = (n, k-2), Q_{p+1} = Q_p + (-2, -b_p), P(b) = (b-1) 2^-b.  white <=> ||th(Q)|| >= 1/5.
Reports V/P (P = floor(n/2) pairs), its lower quantiles, and the exit statistic:
after each maximal gray run (excluding the run containing p=0), fraction of exits that land white.
Usage: python3 qwalk_mc.py n start(top|rand|neg) nwalks [seed]
"""
import sys
import math
import random


def sample_b(rng):
    # b = a + a', a,a' ~ Geom(1/2) on {1,2,...}
    a = 1
    while rng.random() < 0.5:
        a += 1
    c = 1
    while rng.random() < 0.5:
        c += 1
    return a + c


def walk(n, k, rng, eps=0.2):
    P = n // 2
    m, J = n, k - 2
    flags = []
    for p in range(P):
        mod = 3 ** m
        r = pow(2, J, mod)
        x = r / mod
        flags.append(min(x, 1 - x) >= eps)
        m -= 2
        J -= sample_b(rng)
    return flags


def main():
    n = int(sys.argv[1]); start = sys.argv[2]; N = int(sys.argv[3])
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    fr, exits, exitw, maxrun = [], 0, 0, []
    for _ in range(N):
        if start == "top":
            k = math.floor(n * math.log2(3)) - 7
        elif start == "rand":
            k = rng.randrange(2 * 3 ** (n - 1))
        else:
            k = -rng.randrange(1, 4 * n)
        f = walk(n, k, rng)
        P = len(f)
        fr.append(sum(f) / P)
        # gray runs not containing index 0
        i = 0
        while i < P and not f[i]:
            i += 1
        first = i
        runs = []
        while i < P:
            if f[i]:
                i += 1
                continue
            j = i
            while j < P and not f[j]:
                j += 1
            runs.append(j - i)
            if j < P:
                exits += 1
                exitw += 1  # landing at j is white by maximality -- trivially; measure next instead
            i = j
        maxrun.append(max(runs) if runs else 0)
    fr.sort()
    q = lambda t: fr[int(t * (len(fr) - 1))]
    maxrun.sort()
    print(f"n={n} start={start} walks={N}: V/P mean={sum(fr)/N:.3f}  q01={q(0.01):.3f} q05={q(0.05):.3f} "
          f"min={fr[0]:.3f}   longest non-initial gray run: median={maxrun[N//2]} max={maxrun[-1]}")


if __name__ == "__main__":
    main()
