from sens5 import *
grid = sorted(B4)
for m in (134, 135, 136):
    for c in (1, 2, 3, 4, 5, 6, 8, 10, 12):
        hit = None
        for Nmax in grid:
            if hard_trace(checker(hyp(c, Nmax)), m)[0] < 1:
                hit = Nmax; break
        print(f"m={m} c={c}: smallest Nmax closing hard shape = {hit}", flush=True)
