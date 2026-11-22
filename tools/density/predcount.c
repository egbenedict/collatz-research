/* Count #{n <= X : z is in the Collatz orbit of n} by DFS on the inverse tree of C(n) = n/2 | 3n+1.
 * Predecessors of m: 2m, and (m-1)/3 when m = 4 mod 6 (and (m-1)/3 > 1).
 * Nodes above CAP are not expanded (undercount: orbits of n <= X that climb above CAP before z).
 * Usage: ./predcount X CAPFACTOR z1 z2 ...      prints z, count, count/X, z*count/X */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static uint64_t X; static u128 CAP; static uint64_t cnt;
static u128 *stk; static size_t cap_stk = 1<<20;
static void run(uint64_t z) {
    size_t sp = 0; cnt = 0; stk[sp++] = z;
    while (sp) {
        u128 m = stk[--sp];
        if (m <= X) cnt++;
        if (sp + 2 >= cap_stk) { cap_stk *= 2; stk = realloc(stk, cap_stk * sizeof(u128)); }
        if (2*m <= CAP) stk[sp++] = 2*m;
        if (m % 6 == 4) { u128 p = (m - 1) / 3; if (p > 1) stk[sp++] = p; }
    }
}
int main(int argc, char **argv) {
    X = strtoull(argv[1], 0, 10); CAP = (u128)X * strtoull(argv[2], 0, 10);
    stk = malloc(cap_stk * sizeof(u128));
    for (int i = 3; i < argc; i++) {
        uint64_t z = strtoull(argv[i], 0, 10); run(z);
        printf("%llu %llu %.6e %.4f\n", (unsigned long long)z, (unsigned long long)cnt, (double)cnt / X, (double)z * cnt / X);
        fflush(stdout);
    }
}
