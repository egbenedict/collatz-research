/*
 * b3mitm.c -- meet-in-the-middle search for odd n < 2^N whose first three rises total >= T.
 * Derived from b4mitm.c (same arithmetic helpers), one level shorter:
 *   parity vector 1^k0 0^l0 1^k1 0^l1 1^inf has 2-adic integer
 *     nu = A + 2^l0 * R,   A = -1 + u(k0)(1-2^l0),   R = u(K1)(1-2^l1),   K1 = k0+k1,
 *   u(K) = 2^K 3^-K in Z_2. The third rise is at most z - P, P = K1 + l0 + l1, z = lowest set bit of
 *   nu at position >= N. Total = z - L (L = l0+l1); total >= T needs bits [N, T+2) of nu zero, so
 *   w = T + 2 - N; the w-bit window of A is joined against that of 2^l0 R (carry 0 or 1).
 * Bounds: k0 <= N; l0 limited by prefix feasibility; K1 <= KMAX where KMAX >= B2(N);
 *   l1 < N + (delta-1) K1 (same argument as b4mitm for l2).
 * Arithmetic exact mod 2^BITS (BITS = 64*LIMBS, -DLIMBS=...); a hit needs z < BITS (else UNRESOLVED).
 * Usage: b3mitm N T KMAX
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef LIMBS
#define LIMBS 8           /* -DLIMBS=16 for 1024-bit arithmetic (needed when T + 3 > 510) */
#endif
#define BITS (64 * LIMBS)
typedef struct { uint64_t w[LIMBS]; } big;

static void add(big *r, const big *a, const big *b) {
    unsigned __int128 c = 0;
    for (int i = 0; i < LIMBS; i++) { c += (unsigned __int128)a->w[i] + b->w[i]; r->w[i] = (uint64_t)c; c >>= 64; }
}
static void sub(big *r, const big *a, const big *b) {
    unsigned __int128 br = 0;
    for (int i = 0; i < LIMBS; i++) {
        unsigned __int128 d = (unsigned __int128)a->w[i] - b->w[i] - br;
        r->w[i] = (uint64_t)d; br = (d >> 64) ? 1 : 0;
    }
}
static void shl(big *r, const big *a, int s) {
    big t; memset(&t, 0, sizeof t);
    if (s < BITS) {
        int q = s / 64, b = s % 64;
        for (int i = LIMBS - 1; i >= q; i--) {
            uint64_t v = a->w[i - q] << b;
            if (b && i - q - 1 >= 0) v |= a->w[i - q - 1] >> (64 - b);
            t.w[i] = v;
        }
    }
    *r = t;
}
static void mul(big *r, const big *a, const big *b) {
    uint64_t t[LIMBS] = {0};
    for (int i = 0; i < LIMBS; i++) {
        unsigned __int128 c = 0;
        for (int j = 0; i + j < LIMBS; j++) {
            c += (unsigned __int128)a->w[i] * b->w[j] + t[i + j];
            t[i + j] = (uint64_t)c; c >>= 64;
        }
    }
    memcpy(r->w, t, sizeof t);
}
/* bits [pos, pos+len) as an integer, len <= 32, pos+len <= BITS */
static inline uint32_t window(const big *a, int pos, int len) {
    int q = pos >> 6, s = pos & 63;
    uint64_t v = a->w[q] >> s;
    if (s && q + 1 < LIMBS) v |= a->w[q + 1] << (64 - s);
    return (uint32_t)(v & ((1ULL << len) - 1));
}
/* lowest set bit at position >= N, or -1 if none below BITS */
static int lowbit_from(const big *a, int N) {
    for (int p = N; p < BITS; ) {
        int q = p >> 6;
        uint64_t v = a->w[q] >> (p & 63);
        if (v) return p + __builtin_ctzll(v);
        p = (q + 1) * 64;
    }
    return -1;
}
/* are bits [lo, hi) all zero? (hi <= BITS) */
static int zero_range(const big *a, int lo, int hi) {
    if (hi <= lo) return 1;
    int z = lowbit_from(a, lo);
    return z < 0 || z >= hi;
}
/* window of 2^s * R at [N, N+w): bits [N-s, N-s+w) of R if s <= N, else (R << (s-N)) low bits */
static inline uint32_t shifted_window(const big *R, int s, int N, int w) {
    if (s <= N) return window(R, N - s, w);
    big t; shl(&t, R, s - N); return window(&t, 0, w);
}

#define KU 4200
static big U[KU];

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s N T KMAX\n", argv[0]); return 2; }
    int N = atoi(argv[1]), T = atoi(argv[2]), KMAX = atoi(argv[3]);
    int w = T + 2 - N;
    if (w < 1 || w > 30 || T + 2 > BITS - 2 || N + 64 > BITS || KMAX >= KU) { fprintf(stderr, "bad parameters (w=%d)\n", w); return 2; }
    int KB = w < 20 ? w : 20;
    uint32_t mask = (1u << w) - 1;
    big inv3; for (int i = 0; i < LIMBS; i++) inv3.w[i] = 0xAAAAAAAAAAAAAAAAULL; inv3.w[0] = 0xAAAAAAAAAAAAAAABULL;
    { big three; memset(&three, 0, sizeof three); three.w[0] = 3; big chk; mul(&chk, &inv3, &three);
      for (int i = 0; i < LIMBS; i++) if (chk.w[i] != (i == 0)) { fprintf(stderr, "inv3 wrong\n"); return 1; } }
    big p; memset(&p, 0, sizeof p); p.w[0] = 1;
    for (int K = 0; K < KU; K++) { shl(&U[K], &p, K); mul(&p, &p, &inv3); }
    big M1; for (int i = 0; i < LIMBS; i++) M1.w[i] = ~0ULL;

    long cap = 0;
    for (int K1 = 2; K1 <= KMAX; K1++) cap += N + (585 * K1) / 1000 + 3;
    big *R = malloc(sizeof(big) * cap); int *RK = malloc(sizeof(int) * cap), *RL = malloc(sizeof(int) * cap);
    if (!R || !RK || !RL) { fprintf(stderr, "out of memory (cap=%ld)\n", cap); return 1; }
    int nR = 0;
    for (int K1 = 2; K1 <= KMAX; K1++) {
        int lmax = N + (585 * K1) / 1000 + 2;
        for (int l1 = 1; l1 <= lmax; l1++) { big t; shl(&t, &U[K1], l1); sub(&R[nR], &U[K1], &t); RK[nR] = K1; RL[nR] = l1; nR++; }
    }
    uint32_t *key = malloc(sizeof(uint32_t) * nR); int *order = malloc(sizeof(int) * nR);
    int nb = 1 << KB; int *start = malloc(sizeof(int) * (nb + 1)); int *fill = malloc(sizeof(int) * nb);

    long long nleft = 0, nmatch = 0, nhits = 0, nunres = 0;
    int best = -1, bk[4] = {0}, LMAX0 = 4 * N + 100;
    for (int l0 = 1; l0 <= LMAX0; l0++) {
        int built = 0;
        for (int k0 = 1; k0 <= N; k0++) {
            big t, A; shl(&t, &U[k0], l0); sub(&A, &U[k0], &t); add(&A, &A, &M1);
            int P0 = k0 + l0;
            if (P0 > N && !zero_range(&A, N, P0)) continue;   /* prefix 1^k0 0^l0 infeasible (not monotone in k0: no break) */
            if (k0 + 1 > KMAX) break;
            nleft++;
            if (!built) {
                memset(start, 0, sizeof(int) * (nb + 1));
                for (int i = 0; i < nR; i++) { key[i] = shifted_window(&R[i], l0, N, w); start[(key[i] & (nb - 1)) + 1]++; }
                for (int b = 0; b < nb; b++) start[b + 1] += start[b];
                memcpy(fill, start, sizeof(int) * nb);
                for (int i = 0; i < nR; i++) order[fill[key[i] & (nb - 1)]++] = i;
                built = 1;
            }
            uint32_t a_hi = window(&A, N, w);
            for (int c = 0; c <= 1; c++) {
                uint32_t need = (uint32_t)(-(int64_t)a_hi - c) & mask;
                int b = need & (nb - 1);
                for (int j = start[b]; j < start[b + 1]; j++) {
                    int i = order[j];
                    if (key[i] != need || RK[i] <= k0) continue;
                    nmatch++;
                    big S, nu; shl(&S, &R[i], l0); add(&nu, &A, &S);
                    int z = lowbit_from(&nu, N);
                    int L = l0 + RL[i], P = RK[i] + L;
                    if (z < 0) { printf("UNRESOLVED k0=%d l0=%d k1=%d l1=%d (z >= %d)\n", k0, l0, RK[i] - k0, RL[i], BITS); nunres++; continue; }
                    if (z <= P) continue;
                    int tot = z - L;
                    if (tot >= T) {
                        nhits++;
                        printf("HIT total=%d k0=%d l0=%d k1=%d l1=%d k2max=%d\n", tot, k0, l0, RK[i] - k0, RL[i], z - P);
                        if (tot > best) { best = tot; bk[0] = k0; bk[1] = l0; bk[2] = RK[i] - k0; bk[3] = RL[i]; }
                    }
                }
            }
        }
    }
    printf("SUMMARY N=%d T=%d KMAX=%d w=%d bits=%d right=%d left=%lld matches=%lld hits=%lld unresolved=%lld max_total=%d",
           N, T, KMAX, w, BITS, nR, nleft, nmatch, nhits, nunres, best);
    if (best >= 0) printf(" at (%d,%d,%d,%d)", bk[0], bk[1], bk[2], bk[3]);
    printf("\n");
    return nunres ? 3 : 0;
}
