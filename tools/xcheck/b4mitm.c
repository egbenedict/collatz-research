/*
 * b4mitm.c -- meet-in-the-middle search for odd n < 2^N whose first four rises total >= T.
 *
 * Independent implementation for the cross-check (own derivation; see rises.py):
 *   For runs (k0,l0,k1,l1,k2,l2), the 2-adic integer with parity vector
 *   1^k0 0^l0 1^k1 0^l1 1^k2 0^l2 1^inf is
 *     nu = A + 2^(l0+l1) * R,
 *     A  = -1 + u(k0)(1-2^l0) + 2^l0 u(K1)(1-2^l1),   K1 = k0+k1,
 *     R  = u(K2)(1-2^l2),                               K2 = K1+k2,
 *   with u(K) = 2^K 3^-K in Z_2. The fourth rise can be at most z - P, where P = K2+L,
 *   L = l0+l1+l2, and z = lowest set bit of nu at position >= N (Terras + least representative).
 *   The total is z - L, so total >= T needs z >= T + L >= T + 3, i.e. bits [N, T+3) of nu are zero.
 *   With w = T + 3 - N, matching the w-bit window of A against that of 2^L2 R (carry 0 or 1)
 *   is a necessary condition, so each A is joined against a table of R keyed by its window.
 *
 * Bounds (rigorous):
 *   k0 <= N (n+1 = a 2^k0 <= 2^N);
 *   l0, k1, l1 are limited by prefix feasibility (bits [N, P) of the 2-adic integer are zero);
 *   K2 <= KMAX, where KMAX >= B3(N) (the first three rises of n < 2^N);
 *   l2 <  N + (delta-1) K2, because the next minimum is >= 1, so 2^l2 < a2 3^k2 <= 2^(x2 + (delta-1)k2),
 *         and x2 <= N + (delta-1)(k0+k1).
 *   L2 = l0 + l1 is looped up to LMAX2 = 4N + 100. Feasibility, not this cap, limits it.
 *
 * Arithmetic is exact mod 2^512 (8 x 64-bit limbs). Every bit used is below 512, and a hit is
 * reported only when z < 512 (checked).
 *
 * Usage: b4mitm N T KMAX        (prints every hit with total >= T, and the maximum found)
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

#define KU 1100
static big U[KU];          /* U[K] = 2^K * 3^-K mod 2^512 (0 for K >= 512) */

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s N T KMAX\n", argv[0]); return 2; }
    int N = atoi(argv[1]), T = atoi(argv[2]), KMAX = atoi(argv[3]);
    int w = T + 3 - N;
    if (w < 1 || w > 30 || T + 3 > BITS - 2 || N + 64 > BITS) { fprintf(stderr, "bad parameters (w=%d)\n", w); return 2; }
    int KB = w < 20 ? w : 20;                       /* bucket on low KB bits of the key */
    uint32_t mask = (w == 32) ? 0xFFFFFFFFu : ((1u << w) - 1);

    /* inverse of 3 mod 2^512: 0xAAAA...AAAB */
    big inv3; for (int i = 0; i < LIMBS; i++) inv3.w[i] = 0xAAAAAAAAAAAAAAAAULL; inv3.w[0] = 0xAAAAAAAAAAAAAAABULL;
    { big three; memset(&three, 0, sizeof three); three.w[0] = 3; big chk; mul(&chk, &inv3, &three);
      for (int i = 0; i < LIMBS; i++) if (chk.w[i] != (i == 0)) { fprintf(stderr, "inv3 wrong\n"); return 1; } }
    big p; memset(&p, 0, sizeof p); p.w[0] = 1;
    for (int K = 0; K < KU; K++) { shl(&U[K], &p, K); mul(&p, &p, &inv3); }

    big M1; for (int i = 0; i < LIMBS; i++) M1.w[i] = ~0ULL;            /* -1 */

    /* right side: all (K2, l2) */
    int nR = 0, cap = 0;
    for (int K2 = 3; K2 <= KMAX; K2++) cap += N + (585 * K2) / 1000 + 3;
    big *R = malloc(sizeof(big) * cap);
    int *RK = malloc(sizeof(int) * cap), *RL = malloc(sizeof(int) * cap);
    for (int K2 = 3; K2 <= KMAX; K2++) {
        int lmax = N + (585 * K2) / 1000 + 2;
        for (int l2 = 1; l2 <= lmax; l2++) {
            big t; shl(&t, &U[K2], l2); sub(&R[nR], &U[K2], &t);
            RK[nR] = K2; RL[nR] = l2; nR++;
        }
    }
    uint32_t *key = malloc(sizeof(uint32_t) * nR);
    int *order = malloc(sizeof(int) * nR);
    int nb = 1 << KB;
    int *start = malloc(sizeof(int) * (nb + 1));

    long long nleft = 0, nmatch = 0, nhits = 0;
    int best = -1, LMAX2 = 4 * N + 100;
    int bk[6] = {0};
    for (int L2 = 2; L2 <= LMAX2; L2++) {
        int built = 0;
        for (int k0 = 1; k0 <= N; k0++) {
            for (int l0 = 1; l0 <= L2 - 1; l0++) {
                big t, nu0; shl(&t, &U[k0], l0); sub(&nu0, &U[k0], &t); add(&nu0, &nu0, &M1);
                int P0 = k0 + l0;
                if (P0 > N && !zero_range(&nu0, N, P0)) break;           /* bare prefix 1^k0 0^l0 infeasible */
                int l1 = L2 - l0;
                int z0 = lowbit_from(&nu0, N);
                if (z0 < 0) { fprintf(stderr, "overflow z0\n"); return 1; }
                for (int k1 = 1; P0 + k1 <= z0; k1++) {
                    int K1 = k0 + k1;
                    if (K1 >= KMAX) break;                                /* K2 = K1 + k2 must stay <= KMAX */
                    big t1, t2, A;
                    shl(&t1, &U[K1], l1); sub(&t1, &U[K1], &t1); shl(&t2, &t1, l0); add(&A, &nu0, &t2);
                    int PL = P0 + k1 + l1;
                    if (PL > N && !zero_range(&A, N, PL)) continue;      /* prefix ...1^k1 0^l1 infeasible */
                    nleft++;
                    if (!built) {                                         /* counting sort of R by key at shift L2 */
                        memset(start, 0, sizeof(int) * (nb + 1));
                        for (int i = 0; i < nR; i++) { key[i] = shifted_window(&R[i], L2, N, w); start[(key[i] & (nb - 1)) + 1]++; }
                        for (int b = 0; b < nb; b++) start[b + 1] += start[b];
                        int *fill = malloc(sizeof(int) * nb); memcpy(fill, start, sizeof(int) * nb);
                        for (int i = 0; i < nR; i++) order[fill[key[i] & (nb - 1)]++] = i;
                        free(fill); built = 1;
                    }
                    uint32_t a_hi = window(&A, N, w);
                    for (int c = 0; c <= 1; c++) {
                        uint32_t need = (uint32_t)(-(int64_t)a_hi - c) & mask;
                        int b = need & (nb - 1);
                        for (int j = start[b]; j < start[b + 1]; j++) {
                            int i = order[j];
                            if (key[i] != need || RK[i] <= K1) continue;
                            nmatch++;
                            big S, nu; shl(&S, &R[i], L2); add(&nu, &A, &S);
                            int z = lowbit_from(&nu, N);
                            int L = L2 + RL[i], P = RK[i] + L;
                            if (z < 0) { printf("UNRESOLVED k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d (z >= 512)\n", k0, l0, k1, l1, RK[i] - K1, RL[i]); nhits++; continue; }
                            if (z <= P) continue;                         /* fourth rise must be >= 1 */
                            int tot = z - L;
                            if (tot >= T) {
                                nhits++;
                                printf("HIT total=%d k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d k3max=%d\n", tot, k0, l0, k1, l1, RK[i] - K1, RL[i], z - P);
                                if (tot > best) { best = tot; bk[0]=k0; bk[1]=l0; bk[2]=k1; bk[3]=l1; bk[4]=RK[i]-K1; bk[5]=RL[i]; }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("SUMMARY N=%d T=%d KMAX=%d w=%d right=%d left=%lld matches=%lld hits=%lld max_total=%d",
           N, T, KMAX, w, nR, nleft, nmatch, nhits, best);
    if (best >= 0) printf(" at (%d,%d,%d,%d,%d,%d)", bk[0], bk[1], bk[2], bk[3], bk[4], bk[5]);
    printf("\n");
    return 0;
}
