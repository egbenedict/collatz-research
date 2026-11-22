/*
 * b5mitm.c -- meet-in-the-middle search for odd n < 2^N whose first FIVE rises total >= T.
 * Extension of ../b4mitm.c (copied, not modified in place). Same derivation (see ../rises.py):
 *
 *   runs (k0,l0,k1,l1,k2,l2,k3,l3), K1=k0+k1, K2=K1+k2, K3=K2+k3, L1=l0+l1, L3=L1+l2, L=L3+l3.
 *   nu = -1 + u(k0)(1-2^l0) + 2^l0 u(K1)(1-2^l1) + 2^L1 u(K2)(1-2^l2) + 2^L3 u(K3)(1-2^l3),
 *   u(K) = 2^K 3^-K in Z_2.  Left  A = first three pairs (k0..l2),  right R = u(K3)(1-2^l3),
 *   nu = A + 2^L3 R.  z = lowest set bit of nu at position >= N (Terras + least representative),
 *   fifth rise <= z - P (P = K3 + L), total = z - L. total >= T needs z >= T + L >= T + 4, so
 *   bits [N, T+4) of nu are zero: with w = T + 4 - N, window(A) + window(2^L3 R) + carry == 0 mod 2^w,
 *   carry in {0,1}. Each match is then tested exactly mod 2^512.
 *
 * Bounds (rigorous, as in b4mitm):
 *   k0 <= N; l0, k1, l1, k2, l2 limited by prefix feasibility (bits [N,P) of the prefix class are 0);
 *   K3 <= KMAX where KMAX >= B4(N) (first four rises of n < 2^N);
 *   l3 <  N + (delta-1) K3  (next minimum >= 1), looped to N + floor(0.585 K3) + 2.
 *
 * Loop order: L3 outermost, so the right table keyed at shift L3 is built once per L3.
 * Parallel: run NW processes with ID = 0..NW-1; process ID handles L3 = ID mod NW.
 *
 * Usage: b5mitm N T KMAX [ID NW]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LIMBS 8
#define BITS 512
typedef struct { uint64_t w[LIMBS]; } big;

static inline void add(big *r, const big *a, const big *b) {
    unsigned __int128 c = 0;
    for (int i = 0; i < LIMBS; i++) { c += (unsigned __int128)a->w[i] + b->w[i]; r->w[i] = (uint64_t)c; c >>= 64; }
}
static inline void sub(big *r, const big *a, const big *b) {
    unsigned __int128 br = 0;
    for (int i = 0; i < LIMBS; i++) {
        unsigned __int128 d = (unsigned __int128)a->w[i] - b->w[i] - br;
        r->w[i] = (uint64_t)d; br = (d >> 64) ? 1 : 0;
    }
}
static inline void shl(big *r, const big *a, int s) {
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
static inline uint32_t window(const big *a, int pos, int len) {
    int q = pos >> 6, s = pos & 63;
    uint64_t v = a->w[q] >> s;
    if (s && q + 1 < LIMBS) v |= a->w[q + 1] << (64 - s);
    return (uint32_t)(v & ((1ULL << len) - 1));
}
static inline int lowbit_from(const big *a, int N) {
    for (int p = N; p < BITS; ) {
        int q = p >> 6;
        uint64_t v = a->w[q] >> (p & 63);
        if (v) return p + __builtin_ctzll(v);
        p = (q + 1) * 64;
    }
    return -1;
}
static inline int zero_range(const big *a, int lo, int hi) {
    if (hi <= lo) return 1;
    int z = lowbit_from(a, lo);
    return z < 0 || z >= hi;
}
static inline uint32_t shifted_window(const big *R, int s, int N, int w) {
    if (s <= N) return window(R, N - s, w);
    big t; shl(&t, R, s - N); return window(&t, 0, w);
}
/* u(K) << s  -  u(K) << (s + l)   i.e. 2^s u(K) (1 - 2^l) */
static inline void term(big *r, const big *uK, int s, int l) {
    big a, b; shl(&a, uK, s); shl(&b, uK, s + l); sub(r, &a, &b);
}

#define KU 1100
static big U[KU];

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s N T KMAX [ID NW]\n", argv[0]); return 2; }
    int N = atoi(argv[1]), T = atoi(argv[2]), KMAX = atoi(argv[3]);
    int ID = argc > 5 ? atoi(argv[4]) : 0, NW = argc > 5 ? atoi(argv[5]) : 1;
    int w = T + 4 - N;
    if (w < 1 || w > 30 || T + 4 > BITS - 2 || N + 64 > BITS || KMAX >= KU) { fprintf(stderr, "bad parameters (w=%d)\n", w); return 2; }
    int KB = w < 16 ? w : 16;
    uint32_t mask = (1u << w) - 1;

    big inv3; for (int i = 0; i < LIMBS; i++) inv3.w[i] = 0xAAAAAAAAAAAAAAAAULL; inv3.w[0] = 0xAAAAAAAAAAAAAAABULL;
    { big three; memset(&three, 0, sizeof three); three.w[0] = 3; big chk; mul(&chk, &inv3, &three);
      for (int i = 0; i < LIMBS; i++) if (chk.w[i] != (i == 0)) { fprintf(stderr, "inv3 wrong\n"); return 1; } }
    big p; memset(&p, 0, sizeof p); p.w[0] = 1;
    for (int K = 0; K < KU; K++) { shl(&U[K], &p, K); mul(&p, &p, &inv3); }
    big M1; for (int i = 0; i < LIMBS; i++) M1.w[i] = ~0ULL;

    /* right side: all (K3, l3), K3 >= 4 */
    int nR = 0, cap = 0;
    for (int K3 = 4; K3 <= KMAX; K3++) cap += N + (585 * K3) / 1000 + 3;
    big *R = malloc(sizeof(big) * cap);
    int *RK = malloc(sizeof(int) * cap), *RL = malloc(sizeof(int) * cap);
    for (int K3 = 4; K3 <= KMAX; K3++) {
        int lmax = N + (585 * K3) / 1000 + 2;
        for (int l3 = 1; l3 <= lmax; l3++) { term(&R[nR], &U[K3], 0, l3); RK[nR] = K3; RL[nR] = l3; nR++; }
    }
    uint32_t *key = malloc(sizeof(uint32_t) * nR);
    int *order = malloc(sizeof(int) * nR);
    uint32_t *okey = malloc(sizeof(uint32_t) * nR);
    int nb = 1 << KB;
    int *start = malloc(sizeof(int) * (nb + 1));

    long long nleft = 0, nmatch = 0, nhits = 0;
    int best = -1, LMAX3 = 6 * N + 200;
    int bk[8] = {0};
    for (int L3 = 3; L3 <= LMAX3; L3++) {
        if (L3 % NW != ID) continue;
        int built = 0;
        for (int k0 = 1; k0 <= N; k0++) {
            for (int l0 = 1; l0 <= L3 - 2; l0++) {
                big nu0; term(&nu0, &U[k0], 0, l0); add(&nu0, &nu0, &M1);
                int P0 = k0 + l0;
                if (P0 > N && !zero_range(&nu0, N, P0)) break;
                int z0 = lowbit_from(&nu0, N);
                if (z0 < 0) { fprintf(stderr, "overflow z0\n"); return 1; }
                for (int k1 = 1; P0 + k1 <= z0; k1++) {
                    int K1 = k0 + k1;
                    if (K1 + 1 >= KMAX) break;                      /* need K3 = K1 + k2 + k3 <= KMAX */
                    for (int l1 = 1; l0 + l1 <= L3 - 1; l1++) {
                        int L1 = l0 + l1, P1 = P0 + k1 + l1;
                        big A1, t; term(&t, &U[K1], l0, l1); add(&A1, &nu0, &t);
                        if (P1 > N && !zero_range(&A1, N, P1)) break;
                        int z1 = lowbit_from(&A1, N);
                        if (z1 < 0) { fprintf(stderr, "overflow z1\n"); return 1; }
                        int l2 = L3 - L1;
                        for (int k2 = 1; P1 + k2 <= z1; k2++) {
                            int K2 = K1 + k2;
                            if (K2 >= KMAX) break;
                            int P2 = P1 + k2 + l2;
                            big A, t2; term(&t2, &U[K2], L1, l2); add(&A, &A1, &t2);
                            if (P2 > N && !zero_range(&A, N, P2)) continue;
                            nleft++;
                            if (!built) {
                                memset(start, 0, sizeof(int) * (nb + 1));
                                for (int i = 0; i < nR; i++) { key[i] = shifted_window(&R[i], L3, N, w); start[(key[i] & (nb - 1)) + 1]++; }
                                for (int b = 0; b < nb; b++) start[b + 1] += start[b];
                                int *fill = malloc(sizeof(int) * nb); memcpy(fill, start, sizeof(int) * nb);
                                for (int i = 0; i < nR; i++) { int j = fill[key[i] & (nb - 1)]++; order[j] = i; okey[j] = key[i]; }
                                free(fill); built = 1;
                            }
                            uint32_t a_hi = window(&A, N, w);
                            for (int c = 0; c <= 1; c++) {
                                uint32_t need = (uint32_t)(-(int64_t)a_hi - c) & mask;
                                int b = need & (nb - 1);
                                for (int j = start[b]; j < start[b + 1]; j++) {
                                    if (okey[j] != need) continue;
                                    int i = order[j];
                                    if (RK[i] <= K2) continue;
                                    nmatch++;
                                    big S, nu; shl(&S, &R[i], L3); add(&nu, &A, &S);
                                    int z = lowbit_from(&nu, N);
                                    int L = L3 + RL[i], P = RK[i] + L;
                                    if (z < 0) { printf("UNRESOLVED k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d k3=%d l3=%d (z >= 512)\n", k0, l0, k1, l1, k2, l2, RK[i] - K2, RL[i]); nhits++; continue; }
                                    if (z <= P) continue;
                                    int tot = z - L;
                                    if (tot >= T) {
                                        nhits++;
                                        printf("HIT total=%d k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d k3=%d l3=%d k4max=%d\n", tot, k0, l0, k1, l1, k2, l2, RK[i] - K2, RL[i], z - P);
                                        if (tot > best) { best = tot; int v[8] = {k0, l0, k1, l1, k2, l2, RK[i] - K2, RL[i]}; memcpy(bk, v, sizeof v); }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("SUMMARY N=%d T=%d KMAX=%d w=%d ID=%d NW=%d right=%d left=%lld matches=%lld hits=%lld max_total=%d",
           N, T, KMAX, w, ID, NW, nR, nleft, nmatch, nhits, best);
    if (best >= 0) printf(" at (%d,%d,%d,%d,%d,%d,%d,%d)", bk[0], bk[1], bk[2], bk[3], bk[4], bk[5], bk[6], bk[7]);
    printf("\n");
    return 0;
}
