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
 *   carry in {0,1}. Each match is then tested exactly mod 2^BITS (BITS = 64 LIMBS; 512 by default).
 *
 * Bounds (rigorous, as in b4mitm):
 *   k0 <= N; l0, k1, l1, k2, l2 limited by prefix feasibility (bits [N,P) of the prefix class are 0);
 *   K3 <= KMAX where KMAX >= B4(N) (first four rises of n < 2^N);
 *   l3 <  N + (delta-1) K3  (next minimum >= 1), looped to N + floor(0.585 K3) + 2.
 *
 * Implementation (v1): right tables keyed at every shift L3 are built up front (L3 < N + 0.585 KMAX + 3,
 *   since y_{i+1} < y_i + (delta-1)k_i - l_i + 1 and y_3 >= 1; a feasible left with larger L3 aborts).
 *   Left loops nest naturally (k0,l0,k1,l1,k2,l2). For fixed (k0..l1,k2), V = A1 + 2^L1 u(K2) has the
 *   final bits [0, P1+k2+l2) of A for every l2, so the l2 loop runs to max(N, z_N(V)) - P1 - k2 exactly.
 *   Work items (k0,l0) are shared by NTHREADS pthreads through an atomic counter.
 *
 * Window: only the low WK = min(w, 30) bits of the w-bit window are used as the join key. A true hit has
 *   all w bits zero, so in particular the low WK bits; this is a weaker filter, never a lost hit (every
 *   match is still tested exactly mod 2^BITS).
 * Sharding: optional SHARD NSHARD keep only work items (k0,l0) with index % NSHARD == SHARD. The union of
 *   all shards is exactly the unsharded item list.
 *
 * Usage: b5mitm N T KMAX [NTHREADS [SHARD NSHARD]]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef LIMBS
#define LIMBS 8           /* 8 limbs = 512 bits; N >= ~300 needs -DLIMBS=16 (shift cap LCAP < BITS) */
#endif
#define BITS (64 * LIMBS)
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


#include <pthread.h>
#include <stdatomic.h>
#ifndef WKCAP
#define WKCAP 30
#endif
#define KU 1100
static big U[KU];
static int N, T, KMAX, w, WK, KB, nb, nR, LCAP;
static uint32_t mask;
static big *R; static int *RK, *RL;
static int **TSTART; static uint32_t **TKEY; static int **TORD;   /* per shift L3 */
static int nitems; static int (*items)[2]; static atomic_int next_item;
static pthread_mutex_t outmx = PTHREAD_MUTEX_INITIALIZER;
static big M1;

typedef struct { long long nleft, nmatch, nhits, bad; int best; int bk[8]; } stats;

static void build_table(int L3) {
    uint32_t *key = malloc(sizeof(uint32_t) * nR);
    int *start = calloc(nb + 1, sizeof(int)), *ord = malloc(sizeof(int) * nR);
    uint32_t *okey = malloc(sizeof(uint32_t) * nR);
    for (int i = 0; i < nR; i++) { key[i] = shifted_window(&R[i], L3, N, WK); start[(key[i] & (nb - 1)) + 1]++; }
    for (int b = 0; b < nb; b++) start[b + 1] += start[b];
    int *fill = malloc(sizeof(int) * nb); memcpy(fill, start, sizeof(int) * nb);
    for (int i = 0; i < nR; i++) { int j = fill[key[i] & (nb - 1)]++; ord[j] = i; okey[j] = key[i]; }
    free(fill); free(key);
    TSTART[L3] = start; TKEY[L3] = okey; TORD[L3] = ord;
}

static void *worker(void *arg) {
    stats *st = arg;
    for (;;) {
        int it = atomic_fetch_add(&next_item, 1);
        if (it >= nitems) break;
        int k0 = items[it][0], l0 = items[it][1];
        big nu0; term(&nu0, &U[k0], 0, l0); add(&nu0, &nu0, &M1);
        int P0 = k0 + l0;
        if (P0 > N && !zero_range(&nu0, N, P0)) continue;
        int z0 = lowbit_from(&nu0, N);
        if (z0 < 0) { st->bad++; continue; }
        for (int k1 = 1; P0 + k1 <= z0; k1++) {
            int K1 = k0 + k1;
            if (K1 + 2 > KMAX) break;                       /* K3 = K1 + k2 + k3 <= KMAX */
            for (int l1 = 1; ; l1++) {
                int L1 = l0 + l1, P1 = P0 + k1 + l1;
                big A1, t; term(&t, &U[K1], l0, l1); add(&A1, &nu0, &t);
                if (P1 > N && !zero_range(&A1, N, P1)) break;
                int z1 = lowbit_from(&A1, N);
                if (z1 < 0) { st->bad++; break; }
                for (int k2 = 1; P1 + k2 <= z1; k2++) {
                    int K2 = K1 + k2;
                    if (K2 + 1 > KMAX) break;
                    big V, s1; shl(&s1, &U[K2], L1); add(&V, &A1, &s1);
                    int zV = lowbit_from(&V, N);
                    if (zV < 0) { st->bad++; continue; }
                    int top = (zV > N ? zV : N);              /* feasible iff P1 + k2 + l2 <= top */
                    for (int l2 = 1; P1 + k2 + l2 <= top; l2++) {
                        int L3 = L1 + l2;
                        if (L3 >= LCAP) { st->bad++; break; }
                        big A, s2; shl(&s2, &U[K2], L3); sub(&A, &V, &s2);
                        st->nleft++;
                        uint32_t a_hi = window(&A, N, WK);
                        int *start = TSTART[L3]; uint32_t *okey = TKEY[L3]; int *ord = TORD[L3];
                        for (int c = 0; c <= 1; c++) {
                            uint32_t need = (uint32_t)(-(int64_t)a_hi - c) & mask;
                            int b = need & (nb - 1);
                            for (int j = start[b]; j < start[b + 1]; j++) {
                                if (okey[j] != need) continue;
                                int i = ord[j];
                                if (RK[i] <= K2) continue;
                                st->nmatch++;
                                big S, nu; shl(&S, &R[i], L3); add(&nu, &A, &S);
                                int z = lowbit_from(&nu, N);
                                int L = L3 + RL[i], P = RK[i] + L;
                                if (z < 0) {
                                    pthread_mutex_lock(&outmx);
                                    printf("UNRESOLVED k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d k3=%d l3=%d (z >= %d)\n", k0, l0, k1, l1, k2, l2, RK[i] - K2, RL[i], BITS);
                                    pthread_mutex_unlock(&outmx);
                                    st->nhits++; continue;
                                }
                                if (z <= P) continue;
                                int tot = z - L;
                                if (tot >= T) {
                                    st->nhits++;
                                    pthread_mutex_lock(&outmx);
                                    printf("HIT total=%d k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d k3=%d l3=%d k4max=%d\n", tot, k0, l0, k1, l1, k2, l2, RK[i] - K2, RL[i], z - P);
                                    fflush(stdout);
                                    pthread_mutex_unlock(&outmx);
                                    if (tot > st->best) { st->best = tot; int v[8] = {k0, l0, k1, l1, k2, l2, RK[i] - K2, RL[i]}; memcpy(st->bk, v, sizeof v); }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s N T KMAX [NTHREADS [SHARD NSHARD]]\n", argv[0]); return 2; }
    N = atoi(argv[1]); T = atoi(argv[2]); KMAX = atoi(argv[3]);
    int NT = argc > 4 ? atoi(argv[4]) : 1;
    int SH = argc > 6 ? atoi(argv[5]) : 0, NS = argc > 6 ? atoi(argv[6]) : 1;
    if (NS < 1 || SH < 0 || SH >= NS || NT < 1 || NT > 64) { fprintf(stderr, "bad shard/threads\n"); return 2; }
    w = T + 4 - N;
    LCAP = N + (585 * KMAX) / 1000 + 4;
    if (w < 1 || T + 4 > BITS - 2 || N + 64 > BITS || KMAX >= KU || LCAP >= BITS) { fprintf(stderr, "bad parameters (w=%d)\n", w); return 2; }
    WK = w < WKCAP ? w : WKCAP;
    KB = WK < 16 ? WK : 16; nb = 1 << KB; mask = (1u << WK) - 1;

    big inv3; for (int i = 0; i < LIMBS; i++) inv3.w[i] = 0xAAAAAAAAAAAAAAAAULL; inv3.w[0] = 0xAAAAAAAAAAAAAAABULL;
    { big three; memset(&three, 0, sizeof three); three.w[0] = 3; big chk; mul(&chk, &inv3, &three);
      for (int i = 0; i < LIMBS; i++) if (chk.w[i] != (i == 0)) { fprintf(stderr, "inv3 wrong\n"); return 1; } }
    big p; memset(&p, 0, sizeof p); p.w[0] = 1;
    for (int K = 0; K < KU; K++) { shl(&U[K], &p, K); mul(&p, &p, &inv3); }
    for (int i = 0; i < LIMBS; i++) M1.w[i] = ~0ULL;

    int cap = 0;
    for (int K3 = 4; K3 <= KMAX; K3++) cap += N + (585 * K3) / 1000 + 3;
    R = malloc(sizeof(big) * cap); RK = malloc(sizeof(int) * cap); RL = malloc(sizeof(int) * cap);
    nR = 0;
    for (int K3 = 4; K3 <= KMAX; K3++) {
        int lmax = N + (585 * K3) / 1000 + 2;
        for (int l3 = 1; l3 <= lmax; l3++) { term(&R[nR], &U[K3], 0, l3); RK[nR] = K3; RL[nR] = l3; nR++; }
    }
    TSTART = calloc(LCAP, sizeof(int *)); TKEY = calloc(LCAP, sizeof(uint32_t *)); TORD = calloc(LCAP, sizeof(int *));
    for (int L3 = 3; L3 < LCAP; L3++) build_table(L3);

    items = malloc(sizeof(int[2]) * (N + 1) * (LCAP + 1)); nitems = 0; long idx = 0;
    for (int k0 = 1; k0 <= N; k0++) for (int l0 = 1; l0 < LCAP; l0++) { if (idx++ % NS == SH) { items[nitems][0] = k0; items[nitems][1] = l0; nitems++; } }
    atomic_store(&next_item, 0);
    pthread_t th[64]; stats st[64]; memset(st, 0, sizeof st);
    for (int i = 0; i < NT; i++) { st[i].best = -1; pthread_create(&th[i], NULL, worker, &st[i]); }
    stats tot; memset(&tot, 0, sizeof tot); tot.best = -1;
    for (int i = 0; i < NT; i++) {
        pthread_join(th[i], NULL);
        tot.nleft += st[i].nleft; tot.nmatch += st[i].nmatch; tot.nhits += st[i].nhits; tot.bad += st[i].bad;
        if (st[i].best > tot.best) { tot.best = st[i].best; memcpy(tot.bk, st[i].bk, sizeof tot.bk); }
    }
    printf("SUMMARY N=%d T=%d KMAX=%d w=%d wk=%d shard=%d/%d threads=%d right=%d left=%lld matches=%lld hits=%lld bad=%lld max_total=%d",
           N, T, KMAX, w, WK, SH, NS, NT, nR, tot.nleft, tot.nmatch, tot.nhits, tot.bad, tot.best);
    if (tot.best >= 0) printf(" at (%d,%d,%d,%d,%d,%d,%d,%d)", tot.bk[0], tot.bk[1], tot.bk[2], tot.bk[3], tot.bk[4], tot.bk[5], tot.bk[6], tot.bk[7]);
    printf("\n");
    return tot.bad ? 3 : 0;
}
