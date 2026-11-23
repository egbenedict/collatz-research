/*
 * b6mitm.c -- six-rise version of b5mitm.c (same arithmetic, sharding, window and item logic).
 * Left: (k0,l0,k1,l1,k2,l2,k3,l3) with exact prefix feasibility at every level (the l3 loop uses the
 * V-trick of b5mitm: V = A2 + 2^L2 u(K3) fixes the deciding bits). Right: R(K4,l4) = u(K4)(1-2^l4) at
 * shift L4 = l0+..+l3. Total = z - L (L = L4 + l4); total >= T needs bits [N, T+5) of nu zero, so the
 * join window is w = T + 5 - N. Caps: K4 <= KMAX where KMAX >= B5(N); l4 < N + (delta-1)K4;
 * L4 < LCAP = N + 0.585 KMAX + 5 (hard check: exceeding it counts as bad and fails the run).
 * Usage: b6mitm N T KMAX [NTHREADS [SHARD NSHARD]]
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
static int **TSTART; static uint32_t **TKEY; static int **TORD;   /* per shift L4 */
static int nitems; static int (*items)[2]; static atomic_int next_item;
static pthread_mutex_t outmx = PTHREAD_MUTEX_INITIALIZER;
static big M1;

typedef struct { long long nleft, nmatch, nhits, bad; int best; int bk[10]; } stats;

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
            if (K1 + 3 > KMAX) break;
            for (int l1 = 1; ; l1++) {
                int L1 = l0 + l1, P1 = P0 + k1 + l1;
                big A1, t; term(&t, &U[K1], l0, l1); add(&A1, &nu0, &t);
                if (P1 > N && !zero_range(&A1, N, P1)) break;
                int z1 = lowbit_from(&A1, N);
                if (z1 < 0) { st->bad++; break; }
                for (int k2 = 1; P1 + k2 <= z1; k2++) {
                    int K2 = K1 + k2;
                    if (K2 + 2 > KMAX) break;
                    for (int l2 = 1; ; l2++) {
                        int L2 = L1 + l2, P2 = P1 + k2 + l2;
                        big A2, t2; term(&t2, &U[K2], L1, l2); add(&A2, &A1, &t2);
                        if (P2 > N && !zero_range(&A2, N, P2)) break;
                        int z2 = lowbit_from(&A2, N);
                        if (z2 < 0) { st->bad++; break; }
                        for (int k3 = 1; P2 + k3 <= z2; k3++) {
                            int K3 = K2 + k3;
                            if (K3 + 1 > KMAX) break;
                            big V, s1; shl(&s1, &U[K3], L2); add(&V, &A2, &s1);
                            int zV = lowbit_from(&V, N);
                            if (zV < 0) { st->bad++; continue; }
                            int top = (zV > N ? zV : N);
                            for (int l3 = 1; P2 + k3 + l3 <= top; l3++) {
                                int L4 = L2 + l3;
                                if (L4 >= LCAP) { st->bad++; break; }
                                big A, s2; shl(&s2, &U[K3], L4); sub(&A, &V, &s2);
                                st->nleft++;
                                uint32_t a_hi = window(&A, N, WK);
                                int *start = TSTART[L4]; uint32_t *okey = TKEY[L4]; int *ord = TORD[L4];
                                for (int c = 0; c <= 1; c++) {
                                    uint32_t need = (uint32_t)(-(int64_t)a_hi - c) & mask;
                                    int b = need & (nb - 1);
                                    for (int j = start[b]; j < start[b + 1]; j++) {
                                        if (okey[j] != need) continue;
                                        int i = ord[j];
                                        if (RK[i] <= K3) continue;
                                        st->nmatch++;
                                        big S, nu; shl(&S, &R[i], L4); add(&nu, &A, &S);
                                        int z = lowbit_from(&nu, N);
                                        int L = L4 + RL[i], P = RK[i] + L;
                                        if (z < 0) {
                                            pthread_mutex_lock(&outmx);
                                            printf("UNRESOLVED k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d k3=%d l3=%d k4=%d l4=%d (z >= %d)\n", k0, l0, k1, l1, k2, l2, k3, l3, RK[i] - K3, RL[i], BITS);
                                            pthread_mutex_unlock(&outmx);
                                            st->nhits++; st->bad++; continue;
                                        }
                                        if (z <= P) continue;
                                        int tot = z - L;
                                        if (tot >= T) {
                                            st->nhits++;
                                            pthread_mutex_lock(&outmx);
                                            printf("HIT total=%d k0=%d l0=%d k1=%d l1=%d k2=%d l2=%d k3=%d l3=%d k4=%d l4=%d k5max=%d\n", tot, k0, l0, k1, l1, k2, l2, k3, l3, RK[i] - K3, RL[i], z - P);
                                            fflush(stdout);
                                            pthread_mutex_unlock(&outmx);
                                            if (tot > st->best) { st->best = tot; int v[10] = {k0, l0, k1, l1, k2, l2, k3, l3, RK[i] - K3, RL[i]}; memcpy(st->bk, v, sizeof v); }
                                        }
                                    }
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
    if (argc < 4) { fprintf(stderr, "usage: %s N T KMAX [NTHREADS [SHARD NSHARD]] (six rises)\n", argv[0]); return 2; }
    N = atoi(argv[1]); T = atoi(argv[2]); KMAX = atoi(argv[3]);
    int NT = argc > 4 ? atoi(argv[4]) : 1;
    int SH = argc > 6 ? atoi(argv[5]) : 0, NS = argc > 6 ? atoi(argv[6]) : 1;
    if (NS < 1 || SH < 0 || SH >= NS || NT < 1 || NT > 64) { fprintf(stderr, "bad shard/threads\n"); return 2; }
    w = T + 5 - N;
    LCAP = N + (585 * KMAX) / 1000 + 5;
    if (w < 1 || T + 5 > BITS - 2 || N + 64 > BITS || KMAX >= KU || LCAP >= BITS) { fprintf(stderr, "bad parameters (w=%d)\n", w); return 2; }
    WK = w < WKCAP ? w : WKCAP;
    KB = WK < 16 ? WK : 16; nb = 1 << KB; mask = (1u << WK) - 1;

    big inv3; for (int i = 0; i < LIMBS; i++) inv3.w[i] = 0xAAAAAAAAAAAAAAAAULL; inv3.w[0] = 0xAAAAAAAAAAAAAAABULL;
    { big three; memset(&three, 0, sizeof three); three.w[0] = 3; big chk; mul(&chk, &inv3, &three);
      for (int i = 0; i < LIMBS; i++) if (chk.w[i] != (i == 0)) { fprintf(stderr, "inv3 wrong\n"); return 1; } }
    big p; memset(&p, 0, sizeof p); p.w[0] = 1;
    for (int K = 0; K < KU; K++) { shl(&U[K], &p, K); mul(&p, &p, &inv3); }
    for (int i = 0; i < LIMBS; i++) M1.w[i] = ~0ULL;

    int cap = 0;
    for (int K4 = 5; K4 <= KMAX; K4++) cap += N + (585 * K4) / 1000 + 3;
    R = malloc(sizeof(big) * cap); RK = malloc(sizeof(int) * cap); RL = malloc(sizeof(int) * cap);
    nR = 0;
    for (int K4 = 5; K4 <= KMAX; K4++) {
        int lmax = N + (585 * K4) / 1000 + 2;
        for (int l4 = 1; l4 <= lmax; l4++) { term(&R[nR], &U[K4], 0, l4); RK[nR] = K4; RL[nR] = l4; nR++; }
    }
    TSTART = calloc(LCAP, sizeof(int *)); TKEY = calloc(LCAP, sizeof(uint32_t *)); TORD = calloc(LCAP, sizeof(int *));
    for (int L4 = 4; L4 < LCAP; L4++) build_table(L4);

    items = malloc(sizeof(int[2]) * (N + 1) * (LCAP + 1)); nitems = 0; long idx = 0;
    /* Optional restriction (for cheap lower-bound / falsifier searches, NOT for exact values):
       env K0MIN, K0MAX, L0MAX keep only items with K0MIN <= k0 <= K0MAX and l0 <= L0MAX. */
    int k0min = getenv("K0MIN") ? atoi(getenv("K0MIN")) : 1, k0max = getenv("K0MAX") ? atoi(getenv("K0MAX")) : N;
    int l0max = getenv("L0MAX") ? atoi(getenv("L0MAX")) : LCAP;
    int restricted = k0min > 1 || k0max < N || l0max < LCAP;
    for (int k0 = 1; k0 <= N; k0++) for (int l0 = 1; l0 < LCAP; l0++) {
        if (idx++ % NS != SH) continue;
        if (k0 < k0min || k0 > k0max || l0 > l0max) continue;
        items[nitems][0] = k0; items[nitems][1] = l0; nitems++; }
    if (restricted) fprintf(stderr, "RESTRICTED items k0 in [%d,%d], l0 <= %d: %d items (partial search, lower bound only)\n", k0min, k0max, l0max, nitems);
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
    if (tot.best >= 0) printf(" at (%d,%d,%d,%d,%d,%d,%d,%d,%d,%d)", tot.bk[0], tot.bk[1], tot.bk[2], tot.bk[3], tot.bk[4], tot.bk[5], tot.bk[6], tot.bk[7], tot.bk[8], tot.bk[9]);
    printf("\n");
    return tot.bad ? 3 : 0;
}
