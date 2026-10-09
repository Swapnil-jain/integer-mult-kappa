/* FAST PASS mod 2^21-9 in double precision (zero candidates are re-checked mod 2^61-1 by gate.c).
   All-edges C1 checker for the opposite-bank compile of witness 2.
   For ONE fixed basis X (integer entries, seed-determined), for every charged residual idempotent p of rank r:
   the leading principal minors 1..r of X^-1 p X are nonzero.  Arithmetic mod P = 2^61-1.  Soundness: X is an
   integer matrix invertible mod P, every frame Gram matrix is invertible mod P (checked) and 9 is a unit, so
   X^-1 p X has entries in Z_(P); a minor nonzero mod P is nonzero over Q.
   Residual forms (m = h*h, index (x,y) -> x*h+y, P_b = t_b g_b^T, g_b = G t_b/(t_b^T G t_b), Q = I - P):
     side/data step pi = P_B - P_A :  s1 pi (x) P_b,  s2 P_b (x) pi                rank r = dim B - dim A
     aux edge sigma               :  s1 I - Q_s (x) P_b, s2 I - P_b (x) Q_s       rank m - (h - dim sigma)
     data entrance (a,b)          :  Q_a (x) Q_b = I - (P_a (x) I + Q_a (x) P_b)   rank m - 2h + 1
     copy correction (a,b)        :  P_a (x) P_b                                   rank 1
   Low rank p = U V^T: leading minors of the r x r matrix U'[:r] V'[:r]^T (U' = X^-1 U, V' = X^T V).
   High rank p = I - U V^T (rank s complement): minor k = det(I_s - sum_{i<k} V'_i U'_i^T), tracked by
   Sherman-Morrison; every determinant-lemma factor 1 - u^T A^-1 v must be nonzero.
   Usage: gate data.bin MODE [args]   MODE in: small aux copy entr A0 A1 | dense | neg */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
typedef uint64_t u64; typedef unsigned __int128 u128;
static const u64 P = 2097143ULL;   /* 2^21 - 9: fast double-precision pass; zeros are re-checked mod 2^61-1 */
static inline u64 fold(u128 x) { return (x >> 64) ? (u64)(x % P) : (u64)((u64)x % P); }
static inline u64 mul(u64 a, u64 b) { return fold((u128)a * b); }
static inline u64 add(u64 a, u64 b) { u64 r = a + b; return r >= P ? r - P : r; }
static inline u64 sub(u64 a, u64 b) { return a >= b ? a - b : a + P - b; }
static u64 pw(u64 a, u64 e) { u64 r = 1; while (e) { if (e & 1) r = mul(r, a); a = mul(a, a); e >>= 1; } return r; }
static inline u64 inv(u64 a) { return pw(a, P - 2); }
static u64 rs = 0x9E3779B97F4A7C15ULL;
static u64 rnd(void) { u64 z = (rs += 0x9E3779B97F4A7C15ULL); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
static double T0; static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
#define LOG(...) do { fprintf(stderr, "[%7.1fs] ", now() - T0); fprintf(stderr, __VA_ARGS__); fprintf(stderr, "\n"); fflush(stderr); } while (0)

static int h, v, m, nfr, nst, nsig;
static int (*trip)[3]; static int *fdim; static u64 **frow; static int (*st)[4]; static int *sig;
static u64 *X, *Xi;                      /* m x m */
static u64 *Z1, *Z2, *W1, *W2;           /* [i][b][a]  m x v x h */
#define IDX(i, b, a) (((size_t)(i) * v + (b)) * h + (a))
static u64 inv9, inv3, inv2;
static u64 *tv, *gv;                     /* v x h */

static void load(const char *fn) {
    FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(2); }
    int64_t hd[5]; if (fread(hd, 8, 5, f) != 5) exit(2);
    h = hd[0]; v = hd[1]; nfr = hd[2]; nst = hd[3]; nsig = hd[4]; m = h * h;
    trip = malloc(sizeof(int[3]) * v);
    for (int b = 0; b < v; b++) { int64_t t[3]; if (fread(t, 8, 3, f) != 3) exit(2); for (int j = 0; j < 3; j++) trip[b][j] = t[j]; }
    fdim = malloc(sizeof(int) * nfr); frow = malloc(sizeof(u64 *) * nfr);
    for (int i = 0; i < nfr; i++) { int64_t d; if (fread(&d, 8, 1, f) != 1) exit(2); fdim[i] = d; frow[i] = malloc(8 * (d * h + 1)); { uint64_t *tmp = malloc(8 * (d * h + 1)); if (d && fread(tmp, 8, d * h, f) != (size_t)(d * h)) exit(2); const uint64_t P61 = (1ULL << 61) - 1; for (int j = 0; j < d * h; j++) { int64_t x = tmp[j] > P61 / 2 ? -(int64_t)(P61 - tmp[j]) : (int64_t)tmp[j]; int64_t y = x % (int64_t)P; frow[i][j] = (u64)(y < 0 ? y + P : y); } free(tmp); } }
    st = malloc(sizeof(int[4]) * nst);
    for (int i = 0; i < nst; i++) { int64_t t[4]; if (fread(t, 8, 4, f) != 4) exit(2); for (int j = 0; j < 4; j++) st[i][j] = t[j]; }
    sig = malloc(sizeof(int) * nsig);
    for (int i = 0; i < nsig; i++) { int64_t t; if (fread(&t, 8, 1, f) != 1) exit(2); sig[i] = t; }
    fclose(f);
    inv9 = inv(9); inv3 = inv(3); inv2 = inv(2);
    tv = calloc((size_t)v * h, 8); gv = calloc((size_t)v * h, 8);
    for (int b = 0; b < v; b++) {
        for (int j = 0; j < 3; j++) tv[b * h + trip[b][j]] = 1;
        /* G t = t - (3/9) 1, t^T G t = 3 - 9/9 = 2, so g = (t - 1/3)/2 */
        for (int a = 0; a < h; a++) gv[b * h + a] = mul(sub(tv[b * h + a], inv3), inv2);
    }
    LOG("loaded h=%d v=%d frames=%d steps=%d sigmas=%d", h, v, nfr, nst, nsig);
}

/* in-place inverse of n x n (row-major) mod P; returns 0 if singular */
static int matinv(u64 *A, int n) {
    u64 *M = malloc(8 * (size_t)n * 2 * n);
    for (int i = 0; i < n; i++) for (int j = 0; j < 2 * n; j++) M[(size_t)i * 2 * n + j] = j < n ? A[(size_t)i * n + j] : (j - n == i);
    for (int c = 0; c < n; c++) {
        int p = -1; for (int r = c; r < n; r++) if (M[(size_t)r * 2 * n + c]) { p = r; break; }
        if (p < 0) { free(M); return 0; }
        if (p != c) for (int j = 0; j < 2 * n; j++) { u64 t = M[(size_t)p * 2 * n + j]; M[(size_t)p * 2 * n + j] = M[(size_t)c * 2 * n + j]; M[(size_t)c * 2 * n + j] = t; }
        u64 iv = inv(M[(size_t)c * 2 * n + c]); u64 *rc = M + (size_t)c * 2 * n;
        for (int j = 0; j < 2 * n; j++) rc[j] = mul(rc[j], iv);
        for (int r = 0; r < n; r++) if (r != c && M[(size_t)r * 2 * n + c]) {
            u64 fct = M[(size_t)r * 2 * n + c]; u64 *rr = M + (size_t)r * 2 * n;
            for (int j = 0; j < 2 * n; j++) rr[j] = sub(rr[j], mul(fct, rc[j]));
        }
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[(size_t)i * n + j] = M[(size_t)i * 2 * n + n + j];
    free(M); return 1;
}

/* G-orthogonal projector onto the row span of frame F (column convention): P = B^T (B G B^T)^-1 B G */
static int gram_fail = 0;
static void proj(int F, u64 *Pm) {
    int d = fdim[F]; memset(Pm, 0, 8 * h * h);
    if (d == 0) return;
    if (d == h) { for (int i = 0; i < h; i++) Pm[i * h + i] = 1; return; }
    u64 *B = frow[F], BG[23 * 23], Gr[23 * 23], T[23 * 23];
    for (int r = 0; r < d; r++) { u64 s = 0; for (int j = 0; j < h; j++) s = add(s, B[r * h + j]); u64 s9 = mul(s, inv9); for (int j = 0; j < h; j++) BG[r * h + j] = sub(B[r * h + j], s9); }
    for (int a = 0; a < d; a++) for (int b = 0; b < d; b++) { u128 s = 0; for (int j = 0; j < h; j++) s += (u128)BG[a * h + j] * B[b * h + j]; Gr[a * d + b] = fold(s); }
    if (!matinv(Gr, d)) { gram_fail++; return; }
    for (int a = 0; a < d; a++) for (int j = 0; j < h; j++) { u128 s = 0; for (int b = 0; b < d; b++) s += (u128)Gr[a * d + b] * BG[b * h + j]; T[a * h + j] = fold(s); }
    for (int i = 0; i < h; i++) for (int j = 0; j < h; j++) { u128 s = 0; for (int a = 0; a < d; a++) s += (u128)B[a * h + i] * T[a * h + j]; Pm[i * h + j] = fold(s); }
}

/* rank factorization of an h x h idempotent-candidate M of expected rank r: M = A Bt (A h x r, Bt r x h).
   Returns the rank found mod P; checks A Bt == M and M^2 == M. */
static int idem_fail = 0, fact_fail = 0;
static int rankfac(const u64 *M, int r, u64 *A, u64 *Bt) {
    u64 E[23 * 23]; memcpy(E, M, 8 * h * h); int rowp[23], colp[23], rk = 0; int used[23] = {0};
    for (int c = 0; c < h && rk < h; c++) {
        int p = -1; for (int i = 0; i < h; i++) if (!used[i] && E[i * h + c]) { p = i; break; }
        if (p < 0) continue;
        used[p] = 1; rowp[rk] = p; colp[rk] = c; rk++;
        u64 iv = inv(E[p * h + c]);
        for (int i = 0; i < h; i++) if (!used[i] && E[i * h + c]) { u64 fct = mul(E[i * h + c], iv); for (int j = 0; j < h; j++) E[i * h + j] = sub(E[i * h + j], mul(fct, E[p * h + j])); }
    }
    { u64 S[23 * 23]; for (int i = 0; i < h; i++) for (int j = 0; j < h; j++) { u128 s = 0; for (int k = 0; k < h; k++) s += (u128)M[i * h + k] * M[k * h + j]; S[i * h + j] = fold(s); }
      if (memcmp(S, M, 8 * h * h)) idem_fail++; }
    if (rk != r) return rk;
    u64 K[23 * 23];
    for (int i = 0; i < h; i++) for (int c = 0; c < r; c++) A[i * r + c] = M[i * h + colp[c]];
    for (int a = 0; a < r; a++) for (int b = 0; b < r; b++) K[a * r + b] = M[rowp[a] * h + colp[b]];
    if (!matinv(K, r)) { fact_fail++; return -1; }
    for (int a = 0; a < r; a++) for (int j = 0; j < h; j++) { u128 s = 0; for (int b = 0; b < r; b++) s += (u128)K[a * r + b] * M[rowp[b] * h + j]; Bt[a * h + j] = fold(s); }
    for (int i = 0; i < h; i++) for (int j = 0; j < h; j++) { u128 s = 0; for (int c = 0; c < r; c++) s += (u128)A[i * r + c] * Bt[c * h + j]; if (fold(s) != M[i * h + j]) { fact_fail++; return -1; } }
    return rk;
}

static void make_basis(int kind) {
    X = malloc(8 * (size_t)m * m); Xi = malloc(8 * (size_t)m * m);
    for (size_t i = 0; i < (size_t)m * m; i++) { uint64_t xr = rnd() >> 4; X[i] = (u64)(xr % P); }   /* integers in [0, 2^60): the SAME integer matrix for every prime. Entries up to 2e6 were too small: linear forms V' = X^T V with small coefficients vanished exactly over Q on 11 stage-2 edges (b = 312) */
    if (kind == 1) { memset(X, 0, 8 * (size_t)m * m); for (int i = 0; i < m; i++) X[(size_t)i * m + i] = 1; }   /* NEG: X = I */
    memcpy(Xi, X, 8 * (size_t)m * m);
    if (!matinv(Xi, m)) { LOG("X singular mod P"); exit(3); }
    for (int i = 0; i < m; i += 37) for (int j = 0; j < m; j += 41) { u128 s = 0; u64 acc = 0; for (int k = 0; k < m; k++) { s += (u128)X[(size_t)i * m + k] * Xi[(size_t)k * m + j]; if ((k & 31) == 31) { acc = add(acc, fold(s)); s = 0; } } acc = add(acc, fold(s)); if (acc != (u64)(i == j)) { LOG("X Xi != I"); exit(3); } }
    LOG("basis X ready (kind %d), X Xi = I spot-checked", kind);
}
static void make_tables(int rows) {
    size_t n = (size_t)rows * v * h;
    Z1 = malloc(8 * n); Z2 = malloc(8 * n); W1 = malloc(8 * n); W2 = malloc(8 * n);
    u64 *Y = malloc(8 * m), *Xc = malloc(8 * m), rsY[23], csY[23], rsX[23], csX[23];
    for (int i = 0; i < rows; i++) {
        for (int k = 0; k < m; k++) { Y[k] = Xi[(size_t)i * m + k]; Xc[k] = X[(size_t)k * m + i]; }   /* Y[al*h+be] */
        for (int a = 0; a < h; a++) { u64 s = 0, t = 0, s2 = 0, t2 = 0; for (int c = 0; c < h; c++) { s = add(s, Xc[a * h + c]); t = add(t, Xc[c * h + a]); s2 = add(s2, Y[a * h + c]); t2 = add(t2, Y[c * h + a]); } rsX[a] = s; csX[a] = t; rsY[a] = s2; csY[a] = t2; }
        for (int b = 0; b < v; b++) {
            int *T = trip[b];
            for (int a = 0; a < h; a++) {
                Z1[IDX(i, b, a)] = add(add(Y[a * h + T[0]], Y[a * h + T[1]]), Y[a * h + T[2]]);            /* (Y t)[a]   */
                Z2[IDX(i, b, a)] = add(add(Y[T[0] * h + a], Y[T[1] * h + a]), Y[T[2] * h + a]);            /* (t^T Y)[a] */
                u64 x1 = add(add(Xc[a * h + T[0]], Xc[a * h + T[1]]), Xc[a * h + T[2]]);
                u64 x2 = add(add(Xc[T[0] * h + a], Xc[T[1] * h + a]), Xc[T[2] * h + a]);
                W1[IDX(i, b, a)] = mul(sub(x1, mul(rsX[a], inv3)), inv2);                                   /* (Xc g)[a]   */
                W2[IDX(i, b, a)] = mul(sub(x2, mul(csX[a], inv3)), inv2);                                   /* (g^T Xc)[a] */
            }
        }
    }
    (void)rsY; (void)csY; free(Y); free(Xc);
    LOG("tables for %d rows ready (%.0f MB)", rows, 4.0 * 8 * n / 1e6);
}

/* LU without pivoting of r x r M: 1 if all r pivots nonzero */
static int lead_ok(u64 *M, int r) {
    for (int k = 0; k < r; k++) {
        if (!M[k * r + k]) return 0;
        u64 iv = inv(M[k * r + k]);
        for (int i = k + 1; i < r; i++) { u64 fct = mul(M[i * r + k], iv); if (fct) for (int j = k + 1; j < r; j++) M[i * r + j] = sub(M[i * r + j], mul(fct, M[k * r + j])); }
    }
    return 1;
}

static long long nfail = 0; static FILE *FL;
static void fail(const char *cls, long long a, long long b, int stage, int k) { nfail++; fprintf(FL, "%s %lld %lld s%d minor %d\n", cls, a, b, stage, k); }

/* ---- small: every distinct h-step x every invocation b x both stages */
static void run_small(int s0, int s1) {
    u64 PA[529], PB[529], pi[529], A[529], Bt[529], U[23 * 23], V[23 * 23], M[23 * 23];
    long long edges = 0, rankbad = 0; int lastF = -1;
    for (int k = s0; k < s1 && k < nst; k++) {
        int fa = st[k][0], fb = st[k][1], r = st[k][2];
        if (fa != lastF) { proj(fa, PA); lastF = fa; }
        proj(fb, PB);
        for (int i = 0; i < h * h; i++) pi[i] = sub(PB[i], PA[i]);
        int rk = rankfac(pi, r, A, Bt);
        if (rk != r) { rankbad++; fail("rank", k, rk, 0, r); continue; }
        for (int stage = 1; stage <= 2; stage++) {
            u64 *ZT = stage == 1 ? Z1 : Z2, *WT = stage == 1 ? W1 : W2;
            for (int b = 0; b < v; b++) {
                for (int i = 0; i < r; i++) {
                    const u64 *z = ZT + IDX(i, b, 0), *w = WT + IDX(i, b, 0);
                    for (int c = 0; c < r; c++) { u128 s = 0, t = 0; for (int a = 0; a < h; a++) { s += (u128)z[a] * A[a * r + c]; t += (u128)w[a] * Bt[c * h + a]; } U[i * r + c] = fold(s); V[i * r + c] = fold(t); }
                }
                for (int i = 0; i < r; i++) for (int j = 0; j < r; j++) { u128 s = 0; for (int c = 0; c < r; c++) s += (u128)U[i * r + c] * V[j * r + c]; M[i * r + j] = fold(s); }
                if (!lead_ok(M, r)) fail("step", k, b, stage, r);
                edges++;
            }
        }
        if ((k - s0) % 10000 == 9999) LOG("small: %d/%d steps, edges %lld, failures %lld", k + 1, s1, edges, nfail);
    }
    LOG("small DONE steps [%d,%d): edges %lld, rank mismatches %lld, idempotency failures %d, factor failures %d, gram failures %d, failures %lld",
        s0, s1 < nst ? s1 : nst, edges, rankbad, idem_fail, fact_fail, gram_fail, nfail);
}

/* high-rank checker: p = I - U V^T, U', V' rows supplied by a callback, s columns; minors 1..m-s */
typedef void (*rowfn)(int i, u64 *u, u64 *vv, void *ctx);
static int sm_check(int s, rowfn f, void *ctx, int *where) {
    static u64 Ai[46 * 46]; u64 u[46], vv[46], y[46], z[46];
    memset(Ai, 0, 8 * s * s); for (int i = 0; i < s; i++) Ai[i * s + i] = 1;
    for (int k = 0; k < m - s; k++) {
        f(k, u, vv, ctx);
        for (int i = 0; i < s; i++) { u128 a = 0, b = 0; for (int j = 0; j < s; j++) { a += (u128)Ai[i * s + j] * vv[j]; b += (u128)u[j] * Ai[j * s + i]; } y[i] = fold(a); z[i] = fold(b); }
        u128 d = 0; for (int i = 0; i < s; i++) d += (u128)u[i] * y[i];
        u64 del = sub(1, fold(d));
        if (!del) { *where = k + 1; return 0; }
        u64 di = inv(del);
        for (int i = 0; i < s; i++) { u64 yi = mul(y[i], di); u64 *row = Ai + i * s; for (int j = 0; j < s; j++) row[j] = fold((u128)yi * z[j] + row[j]); }
    }
    return 1;
}
struct auxctx { u64 *C, *Dt; int s, b, stage; };
static void aux_row(int i, u64 *u, u64 *vv, void *cx) {
    struct auxctx *c = cx; const u64 *z = (c->stage == 1 ? Z1 : Z2) + IDX(i, c->b, 0), *w = (c->stage == 1 ? W1 : W2) + IDX(i, c->b, 0);
    for (int q = 0; q < c->s; q++) { u128 a = 0, b = 0; for (int al = 0; al < h; al++) { a += (u128)z[al] * c->C[al * c->s + q]; b += (u128)w[al] * c->Dt[q * h + al]; } u[q] = fold(a); vv[q] = fold(b); }
}
static void run_aux(int g0, int g1) {
    u64 Ps[529], Qs[529], C[529], Dt[529]; long long edges = 0;
    for (int g = g0; g < g1 && g < nsig; g++) {
        int F = sig[g], s = h - fdim[F];
        proj(F, Ps); for (int i = 0; i < h; i++) for (int j = 0; j < h; j++) Qs[i * h + j] = sub(i == j, Ps[i * h + j]);
        if (rankfac(Qs, s, C, Dt) != s) { fail("auxrank", g, 0, 0, s); continue; }
        for (int stage = 1; stage <= 2; stage++) for (int b = 0; b < v; b++) {
            struct auxctx cx = { C, Dt, s, b, stage }; int wh = 0;
            if (!sm_check(s, aux_row, &cx, &wh)) fail("aux", g, b, stage, wh);
            edges++;
        }
        if ((g - g0) % 100 == 99) LOG("aux: %d/%d sigmas, edges %lld, failures %lld", g + 1, g1 < nsig ? g1 : nsig, edges, nfail);
    }
    LOG("aux DONE sigmas [%d,%d): edges %lld, idempotency failures %d, gram failures %d, failures %lld", g0, g1 < nsig ? g1 : nsig, edges, idem_fail, gram_fail, nfail);
}
struct entctx { int a, b; };
static void ent_row(int i, u64 *u, u64 *vv, void *cx) {
    struct entctx *c = cx; int a = c->a, b = c->b;
    const u64 *z2 = Z2 + IDX(i, a, 0), *w2 = W2 + IDX(i, a, 0), *z1 = Z1 + IDX(i, b, 0), *w1 = W1 + IDX(i, b, 0);
    const int *T = trip[a]; u64 tz = add(add(z1[T[0]], z1[T[1]]), z1[T[2]]);                 /* t_a . z1 */
    const u64 *ga = gv + (size_t)a * h;
    for (int q = 0; q < h; q++) { u[q] = z2[q]; vv[q] = w2[q]; u[h + q] = sub(z1[q], mul(ga[q], tz)); vv[h + q] = w1[q]; }   /* P_a (x) I | Q_a (x) P_b */
}
static void run_entr(int a0, int a1) {
    long long edges = 0;
    for (int a = a0; a < a1 && a < v; a++) {
        for (int b = 0; b < v; b++) { struct entctx cx = { a, b }; int wh = 0; if (!sm_check(2 * h, ent_row, &cx, &wh)) fail("entrance", a, b, 2, wh); edges++; }
        if ((a - a0) % 20 == 19) LOG("entrance: a %d/%d, pairs %lld, failures %lld", a + 1, a1 < v ? a1 : v, edges, nfail);
    }
    LOG("entrance DONE a in [%d,%d): pairs %lld, failures %lld", a0, a1 < v ? a1 : v, edges, nfail);
}

static double *D1, *D2, *E1, *E2, *gvd;   /* double copies of Z1, Z2, W1, W2, g */
static const double PDd = 2097143.0, PINV = 1.0 / 2097143.0;
static inline double rd(double x) { double r = x - PDd * __builtin_floor(x * PINV); if (r < 0) r += PDd; if (r >= PDd) r -= PDd; return r; }
static void make_dtables(void) {
    size_t n = (size_t)m * v * h;   /* convert in place (same 8-byte width) to stay under the RAM cap */
    D1 = (double *)Z1; D2 = (double *)Z2; E1 = (double *)W1; E2 = (double *)W2;
    for (size_t i = 0; i < n; i++) { u64 a = Z1[i], b = Z2[i], c = W1[i], d = W2[i]; D1[i] = (double)a; D2[i] = (double)b; E1[i] = (double)c; E2[i] = (double)d; }
    gvd = malloc(8 * (size_t)v * h); for (size_t i = 0; i < (size_t)v * h; i++) gvd[i] = gv[i];
    LOG("double tables ready");
}
static inline double dinv(double a) { int64_t r0 = 2097143, r1 = (int64_t)a, t0 = 0, t1 = 1; while (r1) { int64_t q = r0 / r1, t; t = r0 - q * r1; r0 = r1; r1 = t; t = t0 - q * t1; t0 = t1; t1 = t; } if (t0 < 0) t0 += 2097143; return (double)t0; }
typedef void (*drowfn)(int i, double *u, double *vv, void *ctx);
static int smd(int s, drowfn f, void *ctx, int *where) {
    static double Ai[46 * 46]; double u[46], vv[46], y[46], z[46];
    memset(Ai, 0, 8 * s * s); for (int i = 0; i < s; i++) Ai[i * s + i] = 1;
    for (int k = 0; k < m - s; k++) {
        f(k, u, vv, ctx);
        for (int i = 0; i < s; i++) { double a = 0; const double *row = Ai + i * s; for (int j = 0; j < s; j++) a += row[j] * vv[j]; y[i] = rd(a); }
        for (int i = 0; i < s; i++) z[i] = 0;
        for (int j = 0; j < s; j++) { double uj = u[j]; const double *row = Ai + j * s; for (int i = 0; i < s; i++) z[i] += uj * row[i]; }
        for (int i = 0; i < s; i++) z[i] = rd(z[i]);
        double d = 0; for (int i = 0; i < s; i++) d += u[i] * y[i];
        double del = rd(1.0 - rd(d) + PDd);
        if (del == 0) { *where = k + 1; return 0; }
        double di = dinv(del);
        for (int i = 0; i < s; i++) { double yi = rd(y[i] * di); double *row = Ai + i * s; for (int j = 0; j < s; j++) row[j] = rd(row[j] + yi * z[j]); }
    }
    return 1;
}
struct auxd { double *C, *Dt; int s, b, stage; };
static void auxd_row(int i, double *u, double *vv, void *cx) {
    struct auxd *c = cx; const double *z = (c->stage == 1 ? D1 : D2) + IDX(i, c->b, 0), *w = (c->stage == 1 ? E1 : E2) + IDX(i, c->b, 0);
    for (int q = 0; q < c->s; q++) { double a = 0, b = 0; for (int al = 0; al < h; al++) { a += z[al] * c->C[al * c->s + q]; b += w[al] * c->Dt[q * h + al]; } u[q] = rd(a); vv[q] = rd(b); }
}
static void run_auxd(int g0, int g1) {
    u64 Ps[529], Qs[529], C[529], Dt[529]; double Cd[529], Dd[529]; long long edges = 0;
    for (int g = g0; g < g1 && g < nsig; g++) {
        int F = sig[g], s = h - fdim[F];
        proj(F, Ps); for (int i = 0; i < h; i++) for (int j = 0; j < h; j++) Qs[i * h + j] = sub(i == j, Ps[i * h + j]);
        if (rankfac(Qs, s, C, Dt) != s) { fail("auxrank", g, 0, 0, s); continue; }
        for (int i = 0; i < h * s; i++) { Cd[i] = C[i]; Dd[i] = Dt[i]; }
        for (int stage = 1; stage <= 2; stage++) for (int b = 0; b < v; b++) {
            struct auxd cx = { Cd, Dd, s, b, stage }; int wh = 0;
            if (!smd(s, auxd_row, &cx, &wh)) fail("aux", g, b, stage, wh);
            edges++;
        }
        if ((g - g0) % 200 == 199) LOG("aux: %d/%d sigmas, edges %lld, candidates %lld", g + 1, g1 < nsig ? g1 : nsig, edges, nfail);
    }
    LOG("aux DONE sigmas [%d,%d): edges %lld, idempotency failures %d, gram failures %d, zero candidates %lld", g0, g1 < nsig ? g1 : nsig, edges, idem_fail, gram_fail, nfail);
}
struct entd { int a, b; };
static void entd_row(int i, double *u, double *vv, void *cx) {
    struct entd *c = cx; int a = c->a, b = c->b;
    const double *z2 = D2 + IDX(i, a, 0), *w2 = E2 + IDX(i, a, 0), *z1 = D1 + IDX(i, b, 0), *w1 = E1 + IDX(i, b, 0);
    const int *T = trip[a]; double tz = rd(z1[T[0]] + z1[T[1]] + z1[T[2]]); const double *ga = gvd + (size_t)a * h;
    for (int q = 0; q < h; q++) { u[q] = z2[q]; vv[q] = w2[q]; u[h + q] = rd(z1[q] - ga[q] * tz + PDd * PDd); vv[h + q] = w1[q]; }
}
static void run_entrd(int a0, int a1) {
    long long edges = 0;
    for (int a = a0; a < a1 && a < v; a++) {
        for (int b = 0; b < v; b++) { struct entd cx = { a, b }; int wh = 0; if (!smd(2 * h, entd_row, &cx, &wh)) fail("entrance", a, b, 2, wh); edges++; }
        if ((a - a0) % 50 == 49) LOG("entrance: a %d/%d, pairs %lld, candidates %lld", a + 1, a1 < v ? a1 : v, edges, nfail);
    }
    LOG("entrance DONE a in [%d,%d): pairs %lld, zero candidates %lld", a0, a1 < v ? a1 : v, edges, nfail);
}
static void run_copy(void) {
    long long edges = 0;
    for (int a = 0; a < v; a++) for (int b = 0; b < v; b++) {
        /* P_a (x) P_b = (t_a (x) t_b)(g_a (x) g_b)^T; minor 1 = (Xi[0].(t_a(x)t_b)) ((g_a(x)g_b).X[:,0]) */
        const int *T = trip[a]; const u64 *z = Z1 + IDX(0, b, 0), *w = W1 + IDX(0, b, 0), *ga = gv + (size_t)a * h;
        u64 uu = add(add(z[T[0]], z[T[1]]), z[T[2]]); u128 s = 0; for (int q = 0; q < h; q++) s += (u128)ga[q] * w[q];
        if (!mul(uu, fold(s))) fail("copy", a, b, 0, 1);
        edges++;
    }
    LOG("copy DONE pairs %lld, failures %lld", edges, nfail);
}

/* ---- dense cross-check: build p explicitly (m x m) for sampled edges of every class, conjugate, and run plain
   elimination without pivoting; must agree with the structured result.  Also the rank-r tail must vanish. */
static void kron(const u64 *A, const u64 *B, u64 *K) { for (int i = 0; i < h; i++) for (int j = 0; j < h; j++) for (int k = 0; k < h; k++) for (int l = 0; l < h; l++) K[(size_t)(i * h + k) * m + j * h + l] = mul(A[i * h + j], B[k * h + l]); }
static void mm(const u64 *A, const u64 *B, u64 *C, int n) { for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { u128 s = 0; u64 acc = 0; for (int k = 0; k < n; k++) { s += (u128)A[(size_t)i * n + k] * B[(size_t)k * n + j]; if ((k & 31) == 31) { acc = add(acc, fold(s)); s = 0; } } C[(size_t)i * n + j] = add(acc, fold(s)); } }
static int dense_profile(u64 *p, int r) {   /* returns number of nonzero leading pivots before first zero, and checks rank r */
    u64 *T = malloc(8 * (size_t)m * m), *Mx = malloc(8 * (size_t)m * m);
    mm(Xi, p, T, m); mm(T, X, Mx, m);
    int k;
    for (k = 0; k < m; k++) {
        u64 pv = Mx[(size_t)k * m + k]; if (!pv) break;
        u64 iv = inv(pv);
        for (int i = k + 1; i < m; i++) { u64 fct = mul(Mx[(size_t)i * m + k], iv); if (fct) for (int j = k; j < m; j++) Mx[(size_t)i * m + j] = sub(Mx[(size_t)i * m + j], mul(fct, Mx[(size_t)k * m + j])); }
    }
    int tailzero = 1; if (k == r) for (size_t i = r; i < (size_t)m && tailzero; i++) for (size_t j = r; j < (size_t)m; j++) if (Mx[i * m + j]) { tailzero = 0; break; }
    free(T); free(Mx);
    return k == r && tailzero ? r : -1 - k;
}
static void Pline(int b, u64 *Pm) { for (int i = 0; i < h; i++) for (int j = 0; j < h; j++) Pm[i * h + j] = mul(tv[b * h + i], gv[b * h + j]); }
static void run_dense(int nper) {
    u64 *p = malloc(8 * (size_t)m * m); u64 PA[529], PB[529], pi[529], Pb[529], Pa[529], Qa[529], Qb[529], I[529];
    memset(I, 0, sizeof I); for (int i = 0; i < h; i++) I[i * h + i] = 1;
    int ok = 0, tot = 0;
    for (int t = 0; t < nper; t++) {
        int k = rnd() % nst, b = rnd() % v, a = rnd() % v, g = rnd() % nsig;
        proj(st[k][0], PA); proj(st[k][1], PB); for (int i = 0; i < 529; i++) pi[i] = sub(PB[i], PA[i]);
        Pline(b, Pb); Pline(a, Pa); for (int i = 0; i < 529; i++) { Qa[i] = sub(I[i], Pa[i]); Qb[i] = sub(I[i], Pb[i]); }
        int r = st[k][2], res;
        kron(pi, Pb, p); res = dense_profile(p, r); tot++; ok += res == r; LOG("dense step %d s1 b=%d r=%d -> %d", k, b, r, res);
        kron(Pb, pi, p); res = dense_profile(p, r); tot++; ok += res == r; LOG("dense step %d s2 b=%d r=%d -> %d", k, b, r, res);
        u64 Ps[529], Qs[529]; int F = sig[g], s = h - fdim[F]; proj(F, Ps); for (int i = 0; i < 529; i++) Qs[i] = sub(I[i], Ps[i]);
        kron(Qs, Pb, p); for (size_t i = 0; i < (size_t)m * m; i++) p[i] = sub((i / m) == (i % m), p[i]); res = dense_profile(p, m - s); tot++; ok += res == m - s; LOG("dense aux sigma %d s1 b=%d rank %d -> %d", g, b, m - s, res);
        kron(Pb, Qs, p); for (size_t i = 0; i < (size_t)m * m; i++) p[i] = sub((i / m) == (i % m), p[i]); res = dense_profile(p, m - s); tot++; ok += res == m - s; LOG("dense aux sigma %d s2 b=%d rank %d -> %d", g, b, m - s, res);
        kron(Qa, Qb, p); res = dense_profile(p, m - 2 * h + 1); tot++; ok += res == m - 2 * h + 1; LOG("dense entrance (%d,%d) rank %d -> %d", a, b, m - 2 * h + 1, res);
        kron(Pa, Pb, p); res = dense_profile(p, 1); tot++; ok += res == 1; LOG("dense copy (%d,%d) -> %d", a, b, res);
    }
    LOG("dense DONE %d/%d agree (full rank-r profile incl. vanishing tail)", ok, tot);
}

int main(int argc, char **argv) {
    T0 = now(); if (argc < 3) return 1;
    load(argv[1]); const char *mode = argv[2];
    FL = stdout;
    if (!strcmp(mode, "dense")) { make_basis(0); run_dense(argc > 3 ? atoi(argv[3]) : 2); return 0; }
    int negI = argc > 5 && !strcmp(argv[5], "negI");
    make_basis(negI ? 1 : 0);
    int rows = !strcmp(mode, "small") ? 22 : !strcmp(mode, "copy") ? 1 : m;
    make_tables(rows);
    if (negI) LOG("NEGATIVE CONTROL: X = I");
    int lo = argc > 3 ? atoi(argv[3]) : 0, hi = argc > 4 ? atoi(argv[4]) : 1 << 30;
    if (!strcmp(mode, "small")) run_small(lo, hi);
    else if (!strcmp(mode, "aux")) { make_dtables(); run_auxd(lo, hi); }
    else if (!strcmp(mode, "entr")) { make_dtables(); run_entrd(lo, hi); }
    else if (!strcmp(mode, "copy")) run_copy();
    printf("RESULT %s [%d,%d) failures %lld\n", mode, lo, hi, nfail);
    return 0;
}
