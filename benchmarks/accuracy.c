/* benchmarks/accuracy.c -- bench_fft --accuracy: the accuracy of pffft against
 * an exact transform, and the equivalences that must hold bit for bit within
 * one build of pffft.
 *
 * For every transform kind, size and precision, each pffft (the one linked
 * into bench_fft, or each --lib) gets the same random inputs, uniform in
 * [-1, 1) and rounded to the precision: as many as make at least
 * --seeds K * MIN_POINTS points, at least one. Over all its inputs a case is
 * measured by
 *
 *   fwd  the relative RMS error of the ordered forward transform of a random
 *        signal x: sqrt(sum |y - ref|^2 / sum |ref|^2), ref the exact DFT;
 *   bwd  the same for the ordered backward transform of a random spectrum
 *        (any values in the ordered layout are a valid spectrum), against
 *        the exact inverse DFT (unscaled);
 *   rt   the relative RMS error of the round trip, backward(forward(x)) / n
 *        against x: for information only, as errors of the twiddles largely
 *        cancel in it (the backward pass uses the forward one's, conjugated);
 *
 * and fwd_max, bwd_max, rt_max, the largest of these over the single inputs.
 * Pooled over this many points the errors of two builds compare to about a
 * percent at every size. With bits_only (bench_fft --bits-only) the
 * reference is not computed and fwd and bwd are left empty: what remains, the
 * checks and the hashes, is all that a change that moves data but computes
 * nothing new needs, and it is many times faster.
 *
 * The reference is a mixed-radix FFT in double-double arithmetic (about 106
 * bits), with twiddles from an exact integer reduction of the angle to the
 * first octant and sine and cosine from their series: its error is some
 * 1e-30, far below pffft's in both precisions. (Another FFT in double
 * precision, such as PocketFFT, would not do for f64: its error is the size
 * of pffft's and correlated with it, as both descend from FFTPACK, so the
 * distance to it rewards resembling it rather than being accurate.)
 *
 * Within one build these must hold bit for bit (memcmp), else a check fails:
 *
 *   ip-fwd, ip-bwd    the ordered transform in place equals the one out of
 *                     place, forward and backward;
 *   ip-un             the same for the unordered forward transform;
 *   reorder-fwd       the ordered forward transform equals pffft_transform
 *                     followed by pffft_zreorder forward;
 *   reorder-bwd       the ordered backward transform equals pffft_zreorder
 *                     backward followed by pffft_transform;
 *   fwd-limit, bwd-limit, rt-limit  the errors of single inputs stay below
 *                     LIMIT_F32 or LIMIT_F64, about ten times what pffft
 *                     gives.
 *
 * Each row of the CSV (one per library and case) also carries a hash of the
 * ordered forward results, of the ordered backward results and of the
 * unordered forward results of all inputs, so that a script can tell whether
 * two builds (split and unsplit, two revisions) give the same bits:
 * benchmarks/accuracy.py does.
 *
 * The libraries are called through pointers that spell setups and arrays as
 * void * and enums as int, as the engines of engines.c do: a --lib may be of
 * any revision, and no header describes it.
 */
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "bench.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define LIMIT_F32 2e-6
#define LIMIT_F64 5e-15
#define MIN_POINTS 4096

typedef void *(*new_setup_fn)(int n, int transform);
typedef int (*is_valid_size_fn)(int n, int transform);
typedef void (*destroy_setup_fn)(void *setup);
typedef void (*transform_fn)(const void *setup, const void *in, void *out, void *work, int direction);
typedef void (*zreorder_fn)(const void *setup, const void *in, void *out, int direction);

/* The functions of one pffft in one precision; new_setup is NULL if the
 * library lacks the precision. */
typedef struct api {
  new_setup_fn new_setup;
  is_valid_size_fn is_valid_size; /* NULL in old revisions: new_setup refuses */
  destroy_setup_fn destroy_setup;
  transform_fn transform, transform_ordered;
  zreorder_fn zreorder;
} api;

/* The pffft linked into bench_fft. Declared here rather than taken from
 * pffft.h, in the spelling of the pointers above. */
#ifdef BENCH_HAVE_F32
void *pffft_new_setup(int n, int transform);
int pffft_is_valid_size(int n, int transform);
void pffft_destroy_setup(void *setup);
void pffft_transform(const void *setup, const void *in, void *out, void *work, int direction);
void pffft_transform_ordered(const void *setup, const void *in, void *out, void *work, int direction);
void pffft_zreorder(const void *setup, const void *in, void *out, int direction);
#endif
#ifdef BENCH_HAVE_F64
void *pffftd_new_setup(int n, int transform);
int pffftd_is_valid_size(int n, int transform);
void pffftd_destroy_setup(void *setup);
void pffftd_transform(const void *setup, const void *in, void *out, void *work, int direction);
void pffftd_transform_ordered(const void *setup, const void *in, void *out, void *work, int direction);
void pffftd_zreorder(const void *setup, const void *in, void *out, int direction);
#endif

static void builtin_api(int f64, api *a) {
  memset(a, 0, sizeof *a);
#ifdef BENCH_HAVE_F32
  if (!f64) {
    a->new_setup = pffft_new_setup;
    a->is_valid_size = pffft_is_valid_size;
    a->destroy_setup = pffft_destroy_setup;
    a->transform = pffft_transform;
    a->transform_ordered = pffft_transform_ordered;
    a->zreorder = pffft_zreorder;
  }
#endif
#ifdef BENCH_HAVE_F64
  if (f64) {
    a->new_setup = pffftd_new_setup;
    a->is_valid_size = pffftd_is_valid_size;
    a->destroy_setup = pffftd_destroy_setup;
    a->transform = pffftd_transform;
    a->transform_ordered = pffftd_transform_ordered;
    a->zreorder = pffftd_zreorder;
  }
#endif
}

/* Sets *fn, a function pointer of any type, to a function of the library, or
 * to NULL if it lacks it (ISO C does not convert an object pointer to a
 * function pointer, so the bits are copied). */
static void lib_function(void *lib, const char *prefix, const char *name, void *fn) {
  char symbol[64];
  void *address;
  snprintf(symbol, sizeof symbol, "%s_%s", prefix, name);
  address = bench_dl_sym(lib, symbol);
  memcpy(fn, &address, sizeof address);
}

static void lib_api(void *lib, int f64, api *a) {
  const char *prefix = f64 ? "pffftd" : "pffft";
  memset(a, 0, sizeof *a);
  lib_function(lib, prefix, "new_setup", &a->new_setup);
  lib_function(lib, prefix, "is_valid_size", &a->is_valid_size);
  lib_function(lib, prefix, "destroy_setup", &a->destroy_setup);
  lib_function(lib, prefix, "transform", &a->transform);
  lib_function(lib, prefix, "transform_ordered", &a->transform_ordered);
  lib_function(lib, prefix, "zreorder", &a->zreorder);
  if (!a->destroy_setup || !a->transform || !a->transform_ordered || !a->zreorder) a->new_setup = NULL;
}

/* splitmix64 */
static uint64_t next_random(uint64_t *state) {
  uint64_t z = (*state += 0x9E3779B97F4A7C15u);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
  return z ^ (z >> 31);
}

/* A 64-bit hash of a byte array whose length is a multiple of 8. */
static uint64_t hash_bytes(const void *p, size_t bytes) {
  const unsigned char *b = (const unsigned char *)p;
  uint64_t h = 0x84222325CBF29CE4u ^ bytes, w;
  size_t i;
  for (i = 0; i + 8 <= bytes; i += 8) {
    memcpy(&w, b + i, 8);
    h = (h ^ w) * 0x100000001B3u;
    h ^= h >> 29;
  }
  return h;
}

static void put(void *a, size_t i, double v, int f64) {
  if (f64) ((double *)a)[i] = v;
  else ((float *)a)[i] = (float)v;
}


/* ------------------------------------------------------------------------
 * The reference: double-double arithmetic (a value is hi + lo, |lo| <= half
 * an ulp of hi), after Dekker and Knuth; Bailey's QD library has the same
 * operations. The error-free products need a fused multiply-add where the
 * target has one, and Dekker's splitting elsewhere.
 * ---------------------------------------------------------------------- */

typedef struct dd {
  double hi, lo;
} dd;

typedef struct cdd {
  dd re, im;
} cdd;

static dd dd_make(double hi, double lo) {
  dd r;
  r.hi = hi;
  r.lo = lo;
  return r;
}

/* a + b exactly, |a| >= |b| */
static dd quick_two_sum(double a, double b) {
  double s = a + b;
  return dd_make(s, b - (s - a));
}

/* a + b exactly */
static dd two_sum(double a, double b) {
  double s = a + b, bb = s - a;
  return dd_make(s, (a - (s - bb)) + (b - bb));
}

/* a * b exactly */
static dd two_prod(double a, double b) {
  double p = a * b;
#if defined(FP_FAST_FMA)
  return dd_make(p, fma(a, b, -p));
#else
  const double split = 134217729.0; /* 2^27 + 1 */
  double t = split * a, ah = t - (t - a), al = a - ah;
  double u = split * b, bh = u - (u - b), bl = b - bh;
  return dd_make(p, ((ah * bh - p) + ah * bl + al * bh) + al * bl);
#endif
}

static dd dd_add(dd a, dd b) {
  dd s = two_sum(a.hi, b.hi), t = two_sum(a.lo, b.lo);
  s.lo += t.hi;
  s = quick_two_sum(s.hi, s.lo);
  s.lo += t.lo;
  return quick_two_sum(s.hi, s.lo);
}

static dd dd_neg(dd a) {
  return dd_make(-a.hi, -a.lo);
}

static dd dd_sub(dd a, dd b) {
  return dd_add(a, dd_neg(b));
}

static dd dd_mul(dd a, dd b) {
  dd p = two_prod(a.hi, b.hi);
  p.lo += a.hi * b.lo + a.lo * b.hi;
  return quick_two_sum(p.hi, p.lo);
}

static dd dd_div_d(dd a, double b) {
  double q1 = a.hi / b;
  dd p = two_prod(q1, b), s = two_sum(a.hi, -p.hi);
  s.lo += a.lo - p.lo;
  return quick_two_sum(q1, (s.hi + s.lo) / b);
}

static cdd cdd_add(cdd a, cdd b) {
  cdd r;
  r.re = dd_add(a.re, b.re);
  r.im = dd_add(a.im, b.im);
  return r;
}

static cdd cdd_sub(cdd a, cdd b) {
  cdd r;
  r.re = dd_sub(a.re, b.re);
  r.im = dd_sub(a.im, b.im);
  return r;
}

static cdd cdd_mul(cdd a, cdd b) {
  cdd r;
  r.re = dd_sub(dd_mul(a.re, b.re), dd_mul(a.im, b.im));
  r.im = dd_add(dd_mul(a.re, b.im), dd_mul(a.im, b.re));
  return r;
}

/* a * (sign * i) */
static cdd cdd_mul_i(cdd a, int sign) {
  cdd r;
  r.re = sign > 0 ? dd_neg(a.im) : a.im;
  r.im = sign > 0 ? a.re : dd_neg(a.re);
  return r;
}

/* cos and sin of x, 0 <= x < 7, by their series */
static void dd_sincos(dd x, dd *c, dd *s) {
  dd x2 = dd_mul(x, x), ts = x, tc = dd_make(1, 0);
  int k;
  *s = ts;
  *c = tc;
  for (k = 1; k < 200; k += 2) {
    ts = dd_div_d(dd_mul(ts, x2), -(double)(k + 1) * (k + 2));
    tc = dd_div_d(dd_mul(tc, x2), -(double)k * (k + 1));
    *s = dd_add(*s, ts);
    *c = dd_add(*c, tc);
    if (fabs(ts.hi) < 1e-36 && fabs(tc.hi) < 1e-36) break;
  }
}

/* 2 pi j / n in double-double, for 0 <= j <= n */
static dd dd_angle(size_t j, size_t n) {
  static const dd two_pi = {6.283185307179586232, 2.449293598294706e-16};
  return dd_div_d(dd_mul(two_pi, dd_make((double)j, 0)), (double)n);
}

/* The twiddles exp(+-2 pi i j / N). If 8 divides N, tab holds cos and sin of
 * 2 pi i / N for i = 0 .. N/8, the first octant, from which the others follow
 * by symmetry (an exact reduction of the angle); else i = 0 .. N - 1. Each
 * entry is the product of a coarse and a fine one (angles a * B and b), both
 * from the series. */
typedef struct twiddles {
  size_t N, M; /* M = N/8 if octants, else 0 */
  cdd *tab;
} twiddles;

static int twiddles_init(twiddles *t, size_t N) {
  size_t L = (N % 8 == 0) ? N / 8 : N - 1, B = 1, i;
  cdd *coarse, *fine;
  while (B * B < L + 1) B++;
  t->N = N;
  t->M = (N % 8 == 0) ? N / 8 : 0;
  t->tab = (cdd *)malloc((L + 1) * sizeof(cdd));
  coarse = (cdd *)malloc((L / B + 1) * sizeof(cdd));
  fine = (cdd *)malloc(B * sizeof(cdd));
  if (!t->tab || !coarse || !fine) {
    free(t->tab);
    free(coarse);
    free(fine);
    return 1;
  }
  for (i = 0; i <= L / B; i++) dd_sincos(dd_angle(i * B, N), &coarse[i].re, &coarse[i].im);
  for (i = 0; i < B; i++) dd_sincos(dd_angle(i, N), &fine[i].re, &fine[i].im);
  for (i = 0; i <= L; i++) t->tab[i] = (i % B) ? cdd_mul(coarse[i / B], fine[i % B]) : coarse[i / B];
  free(coarse);
  free(fine);
  return 0;
}

/* exp(sign * 2 pi i j / N) */
static cdd twiddle(const twiddles *t, size_t j, int sign) {
  cdd w;
  j %= t->N;
  if (t->M) {
    size_t M = t->M, o = j / M, r = j % M;
    cdd a = t->tab[(o & 1) ? M - r : r];
    dd c, s;
    switch (o & 3) {
    case 0: c = a.re; s = a.im; break;                 /* theta = phi */
    case 1: c = a.im; s = a.re; break;                 /* pi/2 - psi */
    case 2: c = dd_neg(a.im); s = a.re; break;         /* pi/2 + phi */
    default: c = dd_neg(a.re); s = a.im; break;        /* pi - psi */
    }
    if (o & 4) {                                      /* pi + ... */
      c = dd_neg(c);
      s = dd_neg(s);
    }
    w.re = c;
    w.im = s;
  } else {
    w = t->tab[j];
  }
  if (sign < 0) w.im = dd_neg(w.im);
  return w;
}

/* The input of the reference: element i of the sequence transformed, read
 * from the doubles x in one of these layouts */
enum { IN_REAL, IN_CPLX, IN_HERMITIAN };

static cdd load(const double *x, int layout, size_t n, size_t i) {
  cdd v;
  v.im = dd_make(0, 0);
  if (layout == IN_REAL) {
    v.re = dd_make(x[i], 0);
  } else if (layout == IN_CPLX) {
    v.re = dd_make(x[2 * i], 0);
    v.im = dd_make(x[2 * i + 1], 0);
  } else { /* pffft's ordered real spectrum: r0, r(n/2), r1, i1, ... */
    if (i == 0) v.re = dd_make(x[0], 0);
    else if (2 * i == n) v.re = dd_make(x[1], 0);
    else if (2 * i < n) {
      v.re = dd_make(x[2 * i], 0);
      v.im = dd_make(x[2 * i + 1], 0);
    } else {
      v.re = dd_make(x[2 * (n - i)], 0);
      v.im = dd_make(-x[2 * (n - i) + 1], 0);
    }
  }
  return v;
}

/* out[k] = sum_j in(off + j * stride) * exp(sign 2 pi i j k / n), n = N / s:
 * decimation in time, radix 4, 2, 3 and 5. Returns nonzero for another
 * prime factor. */
static int dd_fft(const double *x, int layout, size_t off, size_t stride, cdd *out, size_t n, size_t s,
                  const twiddles *t, int sign) {
  size_t p, m, k, r, q;
  if (n == 1) {
    out[0] = load(x, layout, t->N, off);
    return 0;
  }
  p = n % 4 == 0 ? 4 : n % 2 == 0 ? 2 : n % 3 == 0 ? 3 : n % 5 == 0 ? 5 : 0;
  if (!p) return 1;
  m = n / p;
  for (r = 0; r < p; r++)
    if (dd_fft(x, layout, off + r * stride, stride * p, out + r * m, m, s * p, t, sign)) return 1;
  for (k = 0; k < m; k++) {
    cdd a[5], y[5];
    a[0] = out[k];
    for (r = 1; r < p; r++) a[r] = cdd_mul(out[r * m + k], twiddle(t, r * k * s, sign));
    if (p == 2) {
      y[0] = cdd_add(a[0], a[1]);
      y[1] = cdd_sub(a[0], a[1]);
    } else if (p == 4) {
      cdd s02 = cdd_add(a[0], a[2]), d02 = cdd_sub(a[0], a[2]);
      cdd s13 = cdd_add(a[1], a[3]), d13 = cdd_mul_i(cdd_sub(a[1], a[3]), sign);
      y[0] = cdd_add(s02, s13);
      y[2] = cdd_sub(s02, s13);
      y[1] = cdd_add(d02, d13);
      y[3] = cdd_sub(d02, d13);
    } else {
      for (q = 0; q < p; q++) {
        y[q] = a[0];
        for (r = 1; r < p; r++) y[q] = cdd_add(y[q], cdd_mul(a[r], twiddle(t, (r * q % p) * (t->N / p), sign)));
      }
    }
    for (q = 0; q < p; q++) out[k + q * m] = y[q];
  }
  return 0;
}

/* ------------------------------------------------------------------------
 * The measurement
 * ---------------------------------------------------------------------- */

/* What one library gave over the inputs of a case. */
typedef struct tally {
  double err[3], ref[3]; /* sums of squares, for fwd, bwd, rt */
  double max[3];         /* the largest relative RMS error of one input */
  uint64_t hfwd, hbwd, hun;
  unsigned checks; /* failed checks, by bit */
  int inputs;
} tally;

enum { FWD, BWD, RT };
static const char *const check_names[] = {"ip-fwd", "reorder-fwd", "ip-un", "ip-bwd", "reorder-bwd",
                                          "fwd-limit", "bwd-limit", "rt-limit"};
enum { IP_FWD = 1, REORDER_FWD = 2, IP_UN = 4, IP_BWD = 8, REORDER_BWD = 16, FWD_LIMIT = 32, BWD_LIMIT = 64,
       RT_LIMIT = 128 };
#define NCHECKS 8

static uint64_t mix(uint64_t h, uint64_t v) {
  return (h ^ v) * 0x9E3779B97F4A7C15u + 0x632BE59BD9B4E019u;
}

/* element i of the exact result in pffft's ordered layout */
static dd ref_at(const cdd *r, int layout, size_t n, size_t i) {
  if (layout == IN_CPLX) return (i & 1) ? r[i / 2].im : r[i / 2].re; /* a complex result */
  if (layout == IN_REAL) return r[i].re;                             /* a real signal */
  if (i == 0) return r[0].re;                                        /* a real spectrum */
  if (i == 1) return r[n / 2].re;
  return (i & 1) ? r[i / 2].im : r[i / 2].re;
}

/* adds the error of got against the exact r to t, as measure m */
static void add_error(tally *t, int m, const cdd *r, int layout, const void *got, size_t len, size_t n,
                      int f64, double limit) {
  double e = 0, s = 0, rel;
  size_t i;
  for (i = 0; i < len; i++) {
    dd v = ref_at(r, layout, n, i);
    double g = f64 ? ((const double *)got)[i] : (double)((const float *)got)[i];
    double d = (g - v.hi) - v.lo;
    e += d * d;
    s += v.hi * v.hi;
  }
  t->err[m] += e;
  t->ref[m] += s;
  rel = sqrt(e / (s > 0 ? s : 1));
  if (rel > t->max[m]) t->max[m] = rel;
  if (rel > limit) t->checks |= (m == FWD ? FWD_LIMIT : m == BWD ? BWD_LIMIT : RT_LIMIT);
}

/* the arrays of one case: x, spec and the exact results r in doubles, the
 * others of the precision, all of len scalars (r: n complex) */
typedef struct arrays {
  double *x, *spec;
  cdd *r;
  void *in, *out, *back, *tmp, *tmp2, *work;
} arrays;

static void free_arrays(arrays *a) {
  bench_aligned_free(a->x);
  bench_aligned_free(a->spec);
  bench_aligned_free(a->r);
  bench_aligned_free(a->in);
  bench_aligned_free(a->out);
  bench_aligned_free(a->back);
  bench_aligned_free(a->tmp);
  bench_aligned_free(a->tmp2);
  bench_aligned_free(a->work);
  memset(a, 0, sizeof *a);
}

static int alloc_arrays(arrays *a, size_t len, size_t n, int f64) {
  size_t bytes = len * (f64 ? sizeof(double) : sizeof(float));
  memset(a, 0, sizeof *a);
  a->x = bench_aligned_alloc(len * sizeof(double));
  a->spec = bench_aligned_alloc(len * sizeof(double));
  a->r = n ? bench_aligned_alloc(n * sizeof(cdd)) : NULL; /* n = 0: no exact results */
  a->in = bench_aligned_alloc(bytes);
  a->out = bench_aligned_alloc(bytes);
  a->back = bench_aligned_alloc(bytes);
  a->tmp = bench_aligned_alloc(bytes);
  a->tmp2 = bench_aligned_alloc(bytes);
  a->work = bench_aligned_alloc(bytes);
  if (a->x && a->spec && (a->r || !n) && a->in && a->out && a->back && a->tmp && a->tmp2 && a->work) return 0;
  free_arrays(a);
  return 1;
}

/* random values uniform in [-1, 1), rounded to the precision */
static void random_fill(double *v, size_t len, int f64, uint64_t *state) {
  size_t i;
  for (i = 0; i < len; i++) {
    double u = (double)(int64_t)(next_random(state) >> 11) * (2.0 / 9007199254740992.0) - 1.0;
    v[i] = f64 ? u : (double)(float)u;
  }
}

/* Forward: the ordered transform of x (a->r holds its exact DFT), in and out
 * of place, unordered and reordered, and the round trip. */
static void measure_forward(tally *t, const api *f, const arrays *a, int n, int cplx, int f64, int exact) {
  const int fwd = 0, bwd = 1; /* PFFFT_FORWARD, PFFFT_BACKWARD */
  size_t len = (size_t)(cplx ? 2 * n : n), bytes = len * (f64 ? sizeof(double) : sizeof(float)), i;
  double limit = f64 ? LIMIT_F64 : LIMIT_F32;
  void *s = f->new_setup(n, cplx);
  if (!s) return;
  for (i = 0; i < len; i++) put(a->in, i, a->x[i], f64);

  f->transform_ordered(s, a->in, a->out, a->work, fwd);
  if (exact) add_error(t, FWD, a->r, cplx ? IN_CPLX : IN_HERMITIAN, a->out, len, (size_t)n, f64, limit);
  t->hfwd = mix(t->hfwd, hash_bytes(a->out, bytes));
  memcpy(a->tmp, a->in, bytes);
  f->transform_ordered(s, a->tmp, a->tmp, a->work, fwd);
  if (memcmp(a->tmp, a->out, bytes)) t->checks |= IP_FWD;

  f->transform(s, a->in, a->tmp, a->work, fwd);
  t->hun = mix(t->hun, hash_bytes(a->tmp, bytes));
  f->zreorder(s, a->tmp, a->tmp2, fwd);
  if (memcmp(a->tmp2, a->out, bytes)) t->checks |= REORDER_FWD;
  memcpy(a->tmp2, a->in, bytes);
  f->transform(s, a->tmp2, a->tmp2, a->work, fwd);
  if (memcmp(a->tmp2, a->tmp, bytes)) t->checks |= IP_UN;

  /* the round trip, against x itself */
  f->transform_ordered(s, a->out, a->back, a->work, bwd);
  {
    double e = 0, r = 0, rel;
    for (i = 0; i < len; i++) {
      double g = (f64 ? ((const double *)a->back)[i] : (double)((const float *)a->back)[i]) / n;
      e += (g - a->x[i]) * (g - a->x[i]);
      r += a->x[i] * a->x[i];
    }
    t->err[RT] += e;
    t->ref[RT] += r;
    rel = sqrt(e / (r > 0 ? r : 1));
    if (rel > t->max[RT]) t->max[RT] = rel;
    if (rel > limit) t->checks |= RT_LIMIT;
  }
  f->destroy_setup(s);
  t->inputs++;
}

/* Backward: the ordered transform of the spectrum spec (a->r holds its exact
 * inverse DFT), in and out of place, and reordered first. */
static void measure_backward(tally *t, const api *f, const arrays *a, int n, int cplx, int f64, int exact) {
  const int bwd = 1;
  size_t len = (size_t)(cplx ? 2 * n : n), bytes = len * (f64 ? sizeof(double) : sizeof(float)), i;
  double limit = f64 ? LIMIT_F64 : LIMIT_F32;
  void *s = f->new_setup(n, cplx);
  if (!s) return;
  for (i = 0; i < len; i++) put(a->in, i, a->spec[i], f64);

  f->transform_ordered(s, a->in, a->back, a->work, bwd);
  if (exact) add_error(t, BWD, a->r, cplx ? IN_CPLX : IN_REAL, a->back, len, (size_t)n, f64, limit);
  t->hbwd = mix(t->hbwd, hash_bytes(a->back, bytes));
  memcpy(a->tmp2, a->in, bytes);
  f->transform_ordered(s, a->tmp2, a->tmp2, a->work, bwd);
  if (memcmp(a->tmp2, a->back, bytes)) t->checks |= IP_BWD;

  f->zreorder(s, a->in, a->tmp2, bwd);
  f->transform(s, a->tmp2, a->tmp, a->work, bwd);
  if (memcmp(a->tmp, a->back, bytes)) t->checks |= REORDER_BWD;
  f->destroy_setup(s);
}

/* Writes the CSV row of a library and case; returns its failed checks. */
static int report(FILE *csv, const char *label, const tally *t, int n, int cplx, int f64) {
  static const char *const prec_names[2] = {"f32", "f64"}, *const xform_names[2] = {"real", "cplx"};
  char checks[160] = "";
  int failed = 0, c;
  for (c = 0; c < NCHECKS; c++) {
    if (!(t->checks & (1u << c))) continue;
    if (checks[0]) strcat(checks, ";");
    strcat(checks, check_names[c]);
    failed++;
  }
  if (csv) {
    char fwd[24] = "", bwd[24] = "", fwd_max[24] = "", bwd_max[24] = ""; /* empty with bits_only */
    if (t->ref[FWD] > 0) {
      snprintf(fwd, sizeof fwd, "%.4e", sqrt(t->err[FWD] / t->ref[FWD]));
      snprintf(fwd_max, sizeof fwd_max, "%.4e", t->max[FWD]);
    }
    if (t->ref[BWD] > 0) {
      snprintf(bwd, sizeof bwd, "%.4e", sqrt(t->err[BWD] / t->ref[BWD]));
      snprintf(bwd_max, sizeof bwd_max, "%.4e", t->max[BWD]);
    }
    fprintf(csv, "%s,%s,%s,%d,%d,%s,%s,%.4e,%s,%s,%.4e,%s,%016llx,%016llx,%016llx\n", label, prec_names[f64],
            xform_names[cplx], n, t->inputs, fwd, bwd, sqrt(t->err[RT] / t->ref[RT]), fwd_max, bwd_max, t->max[RT],
            checks, (unsigned long long)t->hfwd, (unsigned long long)t->hbwd, (unsigned long long)t->hun);
  }
  if (failed)
    fprintf(stderr, "bench_fft --accuracy: %s %s %s n=%d fails: %s (largest fwd %.3g, bwd %.3g, rt %.3g)\n", label,
            prec_names[f64], xform_names[cplx], n, checks, t->max[FWD], t->max[BWD], t->max[RT]);
  return failed;
}

void bench_accuracy_columns(FILE *csv) {
  fprintf(csv, "engine,prec,xform,n,inputs,fwd,bwd,rt,fwd_max,bwd_max,rt_max,checks,hash_fwd,hash_bwd,hash_un\n");
}

int bench_accuracy(FILE *csv, const bench_accuracy_lib *libs, int nlibs, const int *sizes, int nsizes,
                   const int *prec, const int *xform, int seeds, unsigned long seed, int bits_only) {
  static const char *const prec_names[2] = {"f32", "f64"}, *const xform_names[2] = {"real", "cplx"};
  int failed = 0, measured = 0;
  for (int f64 = 0; f64 < 2; f64++) {
    api apis[BENCH_MAX_ENGINES];
    if (!prec[f64]) continue;
    for (int l = 0; l < nlibs && l < BENCH_MAX_ENGINES; l++) {
      if (libs[l].lib) lib_api(libs[l].lib, f64, &apis[l]);
      else builtin_api(f64, &apis[l]);
    }
    for (int cplx = 0; cplx < 2; cplx++) {
      if (!xform[cplx]) continue;
      for (int i = 0; i < nsizes; i++) {
        const int n = sizes[i];
        const size_t len = (size_t)(cplx ? 2 * n : n);
        int any = 0, inputs;
        tally t[BENCH_MAX_ENGINES];
        twiddles tw;
        arrays a;
        for (int l = 0; l < nlibs; l++) {
          const api *f = &apis[l];
          if (f->new_setup && (!f->is_valid_size || f->is_valid_size(n, cplx))) any = 1;
        }
        if (!any) continue;
        tw.tab = NULL;
        if (alloc_arrays(&a, len, bits_only ? 0 : (size_t)n, f64) || (!bits_only && twiddles_init(&tw, (size_t)n))) {
          fprintf(stderr, "bench_fft --accuracy: out of memory at %s %s n=%d\n", prec_names[f64],
                  xform_names[cplx], n);
          free_arrays(&a);
          return failed + 1;
        }
        fprintf(stderr, "%s %s n=%d\n", prec_names[f64], xform_names[cplx], n);
        inputs = (int)(((size_t)seeds * MIN_POINTS + (size_t)n - 1) / (size_t)n); /* at least 1 */
        memset(t, 0, sizeof t);
        for (int k = 0; k < inputs; k++) {
          /* the same inputs for every library, whatever else is measured */
          uint64_t state = ((uint64_t)seed << 32) ^ ((uint64_t)n << 8) ^ (uint64_t)(f64 * 2 + cplx), s;
          state ^= (uint64_t)k * 0xD1B54A32D192ED03u;
          s = next_random(&state);
          random_fill(a.x, len, f64, &s);
          random_fill(a.spec, len, f64, &s);
          if (!bits_only && dd_fft(a.x, cplx ? IN_CPLX : IN_REAL, 0, 1, a.r, (size_t)n, 1, &tw, -1)) {
            fprintf(stderr, "bench_fft --accuracy: n=%d has a prime factor above 5\n", n);
            break;
          }
          for (int l = 0; l < nlibs; l++) {
            const api *f = &apis[l];
            if (f->new_setup && (!f->is_valid_size || f->is_valid_size(n, cplx)))
              measure_forward(&t[l], f, &a, n, cplx, f64, !bits_only);
          }
          if (!bits_only) dd_fft(a.spec, cplx ? IN_CPLX : IN_HERMITIAN, 0, 1, a.r, (size_t)n, 1, &tw, +1);
          for (int l = 0; l < nlibs; l++) {
            const api *f = &apis[l];
            if (f->new_setup && (!f->is_valid_size || f->is_valid_size(n, cplx)))
              measure_backward(&t[l], f, &a, n, cplx, f64, !bits_only);
          }
        }
        for (int l = 0; l < nlibs; l++) {
          if (!t[l].inputs) continue;
          failed += report(csv, libs[l].label, &t[l], n, cplx, f64);
          measured++;
        }
        free(tw.tab);
        free_arrays(&a);
      }
    }
  }
  if (!measured) {
    fprintf(stderr, "bench_fft --accuracy: no library does any of the selected sizes\n");
    return 1;
  }
  return failed;
}
