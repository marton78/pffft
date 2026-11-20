/* benchmarks/engines.c -- the FFT libraries that bench_fft measures.
 *
 * An engine wraps one library, or one way of calling it, in the interface of
 * bench.h. This file is compiled once per precision, with BENCH_PREC_F32 or
 * with BENCH_PREC_F64, and both objects are linked into one program. So all
 * of it is static, except the two functions at the end, whose names carry
 * the precision, and FFTPACK, whose functions engines_fftpack.h renames.
 *
 * pffft is always there. The build announces every other library with a
 * definition: HAVE_FFTPACK, HAVE_FFTW, HAVE_KISS_FFT, HAVE_POCKET_FFT,
 * HAVE_GREEN_FFTS, HAVE_VECLIB (Apple vDSP), HAVE_MKL. KissFFT and Green exist
 * for float only and drop out of the double build by themselves.
 *
 * What keeps the comparison fair:
 *   - Every engine transforms the arrays of the driver, so that all engines
 *     of a case work in the same memory: out of place from buf[0] to buf[1]
 *     and on to buf[2], with buf[3] as work area, or in place in buf[0].
 *   - run() holds the calls of the library and nothing else: no copying, no
 *     scaling.
 *   - A size that a library cannot do is refused, never rounded up.
 */
#include "bench.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(BENCH_PREC_F32) && !defined(BENCH_PREC_F64)
#include "pffft/pffft.h"
typedef float scalar;
typedef PFFFT_Setup pf_setup;
#define PF(name) pffft_##name
#define PF_SYMBOL(name) "pffft_" name
#define FFTW(name) fftwf_##name
#define POCKET(head, tail) head##ffts##tail
#define VDSP(name) name
#define MKL_PRECISION DFTI_SINGLE
#define BENCH_ENGINES bench_engines_f32
#define BENCH_LIB_ENGINES bench_lib_engines_f32
#elif defined(BENCH_PREC_F64) && !defined(BENCH_PREC_F32)
#include "pffft/pffft_double.h"
typedef double scalar;
typedef PFFFTD_Setup pf_setup;
#define PF(name) pffftd_##name
#define PF_SYMBOL(name) "pffftd_" name
#define FFTW(name) fftw_##name
#define POCKET(head, tail) head##fft##tail
#define VDSP(name) name##D
#define MKL_PRECISION DFTI_DOUBLE
#define BENCH_ENGINES bench_engines_f64
#define BENCH_LIB_ENGINES bench_lib_engines_f64
#else
#error "compile engines.c with exactly one of BENCH_PREC_F32 and BENCH_PREC_F64"
#endif

/* the array in which a pair leaves its result */
enum { IN_PLACE = 0, OUT_OF_PLACE = 2 };

/* The last step of every create(): tells the driver where a pair reads and
 * writes, and by how much it scales. */
static int instance_ready(bench_instance *inst, void *state, const bench_bufs *bufs,
                          int out, int n, int cplx, double gain)
{
  inst->state = state;
  inst->in = bufs->buf[0];
  inst->out = bufs->buf[out];
  inst->len = (size_t)n * (cplx ? 2 : 1);
  inst->gain = gain;
  return 0;
}

#if (defined(HAVE_GREEN_FFTS) && defined(BENCH_PREC_F32)) || defined(HAVE_VECLIB)
/* log2(n) if n is a power of two, otherwise -1 */
static int exact_log2(int n)
{
  int k = 0;
  if (n < 1 || (n & (n - 1)) != 0) return -1;
  while ((1 << k) < n) k++;
  return k;
}
#endif

/* ------------------------------------------------------------------------
 * pffft, linked in from this source tree: pffft_transform_ordered(), which
 * gives the spectrum in the same order as the other libraries.
 * ---------------------------------------------------------------------- */

typedef struct state_pffft {
  pf_setup *setup;
  bench_bufs bufs;
} state_pffft;

static int create_pffft(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                        bench_instance *inst)
{
  const pffft_transform_t kind = cplx ? PFFFT_COMPLEX : PFFFT_REAL;
  state_pffft *s;
  (void)e;
  if (!PF(is_valid_size)(n, kind)) return 1;
  s = malloc(sizeof *s);
  if (!s) return 1;
  s->setup = PF(new_setup)(n, kind);
  if (!s->setup) {
    free(s);
    return 1;
  }
  s->bufs = *bufs;
  return instance_ready(inst, s, bufs, OUT_OF_PLACE, n, cplx, (double)n);
}

static void run_pffft(const bench_instance *inst, int pairs)
{
  const state_pffft *s = inst->state;
  scalar *x = s->bufs.buf[0], *y = s->bufs.buf[1], *z = s->bufs.buf[2], *work = s->bufs.buf[3];
  int i;
  for (i = 0; i < pairs; i++) {
    PF(transform_ordered)(s->setup, x, y, work, PFFFT_FORWARD);
    PF(transform_ordered)(s->setup, y, z, work, PFFFT_BACKWARD);
  }
}

static void destroy_pffft(bench_instance *inst)
{
  state_pffft *s = inst->state;
  PF(destroy_setup)(s->setup);
  free(s);
}

/* ------------------------------------------------------------------------
 * pffft from a shared library given with --lib: "pffft@label", which calls
 * its pffft_transform_ordered(). The library may be of any revision, so no header describes
 * it: its functions are called through pointers that spell setups and arrays
 * as void * and enums as int.
 * ---------------------------------------------------------------------- */

typedef void *(*lib_new_setup_fn)(int n, int transform);
typedef int (*lib_is_valid_size_fn)(int n, int transform);
typedef void (*lib_destroy_setup_fn)(void *setup);
typedef void (*lib_transform_fn)(const void *setup, const void *in, void *out, void *work,
                                 int direction);

typedef struct state_lib {
  void *setup;
  lib_transform_fn transform;
  lib_destroy_setup_fn destroy_setup;
  bench_bufs bufs;
} state_lib;

/* Sets *fn, a function pointer of any type, to a function of the library, or
 * to NULL if the library lacks it. The bits are copied because ISO C does not
 * convert the object pointer of a symbol lookup to a function pointer. */
static void lib_function(void *lib, const char *symbol, void *fn)
{
  void *address = bench_dl_sym(lib, symbol);
  memcpy(fn, &address, sizeof address);
}

static int create_lib(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                      bench_instance *inst)
{
  const int kind = cplx ? 1 : 0; /* PFFFT_REAL is 0, PFFFT_COMPLEX is 1 */
  lib_new_setup_fn new_setup;
  lib_is_valid_size_fn is_valid_size;
  state_lib *s;
  lib_function(e->lib, PF_SYMBOL("new_setup"), &new_setup);
  lib_function(e->lib, PF_SYMBOL("is_valid_size"), &is_valid_size);
  /* revisions without is_valid_size refuse a size in new_setup */
  if (is_valid_size && !is_valid_size(n, kind)) return 1;
  s = malloc(sizeof *s);
  if (!s) return 1;
  s->setup = new_setup(n, kind);
  if (!s->setup) {
    free(s);
    return 1;
  }
  lib_function(e->lib, PF_SYMBOL("transform_ordered"), &s->transform);
  lib_function(e->lib, PF_SYMBOL("destroy_setup"), &s->destroy_setup);
  s->bufs = *bufs;
  return instance_ready(inst, s, bufs, OUT_OF_PLACE, n, cplx, (double)n);
}

static void run_lib(const bench_instance *inst, int pairs)
{
  const state_lib *s = inst->state;
  void *x = s->bufs.buf[0], *y = s->bufs.buf[1], *z = s->bufs.buf[2], *work = s->bufs.buf[3];
  int i;
  for (i = 0; i < pairs; i++) {
    s->transform(s->setup, x, y, work, 0); /* PFFFT_FORWARD */
    s->transform(s->setup, y, z, work, 1); /* PFFFT_BACKWARD */
  }
}

static void destroy_lib(bench_instance *inst)
{
  state_lib *s = inst->state;
  s->destroy_setup(s->setup);
  free(s);
}

/* ------------------------------------------------------------------------
 * FFTW: "fftw-est" plans with FFTW_ESTIMATE, "fftw-meas" with FFTW_MEASURE.
 * The plans are made for the arrays of the driver. Real transforms are r2c
 * and c2r.
 * ---------------------------------------------------------------------- */
#ifdef HAVE_FFTW
#include <fftw3.h>

typedef struct state_fftw {
  FFTW(plan) forward, backward;
} state_fftw;

/* e->variant: the planner flag */
static int create_fftw(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                       bench_instance *inst)
{
  const unsigned flags = (unsigned)e->variant;
  void *x = bufs->buf[0], *y = bufs->buf[1], *z = bufs->buf[2];
  char *wisdom = NULL;
  state_fftw *s;
  if (flags == FFTW_MEASURE && n > bench_fftw_measure_max) return 1;
  s = malloc(sizeof *s);
  if (!s) return 1;
  /* FFTW keeps what FFTW_MEASURE has found out and plans every later
   * FFTW_ESTIMATE of the same problem from it: once "fftw-meas" has been
   * created, "fftw-est" would get the measured plan too. So the estimate is
   * planned with that knowledge set aside. It is put back afterwards, so that
   * "fftw-meas" need not measure again in every trial. */
  if (flags == FFTW_ESTIMATE) {
    wisdom = FFTW(export_wisdom_to_string)();
    FFTW(forget_wisdom)();
  }
  if (cplx) {
    s->forward = FFTW(plan_dft_1d)(n, x, y, FFTW_FORWARD, flags);
    s->backward = FFTW(plan_dft_1d)(n, y, z, FFTW_BACKWARD, flags);
  } else {
    s->forward = FFTW(plan_dft_r2c_1d)(n, x, y, flags);
    s->backward = FFTW(plan_dft_c2r_1d)(n, y, z, flags);
  }
  if (wisdom) {
    FFTW(import_wisdom_from_string)(wisdom);
    free(wisdom);
  }
  if (!s->forward || !s->backward) {
    if (s->forward) FFTW(destroy_plan)(s->forward);
    if (s->backward) FFTW(destroy_plan)(s->backward);
    free(s);
    return 1;
  }
  return instance_ready(inst, s, bufs, OUT_OF_PLACE, n, cplx, (double)n);
}

static void run_fftw(const bench_instance *inst, int pairs)
{
  const state_fftw *s = inst->state;
  int i;
  for (i = 0; i < pairs; i++) {
    FFTW(execute)(s->forward);
    FFTW(execute)(s->backward);
  }
}

static void destroy_fftw(bench_instance *inst)
{
  state_fftw *s = inst->state;
  FFTW(destroy_plan)(s->forward);
  FFTW(destroy_plan)(s->backward);
  free(s);
}
#endif /* HAVE_FFTW */

/* ------------------------------------------------------------------------
 * KissFFT, built for float only. Out of place: given the same array for
 * input and output it would transform into a buffer of its own and copy.
 * ---------------------------------------------------------------------- */
#if defined(HAVE_KISS_FFT) && defined(BENCH_PREC_F32)
#include <kiss_fft.h>
#include <kiss_fftr.h>

typedef struct state_kiss {
  void *forward, *backward; /* kiss_fft_cfg if cplx, else kiss_fftr_cfg */
  int cplx;
  bench_bufs bufs;
} state_kiss;

static int create_kiss(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                       bench_instance *inst)
{
  state_kiss *s;
  (void)e;
  if (!cplx && (n & 1)) return 1; /* the real transform wants an even size */
  s = malloc(sizeof *s);
  if (!s) return 1;
  if (cplx) {
    s->forward = kiss_fft_alloc(n, 0, NULL, NULL);
    s->backward = kiss_fft_alloc(n, 1, NULL, NULL);
  } else {
    s->forward = kiss_fftr_alloc(n, 0, NULL, NULL);
    s->backward = kiss_fftr_alloc(n, 1, NULL, NULL);
  }
  if (!s->forward || !s->backward) {
    kiss_fft_free(s->forward);
    kiss_fft_free(s->backward);
    free(s);
    return 1;
  }
  s->cplx = cplx;
  s->bufs = *bufs;
  return instance_ready(inst, s, bufs, OUT_OF_PLACE, n, cplx, (double)n);
}

static void run_kiss(const bench_instance *inst, int pairs)
{
  const state_kiss *s = inst->state;
  void *x = s->bufs.buf[0], *y = s->bufs.buf[1], *z = s->bufs.buf[2];
  int i;
  if (s->cplx) {
    for (i = 0; i < pairs; i++) {
      kiss_fft(s->forward, x, y);
      kiss_fft(s->backward, y, z);
    }
  } else {
    for (i = 0; i < pairs; i++) {
      kiss_fftr(s->forward, x, y);
      kiss_fftri(s->backward, y, z);
    }
  }
}

static void destroy_kiss(bench_instance *inst)
{
  state_kiss *s = inst->state;
  kiss_fft_free(s->forward);
  kiss_fft_free(s->backward);
  free(s);
}
#endif /* HAVE_KISS_FFT */

/* ------------------------------------------------------------------------
 * PocketFFT: in place, with scale factor 1 in both directions.
 * ---------------------------------------------------------------------- */
#ifdef HAVE_POCKET_FFT
#ifdef BENCH_PREC_F32
#include <pocketfft_single.h>
#else
#include <pocketfft_double.h>
#endif

typedef struct state_pocket {
  POCKET(r, _plan) rplan; /* one of the two, by the kind of transform */
  POCKET(c, _plan) cplan;
} state_pocket;

static int create_pocket(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                         bench_instance *inst)
{
  state_pocket *s;
  (void)e;
  s = calloc(1, sizeof *s);
  if (!s) return 1;
  if (cplx)
    s->cplan = POCKET(make_c, _plan)((size_t)n);
  else
    s->rplan = POCKET(make_r, _plan)((size_t)n);
  if (!s->cplan && !s->rplan) {
    free(s);
    return 1;
  }
  return instance_ready(inst, s, bufs, IN_PLACE, n, cplx, (double)n);
}

static void run_pocket(const bench_instance *inst, int pairs)
{
  const state_pocket *s = inst->state;
  scalar *x = inst->in;
  int i;
  if (s->cplan) {
    for (i = 0; i < pairs; i++) {
      (void)POCKET(c, _forward)(s->cplan, x, 1);
      (void)POCKET(c, _backward)(s->cplan, x, 1);
    }
  } else {
    for (i = 0; i < pairs; i++) {
      (void)POCKET(r, _forward)(s->rplan, x, 1);
      (void)POCKET(r, _backward)(s->rplan, x, 1);
    }
  }
}

static void destroy_pocket(bench_instance *inst)
{
  state_pocket *s = inst->state;
  if (s->cplan) POCKET(destroy_c, _plan)(s->cplan);
  if (s->rplan) POCKET(destroy_r, _plan)(s->rplan);
  free(s);
}
#endif /* HAVE_POCKET_FFT */

/* ------------------------------------------------------------------------
 * Green's FFT: float only, powers of two only, in place.
 * ---------------------------------------------------------------------- */
#if defined(HAVE_GREEN_FFTS) && defined(BENCH_PREC_F32)
#include "fftext.h"

/* Green keeps its tables inside the library, and fftFree() drops all of them:
 * it may only be called when the last instance goes. */
static int green_instances;

typedef struct state_green {
  long log2n;
  int cplx;
} state_green;

static int create_green(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                        bench_instance *inst)
{
  const int log2n = exact_log2(n);
  state_green *s;
  (void)e;
  if (log2n < 1) return 1;
  s = malloc(sizeof *s);
  if (!s || fftInit(log2n) != 0) {
    free(s);
    return 1;
  }
  green_instances++;
  s->log2n = log2n;
  s->cplx = cplx;
  return instance_ready(inst, s, bufs, IN_PLACE, n, cplx, 1.0);
}

static void run_green(const bench_instance *inst, int pairs)
{
  const state_green *s = inst->state;
  scalar *x = inst->in;
  int i;
  if (s->cplx) {
    for (i = 0; i < pairs; i++) {
      ffts(x, s->log2n, 1);
      iffts(x, s->log2n, 1);
    }
  } else {
    for (i = 0; i < pairs; i++) {
      rffts(x, s->log2n, 1);
      riffts(x, s->log2n, 1);
    }
  }
}

static void destroy_green(bench_instance *inst)
{
  free(inst->state);
  if (--green_instances == 0) fftFree();
}
#endif /* HAVE_GREEN_FFTS */

/* ------------------------------------------------------------------------
 * Apple vDSP (Accelerate): powers of two only, in place. vDSP wants real and
 * imaginary parts in separate arrays: they are the two halves of buf[0].
 * ---------------------------------------------------------------------- */
#ifdef HAVE_VECLIB
#include <Accelerate/Accelerate.h>

#ifdef BENCH_PREC_F32
typedef DSPSplitComplex vdsp_split;
#else
typedef DSPDoubleSplitComplex vdsp_split;
#endif

typedef struct state_vdsp {
  VDSP(FFTSetup) setup;
  vdsp_split halves;
  vDSP_Length log2n;
  int cplx;
} state_vdsp;

static int create_vdsp(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                       bench_instance *inst)
{
  const int log2n = exact_log2(n);
  scalar *x = bufs->buf[0];
  state_vdsp *s;
  (void)e;
  if (log2n < 1) return 1;
  s = malloc(sizeof *s);
  if (!s) return 1;
  s->setup = VDSP(vDSP_create_fftsetup)((vDSP_Length)log2n, kFFTRadix2);
  if (!s->setup) {
    free(s);
    return 1;
  }
  s->halves.realp = x;
  s->halves.imagp = x + (cplx ? n : n / 2);
  s->log2n = (vDSP_Length)log2n;
  s->cplx = cplx;
  /* the real pair scales by 2n: vDSP's real forward transform is twice the
   * mathematical one */
  return instance_ready(inst, s, bufs, IN_PLACE, n, cplx, cplx ? (double)n : 2.0 * n);
}

static void run_vdsp(const bench_instance *inst, int pairs)
{
  const state_vdsp *s = inst->state;
  int i;
  if (s->cplx) {
    for (i = 0; i < pairs; i++) {
      VDSP(vDSP_fft_zip)(s->setup, &s->halves, 1, s->log2n, kFFTDirection_Forward);
      VDSP(vDSP_fft_zip)(s->setup, &s->halves, 1, s->log2n, kFFTDirection_Inverse);
    }
  } else {
    for (i = 0; i < pairs; i++) {
      VDSP(vDSP_fft_zrip)(s->setup, &s->halves, 1, s->log2n, kFFTDirection_Forward);
      VDSP(vDSP_fft_zrip)(s->setup, &s->halves, 1, s->log2n, kFFTDirection_Inverse);
    }
  }
}

static void destroy_vdsp(bench_instance *inst)
{
  state_vdsp *s = inst->state;
  VDSP(vDSP_destroy_fftsetup)(s->setup);
  free(s);
}
#endif /* HAVE_VECLIB */

/* ------------------------------------------------------------------------
 * Intel MKL (DFTI), out of place. The real spectrum is in MKL's default CCS
 * layout, n + 2 scalars.
 * ---------------------------------------------------------------------- */
#ifdef HAVE_MKL
#include <mkl/mkl_dfti.h>

typedef struct state_mkl {
  DFTI_DESCRIPTOR_HANDLE handle;
  bench_bufs bufs;
} state_mkl;

static int create_mkl(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                      bench_instance *inst)
{
  state_mkl *s;
  MKL_LONG status;
  (void)e;
  s = malloc(sizeof *s);
  if (!s) return 1;
  /* the size travels through a variable argument list, where MKL reads a
   * MKL_LONG: the cast is needed */
  status = DftiCreateDescriptor(&s->handle, MKL_PRECISION, cplx ? DFTI_COMPLEX : DFTI_REAL, 1,
                                (MKL_LONG)n);
  if (status != 0) {
    free(s);
    return 1;
  }
  status = DftiSetValue(s->handle, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
  if (status == 0) status = DftiCommitDescriptor(s->handle);
  if (status != 0) {
    DftiFreeDescriptor(&s->handle);
    free(s);
    return 1;
  }
  s->bufs = *bufs;
  return instance_ready(inst, s, bufs, OUT_OF_PLACE, n, cplx, (double)n);
}

static void run_mkl(const bench_instance *inst, int pairs)
{
  const state_mkl *s = inst->state;
  void *x = s->bufs.buf[0], *y = s->bufs.buf[1], *z = s->bufs.buf[2];
  int i;
  for (i = 0; i < pairs; i++) {
    DftiComputeForward(s->handle, x, y);
    DftiComputeBackward(s->handle, y, z);
  }
}

static void destroy_mkl(bench_instance *inst)
{
  state_mkl *s = inst->state;
  DftiFreeDescriptor(&s->handle);
  free(s);
}
#endif /* HAVE_MKL */

/* ------------------------------------------------------------------------
 * FFTPACK, the scalar code that pffft descends from: in place. This section
 * comes last because its header compiles fftpack.c right here, and that
 * brings short names such as real and integer with it.
 * ---------------------------------------------------------------------- */
#ifdef HAVE_FFTPACK
#include "engines_fftpack.h"

typedef struct state_fftpack {
  scalar *wsave;
  int n, cplx;
} state_fftpack;

static int create_fftpack(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                          bench_instance *inst)
{
  state_fftpack *s;
  (void)e;
  s = malloc(sizeof *s);
  if (!s) return 1;
  /* FFTPACK asks for 2n+15 scalars, 4n+15 if complex, and keeps the factors
   * of n in the last 15. Sizes with more than 13 factors, the smallest being
   * 2*3^13, overrun them: the 16 on top are for those. */
  s->wsave = malloc(((size_t)n * (cplx ? 4 : 2) + 15 + 16) * sizeof(scalar));
  if (!s->wsave) {
    free(s);
    return 1;
  }
  if (cplx)
    cffti(n, s->wsave);
  else
    rffti(n, s->wsave);
  s->n = n;
  s->cplx = cplx;
  return instance_ready(inst, s, bufs, IN_PLACE, n, cplx, (double)n);
}

static void run_fftpack(const bench_instance *inst, int pairs)
{
  const state_fftpack *s = inst->state;
  scalar *x = inst->in;
  int i;
  if (s->cplx) {
    for (i = 0; i < pairs; i++) {
      cfftf(s->n, x, s->wsave);
      cfftb(s->n, x, s->wsave);
    }
  } else {
    for (i = 0; i < pairs; i++) {
      rfftf(s->n, x, s->wsave);
      rfftb(s->n, x, s->wsave);
    }
  }
}

static void destroy_fftpack(bench_instance *inst)
{
  state_fftpack *s = inst->state;
  free(s->wsave);
  free(s);
}
#endif /* HAVE_FFTPACK */

/* ------------------------------------------------------------------------
 * The lists
 * ---------------------------------------------------------------------- */

int BENCH_ENGINES(bench_engine *list, int max)
{
  static const bench_engine all[] = {
    {"pffft", create_pffft, run_pffft, destroy_pffft, NULL, 0},
#ifdef HAVE_FFTPACK
    {"fftpack", create_fftpack, run_fftpack, destroy_fftpack, NULL, 0},
#endif
#ifdef HAVE_FFTW
    {"fftw-est", create_fftw, run_fftw, destroy_fftw, NULL, FFTW_ESTIMATE},
    {"fftw-meas", create_fftw, run_fftw, destroy_fftw, NULL, FFTW_MEASURE},
#endif
#if defined(HAVE_KISS_FFT) && defined(BENCH_PREC_F32)
    {"kissfft", create_kiss, run_kiss, destroy_kiss, NULL, 0},
#endif
#ifdef HAVE_POCKET_FFT
    {"pocketfft", create_pocket, run_pocket, destroy_pocket, NULL, 0},
#endif
#if defined(HAVE_GREEN_FFTS) && defined(BENCH_PREC_F32)
    {"green", create_green, run_green, destroy_green, NULL, 0},
#endif
#ifdef HAVE_VECLIB
    {"vdsp", create_vdsp, run_vdsp, destroy_vdsp, NULL, 0},
#endif
#ifdef HAVE_MKL
    {"mkl", create_mkl, run_mkl, destroy_mkl, NULL, 0},
#endif
  };
  int count = (int)(sizeof all / sizeof all[0]), i;
  if (count > max) count = max;
  for (i = 0; i < count; i++) list[i] = all[i];
  return count;
}

int BENCH_LIB_ENGINES(const char *label, const char *path, bench_engine *list, int max,
                      char *err, size_t errlen)
{
  static const char *const needed[3] = {PF_SYMBOL("new_setup"), PF_SYMBOL("destroy_setup"),
                                        PF_SYMBOL("transform_ordered")};
  void *lib;
  int i;
  if (max < 1) {
    snprintf(err, errlen, "no room for the engines of %s", path);
    return -1;
  }
  lib = bench_dl_open(path, err, errlen);
  if (!lib) return -1;
  for (i = 0; i < 3; i++)
    if (!bench_dl_sym(lib, needed[i])) return 0; /* built without this precision */
  snprintf(list->name, sizeof list->name, "pffft@%s", label);
  list->create = create_lib;
  list->run = run_lib;
  list->destroy = destroy_lib;
  list->lib = lib;
  list->variant = 0;
  return 1;
}
