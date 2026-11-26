/* benchmarks/bench.h -- internal interfaces of bench_fft, the pffft FFT
 * benchmark.
 *
 * Three translation units share this header:
 *
 *   bench_fft.c  the driver: command line, measurement loop, statistics, output
 *   platform.c   everything that depends on the operating system: clocks, CPU
 *                cycle counters, core selection, shared-library loading
 *   engines.c    one adapter per FFT library; compiled once per precision
 *
 * How a measurement works is described in benchmarks/README.md.
 */
#ifndef PFFFT_BENCH_H
#define PFFFT_BENCH_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------
 * platform.c
 * ---------------------------------------------------------------------- */

/* A reading of the clocks of the calling thread. Only differences between
 * two readings are meaningful. */
typedef struct bench_snap {
  uint64_t ns;           /* monotonic wall clock, nanoseconds */
  uint64_t cycles;       /* CPU cycles spent on the preferred kind of core */
  uint64_t stray_cycles; /* CPU cycles spent on any other kind of core */
  uint64_t instructions; /* instructions retired */
} bench_snap;

/* pin_cpu argument of bench_platform_init() */
#define BENCH_PIN_AUTO (-1) /* choose the fastest core where pinning is possible */
#define BENCH_PIN_NONE (-2) /* do not pin the thread to a core */

/* Prepare the calling thread for measuring: steer it to the fastest cores and
 * open the CPU cycle counter. Returns 1 if cycles are counted, 0 if only the
 * wall clock is available. Call once, from the thread that measures. */
int bench_platform_init(int pin_cpu);

/* Read the clocks. bench_snap_begin() goes immediately before the timed code
 * and bench_snap_end() immediately after it; the two differ only in the order
 * in which they read, so that the cycle window is as tight as possible.
 * Without a cycle counter the cycle and instruction fields stay 0. */
void bench_snap_begin(bench_snap *s);
void bench_snap_end(bench_snap *s);

/* Write "<prefix>key: value\n" lines that describe the machine: os, arch, cpu,
 * counter, pinned_cpu, l1d_bytes, l2_bytes, l3_bytes (those that are known). */
void bench_platform_describe(FILE *f, const char *prefix);

/* Memory aligned to 64 bytes. */
void *bench_aligned_alloc(size_t bytes);
void bench_aligned_free(void *p);

/* Shared libraries. bench_dl_open() returns NULL on failure and explains why
 * in err. The library's symbols stay private to the returned handle. */
void *bench_dl_open(const char *path, char *err, size_t errlen);
void *bench_dl_sym(void *lib, const char *symbol);

/* ------------------------------------------------------------------------
 * engines.c
 * ---------------------------------------------------------------------- */

/* The work arrays of one case (precision, transform kind, size). All engines
 * of the case share them, so every engine sees the same memory. They hold
 * zeros while timing: a transform of zeros stays zero, which keeps in-place
 * engines from overflowing and lets the engines share the arrays. */
typedef struct bench_bufs {
  void *buf[4];  /* each 64-byte aligned */
  size_t bytes;  /* size of each array: room for at least 2*n + 64 scalars */
} bench_bufs;

/* One engine set up for one case. */
typedef struct bench_instance {
  void *state;  /* private to the engine */
  void *in;     /* where a pair reads its input: one of bench_bufs.buf[] */
  void *out;    /* where a pair leaves its result; equals in for in-place */
  size_t len;   /* scalars of signal in in and out: n if real, 2*n if complex */
  double gain;  /* after one pair, out[i] == gain * (in[i] before the pair) */
} bench_instance;

/* An FFT library, or one way of calling it. A "pair" is one forward transform
 * followed by one backward transform of its result. */
typedef struct bench_engine bench_engine;
struct bench_engine {
  char name[48]; /* "pffft", "fftw-est", "pffft@master", ... */

  /* Set up for n >= 1 points (n real samples, or n complex samples if cplx).
   * Returns 0 on success and nonzero if the engine cannot do this case; it
   * then leaves nothing allocated. May scribble on the arrays: the driver
   * zero-fills them again before timing. */
  int (*create)(const bench_engine *e, int n, int cplx, const bench_bufs *bufs,
                bench_instance *inst);

  /* Run `pairs` pairs. This is the timed code: no allocation, no checks. */
  void (*run)(const bench_instance *inst, int pairs);

  void (*destroy)(bench_instance *inst);

  void *lib;   /* engines from --lib: the loaded library; NULL for built-ins */
  int variant; /* private to the engine */
  int on_request; /* runs only when --engines names it, not by default */
};

#define BENCH_MAX_ENGINES 32

/* The engines compiled into this program, for one precision. Appends at most
 * `max` entries to list and returns how many. Which of the two functions
 * exist is announced by BENCH_HAVE_F32 and BENCH_HAVE_F64, set by the build. */
int bench_engines_f32(bench_engine *list, int max);
int bench_engines_f64(bench_engine *list, int max);

/* The engine "pffft@<label>" of a pffft shared library. Returns how many were
 * appended: 1, 0 if the library was built without this precision, -1 on error
 * with the reason in err. */
int bench_lib_engines_f32(const char *label, const char *path,
                          bench_engine *list, int max, char *err, size_t errlen);
int bench_lib_engines_f64(const char *label, const char *path,
                          bench_engine *list, int max, char *err, size_t errlen);

/* ------------------------------------------------------------------------
 * accuracy.c: bench_fft --accuracy
 * ---------------------------------------------------------------------- */

/* A pffft to check: lib is a library from bench_dl_open(), or NULL for the
 * pffft linked into this program. */
typedef struct bench_accuracy_lib {
  const char *label;
  void *lib;
} bench_accuracy_lib;

/* The column line of the CSV that bench_accuracy() writes. */
void bench_accuracy_columns(FILE *csv);

/* Measures every library at every case selected by prec[2] (f32, f64),
 * xform[2] (real, cplx) and sizes, with inputs derived from seed that make
 * at least seeds * 4096 points per case; writes a CSV row per library and
 * case to csv (if not NULL) and the failures to stderr. With bits_only, no
 * errors against the exact transform: only the checks and the hashes.
 * Returns the number of failed checks, or nonzero if nothing could be
 * measured. accuracy.c describes the measures. */
int bench_accuracy(FILE *csv, const bench_accuracy_lib *libs, int nlibs, const int *sizes, int nsizes,
                   const int *prec, const int *xform, int seeds, unsigned long seed, int bits_only);

/* Set by --fftw-measure-max: the largest size that the "fftw-meas" engine
 * plans with FFTW_MEASURE; above it the engine reports the case as
 * unsupported. Defined in bench_fft.c. */
extern int bench_fftw_measure_max;

#ifdef __cplusplus
}
#endif

#endif /* PFFFT_BENCH_H */
