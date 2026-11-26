/* test_accuracy: this build's pffft against an exact FFT, and the
 * equivalences that must hold bit for bit within it (benchmarks/accuracy.c),
 * on short lengths of every kind of pass. A test program of its own rather
 * than bench_fft --accuracy, so that it runs wherever the tests do, without
 * the benchmarks: in CI on every SIMD backend, emulators included, also
 * with --split (test_accuracy_split). Exits 1 if a check fails.
 *
 *   test_accuracy [--split] [N ...]    (default: the lengths below)
 *
 * --split: every setup splits arrays of more than 1 KiB into blocks
 * (pffft_new_setup_tuned()), so that the split transforms run at every
 * length, as deep as they go.
 */
#include "bench.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* what accuracy.c takes from platform.c in bench_fft */
void *bench_aligned_alloc(size_t bytes) {
  unsigned char *raw = (unsigned char *)malloc(bytes + 64 + sizeof(void *)), *p;
  if (!raw) return NULL;
  p = raw + sizeof(void *);
  p += (64 - (uintptr_t)p % 64) % 64;
  ((void **)p)[-1] = raw;
  return p;
}

void bench_aligned_free(void *p) {
  if (p) free(((void **)p)[-1]);
}

void *bench_dl_sym(void *lib, const char *symbol) {
  (void)lib;
  (void)symbol;
  return NULL;
}

int main(int argc, char **argv) {
  /* powers of two, and lengths with radix 3, 5 and 6 passes and the
     twiddle-free radix-8 and radix-6 last passes; those that are not valid
     for the backend are skipped */
  static const int defaults[] = {16,  32,  48,  64,   96,   128,  160,  192,  256,  320,  384,
                                 480, 512, 640, 768, 1024, 1152, 1440, 2048, 2400, 3072, 3200, 4096};
  int sizes[64], nsizes = 0, prec[2] = {0, 0}, xform[2] = {1, 1}, failed;
  bench_tune split;
  bench_accuracy_lib lib = {"pffft", NULL, NULL};
#ifdef BENCH_HAVE_F32
  prec[0] = 1;
#endif
#ifdef BENCH_HAVE_F64
  prec[1] = 1;
#endif
  memset(&split, 0, sizeof split);
  split.set = BENCH_TUNE_SPLIT_REAL | BENCH_TUNE_SPLIT_COMPLEX;
  split.tuning.split_max_real = split.tuning.split_max_complex = 1024;
  for (int i = 1; i < argc && nsizes < 64; i++) {
    if (!strcmp(argv[i], "--split")) lib.tune = &split;
    else sizes[nsizes++] = atoi(argv[i]);
  }
  if (!nsizes)
    for (size_t i = 0; i < sizeof defaults / sizeof defaults[0]; i++) sizes[nsizes++] = defaults[i];
  bench_accuracy_columns(stdout);
  failed = bench_accuracy(stdout, &lib, 1, sizes, nsizes, prec, xform, 2, 1, 0);
  printf("%s\n", failed ? "FAILED" : "all checks passed");
  return failed ? 1 : 0;
}
