/* benchmarks/engines_ooura.h -- Ooura's FFT as a part of engines.c.
 *
 * The submodule ooura (https://github.com/hayguen/OouraFFT, a fork of Takuya
 * Ooura's package that also builds for float) exports the same function names
 * for either precision, while bench_fft links both precisions into one
 * program. So, as with FFTPACK, the functions get a name per precision here,
 * and the source is compiled as a part of engines.c. Its other functions are
 * static.
 *
 * fftsg.c is the split-radix version, the fastest of the package.
 */
#ifndef PFFFT_BENCH_ENGINES_OOURA_H
#define PFFFT_BENCH_ENGINES_OOURA_H

#include <math.h> /* before the renaming below, which must not reach into it */

#undef OOURA_SINGLE_PREC
#if defined(BENCH_PREC_F32)
#define OOURA_SINGLE_PREC
#define BENCH_OOURA_NAME(function) bench_ooura_f32_##function
#else
#define BENCH_OOURA_NAME(function) bench_ooura_f64_##function
#endif

/* every function that fftsg.c exports, used by the benchmark or not */
#define cdft BENCH_OOURA_NAME(cdft)
#define rdft BENCH_OOURA_NAME(rdft)
#define ddct BENCH_OOURA_NAME(ddct)
#define ddst BENCH_OOURA_NAME(ddst)
#define dfct BENCH_OOURA_NAME(dfct)
#define dfst BENCH_OOURA_NAME(dfst)
#define ooura_prec BENCH_OOURA_NAME(ooura_prec)

#include "../ooura/fftsg.c"

#endif /* PFFFT_BENCH_ENGINES_OOURA_H */
