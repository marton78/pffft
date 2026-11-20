/* benchmarks/engines_fftpack.h -- FFTPACK as a part of engines.c.
 *
 * src/fftpack.c builds for float or for double and exports the same function
 * names either way, while bench_fft links both precisions into one program.
 * So the functions get a name per precision here, and the source is compiled
 * as a part of engines.c instead of being linked as a library.
 *
 * Include this after all other headers: fftpack.c defines short names of its
 * own, such as the types real and integer.
 */
#ifndef PFFFT_BENCH_ENGINES_FFTPACK_H
#define PFFFT_BENCH_ENGINES_FFTPACK_H

#include <math.h> /* before the renaming below, which must not reach into it */

#undef FFTPACK_DOUBLE_PRECISION
#if defined(BENCH_PREC_F64)
#define FFTPACK_DOUBLE_PRECISION
#define BENCH_FFTPACK_NAME(function) bench_fftpack_f64_##function
#else
#define BENCH_FFTPACK_NAME(function) bench_fftpack_f32_##function
#endif

/* every function that fftpack.c exports, used by the benchmark or not */
#define cffti BENCH_FFTPACK_NAME(cffti)
#define cfftf BENCH_FFTPACK_NAME(cfftf)
#define cfftb BENCH_FFTPACK_NAME(cfftb)
#define rffti BENCH_FFTPACK_NAME(rffti)
#define rfftf BENCH_FFTPACK_NAME(rfftf)
#define rfftb BENCH_FFTPACK_NAME(rfftb)
#define cosqi BENCH_FFTPACK_NAME(cosqi)
#define cosqf BENCH_FFTPACK_NAME(cosqf)
#define cosqb BENCH_FFTPACK_NAME(cosqb)
#define costi BENCH_FFTPACK_NAME(costi)
#define cost BENCH_FFTPACK_NAME(cost)
#define sinqi BENCH_FFTPACK_NAME(sinqi)
#define sinqb BENCH_FFTPACK_NAME(sinqb)
#define sinqf BENCH_FFTPACK_NAME(sinqf)
#define sinti BENCH_FFTPACK_NAME(sinti)
#define sint BENCH_FFTPACK_NAME(sint)

#include "../src/fftpack.c"

#endif /* PFFFT_BENCH_ENGINES_FFTPACK_H */
