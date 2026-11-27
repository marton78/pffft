/*
   What pffft.h (single precision) and pffft_double.h (double precision)
   have in common: the export macros, the enums, the tuning of a setup,
   and the helpers that serve both. Both headers include it. Same license
   as pffft.h.
*/

#ifndef PFFFT_TYPES_H
#define PFFFT_TYPES_H

#include <stddef.h> /* for size_t */

#ifndef PFFFT_EXPORT
  #ifdef PFFFT_STATIC_DEFINE
    #define PFFFT_EXPORT
  #elif defined(_WIN32) || defined(__CYGWIN__)
    #ifdef PFFFT_EXPORTS
      #define PFFFT_EXPORT __declspec(dllexport)
    #else
      #define PFFFT_EXPORT __declspec(dllimport)
    #endif
  #elif defined(PFFFT_EXPORTS) && (defined(__GNUC__) || defined(__clang__))
    #define PFFFT_EXPORT __attribute__((visibility("default")))
  #else
    #define PFFFT_EXPORT
  #endif
#endif

#ifndef PFFFT_DEPRECATED
  #if defined(__GNUC__) || defined(__clang__)
    #define PFFFT_DEPRECATED __attribute__((deprecated))
  #elif defined(_MSC_VER)
    #define PFFFT_DEPRECATED __declspec(deprecated)
  #else
    #define PFFFT_DEPRECATED
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

  /* direction of the transform */
  typedef enum { PFFFT_FORWARD, PFFFT_BACKWARD } pffft_direction_t;

  /* type of transform */
  typedef enum { PFFFT_REAL, PFFFT_COMPLEX } pffft_transform_t;

  /*
    The thresholds of a setup that depend on the CPU. pffft_new_setup()
    takes those that pffft_default_tuning() / pffftd_default_tuning() give
    for the CPU it runs on, pffft_new_setup_tuned() those of the caller.
    Sizes are in bytes of one array of the transform: N values for a real
    transform, 2N for a complex one, of float for pffft_, of double for
    pffftd_. (size_t)-1 means never.
  */
  typedef struct pffft_tuning_t {
    /* the bytes of the struct that both the caller and the library know:
       the default functions set it to the smaller of the size they are
       given and their own, and the setup reads only the fields within it,
       so that a program and a library of different versions understand
       each other */
    size_t size;
    /* real transforms of larger arrays are split into blocks that are
       transformed one after the other while they are in the cache */
    size_t split_max_real;
    /* the same for complex transforms */
    size_t split_max_complex;
    /* x86: the passes of a real transform over arrays of this size or more
       prefetch the lines they write, the first pass of a backward real
       transform from half of it on; without effect on other CPUs, and in
       builds with PFFFT_DISABLE_PREFETCH */
    size_t prefetch_min_real;
    /* the same for complex transforms */
    size_t prefetch_min_complex;
    /* real transforms of arrays up to this size, whose length has a factor
       3, end with a radix-6 pass instead of a radix-8 one. The only field
       that changes the results (in the last bits): all others choose
       between ways that give the same bits. */
    size_t radix6_max;
  } pffft_tuning_t;

  /* simple helper to determine next power of 2
     - without inexact/rounding floating point operations
  */
  PFFFT_EXPORT int pffft_next_power_of_two(int N);

  /* simple helper to determine if power of 2 - returns bool */
  PFFFT_EXPORT int pffft_is_power_of_two(int N);

  /*
    buffers with the alignment that the transforms of both precisions need
    (16 bytes for float with SSE, NEON or Altivec, 32 for double with AVX):
    they are aligned to 64 bytes. pffftd_aligned_malloc() and
    pffftd_aligned_free() in pffft_double.h are the same functions.
  */
  PFFFT_EXPORT void *pffft_aligned_malloc(size_t nb_bytes);
  PFFFT_EXPORT void pffft_aligned_free(void *);

#ifdef __cplusplus
}
#endif

#endif /* PFFFT_TYPES_H */
