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
