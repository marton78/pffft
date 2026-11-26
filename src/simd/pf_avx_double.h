/*
   Copyright (c) 2020  Dario Mambro ( dario.mambro@gmail.com )
*/

/* Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )

   Redistribution and use of the Software in source and binary forms,
   with or without modification, is permitted provided that the
   following conditions are met:

   - Neither the names of NCAR's Computational and Information Systems
   Laboratory, the University Corporation for Atmospheric Research,
   nor the names of its sponsors or contributors may be used to
   endorse or promote products derived from this Software without
   specific prior written permission.

   - Redistributions of source code must retain the above copyright
   notices, this list of conditions, and the disclaimer below.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the disclaimer below in the
   documentation and/or other materials provided with the
   distribution.

   THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT
   HOLDERS BE LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
   SOFTWARE.
*/

#ifndef PF_AVX_DBL_H
#define PF_AVX_DBL_H

/*
   vector support macros: the rest of the code is independant of
   AVX -- adding support for other platforms with 4-element
   vectors should be limited to these macros
*/


/*
  AVX support macros
*/
#if !defined(SIMD_SZ) && !defined(PFFFT_SIMD_DISABLE) && defined(__AVX__)
#include <immintrin.h>
typedef __m256d v4sf;

/* 4 doubles by simd vector */
#  define SIMD_SZ 4

typedef union v4sf_union {
  v4sf  v;
  double f[SIMD_SZ];
} v4sf_union;

#  define VARCH "AVX"
#  define VREQUIRES_ALIGN 1
#  define VZERO() _mm256_setzero_pd()
#  define VMUL(a,b) _mm256_mul_pd(a,b)
#  define VADD(a,b) _mm256_add_pd(a,b)
/* see pf_fma_x86.h for why this is not just defined(__FMA__) */
#include "pf_fma_x86.h"
#ifdef PFFFT_X86_HAVE_FMA
#  define VMADD(a,b,c) _mm256_fmadd_pd(a,b,c)
#  define VMSUB(a,b,c) _mm256_fnmadd_pd(a,b,c)
#else
#  define VMADD(a,b,c) _mm256_add_pd(_mm256_mul_pd(a,b), c)
#  define VMSUB(a,b,c) _mm256_sub_pd(c, _mm256_mul_pd(a,b))
#endif
#  define VSUB(a,b) _mm256_sub_pd(a,b)
#  define LD_PS1(p) _mm256_set1_pd(p)
#  define VLOAD_UNALIGNED(ptr)  _mm256_loadu_pd(ptr)
#  define VLOAD_ALIGNED(ptr)    _mm256_load_pd(ptr)

/* INTERLEAVE2 (in1, in2, out1, out2) pseudo code:
out1 = [ in1[0], in2[0], in1[1], in2[1] ]
out2 = [ in1[2], in2[2], in1[3], in2[3] ]
Two shuffles within the 128-bit lanes and two across them. Spelled with
extracts and inserts of the 128-bit halves it took eight, which GCC and MSVC
emit as written (Clang merges some of them). */
#  define INTERLEAVE2(in1, in2, out1, out2) {							\
	__m256d lo__ = _mm256_unpacklo_pd(in1, in2); /* in1[0] in2[0] in1[2] in2[2] */	\
	__m256d hi__ = _mm256_unpackhi_pd(in1, in2); /* in1[1] in2[1] in1[3] in2[3] */	\
	out1 = _mm256_permute2f128_pd(lo__, hi__, 0x20);					\
	out2 = _mm256_permute2f128_pd(lo__, hi__, 0x31);					\
}

/*UNINTERLEAVE2(in1, in2, out1, out2) pseudo code:
out1 = [ in1[0], in1[2], in2[0], in2[2] ]
out2 = [ in1[1], in1[3], in2[1], in2[3] ]
*/
#  define UNINTERLEAVE2(in1, in2, out1, out2) {							\
	__m256d lo__ = _mm256_permute2f128_pd(in1, in2, 0x20); /* in1[0] in1[1] in2[0] in2[1] */	\
	__m256d hi__ = _mm256_permute2f128_pd(in1, in2, 0x31); /* in1[2] in1[3] in2[2] in2[3] */	\
	out1 = _mm256_unpacklo_pd(lo__, hi__);								\
	out2 = _mm256_unpackhi_pd(lo__, hi__);								\
}

#  define VTRANSPOSE4(row0, row1, row2, row3) {				\
        __m256d tmp3, tmp2, tmp1, tmp0;                     \
                                                            \
        tmp0 = _mm256_shuffle_pd((row0),(row1), 0x0);       \
        tmp2 = _mm256_shuffle_pd((row0),(row1), 0xF);       \
        tmp1 = _mm256_shuffle_pd((row2),(row3), 0x0);       \
        tmp3 = _mm256_shuffle_pd((row2),(row3), 0xF);       \
                                                            \
        (row0) = _mm256_permute2f128_pd(tmp0, tmp1, 0x20);	\
        (row1) = _mm256_permute2f128_pd(tmp2, tmp3, 0x20);  \
        (row2) = _mm256_permute2f128_pd(tmp0, tmp1, 0x31);  \
        (row3) = _mm256_permute2f128_pd(tmp2, tmp3, 0x31);  \
    }

/*VSWAPHL(a, b) pseudo code:
return [ b[0], b[1], a[2], a[3] ]
A blend, which is not a shuffle: on Intel it runs on more ports than port 5.
*/
#  define VSWAPHL(a,b) _mm256_blend_pd(a, b, 0x3)

/* a with lane 0 replaced by the scalar f */
#  define VSET_LANE0(a,f) _mm256_blend_pd(a, _mm256_castpd128_pd256(_mm_set_sd(f)), 0x1)

/* The 128-bit halves of a v4sf as pairs of doubles, for the reorder kernels
   (pffft_priv_impl.h): with them INTERLEAVE2 into memory and UNINTERLEAVE2
   from memory need no shuffle across the halves, the work of port 5 on Intel.
   VPAIRS(a, b, lo, hi)  lo = [a0 b0 a2 b2], hi = [a1 b1 a3 b3]
   VSTORE_PAIR_LO(p, v)  [v0 v1] to p, a double pointer aligned to 16 bytes
   VSTORE_PAIR_HI(p, v)  [v2 v3] to p
   VLOAD_PAIRS(p, q)     [p[0] p[1] q[0] q[1]] */
#  define VPAIRS(a, b, lo, hi) { lo = _mm256_unpacklo_pd(a, b); hi = _mm256_unpackhi_pd(a, b); }
#  define VSTORE_PAIR_LO(p, v) _mm_store_pd(p, _mm256_castpd256_pd128(v))
#  define VSTORE_PAIR_HI(p, v) _mm_store_pd(p, _mm256_extractf128_pd(v, 1))
#  define VLOAD_PAIRS(p, q) _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_load_pd(p)), _mm_load_pd(q), 1)

/* reverse/flip all floats */
#  define VREV_S(a)    _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm_permute_pd(_mm256_extractf128_pd(a, 1),1)), _mm_permute_pd(_mm256_castpd256_pd128(a), 1), 1)

/* reverse/flip complex floats */
#  define VREV_C(a)    _mm256_insertf128_pd(_mm256_castpd128_pd256(_mm256_extractf128_pd(a, 1)), _mm256_castpd256_pd128(a), 1)

#  define VALIGNED(ptr) ((((uintptr_t)(ptr)) & 0x1F) == 0)

#endif

#endif /* PF_AVX_DBL_H */

