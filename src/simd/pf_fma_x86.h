/* Copyright (c) 2024  hayguen@web.de
   Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )

   Redistribution and use of the Software in source and binary forms,
   with or without modification, is permitted provided that the
   following conditions are met:

   - Neither the names of PFFFT, nor the names of its
   sponsors or contributors may be used to endorse or promote products
   derived from this Software without specific prior written permission.

   - Redistributions of source code must retain the above copyright
   notices, this list of conditions, and the disclaimer below.

   - Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions, and the disclaimer below in the
   documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
   THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
   DAMAGE.
*/

#ifndef PF_FMA_X86_H
#define PF_FMA_X86_H

/*
  x86 target features, defined in one place for every file that picks an
  SSE code path (pf_sse1_float.h, pf_sse2_double.h, pf_mixer.cpp).

  Does the compiler target SSE, so that the _mm_*_ps intrinsics may be
  used? And SSE2 for the _mm_*_pd ones?

  - x86_64 always has SSE and SSE2. gcc and clang define __x86_64__,
    MSVC defines _M_X64.

  - 32-bit x86 has them only if the compiler targets them: gcc and clang
    define __SSE__ / __SSE2__ (e.g. -msse2, -march=pentium4; not for the
    -march=i686 default of some systems), MSVC sets _M_IX86_FP to 1 for
    /arch:SSE and to 2 for /arch:SSE2 and above.

  - ARM64EC: MSVC defines _M_X64 there too, but the target is ARM, so it
    is excluded.
*/
#if !defined(_M_ARM64EC) && (defined(__x86_64__) || defined(_M_X64) || \
    defined(__SSE__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 1))
#  define PFFFT_X86_HAVE_SSE 1
#endif
#if !defined(_M_ARM64EC) && (defined(__x86_64__) || defined(_M_X64) || \
    defined(__SSE2__) || (defined(_M_IX86_FP) && _M_IX86_FP >= 2))
#  define PFFFT_X86_HAVE_SSE2 1
#endif

/*
  Does this x86 target guarantee FMA3, so that _mm*_fmadd_ps/_pd and
  _mm*_fnmadd_ps/_pd may be emitted?

  Three separate facts are folded in here, which is why it lives in one
  place rather than being repeated in every backend that needs it:

  - gcc and clang define __FMA__ for -mfma, and for -march=haswell and
    later. Note that -mavx2 alone does NOT imply it.

  - MSVC defines no __FMA__ at all; it is absent from the predefined
    macro list. Its /arch:AVX2 (and /arch:AVX512, /arch:AVX10.*) does
    "enable use of Fused Multiply-Add (FMA) instructions", and every CPU
    with AVX2 has FMA3, so __AVX2__ is the equivalent gate there.

  - clang-cl also defines _MSC_VER, but unlike MSVC it rejects the
    intrinsic without the fma target feature, so it must go through
    __FMA__ instead of the MSVC branch.

  https://learn.microsoft.com/en-us/cpp/preprocessor/predefined-macros
  https://learn.microsoft.com/en-us/cpp/build/reference/arch-x64
*/
#if defined(__FMA__) || (defined(_MSC_VER) && !defined(__clang__) && defined(__AVX2__))
#  define PFFFT_X86_HAVE_FMA 1
#endif

#endif /* PF_FMA_X86_H */
