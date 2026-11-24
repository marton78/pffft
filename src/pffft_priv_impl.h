/* Copyright (c) 2013  Julien Pommier ( pommier@modartt.com )
   Copyright (c) 2020  Hayati Ayguen ( h_ayguen@web.de )
   Copyright (c) 2020  Dario Mambro ( dario.mambro@gmail.com )

   Based on original fortran 77 code from FFTPACKv4 from NETLIB
   (http://www.netlib.org/fftpack), authored by Dr Paul Swarztrauber
   of NCAR, in 1985.

   As confirmed by the NCAR fftpack software curators, the following
   FFTPACKv5 license applies to FFTPACKv4 sources. My changes are
   released under the same terms.

   FFTPACK license:

   http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html

   Copyright (c) 2004 the University Corporation for Atmospheric
   Research ("UCAR"). All rights reserved. Developed by NCAR's
   Computational and Information Systems Laboratory, UCAR,
   www.cisl.ucar.edu.

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


   PFFFT : a Pretty Fast FFT.

   This file is largerly based on the original FFTPACK implementation, modified in
   order to take advantage of SIMD instructions of modern CPUs.
*/

/* this file requires architecture specific preprocessor definitions
 * it's only for library internal use
 */


/* define own constants required to turn off g++ extensions .. */
#ifndef M_PI
  #define M_PI    3.14159265358979323846  /* pi */
#endif

#ifndef M_SQRT2
  #define M_SQRT2 1.41421356237309504880  /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
  #define M_SQRT1_2 0.707106781186547524400844362104849039  /* 1/sqrt(2) */
#endif


int FUNC_SIMD_SIZE(void) { return SIMD_SZ; }

int FUNC_MIN_FFT_SIZE(pffft_transform_t transform) {
  /* unfortunately, the fft size must be a multiple of 16 for complex FFTs
     and 32 for real FFTs -- a lot of stuff would need to be rewritten to
     handle other cases (or maybe just switch to a scalar fft, I don't know..) */
  int simdSz = FUNC_SIMD_SIZE();
  if (transform == PFFFT_REAL)
    return ( 2 * simdSz * simdSz );
  else if (transform == PFFFT_COMPLEX)
    return ( simdSz * simdSz );
  else
    return 1;
}

int FUNC_IS_VALID_SIZE(int N, pffft_transform_t cplx) {
  const int N_min = FUNC_MIN_FFT_SIZE(cplx);
  int R = N;
  while (R >= 5*N_min && (R % 5) == 0)  R /= 5;
  while (R >= 3*N_min && (R % 3) == 0)  R /= 3;
  while (R >= 2*N_min && (R % 2) == 0)  R /= 2;
  return (R == N_min) ? 1 : 0;
}

int FUNC_NEAREST_SIZE(int N, pffft_transform_t cplx, int higher) {
  int d;
  const int N_min = FUNC_MIN_FFT_SIZE(cplx);
  if (N < N_min)
    N = N_min;
  d = (higher) ? N_min : -N_min;
  if (d > 0)
    N = N_min * ((N+N_min-1) / N_min);  /* round up */
  else
    N = N_min * (N / N_min);  /* round down */

  for (; ; N += d)
    if (FUNC_IS_VALID_SIZE(N, cplx))
      return N;
}

const char * FUNC_SIMD_ARCH(void) { return VARCH; }


/*
  a + fsign*b and a - fsign*b, for fsign == -1 or +1, as one VSUB or VADD.
  In the always-inline helpers below fsign is a constant, so the choice is
  made at compile time, with any compiler (MSVC, for one, does not fold a
  multiply by a constant -1 into an addition).
*/
#define VADD_FSIGN(a, b) (fsign < 0 ? VSUB(a, b) : VADD(a, b))
#define VSUB_FSIGN(a, b) (fsign < 0 ? VADD(a, b) : VSUB(a, b))

/*
  Prefetching the output of the radix passes (x86).

  Out of L1 these passes are bound by their stores: each iteration writes to
  2 to 8 output streams. On a Xeon 6975P-C, perf counted 40-50% of the cycles
  of the real transform stalled on a full store buffer, with no more L1/L2
  traffic than without the prefetches; the lines of the store misses
  apparently do not arrive fast enough. Prefetching each output stream
  PFFFT_PREFETCH_DIST bytes ahead made the transforms out of L1 up to 40%
  faster on Intel CPUs (Xeon 6975P-C, Xeon Platinum 8573C, Core i7-8700B);
  docs/benchmarks.md has the numbers. A read prefetch measured as good as one
  for writing (PREFETCHW), and every x86 with SSE has it.

  - Passes over less than PFFFT_PREFETCH_MIN bytes mostly stay in L1, where
    the prefetches cost 0.5-5%. All passes of a transform write the same
    number of vectors, so the transform picks the kernels once (PF_CALL):
    f_pf_ps with the prefetches or f_ps without, both from one f_loop_ps
    that takes pf as a constant. Deciding in each kernel instead cost up to
    4% at small n.
  - Not in passf5_ps, where they cost more than they gave (see there).
  - AMD Zen cores did not gain (up to 6% slower with gcc and 14% with MSVC
    at single lengths), so the default is off when the compiler targets Zen:
    -march=znverN, which -march=native picks on compilers that know the CPU.
    Builds for generic x86, for other CPUs or with MSVC cannot tell: the
    README asks to set PFFFT_PREFETCH_DIST to 0 there for AMD CPUs.
  - Never on WebAssembly, where Emscripten's _mm_prefetch does nothing.
  The x86 SIMD headers have already included pf_fma_x86.h and _mm_prefetch.
*/
#ifndef PFFFT_PREFETCH_DIST
   /* extend this list when compilers add later Zen targets */
#  if defined(__znver1__) || defined(__znver2__) || defined(__znver3__) || \
      defined(__znver4__) || defined(__znver5__) || defined(__znver6__)
#    define PFFFT_PREFETCH_DIST 0
#  else
#    define PFFFT_PREFETCH_DIST 128
#  endif
#endif
#ifndef PFFFT_PREFETCH_MIN
#  define PFFFT_PREFETCH_MIN 24576
#endif
#if PFFFT_PREFETCH_DIST > 0 && defined(PFFFT_X86_HAVE_SSE) && !defined(PFFFT_SIMD_DISABLE) && \
    !defined(__EMSCRIPTEN__) && !defined(__wasm__)
   /* Prefetch the line PFFFT_PREFETCH_DIST bytes after (dir = 1) or before
      (dir = -1) p if pf. Near the ends of a pass that address lies outside
      the array: a prefetch never faults, and pointer arithmetic (rather than
      integer arithmetic, which cost clang an instruction per prefetch and up
      to 5% on a Core i7-8700B) lets the compiler fold the offset into the
      instruction. */
#  define PF_PREFETCH(pf, p, dir) do { if (pf) \
     _mm_prefetch((const char *)(p) + (dir) * PFFFT_PREFETCH_DIST, _MM_HINT_T0); } while (0)
   /* whether the passes over nvec vectors (every pass of a transform writes
      all of them) take the kernels with the prefetches */
#  define PF_USE_PREFETCH(nvec) ((size_t)(nvec) * sizeof(v4sf) >= (size_t)PFFFT_PREFETCH_MIN)
   /* call the kernel f##_pf_ps, with the prefetches, if pf, else f##_ps */
#  define PF_CALL(pf, f, args) do { if (pf) f##_pf_ps args; else f##_ps args; } while (0)
#  define PF_ON 1
#else
#  define PF_PREFETCH(pf, p, dir) ((void)(pf))
#  define PF_USE_PREFETCH(nvec) 0
#  define PF_CALL(pf, f, args) f##_ps args
#  define PF_ON 0
#endif

/*
  passf2 and passb2 has been merged here, fsign = -1 for passf2, +1 for passb2
*/
static NEVER_INLINE(void) passf2_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float *wa1, float fsign) {
  int k, i;
  int l1ido = l1*ido;
  if (ido <= 2) {
    for (k=0; k < l1ido; k += ido, ch += ido, cc+= 2*ido) {
      ch[0]         = VADD(cc[0], cc[ido+0]);
      ch[l1ido]     = VSUB(cc[0], cc[ido+0]);
      ch[1]         = VADD(cc[1], cc[ido+1]);
      ch[l1ido + 1] = VSUB(cc[1], cc[ido+1]);
    }
  } else {
    for (k=0; k < l1ido; k += ido, ch += ido, cc += 2*ido) {
      for (i=0; i<ido-1; i+=2) {
        v4sf tr2 = VSUB(cc[i+0], cc[i+ido+0]);
        v4sf ti2 = VSUB(cc[i+1], cc[i+ido+1]);
        v4sf wr = LD_PS1(wa1[i]), wi = VMUL(LD_PS1(fsign), LD_PS1(wa1[i+1]));
        ch[i]   = VADD(cc[i+0], cc[i+ido+0]);
        ch[i+1] = VADD(cc[i+1], cc[i+ido+1]);
        VCPLXMUL(tr2, ti2, wr, wi);
        ch[i+l1ido]   = tr2;
        ch[i+l1ido+1] = ti2;
      }
    }
  }
}

/*
  passf3 and passb3 has been merged here, fsign = -1 for passf3, +1 for passb3

  The radix-3/5 twiddle constants are written as full-precision double literals
  without an 'f' suffix, for when this file is compiled with '#define float double'.
  In the float build, the values round to the same floats as the old suffixed literals.
*/
static ALWAYS_INLINE(void) passf3_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                          const float *wa1, const float *wa2, const float fsign,
                                          const int pf) {
  static const float taur = (float)-0.5;                     /* cos(2*pi/3)   */
  float taui = (float)0.86602540378443864676*fsign;          /* sin(2*pi/3) = sqrt(3)/2 */
  int i, k;
  v4sf tr2, ti2, cr2, ci2, cr3, ci3, dr2, di2, dr3, di3;
  int l1ido = l1*ido;
  float wr1, wi1, wr2, wi2;
  assert(ido > 2);
  for (k=0; k< l1ido; k += ido, cc+= 3*ido, ch +=ido) {
    for (i=0; i<ido-1; i+=2) {
      PF_PREFETCH(pf, ch + i, 1);
      PF_PREFETCH(pf, ch + i + l1ido, 1);
      PF_PREFETCH(pf, ch + i + 2*l1ido, 1);
      tr2 = VADD(cc[i+ido], cc[i+2*ido]);
      cr2 = SVMADD(taur, tr2, cc[i]);
      ch[i]    = VADD(cc[i], tr2);
      ti2 = VADD(cc[i+ido+1], cc[i+2*ido+1]);
      ci2 = SVMADD(taur, ti2, cc[i    +1]);
      ch[i+1]  = VADD(cc[i+1], ti2);
      cr3 = SVMUL(taui, VSUB(cc[i+ido], cc[i+2*ido]));
      ci3 = SVMUL(taui, VSUB(cc[i+ido+1], cc[i+2*ido+1]));
      dr2 = VSUB(cr2, ci3);
      dr3 = VADD(cr2, ci3);
      di2 = VADD(ci2, cr3);
      di3 = VSUB(ci2, cr3);
      wr1=wa1[i]; wi1=fsign*wa1[i+1]; wr2=wa2[i]; wi2=fsign*wa2[i+1];
      VCPLXMUL(dr2, di2, LD_PS1(wr1), LD_PS1(wi1));
      ch[i+l1ido] = dr2;
      ch[i+l1ido + 1] = di2;
      VCPLXMUL(dr3, di3, LD_PS1(wr2), LD_PS1(wi2));
      ch[i+2*l1ido] = dr3;
      ch[i+2*l1ido+1] = di3;
    }
  }
}

/*
  fsign goes to passf3_loop_ps as a constant, as in passf4_ps and
  passf8_ps: fsign*wa[i+1] of the twiddles is then a negation or nothing.
*/
static NEVER_INLINE(void) passf3_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                    const float *wa1, const float *wa2, float fsign) {
  if (fsign < 0) {
    passf3_loop_ps(ido, l1, cc, ch, wa1, wa2, -1, 0);
  } else {
    passf3_loop_ps(ido, l1, cc, ch, wa1, wa2, +1, 0);
  }
} /* passf3 */
#if PF_ON
static NEVER_INLINE(void) passf3_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                       const float *wa1, const float *wa2, float fsign) {
  if (fsign < 0) {
    passf3_loop_ps(ido, l1, cc, ch, wa1, wa2, -1, 1);
  } else {
    passf3_loop_ps(ido, l1, cc, ch, wa1, wa2, +1, 1);
  }
}
#endif

static ALWAYS_INLINE(void) passf4_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                          const float *wa1, const float *wa2, const float *wa3, const float fsign,
                                          const int pf) {
  /* isign == -1 for forward transform and +1 for backward transform */

  int i, k;
  v4sf ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
  int l1ido = l1*ido;
  if (ido == 2) {
    for (k=0; k < l1ido; k += ido, ch += ido, cc += 4*ido) {
      PF_PREFETCH(pf, ch, 1);
      PF_PREFETCH(pf, ch + 1*l1ido, 1);
      PF_PREFETCH(pf, ch + 2*l1ido, 1);
      PF_PREFETCH(pf, ch + 3*l1ido, 1);
      tr1 = VSUB(cc[0], cc[2*ido + 0]);
      tr2 = VADD(cc[0], cc[2*ido + 0]);
      ti1 = VSUB(cc[1], cc[2*ido + 1]);
      ti2 = VADD(cc[1], cc[2*ido + 1]);
      ti4 = VSUB(cc[1*ido + 0], cc[3*ido + 0]);  /* times fsign below */
      tr4 = VSUB(cc[3*ido + 1], cc[1*ido + 1]);
      tr3 = VADD(cc[ido + 0], cc[3*ido + 0]);
      ti3 = VADD(cc[ido + 1], cc[3*ido + 1]);

      ch[0*l1ido + 0] = VADD(tr2, tr3);
      ch[0*l1ido + 1] = VADD(ti2, ti3);
      ch[1*l1ido + 0] = VADD_FSIGN(tr1, tr4);
      ch[1*l1ido + 1] = VADD_FSIGN(ti1, ti4);
      ch[2*l1ido + 0] = VSUB(tr2, tr3);
      ch[2*l1ido + 1] = VSUB(ti2, ti3);
      ch[3*l1ido + 0] = VSUB_FSIGN(tr1, tr4);
      ch[3*l1ido + 1] = VSUB_FSIGN(ti1, ti4);
    }
  } else {
    for (k=0; k < l1ido; k += ido, ch+=ido, cc += 4*ido) {
      for (i=0; i<ido-1; i+=2) {
        float wr1, wi1, wr2, wi2, wr3, wi3;
        PF_PREFETCH(pf, ch + i, 1);
        PF_PREFETCH(pf, ch + i + 1*l1ido, 1);
        PF_PREFETCH(pf, ch + i + 2*l1ido, 1);
        PF_PREFETCH(pf, ch + i + 3*l1ido, 1);
        tr1 = VSUB(cc[i + 0], cc[i + 2*ido + 0]);
        tr2 = VADD(cc[i + 0], cc[i + 2*ido + 0]);
        ti1 = VSUB(cc[i + 1], cc[i + 2*ido + 1]);
        ti2 = VADD(cc[i + 1], cc[i + 2*ido + 1]);
        tr4 = VSUB(cc[i + 3*ido + 1], cc[i + 1*ido + 1]);  /* times fsign below */
        ti4 = VSUB(cc[i + 1*ido + 0], cc[i + 3*ido + 0]);
        tr3 = VADD(cc[i + ido + 0], cc[i + 3*ido + 0]);
        ti3 = VADD(cc[i + ido + 1], cc[i + 3*ido + 1]);

        ch[i] = VADD(tr2, tr3);
        cr3    = VSUB(tr2, tr3);
        ch[i + 1] = VADD(ti2, ti3);
        ci3 = VSUB(ti2, ti3);

        cr2 = VADD_FSIGN(tr1, tr4);
        cr4 = VSUB_FSIGN(tr1, tr4);
        ci2 = VADD_FSIGN(ti1, ti4);
        ci4 = VSUB_FSIGN(ti1, ti4);
        wr1=wa1[i]; wi1=fsign*wa1[i+1];
        VCPLXMUL(cr2, ci2, LD_PS1(wr1), LD_PS1(wi1));
        wr2=wa2[i]; wi2=fsign*wa2[i+1];
        ch[i + l1ido] = cr2;
        ch[i + l1ido + 1] = ci2;

        VCPLXMUL(cr3, ci3, LD_PS1(wr2), LD_PS1(wi2));
        wr3=wa3[i]; wi3=fsign*wa3[i+1];
        ch[i + 2*l1ido] = cr3;
        ch[i + 2*l1ido + 1] = ci3;

        VCPLXMUL(cr4, ci4, LD_PS1(wr3), LD_PS1(wi3));
        ch[i + 3*l1ido] = cr4;
        ch[i + 3*l1ido + 1] = ci4;
      }
    }
  }
}

/*
  fsign goes to passf4_loop_ps as a constant: VADD_FSIGN and VSUB_FSIGN
  are then one VADD or VSUB. gcc and clang also make the twiddles'
  fsign*wa a negation; MSVC does that only for fsign == +1.
*/
static NEVER_INLINE(void) passf4_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                    const float *wa1, const float *wa2, const float *wa3, float fsign) {
  if (fsign < 0) {
    passf4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, -1, 0);
  } else {
    passf4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, +1, 0);
  }
} /* passf4 */
#if PF_ON
static NEVER_INLINE(void) passf4_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                       const float *wa1, const float *wa2, const float *wa3, float fsign) {
  if (fsign < 0) {
    passf4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, -1, 1);
  } else {
    passf4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, +1, 1);
  }
}
#endif

/*
  passf5 and passb5 has been merged here, fsign = -1 for passf5, +1 for passb5
*/
static ALWAYS_INLINE(void) passf5_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                          const float *wa1, const float *wa2,
                                          const float *wa3, const float *wa4, const float fsign) {
  static const float tr11 = (float)0.30901699437494742410;   /* cos(2*pi/5) = (sqrt(5)-1)/4 */
  const float ti11 = (float)0.95105651629515357212*fsign;    /* sin(2*pi/5) */
  static const float tr12 = (float)-0.80901699437494742410;  /* cos(4*pi/5) = -(sqrt(5)+1)/4 */
  const float ti12 = (float)0.58778525229247312917*fsign;    /* sin(4*pi/5) */

  /* Local variables */
  int i, k;
  v4sf ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
    ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;

  float wr1, wi1, wr2, wi2, wr3, wi3, wr4, wi4;

#define cc_ref(a_1,a_2) cc[(a_2-1)*ido + a_1 + 1]
#define ch_ref(a_1,a_3) ch[(a_3-1)*l1*ido + a_1 + 1]

  assert(ido > 2);
  for (k = 0; k < l1; ++k, cc += 5*ido, ch += ido) {
    for (i = 0; i < ido-1; i += 2) {
      ti5 = VSUB(cc_ref(i  , 2), cc_ref(i  , 5));
      ti2 = VADD(cc_ref(i  , 2), cc_ref(i  , 5));
      ti4 = VSUB(cc_ref(i  , 3), cc_ref(i  , 4));
      ti3 = VADD(cc_ref(i  , 3), cc_ref(i  , 4));
      tr5 = VSUB(cc_ref(i-1, 2), cc_ref(i-1, 5));
      tr2 = VADD(cc_ref(i-1, 2), cc_ref(i-1, 5));
      tr4 = VSUB(cc_ref(i-1, 3), cc_ref(i-1, 4));
      tr3 = VADD(cc_ref(i-1, 3), cc_ref(i-1, 4));
      ch_ref(i-1, 1) = VADD(cc_ref(i-1, 1), VADD(tr2, tr3));
      ch_ref(i  , 1) = VADD(cc_ref(i  , 1), VADD(ti2, ti3));
      cr2 = VADD(cc_ref(i-1, 1), SVMADD(tr11, tr2, SVMUL(tr12, tr3)));
      ci2 = VADD(cc_ref(i  , 1), SVMADD(tr11, ti2, SVMUL(tr12, ti3)));
      cr3 = VADD(cc_ref(i-1, 1), SVMADD(tr12, tr2, SVMUL(tr11, tr3)));
      ci3 = VADD(cc_ref(i  , 1), SVMADD(tr12, ti2, SVMUL(tr11, ti3)));
      cr5 = SVMADD(ti11, tr5, SVMUL(ti12, tr4));
      ci5 = SVMADD(ti11, ti5, SVMUL(ti12, ti4));
      cr4 = SVMSUB(ti11, tr4, SVMUL(ti12, tr5));
      ci4 = SVMSUB(ti11, ti4, SVMUL(ti12, ti5));
      dr3 = VSUB(cr3, ci4);
      dr4 = VADD(cr3, ci4);
      di3 = VADD(ci3, cr4);
      di4 = VSUB(ci3, cr4);
      dr5 = VADD(cr2, ci5);
      dr2 = VSUB(cr2, ci5);
      di5 = VSUB(ci2, cr5);
      di2 = VADD(ci2, cr5);
      wr1=wa1[i]; wi1=fsign*wa1[i+1]; wr2=wa2[i]; wi2=fsign*wa2[i+1];
      wr3=wa3[i]; wi3=fsign*wa3[i+1]; wr4=wa4[i]; wi4=fsign*wa4[i+1];
      VCPLXMUL(dr2, di2, LD_PS1(wr1), LD_PS1(wi1));
      ch_ref(i - 1, 2) = dr2;
      ch_ref(i, 2)     = di2;
      VCPLXMUL(dr3, di3, LD_PS1(wr2), LD_PS1(wi2));
      ch_ref(i - 1, 3) = dr3;
      ch_ref(i, 3)     = di3;
      VCPLXMUL(dr4, di4, LD_PS1(wr3), LD_PS1(wi3));
      ch_ref(i - 1, 4) = dr4;
      ch_ref(i, 4)     = di4;
      VCPLXMUL(dr5, di5, LD_PS1(wr4), LD_PS1(wi4));
      ch_ref(i - 1, 5) = dr5;
      ch_ref(i, 5)     = di5;
    }
  }
#undef ch_ref
#undef cc_ref
}

/* fsign as a constant, as in passf3_ps */
static NEVER_INLINE(void) passf5_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                    const float *wa1, const float *wa2,
                                    const float *wa3, const float *wa4, float fsign) {
  /* no prefetches here: with them, single-precision complex lengths with
     radix-5 passes were up to 2.6% slower on a Xeon 6975P-C, about as much as
     double precision gained there */
  if (fsign < 0) {
    passf5_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, wa4, -1);
  } else {
    passf5_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, wa4, +1);
  }
}

/*
  passf8 and passb8 merged here, fsign = -1 for passf8, +1 for passb8

  Radix-8 DFT via two radix-4 DFTs:
    even outputs from DFT-4 of (x0+x4, x1+x5, x2+x6, x3+x7)
    odd  outputs from DFT-4 of (x0-x4, x1-x5, x2-x6, x3-x7)*W8^k

  cffti1_ps puts a factor 8 only last, where its ido is 1 and so the ido
  passed here (idot of cfftf1_ps) is 2: each leg is one complex vector
  and the twiddles are all 1, as in the ido == 2 branch of passf4_ps. So
  this pass takes no twiddles, and the strides are constants.
*/
static ALWAYS_INLINE(void) passf8_loop_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float fsign,
                                          const int pf) {
  static const float K8 = (float)M_SQRT1_2;   /* cos(pi/4) = sin(pi/4) = 1/sqrt(2) */
  const float nK8 = -K8;
  const int l1ido = 2*l1;
  int k;

  for (k = 0; k < l1ido; k += 2, cc += 16, ch += 2) {
    /* written out: MSVC does not unroll a loop of these */
    PF_PREFETCH(pf, ch, 1);
    PF_PREFETCH(pf, ch + 1*l1ido, 1);
    PF_PREFETCH(pf, ch + 2*l1ido, 1);
    PF_PREFETCH(pf, ch + 3*l1ido, 1);
    PF_PREFETCH(pf, ch + 4*l1ido, 1);
    PF_PREFETCH(pf, ch + 5*l1ido, 1);
    PF_PREFETCH(pf, ch + 6*l1ido, 1);
    PF_PREFETCH(pf, ch + 7*l1ido, 1);
    /* Load 8 inputs: cc[2*leg + 0] = real, cc[2*leg + 1] = imag */
    v4sf x0r = cc[ 0], x0i = cc[ 1];
    v4sf x1r = cc[ 2], x1i = cc[ 3];
    v4sf x2r = cc[ 4], x2i = cc[ 5];
    v4sf x3r = cc[ 6], x3i = cc[ 7];
    v4sf x4r = cc[ 8], x4i = cc[ 9];
    v4sf x5r = cc[10], x5i = cc[11];
    v4sf x6r = cc[12], x6i = cc[13];
    v4sf x7r = cc[14], x7i = cc[15];

    /* Step 1: sums and diffs at stride 4 */
    v4sf s0r = VADD(x0r, x4r), s0i = VADD(x0i, x4i);
    v4sf s1r = VADD(x1r, x5r), s1i = VADD(x1i, x5i);
    v4sf s2r = VADD(x2r, x6r), s2i = VADD(x2i, x6i);
    v4sf s3r = VADD(x3r, x7r), s3i = VADD(x3i, x7i);
    v4sf d0r = VSUB(x0r, x4r), d0i = VSUB(x0i, x4i);
    v4sf d1r = VSUB(x1r, x5r), d1i = VSUB(x1i, x5i);
    v4sf d2r = VSUB(x2r, x6r), d2i = VSUB(x2i, x6i);
    v4sf d3r = VSUB(x3r, x7r), d3i = VSUB(x3i, x7i);

    /* Step 2: twiddle diffs by W8^k */
    /* d0 *= W8^0 = 1 (no-op) */
    /* d1 *= W8^1 = K8*(1 + fsign*j) */
    { v4sf a = d1r, b = d1i;
      d1r = VMUL(LD_PS1(K8), VSUB_FSIGN(a, b));
      d1i = VMUL(LD_PS1(K8), VADD_FSIGN(b, a)); }
    /* d2 *= W8^2 = fsign*j: in Step 4 */
    /* d3 *= W8^3 = K8*(-1 + fsign*j) */
    { v4sf a = d3r, b = d3i;
      d3r = VMUL(LD_PS1(nK8), VADD_FSIGN(a, b));
      d3i = VMUL(LD_PS1(nK8), VSUB_FSIGN(b, a)); }

    /* Step 3: DFT-4 on sums -> even outputs (0,2,4,6) */
    {
      v4sf t1r = VSUB(s0r, s2r), t1i = VSUB(s0i, s2i);
      v4sf t2r = VADD(s0r, s2r), t2i = VADD(s0i, s2i);
      v4sf t4i = VSUB(s1r, s3r);  /* times fsign below */
      v4sf t4r = VSUB(s3i, s1i);
      v4sf t3r = VADD(s1r, s3r), t3i = VADD(s1i, s3i);

      ch[0*l1ido + 0] = VADD(t2r, t3r);
      ch[0*l1ido + 1] = VADD(t2i, t3i);
      ch[2*l1ido + 0] = VADD_FSIGN(t1r, t4r);
      ch[2*l1ido + 1] = VADD_FSIGN(t1i, t4i);
      ch[4*l1ido + 0] = VSUB(t2r, t3r);
      ch[4*l1ido + 1] = VSUB(t2i, t3i);
      ch[6*l1ido + 0] = VSUB_FSIGN(t1r, t4r);
      ch[6*l1ido + 1] = VSUB_FSIGN(t1i, t4i);
    }

    /* Step 4: DFT-4 on diffs -> odd outputs (1,3,5,7) */
    {
      /* d2 times fsign*j: d2r' = -fsign*d2i, d2i' = fsign*d2r */
      v4sf t1r = VADD_FSIGN(d0r, d2i), t1i = VSUB_FSIGN(d0i, d2r);
      v4sf t2r = VSUB_FSIGN(d0r, d2i), t2i = VADD_FSIGN(d0i, d2r);
      v4sf t4i = VSUB(d1r, d3r);  /* times fsign below */
      v4sf t4r = VSUB(d3i, d1i);
      v4sf t3r = VADD(d1r, d3r), t3i = VADD(d1i, d3i);

      ch[1*l1ido + 0] = VADD(t2r, t3r);
      ch[1*l1ido + 1] = VADD(t2i, t3i);
      ch[3*l1ido + 0] = VADD_FSIGN(t1r, t4r);
      ch[3*l1ido + 1] = VADD_FSIGN(t1i, t4i);
      ch[5*l1ido + 0] = VSUB(t2r, t3r);
      ch[5*l1ido + 1] = VSUB(t2i, t3i);
      ch[7*l1ido + 0] = VSUB_FSIGN(t1r, t4r);
      ch[7*l1ido + 1] = VSUB_FSIGN(t1i, t4i);
    }
  }
}

/*
  fsign goes to passf8_loop_ps as a constant: VADD_FSIGN and VSUB_FSIGN
  are then one VADD or VSUB, and each multiply by W8^1 or W8^3 takes 2
  vector multiplies instead of 4.
*/
static NEVER_INLINE(void) passf8_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, float fsign) {
  assert(ido == 2);
  (void)ido;
  if (fsign < 0) {
    passf8_loop_ps(l1, cc, ch, -1, 0);
  } else {
    passf8_loop_ps(l1, cc, ch, +1, 0);
  }
} /* passf8 */
#if PF_ON
static NEVER_INLINE(void) passf8_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, float fsign) {
  assert(ido == 2);
  (void)ido;
  if (fsign < 0) {
    passf8_loop_ps(l1, cc, ch, -1, 1);
  } else {
    passf8_loop_ps(l1, cc, ch, +1, 1);
  }
}
#endif

/*
  passf6 and passb6 merged here, fsign = -1 for passf6, +1 for passb6

  cffti1_ps puts a factor 6 only last, in place of the lone 2 and a 3, so
  as for passf8_ps the ido passed here is 2: each leg is one complex
  vector, and no twiddles. As 2 and 3 are coprime, the 6-point DFT needs
  no twiddles inside either (Good-Thomas): with the input index
  n = 3*n1 + 2*n2 mod 6 and the output index m = m1 mod 2, m = m2 mod 3,
  w6^(n*m) = (-1)^(n1*m1) * w3^(n2*m2). So a DFT-2 of each pair (x0, x3),
  (x2, x5), (x4, x1) (n1 = 0, 1 for n2 = 0, 1, 2), then a DFT-3 of the
  sums (m1 = 0) gives y0, y4, y2 and one of the differences (m1 = 1)
  y3, y1, y5 (m2 = 0, 1, 2). That is 32 vector adds, 4 multiply-adds and
  4 multiplies per k, for 6 legs; passf8_ps takes 52 adds and 4
  multiplies for its 8.
*/
static ALWAYS_INLINE(void) passf6_loop_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float fsign,
                                          const int pf) {
  static const float taur = (float)-0.5;                    /* cos(2*pi/3) */
  static const float taui = (float)0.86602540378443864676;  /* sin(2*pi/3) = sqrt(3)/2 */
  const int l1ido = 2*l1;
  int k;

  for (k = 0; k < l1ido; k += 2, cc += 12, ch += 2) {
    PF_PREFETCH(pf, ch, 1);
    PF_PREFETCH(pf, ch + 1*l1ido, 1);
    PF_PREFETCH(pf, ch + 2*l1ido, 1);
    PF_PREFETCH(pf, ch + 3*l1ido, 1);
    PF_PREFETCH(pf, ch + 4*l1ido, 1);
    PF_PREFETCH(pf, ch + 5*l1ido, 1);
    /* Load 6 inputs: cc[2*leg + 0] = real, cc[2*leg + 1] = imag */
    v4sf x0r = cc[ 0], x0i = cc[ 1];
    v4sf x1r = cc[ 2], x1i = cc[ 3];
    v4sf x2r = cc[ 4], x2i = cc[ 5];
    v4sf x3r = cc[ 6], x3i = cc[ 7];
    v4sf x4r = cc[ 8], x4i = cc[ 9];
    v4sf x5r = cc[10], x5i = cc[11];

    /* DFT-2 of (x0, x3), (x2, x5), (x4, x1) */
    v4sf s0r = VADD(x0r, x3r), s0i = VADD(x0i, x3i);
    v4sf s1r = VADD(x2r, x5r), s1i = VADD(x2i, x5i);
    v4sf s2r = VADD(x4r, x1r), s2i = VADD(x4i, x1i);
    v4sf d0r = VSUB(x0r, x3r), d0i = VSUB(x0i, x3i);
    v4sf d1r = VSUB(x2r, x5r), d1i = VSUB(x2i, x5i);
    v4sf d2r = VSUB(x4r, x1r), d2i = VSUB(x4i, x1i);

    /* DFT-3 of the sums -> y0, y4, y2. With c = s0 + taur*(s1 + s2) and
       u = taui*(s1 - s2): Y1 = c + fsign*j*u, Y2 = c - fsign*j*u */
    {
      v4sf tr = VADD(s1r, s2r), ti = VADD(s1i, s2i);
      v4sf cr = SVMADD(taur, tr, s0r), ci = SVMADD(taur, ti, s0i);
      v4sf ur = SVMUL(taui, VSUB(s1r, s2r)), ui = SVMUL(taui, VSUB(s1i, s2i));
      ch[0*l1ido + 0] = VADD(s0r, tr);
      ch[0*l1ido + 1] = VADD(s0i, ti);
      ch[4*l1ido + 0] = VSUB_FSIGN(cr, ui);
      ch[4*l1ido + 1] = VADD_FSIGN(ci, ur);
      ch[2*l1ido + 0] = VADD_FSIGN(cr, ui);
      ch[2*l1ido + 1] = VSUB_FSIGN(ci, ur);
    }

    /* DFT-3 of the differences -> y3, y1, y5 */
    {
      v4sf tr = VADD(d1r, d2r), ti = VADD(d1i, d2i);
      v4sf cr = SVMADD(taur, tr, d0r), ci = SVMADD(taur, ti, d0i);
      v4sf ur = SVMUL(taui, VSUB(d1r, d2r)), ui = SVMUL(taui, VSUB(d1i, d2i));
      ch[3*l1ido + 0] = VADD(d0r, tr);
      ch[3*l1ido + 1] = VADD(d0i, ti);
      ch[1*l1ido + 0] = VSUB_FSIGN(cr, ui);
      ch[1*l1ido + 1] = VADD_FSIGN(ci, ur);
      ch[5*l1ido + 0] = VADD_FSIGN(cr, ui);
      ch[5*l1ido + 1] = VSUB_FSIGN(ci, ur);
    }
  }
}

/* fsign as a constant, as in passf8_ps */
static NEVER_INLINE(void) passf6_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, float fsign) {
  assert(ido == 2);
  (void)ido;
  if (fsign < 0) {
    passf6_loop_ps(l1, cc, ch, -1, 0);
  } else {
    passf6_loop_ps(l1, cc, ch, +1, 0);
  }
} /* passf6 */
#if PF_ON
static NEVER_INLINE(void) passf6_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, float fsign) {
  assert(ido == 2);
  (void)ido;
  if (fsign < 0) {
    passf6_loop_ps(l1, cc, ch, -1, 1);
  } else {
    passf6_loop_ps(l1, cc, ch, +1, 1);
  }
}
#endif

#undef VADD_FSIGN
#undef VSUB_FSIGN


static ALWAYS_INLINE(void) radf2_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float *wa1,
                                         const int pf) {
  static const float minus_one = -1.f;
  int i, k, l1ido = l1*ido;
  for (k=0; k < l1ido; k += ido) {
    v4sf a = cc[k], b = cc[k + l1ido];
    ch[2*k] = VADD(a, b);
    ch[2*(k+ido)-1] = VSUB(a, b);
  }
  if (ido < 2) return;
  if (ido != 2) {
    for (k=0; k < l1ido; k += ido) {
      for (i=2; i<ido; i+=2) {
        v4sf tr2, ti2, br, bi;
        PF_PREFETCH(pf, ch + i - 1 + 2*k, 1);
        PF_PREFETCH(pf, ch + 2*(k+ido) - i - 1, -1);
        tr2 = cc[i - 1 + k + l1ido]; ti2 = cc[i + k + l1ido];
        br = cc[i - 1 + k]; bi = cc[i + k];
        VCPLXMULCONJ(tr2, ti2, LD_PS1(wa1[i - 2]), LD_PS1(wa1[i - 1]));
        ch[i + 2*k] = VADD(bi, ti2);
        ch[2*(k+ido) - i] = VSUB(ti2, bi);
        ch[i - 1 + 2*k] = VADD(br, tr2);
        ch[2*(k+ido) - i -1] = VSUB(br, tr2);
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k=0; k < l1ido; k += ido) {
    ch[2*k + ido] = SVMUL(minus_one, cc[ido-1 + k + l1ido]);
    ch[2*k + ido-1] = cc[k + ido-1];
  }
}

static NEVER_INLINE(void) radf2_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float *wa1) {
  radf2_loop_ps(ido, l1, cc, ch, wa1, 0);
} /* radf2 */
#if PF_ON
static NEVER_INLINE(void) radf2_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float *wa1) {
  radf2_loop_ps(ido, l1, cc, ch, wa1, 1);
}
#endif


static ALWAYS_INLINE(void) radb2_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float *wa1,
                                         const int pf) {
  static const float minus_two=-2;
  int i, k, l1ido = l1*ido;
  v4sf a,b,c,d, tr2, ti2;
  for (k=0; k < l1ido; k += ido) {
    a = cc[2*k]; b = cc[2*(k+ido) - 1];
    ch[k] = VADD(a, b);
    ch[k + l1ido] =VSUB(a, b);
  }
  if (ido < 2) return;
  if (ido != 2) {
    for (k = 0; k < l1ido; k += ido) {
      for (i = 2; i < ido; i += 2) {
        PF_PREFETCH(pf, ch + i - 1 + k, 1);
        PF_PREFETCH(pf, ch + i - 1 + k + l1ido, 1);
        a = cc[i-1 + 2*k]; b = cc[2*(k + ido) - i - 1];
        c = cc[i+0 + 2*k]; d = cc[2*(k + ido) - i + 0];
        ch[i-1 + k] = VADD(a, b);
        tr2 = VSUB(a, b);
        ch[i+0 + k] = VSUB(c, d);
        ti2 = VADD(c, d);
        VCPLXMUL(tr2, ti2, LD_PS1(wa1[i - 2]), LD_PS1(wa1[i - 1]));
        ch[i-1 + k + l1ido] = tr2;
        ch[i+0 + k + l1ido] = ti2;
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k = 0; k < l1ido; k += ido) {
    a = cc[2*k + ido-1]; b = cc[2*k + ido];
    ch[k + ido-1] = VADD(a,a);
    ch[k + ido-1 + l1ido] = SVMUL(minus_two, b);
  }
}

static NEVER_INLINE(void) radb2_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float *wa1) {
  radb2_loop_ps(ido, l1, cc, ch, wa1, 0);
} /* radb2 */
#if PF_ON
static NEVER_INLINE(void) radb2_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const float *wa1) {
  radb2_loop_ps(ido, l1, cc, ch, wa1, 1);
}
#endif

/*
  The odd-radix passes, radf3/radb3 and radf5/radb5, take an even ido
  too, so that rffti1_ps can put a radix-8 or radix-6 pass last, after
  the 3s and 5s. The i-loop handles the complex pairs (i-1, i) of each leg; with ido
  even, that leaves the last element, ido-1: the leg's Nyquist term, a
  real x_j at frequency ido/2. Its twiddle is exp(-i*pi*j/p), the table's
  exp(-2*i*pi*j*f/(p*ido)) at f = ido/2, so it needs no table entry. For
  a radix p pass, the p of them make the output frequencies (2m+1)*ido/2:
    X_m = sum_j x_j exp(-i*pi*j*(2m+1)/p),  m = 0 .. (p-1)/2,
  the rest being their conjugates. Re X_m goes to ch(ido-1, 2m) and
  Im X_m to ch(0, 2m+1), the half-complex places of frequency
  (2m+1)*ido/2 in the output block of p*ido, and the real
  X_(p-1)/2 = x0 - x1 + x2 - ... is the block's Nyquist term, at
  ch(ido-1, p-1). Backward, unnormalized like the rest,
    x_j = 2 Re sum_{m < (p-1)/2} X_m exp(i*pi*j*(2m+1)/p) + (-1)^j X_(p-1)/2.
  radf4/radb4 do the same in their ido % 2 == 0 tails.
*/
/* the i-loop of radf3_loop_ps over all k; with even (ido even), each k
   ends with its Nyquist terms, while its blocks are still in cache.
   even is a constant, as pf is, so the odd-ido loop compiles as if
   the tail were not there. */
static ALWAYS_INLINE(void) radf3_rows_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                       const float *wa1, const float *wa2,
                                       const int pf, const int even)
{
  static const float taur = (float)-0.5;                     /* cos(2*pi/3)   */
  static const float taui = (float)0.86602540378443864676;   /* sin(2*pi/3) = sqrt(3)/2 */
  static const float minus_taui = (float)-0.86602540378443864676;
  int i, k, ic;
  v4sf ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3, wr1, wi1, wr2, wi2;
  for (k=0; k<l1; k++) {
    for (i=2; i<ido; i+=2) {
      ic = ido - i;
      PF_PREFETCH(pf, ch + i - 1 + 3*k*ido, 1);
      PF_PREFETCH(pf, ch + i - 1 + (3*k + 2)*ido, 1);
      PF_PREFETCH(pf, ch + ic - 1 + (3*k + 1)*ido, -1);
      wr1 = LD_PS1(wa1[i - 2]); wi1 = LD_PS1(wa1[i - 1]);
      dr2 = cc[i - 1 + (k + l1)*ido]; di2 = cc[i + (k + l1)*ido];
      VCPLXMULCONJ(dr2, di2, wr1, wi1);

      wr2 = LD_PS1(wa2[i - 2]); wi2 = LD_PS1(wa2[i - 1]);
      dr3 = cc[i - 1 + (k + l1*2)*ido]; di3 = cc[i + (k + l1*2)*ido];
      VCPLXMULCONJ(dr3, di3, wr2, wi2);

      cr2 = VADD(dr2, dr3);
      ci2 = VADD(di2, di3);
      ch[i - 1 + 3*k*ido] = VADD(cc[i - 1 + k*ido], cr2);
      ch[i + 3*k*ido] = VADD(cc[i + k*ido], ci2);
      tr2 = SVMADD(taur, cr2, cc[i - 1 + k*ido]);
      ti2 = SVMADD(taur, ci2, cc[i + k*ido]);
      tr3 = SVMUL(taui, VSUB(di2, di3));
      ti3 = SVMUL(taui, VSUB(dr3, dr2));
      ch[i - 1 + (3*k + 2)*ido] = VADD(tr2, tr3);
      ch[ic - 1 + (3*k + 1)*ido] = VSUB(tr2, tr3);
      ch[i + (3*k + 2)*ido] = VADD(ti2, ti3);
      ch[ic + (3*k + 1)*ido] = VSUB(ti3, ti2);
    }
    if (even) {
      /* the Nyquist terms x_j (see above): with W = exp(-i*pi/3),
         X0 = x0 + x1 W + x2 W^2 = x0 + (x1 - x2)/2 - i sin(pi/3) (x1 + x2)
         and X1 = x0 - x1 + x2 */
      v4sf x0 = cc[ido-1 + k*ido], x1 = cc[ido-1 + (k + l1)*ido], x2 = cc[ido-1 + (k + 2*l1)*ido];
      v4sf d = VSUB(x2, x1);
      ch[ido-1 + 3*k*ido] = SVMADD(taur, d, x0);                /* Re X0 */
      ch[(3*k + 1)*ido] = SVMUL(minus_taui, VADD(x1, x2));      /* Im X0 */
      ch[ido-1 + (3*k + 2)*ido] = VADD(x0, d);                  /* X1 */
    }
  }
}

static ALWAYS_INLINE(void) radf3_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                     const float *wa1, const float *wa2, const int pf) {
  static const float taur = (float)-0.5;                     /* cos(2*pi/3)   */
  static const float taui = (float)0.86602540378443864676;   /* sin(2*pi/3) = sqrt(3)/2 */
  int k;
  v4sf ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3, wr1, wi1, wr2, wi2;
  for (k=0; k<l1; k++) {
    cr2 = VADD(cc[(k + l1)*ido], cc[(k + 2*l1)*ido]);
    ch[3*k*ido] = VADD(cc[k*ido], cr2);
    ch[(3*k+2)*ido] = SVMUL(taui, VSUB(cc[(k + l1*2)*ido], cc[(k + l1)*ido]));
    ch[ido-1 + (3*k + 1)*ido] = SVMADD(taur, cr2, cc[k*ido]);
  }
  if (ido == 1) return;
  if (ido == 3) {
    /* Fast path: inner loop executes exactly once with i=2, ic=1 */
    wr1 = LD_PS1(wa1[0]); wi1 = LD_PS1(wa1[1]);
    wr2 = LD_PS1(wa2[0]); wi2 = LD_PS1(wa2[1]);
    for (k=0; k<l1; k++) {
      dr2 = cc[1 + (k + l1)*ido]; di2 = cc[2 + (k + l1)*ido];
      VCPLXMULCONJ(dr2, di2, wr1, wi1);
      dr3 = cc[1 + (k + l1*2)*ido]; di3 = cc[2 + (k + l1*2)*ido];
      VCPLXMULCONJ(dr3, di3, wr2, wi2);
      cr2 = VADD(dr2, dr3);
      ci2 = VADD(di2, di3);
      ch[1 + 3*k*ido] = VADD(cc[1 + k*ido], cr2);
      ch[2 + 3*k*ido] = VADD(cc[2 + k*ido], ci2);
      tr2 = SVMADD(taur, cr2, cc[1 + k*ido]);
      ti2 = SVMADD(taur, ci2, cc[2 + k*ido]);
      tr3 = SVMUL(taui, VSUB(di2, di3));
      ti3 = SVMUL(taui, VSUB(dr3, dr2));
      ch[1 + (3*k + 2)*ido] = VADD(tr2, tr3);
      ch[0 + (3*k + 1)*ido] = VSUB(tr2, tr3);
      ch[2 + (3*k + 2)*ido] = VADD(ti2, ti3);
      ch[1 + (3*k + 1)*ido] = VSUB(ti3, ti2);
    }
    return;
  }
  if (ido & 1) {
    radf3_rows_ps(ido, l1, cc, ch, wa1, wa2, pf, 0);
  } else {
    radf3_rows_ps(ido, l1, cc, ch, wa1, wa2, pf, 1);
  }
}

static NEVER_INLINE(void) radf3_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                     const float *wa1, const float *wa2) {
  radf3_loop_ps(ido, l1, cc, ch, wa1, wa2, 0);
} /* radf3 */
#if PF_ON
static NEVER_INLINE(void) radf3_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                     const float *wa1, const float *wa2) {
  radf3_loop_ps(ido, l1, cc, ch, wa1, wa2, 1);
}
#endif


/* the i-loop of radb3_loop_ps over all k; with even (ido even), each k
   ends with its Nyquist terms, while its blocks are still in cache.
   even is a constant, as pf is, so the odd-ido loop compiles as if
   the tail were not there. */
static ALWAYS_INLINE(void) radb3_rows_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                       const float *wa1, const float *wa2,
                                       const int pf, const int even)
{
  static const float taur = (float)-0.5;                     /* cos(2*pi/3)   */
  static const float taui = (float)0.86602540378443864676;   /* sin(2*pi/3) = sqrt(3)/2 */
  static const float minus_taui_2 = (float)-1.7320508075688772935;
  int i, k, ic;
  v4sf ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
  for (k=0; k<l1; k++) {
    for (i=2; i<ido; i+=2) {
      ic = ido - i;
      PF_PREFETCH(pf, ch + i - 1 + k*ido, 1);
      PF_PREFETCH(pf, ch + i - 1 + (k + l1)*ido, 1);
      PF_PREFETCH(pf, ch + i - 1 + (k + 2*l1)*ido, 1);
      tr2 = VADD(cc[i - 1 + (3*k + 2)*ido], cc[ic - 1 + (3*k + 1)*ido]);
      cr2 = SVMADD(taur, tr2, cc[i - 1 + 3*k*ido]);
      ch[i - 1 + k*ido] = VADD(cc[i - 1 + 3*k*ido], tr2);
      ti2 = VSUB(cc[i + (3*k + 2)*ido], cc[ic + (3*k + 1)*ido]);
      ci2 = SVMADD(taur, ti2, cc[i + 3*k*ido]);
      ch[i + k*ido] = VADD(cc[i + 3*k*ido], ti2);
      cr3 = SVMUL(taui, VSUB(cc[i - 1 + (3*k + 2)*ido], cc[ic - 1 + (3*k + 1)*ido]));
      ci3 = SVMUL(taui, VADD(cc[i + (3*k + 2)*ido], cc[ic + (3*k + 1)*ido]));
      dr2 = VSUB(cr2, ci3);
      dr3 = VADD(cr2, ci3);
      di2 = VADD(ci2, cr3);
      di3 = VSUB(ci2, cr3);
      VCPLXMUL(dr2, di2, LD_PS1(wa1[i-2]), LD_PS1(wa1[i-1]));
      ch[i - 1 + (k + l1)*ido] = dr2;
      ch[i + (k + l1)*ido] = di2;
      VCPLXMUL(dr3, di3, LD_PS1(wa2[i-2]), LD_PS1(wa2[i-1]));
      ch[i - 1 + (k + 2*l1)*ido] = dr3;
      ch[i + (k + 2*l1)*ido] = di3;
    }
    if (even) {
      /* the Nyquist terms, from Re X0, Im X0 and X1 (see above radf3_loop_ps):
         x0 = 2 Re X0 + X1, x1, x2 = +-(Re X0 - X1) - sqrt(3) Im X0 */
      v4sf r = cc[ido-1 + 3*k*ido], im = cc[(3*k + 1)*ido], q = cc[ido-1 + (3*k + 2)*ido];
      v4sf d = VSUB(r, q), w = SVMUL(minus_taui_2, im);
      ch[ido-1 + k*ido] = VADD(VADD(r, r), q);
      ch[ido-1 + (k + l1)*ido] = VADD(d, w);
      ch[ido-1 + (k + 2*l1)*ido] = VSUB(w, d);
    }
  }
}

static ALWAYS_INLINE(void) radb3_loop_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf *RESTRICT ch,
                     const float *wa1, const float *wa2, const int pf)
{
  static const float taur = (float)-0.5;                     /* cos(2*pi/3)   */
  static const float taui = (float)0.86602540378443864676;   /* sin(2*pi/3) = sqrt(3)/2 */
  static const float taui_2 = (float)1.7320508075688772935;  /* 2*sin(2*pi/3) = sqrt(3) */
  int k;
  v4sf ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
  for (k=0; k<l1; k++) {
    tr2 = cc[ido-1 + (3*k + 1)*ido]; tr2 = VADD(tr2,tr2);
    cr2 = SVMADD(taur, tr2, cc[3*k*ido]);
    ch[k*ido] = VADD(cc[3*k*ido], tr2);
    ci3 = SVMUL(taui_2, cc[(3*k + 2)*ido]);
    ch[(k + l1)*ido] = VSUB(cr2, ci3);
    ch[(k + 2*l1)*ido] = VADD(cr2, ci3);
  }
  if (ido == 1) return;
  if (ido == 3) {
    /* Fast path: inner loop executes exactly once with i=2, ic=1 */
    v4sf wa1_r = LD_PS1(wa1[0]), wa1_i = LD_PS1(wa1[1]);
    v4sf wa2_r = LD_PS1(wa2[0]), wa2_i = LD_PS1(wa2[1]);
    for (k=0; k<l1; k++) {
      tr2 = VADD(cc[1 + (3*k + 2)*ido], cc[0 + (3*k + 1)*ido]);
      cr2 = SVMADD(taur, tr2, cc[1 + 3*k*ido]);
      ch[1 + k*ido] = VADD(cc[1 + 3*k*ido], tr2);
      ti2 = VSUB(cc[2 + (3*k + 2)*ido], cc[1 + (3*k + 1)*ido]);
      ci2 = SVMADD(taur, ti2, cc[2 + 3*k*ido]);
      ch[2 + k*ido] = VADD(cc[2 + 3*k*ido], ti2);
      cr3 = SVMUL(taui, VSUB(cc[1 + (3*k + 2)*ido], cc[0 + (3*k + 1)*ido]));
      ci3 = SVMUL(taui, VADD(cc[2 + (3*k + 2)*ido], cc[1 + (3*k + 1)*ido]));
      dr2 = VSUB(cr2, ci3);
      dr3 = VADD(cr2, ci3);
      di2 = VADD(ci2, cr3);
      di3 = VSUB(ci2, cr3);
      VCPLXMUL(dr2, di2, wa1_r, wa1_i);
      ch[1 + (k + l1)*ido] = dr2;
      ch[2 + (k + l1)*ido] = di2;
      VCPLXMUL(dr3, di3, wa2_r, wa2_i);
      ch[1 + (k + 2*l1)*ido] = dr3;
      ch[2 + (k + 2*l1)*ido] = di3;
    }
    return;
  }
  if (ido & 1) {
    radb3_rows_ps(ido, l1, cc, ch, wa1, wa2, pf, 0);
  } else {
    radb3_rows_ps(ido, l1, cc, ch, wa1, wa2, pf, 1);
  }
}

static NEVER_INLINE(void) radb3_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf *RESTRICT ch,
                     const float *wa1, const float *wa2) {
  radb3_loop_ps(ido, l1, cc, ch, wa1, wa2, 0);
} /* radb3 */
#if PF_ON
static NEVER_INLINE(void) radb3_pf_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf *RESTRICT ch,
                     const float *wa1, const float *wa2) {
  radb3_loop_ps(ido, l1, cc, ch, wa1, wa2, 1);
}
#endif

static ALWAYS_INLINE(void) radf4_loop_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf * RESTRICT ch,
                                   const float * RESTRICT wa1, const float * RESTRICT wa2, const float * RESTRICT wa3,
                                   const int pf)
{
  static const float minus_hsqt2 = (float)-M_SQRT1_2;  /* -sin(pi/4) = -1/sqrt(2) */
  int i, k, l1ido = l1*ido;
  {
    const v4sf *RESTRICT cc_ = cc, * RESTRICT cc_end = cc + l1ido;
    v4sf * RESTRICT ch_ = ch;
    while (cc < cc_end) {
      /* this loop represents between 25% and 40% of total radf4_ps cost ! */
      v4sf a0 = cc[0], a1 = cc[l1ido];
      v4sf a2 = cc[2*l1ido], a3 = cc[3*l1ido];
      v4sf tr1 = VADD(a1, a3);
      v4sf tr2 = VADD(a0, a2);
      ch[2*ido-1] = VSUB(a0, a2);
      ch[2*ido  ] = VSUB(a3, a1);
      ch[0      ] = VADD(tr1, tr2);
      ch[4*ido-1] = VSUB(tr2, tr1);
      cc += ido; ch += 4*ido;
    }
    cc = cc_; ch = ch_;
  }
  if (ido < 2) return;
  if (ido == 3) {
    /* Fast path: inner loop executes exactly once with i=2, ic=ido-i=1 */
    v4sf wr, wi, cr2, ci2, cr3, ci3, cr4, ci4;
    v4sf tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4;
    for (k = 0; k < l1ido; k += ido) {
      cr2 = cc[1 + k + 1*l1ido];
      ci2 = cc[1 + k + 1*l1ido + 1];
      wr = LD_PS1(wa1[0]); wi = LD_PS1(wa1[1]);
      VCPLXMULCONJ(cr2, ci2, wr, wi);

      cr3 = cc[1 + k + 2*l1ido];
      ci3 = cc[1 + k + 2*l1ido + 1];
      wr = LD_PS1(wa2[0]); wi = LD_PS1(wa2[1]);
      VCPLXMULCONJ(cr3, ci3, wr, wi);

      cr4 = cc[1 + k + 3*l1ido];
      ci4 = cc[1 + k + 3*l1ido + 1];
      wr = LD_PS1(wa3[0]); wi = LD_PS1(wa3[1]);
      VCPLXMULCONJ(cr4, ci4, wr, wi);

      tr1 = VADD(cr2, cr4);
      tr4 = VSUB(cr4, cr2);
      tr2 = VADD(cc[1 + k], cr3);
      tr3 = VSUB(cc[1 + k], cr3);
      ch[1 + 4*k] = VADD(tr1, tr2);
      ch[0 + 4*k + 3*ido] = VSUB(tr2, tr1);
      ti1 = VADD(ci2, ci4);
      ti4 = VSUB(ci2, ci4);
      ch[1 + 4*k + 2*ido] = VADD(ti4, tr3);
      ch[0 + 4*k + 1*ido] = VSUB(tr3, ti4);
      ti2 = VADD(cc[1 + k + 1], ci3);
      ti3 = VSUB(cc[1 + k + 1], ci3);
      ch[2 + 4*k] = VADD(ti1, ti2);
      ch[1 + 4*k + 3*ido] = VSUB(ti1, ti2);
      ch[2 + 4*k + 2*ido] = VADD(tr4, ti3);
      ch[1 + 4*k + 1*ido] = VSUB(tr4, ti3);
    }
    return;
  }
  if (ido != 2) {
    for (k = 0; k < l1ido; k += ido) {
      const v4sf * RESTRICT pc = (v4sf*)(cc + 1 + k);
      for (i=2; i<ido; i += 2, pc += 2) {
        int ic = ido - i;
        v4sf wr, wi, cr2, ci2, cr3, ci3, cr4, ci4;
        v4sf tr1, ti1, tr2, ti2, tr3, ti3, tr4, ti4;

        PF_PREFETCH(pf, ch + i - 1 + 4*k, 1);
        PF_PREFETCH(pf, ch + i - 1 + 4*k + 2*ido, 1);
        PF_PREFETCH(pf, ch + ic - 1 + 4*k + 1*ido, -1);
        PF_PREFETCH(pf, ch + ic - 1 + 4*k + 3*ido, -1);
        cr2 = pc[1*l1ido+0];
        ci2 = pc[1*l1ido+1];
        wr=LD_PS1(wa1[i - 2]);
        wi=LD_PS1(wa1[i - 1]);
        VCPLXMULCONJ(cr2,ci2,wr,wi);

        cr3 = pc[2*l1ido+0];
        ci3 = pc[2*l1ido+1];
        wr = LD_PS1(wa2[i-2]);
        wi = LD_PS1(wa2[i-1]);
        VCPLXMULCONJ(cr3, ci3, wr, wi);

        cr4 = pc[3*l1ido];
        ci4 = pc[3*l1ido+1];
        wr = LD_PS1(wa3[i-2]);
        wi = LD_PS1(wa3[i-1]);
        VCPLXMULCONJ(cr4, ci4, wr, wi);

        /* at this point, on SSE, five of "cr2 cr3 cr4 ci2 ci3 ci4" should be loaded in registers */

        tr1 = VADD(cr2,cr4);
        tr4 = VSUB(cr4,cr2);
        tr2 = VADD(pc[0],cr3);
        tr3 = VSUB(pc[0],cr3);
        ch[i - 1 + 4*k] = VADD(tr1,tr2);
        ch[ic - 1 + 4*k + 3*ido] = VSUB(tr2,tr1); /* at this point tr1 and tr2 can be disposed */
        ti1 = VADD(ci2,ci4);
        ti4 = VSUB(ci2,ci4);
        ch[i - 1 + 4*k + 2*ido] = VADD(ti4,tr3);
        ch[ic - 1 + 4*k + 1*ido] = VSUB(tr3,ti4); /* dispose tr3, ti4 */
        ti2 = VADD(pc[1],ci3);
        ti3 = VSUB(pc[1],ci3);
        ch[i + 4*k] = VADD(ti1, ti2);
        ch[ic + 4*k + 3*ido] = VSUB(ti1, ti2);
        ch[i + 4*k + 2*ido] = VADD(tr4, ti3);
        ch[ic + 4*k + 1*ido] = VSUB(tr4, ti3);
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k=0; k<l1ido; k += ido) {
    v4sf a = cc[ido-1 + k + l1ido], b = cc[ido-1 + k + 3*l1ido];
    v4sf c = cc[ido-1 + k], d = cc[ido-1 + k + 2*l1ido];
    v4sf ti1 = SVMUL(minus_hsqt2, VADD(a, b));
    v4sf tr1 = SVMUL(minus_hsqt2, VSUB(b, a));
    ch[ido-1 + 4*k] = VADD(tr1, c);
    ch[ido-1 + 4*k + 2*ido] = VSUB(c, tr1);
    ch[4*k + 1*ido] = VSUB(ti1, d);
    ch[4*k + 3*ido] = VADD(ti1, d);
  }
}

static NEVER_INLINE(void) radf4_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf * RESTRICT ch,
                                   const float * RESTRICT wa1, const float * RESTRICT wa2, const float * RESTRICT wa3) {
  radf4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, 0);
} /* radf4 */
#if PF_ON
static NEVER_INLINE(void) radf4_pf_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf * RESTRICT ch,
                                      const float * RESTRICT wa1, const float * RESTRICT wa2, const float * RESTRICT wa3) {
  radf4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, 1);
}
#endif


static ALWAYS_INLINE(void) radb4_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                   const float * RESTRICT wa1, const float * RESTRICT wa2, const float *RESTRICT wa3,
                                   const int pf)
{
  static const float minus_sqrt2 = (float)-M_SQRT2;     /* -2*sin(pi/4) = -sqrt(2) */
  static const float two = 2.f;
  int i, k, l1ido = l1*ido;
  v4sf ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
  {
    const v4sf *RESTRICT cc_ = cc, * RESTRICT ch_end = ch + l1ido;
    v4sf *ch_ = ch;
    while (ch < ch_end) {
      v4sf a = cc[0], b = cc[4*ido-1];
      v4sf c = cc[2*ido], d = cc[2*ido-1];
      tr3 = SVMUL(two,d);
      tr2 = VADD(a,b);
      tr1 = VSUB(a,b);
      tr4 = SVMUL(two,c);
      ch[0*l1ido] = VADD(tr2, tr3);
      ch[2*l1ido] = VSUB(tr2, tr3);
      ch[1*l1ido] = VSUB(tr1, tr4);
      ch[3*l1ido] = VADD(tr1, tr4);

      cc += 4*ido; ch += ido;
    }
    cc = cc_; ch = ch_;
  }
  if (ido < 2) return;
  if (ido == 3) {
    /* Fast path: inner loop executes exactly once with i=2, ic=ido-i=1 */
    v4sf wa1_r = LD_PS1(wa1[0]), wa1_i = LD_PS1(wa1[1]);
    v4sf wa2_r = LD_PS1(wa2[0]), wa2_i = LD_PS1(wa2[1]);
    v4sf wa3_r = LD_PS1(wa3[0]), wa3_i = LD_PS1(wa3[1]);
    for (k = 0; k < l1ido; k += ido) {
      tr1 = VSUB(cc[1 + 4*k], cc[9 + 4*k]);
      tr2 = VADD(cc[1 + 4*k], cc[9 + 4*k]);
      ti4 = VSUB(cc[7 + 4*k], cc[3 + 4*k]);
      tr3 = VADD(cc[7 + 4*k], cc[3 + 4*k]);
      ch[k + 1] = VADD(tr2, tr3);
      cr3 = VSUB(tr2, tr3);

      ti3 = VSUB(cc[8 + 4*k], cc[4 + 4*k]);
      tr4 = VADD(cc[8 + 4*k], cc[4 + 4*k]);
      cr2 = VSUB(tr1, tr4);
      cr4 = VADD(tr1, tr4);

      ti1 = VADD(cc[2 + 4*k], cc[10 + 4*k]);
      ti2 = VSUB(cc[2 + 4*k], cc[10 + 4*k]);

      ch[k + 2] = VADD(ti2, ti3);
      ci3 = VSUB(ti2, ti3);
      ci2 = VADD(ti1, ti4);
      ci4 = VSUB(ti1, ti4);
      VCPLXMUL(cr2, ci2, wa1_r, wa1_i);
      ch[k + 1 + l1ido] = cr2;
      ch[k + 2 + l1ido] = ci2;
      VCPLXMUL(cr3, ci3, wa2_r, wa2_i);
      ch[k + 1 + 2*l1ido] = cr3;
      ch[k + 2 + 2*l1ido] = ci3;
      VCPLXMUL(cr4, ci4, wa3_r, wa3_i);
      ch[k + 1 + 3*l1ido] = cr4;
      ch[k + 2 + 3*l1ido] = ci4;
    }
    return;
  }
  if (ido != 2) {
    for (k = 0; k < l1ido; k += ido) {
      const v4sf * RESTRICT pc = (v4sf*)(cc - 1 + 4*k);
      v4sf * RESTRICT ph = (v4sf*)(ch + k + 1);
      for (i = 2; i < ido; i += 2) {
        PF_PREFETCH(pf, ph, 1);
        PF_PREFETCH(pf, ph + l1ido, 1);
        PF_PREFETCH(pf, ph + 2*l1ido, 1);
        PF_PREFETCH(pf, ph + 3*l1ido, 1);

        tr1 = VSUB(pc[i], pc[4*ido - i]);
        tr2 = VADD(pc[i], pc[4*ido - i]);
        ti4 = VSUB(pc[2*ido + i], pc[2*ido - i]);
        tr3 = VADD(pc[2*ido + i], pc[2*ido - i]);
        ph[0] = VADD(tr2, tr3);
        cr3 = VSUB(tr2, tr3);

        ti3 = VSUB(pc[2*ido + i + 1], pc[2*ido - i + 1]);
        tr4 = VADD(pc[2*ido + i + 1], pc[2*ido - i + 1]);
        cr2 = VSUB(tr1, tr4);
        cr4 = VADD(tr1, tr4);

        ti1 = VADD(pc[i + 1], pc[4*ido - i + 1]);
        ti2 = VSUB(pc[i + 1], pc[4*ido - i + 1]);

        ph[1] = VADD(ti2, ti3); ph += l1ido;
        ci3 = VSUB(ti2, ti3);
        ci2 = VADD(ti1, ti4);
        ci4 = VSUB(ti1, ti4);
        VCPLXMUL(cr2, ci2, LD_PS1(wa1[i-2]), LD_PS1(wa1[i-1]));
        ph[0] = cr2;
        ph[1] = ci2; ph += l1ido;
        VCPLXMUL(cr3, ci3, LD_PS1(wa2[i-2]), LD_PS1(wa2[i-1]));
        ph[0] = cr3;
        ph[1] = ci3; ph += l1ido;
        VCPLXMUL(cr4, ci4, LD_PS1(wa3[i-2]), LD_PS1(wa3[i-1]));
        ph[0] = cr4;
        ph[1] = ci4; ph = ph - 3*l1ido + 2;
      }
    }
    if (ido % 2 == 1) return;
  }
  for (k=0; k < l1ido; k+=ido) {
    int i0 = 4*k + ido;
    v4sf c = cc[i0-1], d = cc[i0 + 2*ido-1];
    v4sf a = cc[i0+0], b = cc[i0 + 2*ido+0];
    tr1 = VSUB(c,d);
    tr2 = VADD(c,d);
    ti1 = VADD(b,a);
    ti2 = VSUB(b,a);
    ch[ido-1 + k + 0*l1ido] = VADD(tr2,tr2);
    ch[ido-1 + k + 1*l1ido] = SVMUL(minus_sqrt2, VSUB(ti1, tr1));
    ch[ido-1 + k + 2*l1ido] = VADD(ti2, ti2);
    ch[ido-1 + k + 3*l1ido] = SVMUL(minus_sqrt2, VADD(ti1, tr1));
  }
}

static NEVER_INLINE(void) radb4_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                   const float * RESTRICT wa1, const float * RESTRICT wa2, const float *RESTRICT wa3) {
  radb4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, 0);
} /* radb4 */
#if PF_ON
static NEVER_INLINE(void) radb4_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                      const float * RESTRICT wa1, const float * RESTRICT wa2, const float *RESTRICT wa3) {
  radb4_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, 1);
}
#endif

/* the i-loop of radf5_loop_ps over all k; with even (ido even), each k
   ends with its Nyquist terms, while its blocks are still in cache.
   even is a constant, as pf is, so the odd-ido loop compiles as if
   the tail were not there. */
static ALWAYS_INLINE(void) radf5_rows_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                       const float *wa1, const float *wa2,
                                       const float *wa3, const float *wa4,
                                       const int pf, const int even)
{
  static const float tr11 = (float)0.30901699437494742410;   /* cos(2*pi/5) = (sqrt(5)-1)/4 */
  static const float ti11 = (float)0.95105651629515357212;   /* sin(2*pi/5) */
  static const float tr12 = (float)-0.80901699437494742410;  /* cos(4*pi/5) = -(sqrt(5)+1)/4 */
  static const float ti12 = (float)0.58778525229247312917;   /* sin(4*pi/5) */
  static const float minus_ti11 = (float)-0.95105651629515357212;

  /* System generated locals */
  int cc_offset, ch_offset;

  /* Local variables */
  int i, k, ic;
  v4sf ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, dr4, dr5,
    cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
  int idp2;


#define cc_ref(a_1,a_2,a_3) cc[((a_3)*l1 + (a_2))*ido + a_1]
#define ch_ref(a_1,a_2,a_3) ch[((a_3)*5 + (a_2))*ido + a_1]

  /* Parameter adjustments */
  ch_offset = 1 + ido * 6;
  ch -= ch_offset;
  cc_offset = 1 + ido * (1 + l1);
  cc -= cc_offset;

  idp2 = ido + 2;
  for (k = 1; k <= l1; ++k) {
    for (i = 3; i <= ido; i += 2) {
      ic = idp2 - i;
      PF_PREFETCH(pf, &ch_ref(i - 1, 1, k), 1);
      PF_PREFETCH(pf, &ch_ref(i - 1, 3, k), 1);
      PF_PREFETCH(pf, &ch_ref(i - 1, 5, k), 1);
      PF_PREFETCH(pf, &ch_ref(ic - 1, 2, k), -1);
      PF_PREFETCH(pf, &ch_ref(ic - 1, 4, k), -1);
      dr2 = LD_PS1(wa1[i-3]); di2 = LD_PS1(wa1[i-2]);
      dr3 = LD_PS1(wa2[i-3]); di3 = LD_PS1(wa2[i-2]);
      dr4 = LD_PS1(wa3[i-3]); di4 = LD_PS1(wa3[i-2]);
      dr5 = LD_PS1(wa4[i-3]); di5 = LD_PS1(wa4[i-2]);
      VCPLXMULCONJ(dr2, di2, cc_ref(i-1, k, 2), cc_ref(i, k, 2));
      VCPLXMULCONJ(dr3, di3, cc_ref(i-1, k, 3), cc_ref(i, k, 3));
      VCPLXMULCONJ(dr4, di4, cc_ref(i-1, k, 4), cc_ref(i, k, 4));
      VCPLXMULCONJ(dr5, di5, cc_ref(i-1, k, 5), cc_ref(i, k, 5));
      cr2 = VADD(dr2, dr5);
      ci5 = VSUB(dr5, dr2);
      cr5 = VSUB(di2, di5);
      ci2 = VADD(di2, di5);
      cr3 = VADD(dr3, dr4);
      ci4 = VSUB(dr4, dr3);
      cr4 = VSUB(di3, di4);
      ci3 = VADD(di3, di4);
      ch_ref(i - 1, 1, k) = VADD(cc_ref(i - 1, k, 1), VADD(cr2, cr3));
      ch_ref(i, 1, k) = VSUB(cc_ref(i, k, 1), VADD(ci2, ci3));
      tr2 = VADD(cc_ref(i - 1, k, 1), SVMADD(tr11, cr2, SVMUL(tr12, cr3)));
      ti2 = VSUB(cc_ref(i, k, 1), SVMADD(tr11, ci2, SVMUL(tr12, ci3)));
      tr3 = VADD(cc_ref(i - 1, k, 1), SVMADD(tr12, cr2, SVMUL(tr11, cr3)));
      ti3 = VSUB(cc_ref(i, k, 1), SVMADD(tr12, ci2, SVMUL(tr11, ci3)));
      tr5 = SVMADD(ti11, cr5, SVMUL(ti12, cr4));
      ti5 = SVMADD(ti11, ci5, SVMUL(ti12, ci4));
      tr4 = SVMSUB(ti11, cr4, SVMUL(ti12, cr5));
      ti4 = SVMSUB(ti11, ci4, SVMUL(ti12, ci5));
      ch_ref(i - 1, 3, k) = VSUB(tr2, tr5);
      ch_ref(ic - 1, 2, k) = VADD(tr2, tr5);
      ch_ref(i, 3, k) = VADD(ti2, ti5);
      ch_ref(ic, 2, k) = VSUB(ti5, ti2);
      ch_ref(i - 1, 5, k) = VSUB(tr3, tr4);
      ch_ref(ic - 1, 4, k) = VADD(tr3, tr4);
      ch_ref(i, 5, k) = VADD(ti3, ti4);
      ch_ref(ic, 4, k) = VSUB(ti4, ti3);
    }
    if (even) {
      /* the Nyquist terms x_j (see above radf3_loop_ps): with
         W = exp(-i*pi/5), X0 = sum_j x_j W^j and X1 = sum_j x_j W^(3j) are
           Re X0 = x0 - cos(pi/5) (x4 - x1) + cos(2*pi/5) (x2 - x3),
           Im X0 = -sin(pi/5) (x1 + x4) - sin(2*pi/5) (x2 + x3),
           Re X1 = x0 + cos(2*pi/5) (x4 - x1) - cos(pi/5) (x2 - x3),
           Im X1 = sin(pi/5) (x2 + x3) - sin(2*pi/5) (x1 + x4),
         and X2 = x0 - x1 + x2 - x3 + x4. cos(pi/5) = -tr12, sin(pi/5) = ti12. */
      v4sf x0 = cc_ref(ido, k, 1);
      v4sf a = VSUB(cc_ref(ido, k, 5), cc_ref(ido, k, 2));   /* x4 - x1 */
      v4sf b = VSUB(cc_ref(ido, k, 3), cc_ref(ido, k, 4));   /* x2 - x3 */
      v4sf c = VADD(cc_ref(ido, k, 2), cc_ref(ido, k, 5));   /* x1 + x4 */
      v4sf d = VADD(cc_ref(ido, k, 3), cc_ref(ido, k, 4));   /* x2 + x3 */
      ch_ref(ido, 1, k) = VADD(x0, SVMADD(tr12, a, SVMUL(tr11, b)));
      ch_ref(1, 2, k) = SVMSUB(ti12, c, SVMUL(minus_ti11, d));
      ch_ref(ido, 3, k) = VADD(x0, SVMADD(tr11, a, SVMUL(tr12, b)));
      ch_ref(1, 4, k) = SVMSUB(ti11, c, SVMUL(ti12, d));
      ch_ref(ido, 5, k) = VADD(x0, VADD(a, b));
    }
  }
#undef cc_ref
#undef ch_ref
}

static ALWAYS_INLINE(void) radf5_loop_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                     const float *wa1, const float *wa2, const float *wa3, const float *wa4,
                     const int pf)
{
  static const float tr11 = (float)0.30901699437494742410;   /* cos(2*pi/5) = (sqrt(5)-1)/4 */
  static const float ti11 = (float)0.95105651629515357212;   /* sin(2*pi/5) */
  static const float tr12 = (float)-0.80901699437494742410;  /* cos(4*pi/5) = -(sqrt(5)+1)/4 */
  static const float ti12 = (float)0.58778525229247312917;   /* sin(4*pi/5) */

  /* System generated locals */
  int cc_offset, ch_offset;

  /* Local variables */
  int k;
  v4sf ci4, ci5, cr2, cr3;


#define cc_ref(a_1,a_2,a_3) cc[((a_3)*l1 + (a_2))*ido + a_1]
#define ch_ref(a_1,a_2,a_3) ch[((a_3)*5 + (a_2))*ido + a_1]

  /* Parameter adjustments */
  ch_offset = 1 + ido * 6;
  ch -= ch_offset;
  cc_offset = 1 + ido * (1 + l1);
  cc -= cc_offset;

  /* Function Body */
  for (k = 1; k <= l1; ++k) {
    cr2 = VADD(cc_ref(1, k, 5), cc_ref(1, k, 2));
    ci5 = VSUB(cc_ref(1, k, 5), cc_ref(1, k, 2));
    cr3 = VADD(cc_ref(1, k, 4), cc_ref(1, k, 3));
    ci4 = VSUB(cc_ref(1, k, 4), cc_ref(1, k, 3));
    ch_ref(1, 1, k) = VADD(cc_ref(1, k, 1), VADD(cr2, cr3));
    ch_ref(ido, 2, k) = VADD(cc_ref(1, k, 1), SVMADD(tr11, cr2, SVMUL(tr12, cr3)));
    ch_ref(1, 3, k) = SVMADD(ti11, ci5, SVMUL(ti12, ci4));
    ch_ref(ido, 4, k) = VADD(cc_ref(1, k, 1), SVMADD(tr12, cr2, SVMUL(tr11, cr3)));
    ch_ref(1, 5, k) = SVMSUB(ti11, ci4, SVMUL(ti12, ci5));
    /* printf("pffft: radf5, k=%d ch_ref=%f, ci4=%f\n", k, ch_ref(1, 5, k), ci4); */
  }
  if (ido == 1) {
    return;
  }
  if (ido & 1) {
    radf5_rows_ps(ido, l1, cc + cc_offset, ch + ch_offset, wa1, wa2, wa3, wa4, pf, 0);
  } else {
    radf5_rows_ps(ido, l1, cc + cc_offset, ch + ch_offset, wa1, wa2, wa3, wa4, pf, 1);
  }
#undef cc_ref
#undef ch_ref
}

static NEVER_INLINE(void) radf5_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                     const float *wa1, const float *wa2, const float *wa3, const float *wa4) {
  radf5_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, wa4, 0);
} /* radf5 */
#if PF_ON
static NEVER_INLINE(void) radf5_pf_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                     const float *wa1, const float *wa2, const float *wa3, const float *wa4) {
  radf5_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, wa4, 1);
}
#endif

/* the i-loop of radb5_loop_ps over all k; with even (ido even), each k
   ends with its Nyquist terms, while its blocks are still in cache.
   even is a constant, as pf is, so the odd-ido loop compiles as if
   the tail were not there. */
static ALWAYS_INLINE(void) radb5_rows_ps(int ido, int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch,
                                       const float *wa1, const float *wa2,
                                       const float *wa3, const float *wa4,
                                       const int pf, const int even)
{
  static const float tr11 = (float)0.30901699437494742410;   /* cos(2*pi/5) = (sqrt(5)-1)/4 */
  static const float ti11 = (float)0.95105651629515357212;   /* sin(2*pi/5) */
  static const float tr12 = (float)-0.80901699437494742410;  /* cos(4*pi/5) = -(sqrt(5)+1)/4 */
  static const float ti12 = (float)0.58778525229247312917;   /* sin(4*pi/5) */
  static const float minus_ti11 = (float)-0.95105651629515357212;

  int cc_offset, ch_offset;

  /* Local variables */
  int i, k, ic;
  v4sf ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
    ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
  int idp2;

#define cc_ref(a_1,a_2,a_3) cc[((a_3)*5 + (a_2))*ido + a_1]
#define ch_ref(a_1,a_2,a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]

  /* Parameter adjustments */
  ch_offset = 1 + ido * (1 + l1);
  ch -= ch_offset;
  cc_offset = 1 + ido * 6;
  cc -= cc_offset;

  idp2 = ido + 2;
  for (k = 1; k <= l1; ++k) {
    for (i = 3; i <= ido; i += 2) {
      ic = idp2 - i;
      PF_PREFETCH(pf, &ch_ref(i - 1, k, 1), 1);
      PF_PREFETCH(pf, &ch_ref(i - 1, k, 2), 1);
      PF_PREFETCH(pf, &ch_ref(i - 1, k, 3), 1);
      PF_PREFETCH(pf, &ch_ref(i - 1, k, 4), 1);
      PF_PREFETCH(pf, &ch_ref(i - 1, k, 5), 1);
      ti5 = VADD(cc_ref(i  , 3, k), cc_ref(ic  , 2, k));
      ti2 = VSUB(cc_ref(i  , 3, k), cc_ref(ic  , 2, k));
      ti4 = VADD(cc_ref(i  , 5, k), cc_ref(ic  , 4, k));
      ti3 = VSUB(cc_ref(i  , 5, k), cc_ref(ic  , 4, k));
      tr5 = VSUB(cc_ref(i-1, 3, k), cc_ref(ic-1, 2, k));
      tr2 = VADD(cc_ref(i-1, 3, k), cc_ref(ic-1, 2, k));
      tr4 = VSUB(cc_ref(i-1, 5, k), cc_ref(ic-1, 4, k));
      tr3 = VADD(cc_ref(i-1, 5, k), cc_ref(ic-1, 4, k));
      ch_ref(i - 1, k, 1) = VADD(cc_ref(i-1, 1, k), VADD(tr2, tr3));
      ch_ref(i, k, 1) = VADD(cc_ref(i, 1, k), VADD(ti2, ti3));
      cr2 = VADD(cc_ref(i-1, 1, k), SVMADD(tr11, tr2, SVMUL(tr12, tr3)));
      ci2 = VADD(cc_ref(i  , 1, k), SVMADD(tr11, ti2, SVMUL(tr12, ti3)));
      cr3 = VADD(cc_ref(i-1, 1, k), SVMADD(tr12, tr2, SVMUL(tr11, tr3)));
      ci3 = VADD(cc_ref(i  , 1, k), SVMADD(tr12, ti2, SVMUL(tr11, ti3)));
      cr5 = SVMADD(ti11, tr5, SVMUL(ti12, tr4));
      ci5 = SVMADD(ti11, ti5, SVMUL(ti12, ti4));
      cr4 = SVMSUB(ti11, tr4, SVMUL(ti12, tr5));
      ci4 = SVMSUB(ti11, ti4, SVMUL(ti12, ti5));
      dr3 = VSUB(cr3, ci4);
      dr4 = VADD(cr3, ci4);
      di3 = VADD(ci3, cr4);
      di4 = VSUB(ci3, cr4);
      dr5 = VADD(cr2, ci5);
      dr2 = VSUB(cr2, ci5);
      di5 = VSUB(ci2, cr5);
      di2 = VADD(ci2, cr5);
      VCPLXMUL(dr2, di2, LD_PS1(wa1[i-3]), LD_PS1(wa1[i-2]));
      VCPLXMUL(dr3, di3, LD_PS1(wa2[i-3]), LD_PS1(wa2[i-2]));
      VCPLXMUL(dr4, di4, LD_PS1(wa3[i-3]), LD_PS1(wa3[i-2]));
      VCPLXMUL(dr5, di5, LD_PS1(wa4[i-3]), LD_PS1(wa4[i-2]));

      ch_ref(i-1, k, 2) = dr2; ch_ref(i, k, 2) = di2;
      ch_ref(i-1, k, 3) = dr3; ch_ref(i, k, 3) = di3;
      ch_ref(i-1, k, 4) = dr4; ch_ref(i, k, 4) = di4;
      ch_ref(i-1, k, 5) = dr5; ch_ref(i, k, 5) = di5;
    }
    if (even) {
      /* the Nyquist terms, from X0, X1 and X2 (see above radf3_loop_ps):
           x0 = 2 Re X0 + 2 Re X1 + X2,
           x1, x4 = -+(X2 - 2 cos(pi/5) Re X0 + 2 cos(2*pi/5) Re X1)
                    - 2 sin(pi/5) Im X0 - 2 sin(2*pi/5) Im X1,
           x2, x3 = +-(X2 + 2 cos(2*pi/5) Re X0 - 2 cos(pi/5) Re X1)
                    - 2 sin(2*pi/5) Im X0 + 2 sin(pi/5) Im X1 */
      v4sf q  = cc_ref(ido, 5, k);
      v4sf r0 = VADD(cc_ref(ido, 1, k), cc_ref(ido, 1, k));
      v4sf i0 = VADD(cc_ref(1, 2, k), cc_ref(1, 2, k));
      v4sf r1 = VADD(cc_ref(ido, 3, k), cc_ref(ido, 3, k));
      v4sf i1 = VADD(cc_ref(1, 4, k), cc_ref(1, 4, k));
      v4sf p1 = VADD(q, SVMADD(tr12, r0, SVMUL(tr11, r1)));
      v4sf p2 = VADD(q, SVMADD(tr11, r0, SVMUL(tr12, r1)));
      v4sf s1 = SVMSUB(ti12, i0, SVMUL(minus_ti11, i1));
      v4sf s2 = SVMSUB(ti11, i0, SVMUL(ti12, i1));
      ch_ref(ido, k, 1) = VADD(q, VADD(r0, r1));
      ch_ref(ido, k, 2) = VSUB(s1, p1);
      ch_ref(ido, k, 3) = VADD(s2, p2);
      ch_ref(ido, k, 4) = VSUB(s2, p2);
      ch_ref(ido, k, 5) = VADD(s1, p1);
    }
  }
#undef cc_ref
#undef ch_ref
}

static ALWAYS_INLINE(void) radb5_loop_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf *RESTRICT ch,
                  const float *wa1, const float *wa2, const float *wa3, const float *wa4,
                  const int pf)
{
  static const float tr11 = (float)0.30901699437494742410;   /* cos(2*pi/5) = (sqrt(5)-1)/4 */
  static const float ti11 = (float)0.95105651629515357212;   /* sin(2*pi/5) */
  static const float tr12 = (float)-0.80901699437494742410;  /* cos(4*pi/5) = -(sqrt(5)+1)/4 */
  static const float ti12 = (float)0.58778525229247312917;   /* sin(4*pi/5) */

  int cc_offset, ch_offset;

  /* Local variables */
  int k;
  v4sf ci4, ci5, cr2, cr3, ti4, ti5, tr2, tr3;

#define cc_ref(a_1,a_2,a_3) cc[((a_3)*5 + (a_2))*ido + a_1]
#define ch_ref(a_1,a_2,a_3) ch[((a_3)*l1 + (a_2))*ido + a_1]

  /* Parameter adjustments */
  ch_offset = 1 + ido * (1 + l1);
  ch -= ch_offset;
  cc_offset = 1 + ido * 6;
  cc -= cc_offset;

  /* Function Body */
  for (k = 1; k <= l1; ++k) {
    ti5 = VADD(cc_ref(1, 3, k), cc_ref(1, 3, k));
    ti4 = VADD(cc_ref(1, 5, k), cc_ref(1, 5, k));
    tr2 = VADD(cc_ref(ido, 2, k), cc_ref(ido, 2, k));
    tr3 = VADD(cc_ref(ido, 4, k), cc_ref(ido, 4, k));
    ch_ref(1, k, 1) = VADD(cc_ref(1, 1, k), VADD(tr2, tr3));
    cr2 = VADD(cc_ref(1, 1, k), SVMADD(tr11, tr2, SVMUL(tr12, tr3)));
    cr3 = VADD(cc_ref(1, 1, k), SVMADD(tr12, tr2, SVMUL(tr11, tr3)));
    ci5 = SVMADD(ti11, ti5, SVMUL(ti12, ti4));
    ci4 = SVMSUB(ti11, ti4, SVMUL(ti12, ti5));
    ch_ref(1, k, 2) = VSUB(cr2, ci5);
    ch_ref(1, k, 3) = VSUB(cr3, ci4);
    ch_ref(1, k, 4) = VADD(cr3, ci4);
    ch_ref(1, k, 5) = VADD(cr2, ci5);
  }
  if (ido == 1) {
    return;
  }
  if (ido & 1) {
    radb5_rows_ps(ido, l1, cc + cc_offset, ch + ch_offset, wa1, wa2, wa3, wa4, pf, 0);
  } else {
    radb5_rows_ps(ido, l1, cc + cc_offset, ch + ch_offset, wa1, wa2, wa3, wa4, pf, 1);
  }
#undef cc_ref
#undef ch_ref
}

static NEVER_INLINE(void) radb5_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf *RESTRICT ch,
                  const float *wa1, const float *wa2, const float *wa3, const float *wa4) {
  radb5_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, wa4, 0);
} /* radb5 */
#if PF_ON
static NEVER_INLINE(void) radb5_pf_ps(int ido, int l1, const v4sf *RESTRICT cc, v4sf *RESTRICT ch,
                  const float *wa1, const float *wa2, const float *wa3, const float *wa4) {
  radb5_loop_ps(ido, l1, cc, ch, wa1, wa2, wa3, wa4, 1);
}
#endif

/*
  Real radix-8 passes, forward and backward, only for ido == 1.

  rffti1_ps puts a factor 8 only last, in place of the lone 2 and a 4.
  That factor's pass is the first of rfftf1_ps and the last of rfftb1_ps,
  where ido == 1, so it takes no twiddles: one pass fewer over the data,
  and the radix-2 pass that goes away is the one with the largest ido,
  with twiddles. For each k, the 8 legs are the real vectors cc[k + j*l1]
  (forward input, backward output), and the half-complex block
  X0, Re X1, Im X1, Re X2, Im X2, Re X3, Im X3, X4 at 8*k (forward
  output, backward input), as radf4_ps lays out X0, Re X1, Im X1, X2.
*/
static NEVER_INLINE(void) radf8_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch) {
  static const float minus_hsqt2 = (float)-M_SQRT1_2;  /* -sin(pi/4) = -1/sqrt(2) */
  const v4sf * RESTRICT cc_end = cc + l1;
  while (cc < cc_end) {
    v4sf x0 = cc[0],    x1 = cc[1*l1], x2 = cc[2*l1], x3 = cc[3*l1];
    v4sf x4 = cc[4*l1], x5 = cc[5*l1], x6 = cc[6*l1], x7 = cc[7*l1];
    /* the even outputs are the DFT-4 of the sums at stride 4 */
    v4sf a0 = VADD(x0, x4), a1 = VADD(x1, x5), a2 = VADD(x2, x6), a3 = VADD(x3, x7);
    /* the odd ones that of the differences times W8^j */
    v4sf b0 = VSUB(x0, x4), b1 = VSUB(x1, x5), b2 = VSUB(x2, x6), b3 = VSUB(x3, x7);
    v4sf s02 = VADD(a0, a2), s13 = VADD(a1, a3);
    v4sf t1 = SVMUL(minus_hsqt2, VSUB(b3, b1));
    v4sf t2 = SVMUL(minus_hsqt2, VADD(b1, b3));
    ch[0] = VADD(s02, s13);   /* X0 */
    ch[1] = VADD(b0, t1);     /* Re X1 */
    ch[2] = VSUB(t2, b2);     /* Im X1 */
    ch[3] = VSUB(a0, a2);     /* Re X2 */
    ch[4] = VSUB(a3, a1);     /* Im X2 */
    ch[5] = VSUB(b0, t1);     /* Re X3 */
    ch[6] = VADD(b2, t2);     /* Im X3 */
    ch[7] = VSUB(s02, s13);   /* X4 */
    cc += 1; ch += 8;
  }
}

/* The backward passes write 8 (radb6: 6) output streams, l1 vectors
   apart, which they prefetch on x86 as passf8_ps does (see PF_PREFETCH);
   the forward ones write a single stream, 8 (6) vectors per k, and read
   the 8 (6), which the hardware prefetchers follow. */
static ALWAYS_INLINE(void) radb8_loop_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const int pf) {
  static const float two = 2.f;
  static const float sqrt2 = (float)M_SQRT2;
  v4sf * RESTRICT ch_end = ch + l1;
  while (ch < ch_end) {
    PF_PREFETCH(pf, ch + 0*l1, 1);
    PF_PREFETCH(pf, ch + 1*l1, 1);
    PF_PREFETCH(pf, ch + 2*l1, 1);
    PF_PREFETCH(pf, ch + 3*l1, 1);
    PF_PREFETCH(pf, ch + 4*l1, 1);
    PF_PREFETCH(pf, ch + 5*l1, 1);
    PF_PREFETCH(pf, ch + 6*l1, 1);
    PF_PREFETCH(pf, ch + 7*l1, 1);
    v4sf x0 = cc[0], r1 = cc[1], i1 = cc[2], r2 = cc[3];
    v4sf i2 = cc[4], r3 = cc[5], i3 = cc[6], x4 = cc[7];
    /* x_j = e_j + o_j, x_{j+4} = e_j - o_j: e_j from X0, X2, X4 and
       o_j from X1, X3 */
    v4sf s = VADD(x0, x4), d = VSUB(x0, x4);
    v4sf r2x2 = SVMUL(two, r2), i2x2 = SVMUL(two, i2);
    v4sf e0 = VADD(s, r2x2), e2 = VSUB(s, r2x2);
    v4sf e1 = VSUB(d, i2x2), e3 = VADD(d, i2x2);
    v4sf p = VSUB(r1, r3), q = VADD(i1, i3);
    v4sf o0 = SVMUL(two, VADD(r1, r3));
    v4sf o2 = SVMUL(two, VSUB(i3, i1));
    v4sf o1 = SVMUL(sqrt2, VSUB(p, q));
    v4sf o3 = SVMUL(sqrt2, VADD(p, q));   /* minus o_3 */
    ch[0*l1] = VADD(e0, o0);
    ch[4*l1] = VSUB(e0, o0);
    ch[1*l1] = VADD(e1, o1);
    ch[5*l1] = VSUB(e1, o1);
    ch[2*l1] = VADD(e2, o2);
    ch[6*l1] = VSUB(e2, o2);
    ch[3*l1] = VSUB(e3, o3);
    ch[7*l1] = VADD(e3, o3);
    cc += 8; ch += 1;
  }
}
static NEVER_INLINE(void) radb8_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch) {
  radb8_loop_ps(l1, cc, ch, 0);
}
#if PF_ON
static NEVER_INLINE(void) radb8_pf_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch) {
  radb8_loop_ps(l1, cc, ch, 1);
}
#endif

/*
  Real radix-6 passes, forward and backward, only for ido == 1.

  rffti1_ps puts a factor 6 only last, in place of the lone 2 and a 3,
  where the transform is small (see there). The layout is that of the 8
  above: legs cc[k + j*l1], half-complex block X0, Re X1, Im X1, Re X2,
  Im X2, X3 at 6*k.

  2 and 3 are coprime, so the 6-point DFT is a 2-point DFT on each of the
  pairs (x0, x3), (x2, x5), (x4, x1) followed by a 3-point DFT on the
  sums and one on the differences, with no twiddles in between (Good-Thomas:
  input index 3*j1 + 2*j2, output index 3*k1 + 4*k2, mod 6). The sums give
  X0, X4, X2 and the differences X3, X1, X5, of which only X0..X3 are
  stored. That is 12 adds, 2 multiply-adds and 2 multiplies each way. The
  direct real 6-point DFT with the symmetric sums x_j +- x_{6-j} costs the
  same; the Cooley-Tukey 2x3 split would add twiddles by W6^j.
*/
static NEVER_INLINE(void) radf6_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch) {
  static const float taur = (float)-0.5;                     /* cos(2*pi/3) */
  static const float taui = (float)0.86602540378443864676;   /* sin(2*pi/3) = sqrt(3)/2 */
  const v4sf * RESTRICT cc_end = cc + l1;
  while (cc < cc_end) {
    v4sf x0 = cc[0],    x1 = cc[1*l1], x2 = cc[2*l1];
    v4sf x3 = cc[3*l1], x4 = cc[4*l1], x5 = cc[5*l1];
    /* the 2-point DFTs, in Good-Thomas order */
    v4sf u0 = VADD(x0, x3), u1 = VADD(x2, x5), u2 = VADD(x4, x1);
    v4sf v0 = VSUB(x0, x3), v1 = VSUB(x2, x5), v2 = VSUB(x4, x1);
    /* the 3-point DFTs: of u into X0, X4 = conj(X2); of v into X3, X1 */
    v4sf su = VADD(u1, u2), sv = VADD(v1, v2);
    ch[0] = VADD(u0, su);                      /* X0 */
    ch[1] = SVMADD(taur, sv, v0);              /* Re X1 */
    ch[2] = SVMUL(taui, VSUB(v2, v1));         /* Im X1 */
    ch[3] = SVMADD(taur, su, u0);              /* Re X2 */
    ch[4] = SVMUL(taui, VSUB(u1, u2));         /* Im X2 */
    ch[5] = VADD(v0, sv);                      /* X3 */
    cc += 1; ch += 6;
  }
}

static ALWAYS_INLINE(void) radb6_loop_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch, const int pf) {
  static const float two = 2.f;
  static const float sqrt3 = (float)1.7320508075688772935;  /* 2*sin(2*pi/3) = sqrt(3) */
  v4sf * RESTRICT ch_end = ch + l1;
  while (ch < ch_end) {
    PF_PREFETCH(pf, ch + 0*l1, 1);
    PF_PREFETCH(pf, ch + 1*l1, 1);
    PF_PREFETCH(pf, ch + 2*l1, 1);
    PF_PREFETCH(pf, ch + 3*l1, 1);
    PF_PREFETCH(pf, ch + 4*l1, 1);
    PF_PREFETCH(pf, ch + 5*l1, 1);
    v4sf x0 = cc[0], r1 = cc[1], i1 = cc[2];
    v4sf r2 = cc[3], i2 = cc[4], x3 = cc[5];
    /* the inverse 3-point DFTs, of X0, X4, X2 into u and of X3, X1, X5
       into v, then the 2-point ones into the legs in Good-Thomas order */
    v4sf u0 = SVMADD(two, r2, x0), a = VSUB(x0, r2), b = SVMUL(sqrt3, i2);
    v4sf v0 = SVMADD(two, r1, x3), c = VSUB(x3, r1), d = SVMUL(sqrt3, i1);
    v4sf u1 = VADD(a, b), u2 = VSUB(a, b);
    v4sf v1 = VSUB(c, d), v2 = VADD(c, d);
    ch[0*l1] = VADD(u0, v0);
    ch[3*l1] = VSUB(u0, v0);
    ch[2*l1] = VADD(u1, v1);
    ch[5*l1] = VSUB(u1, v1);
    ch[4*l1] = VADD(u2, v2);
    ch[1*l1] = VSUB(u2, v2);
    cc += 6; ch += 1;
  }
}
static NEVER_INLINE(void) radb6_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch) {
  radb6_loop_ps(l1, cc, ch, 0);
}
#if PF_ON
static NEVER_INLINE(void) radb6_pf_ps(int l1, const v4sf * RESTRICT cc, v4sf * RESTRICT ch) {
  radb6_loop_ps(l1, cc, ch, 1);
}
#endif

static NEVER_INLINE(v4sf *) rfftf1_ps(int n, const v4sf *input_readonly, v4sf *work1, v4sf *work2,
                                      const float *wa, const int *ifac) {
  v4sf *in  = (v4sf*)input_readonly;
  v4sf *out = (in == work2 ? work1 : work2);
  int nf = ifac[1], k1;
  int l2 = n;
  int iw = n-1;
  const int pf = PF_USE_PREFETCH(n); /* the same for every pass */
  (void)pf;
  assert(in != out && work1 != work2);
  for (k1 = 1; k1 <= nf; ++k1) {
    int kh = nf - k1;
    int ip = ifac[kh + 2];
    int l1 = l2 / ip;
    int ido = n / l2;
    iw -= (ip - 1)*ido;
    switch (ip) {
      case 5: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        int ix4 = ix3 + ido;
        PF_CALL(pf, radf5, (ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]));
      } break;
      case 8:
        assert(ido == 1);
        radf8_ps(l1, in, out);
        break;
      case 6:
        assert(ido == 1);
        radf6_ps(l1, in, out);
        break;
      case 4: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        PF_CALL(pf, radf4, (ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3]));
      } break;
      case 3: {
        int ix2 = iw + ido;
        PF_CALL(pf, radf3, (ido, l1, in, out, &wa[iw], &wa[ix2]));
      } break;
      case 2:
        PF_CALL(pf, radf2, (ido, l1, in, out, &wa[iw]));
        break;
      default:
        assert(0);
        break;
    }
    l2 = l1;
    if (out == work2) {
      out = work1; in = work2;
    } else {
      out = work2; in = work1;
    }
  }
  return in; /* this is in fact the output .. */
} /* rfftf1 */

static NEVER_INLINE(v4sf *) rfftb1_ps(int n, const v4sf *input_readonly, v4sf *work1, v4sf *work2,
                                      const float *wa, const int *ifac) {
  v4sf *in  = (v4sf*)input_readonly;
  v4sf *out = (in == work2 ? work1 : work2);
  int nf = ifac[1], k1;
  int l1 = 1;
  int iw = 0;
  const int pf = PF_USE_PREFETCH(n); /* the same for every pass */
  (void)pf;
  assert(in != out);
  for (k1=1; k1<=nf; k1++) {
    int ip = ifac[k1 + 1];
    int l2 = ip*l1;
    int ido = n / l2;
    switch (ip) {
      case 5: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        int ix4 = ix3 + ido;
        PF_CALL(pf, radb5, (ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]));
      } break;
      case 8:
        assert(ido == 1);
        PF_CALL(pf, radb8, (l1, in, out));
        break;
      case 6:
        assert(ido == 1);
        PF_CALL(pf, radb6, (l1, in, out));
        break;
      case 4: {
        int ix2 = iw + ido;
        int ix3 = ix2 + ido;
        PF_CALL(pf, radb4, (ido, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3]));
      } break;
      case 3: {
        int ix2 = iw + ido;
        PF_CALL(pf, radb3, (ido, l1, in, out, &wa[iw], &wa[ix2]));
      } break;
      case 2:
        PF_CALL(pf, radb2, (ido, l1, in, out, &wa[iw]));
        break;
      default:
        assert(0);
        break;
    }
    l1 = l2;
    iw += (ip - 1)*ido;

    if (out == work2) {
      out = work1; in = work2;
    } else {
      out = work2; in = work1;
    }
  }
  return in; /* this is in fact the output .. */
}

#define IFAC_MAX_SIZE 25 /* max number of integer factors for the decomposition, +2 */
static int decompose(int n, int *ifac, const int *ntryh) {
  int nl = n, nf = 0, i, j = 0;
  for (j=0; ntryh[j]; ++j) {
    int ntry = ntryh[j];
    while (nl != 1) {
      int nq = nl / ntry;
      int nr = nl - ntry * nq;
      if (nr == 0) {
        assert(2 + nf < IFAC_MAX_SIZE);
        ifac[2+nf++] = ntry;
        nl = nq;
        if (ntry == 2 && nf != 1) {
          for (i = 2; i <= nf; ++i) {
            int ib = nf - i + 2;
            ifac[ib + 1] = ifac[ib];
          }
          ifac[2] = 2;
        }
      } else break;
    }
  }
  ifac[0] = n;
  ifac[1] = nf;
  return nf;
}



/* the largest transform, in bytes per array of n vectors, that takes a
   radix-6 instead of a radix-8 pass where it can (see rffti1_ps) */
#define RADIX6_MAX_BYTES 32768

static void rffti1_ps(int n, float *wa, int *ifac)
{
  static const int ntryh[] = { 4,2,3,5,0 };
  int k1, j, ii;

  int nf = decompose(n,ifac,ntryh);

  /* decompose() puts the lone factor 2, if there is one, first, then the
     4s, then the 3s and 5s: [2, 4, ..., 3, ..., 5, ...]. Merge the 2 and a
     4 into a radix-8 pass at the end, where ido == 1 and radf8_ps/radb8_ps
     need no twiddles, as cffti1_ps does: [4, ..., 3, ..., 5, ..., 8]. The
     3s and 5s, which {4,2,3,5} puts last so that their ido is odd, then
     get an even ido, which their passes handle (see radf3_loop_ps). The
     twiddle table is filled the same way for an even ido: (ido-1)/2 pairs
     per leg, i.e. none for the element ido-1, which needs none.

     Where there is a 3 and an array of the transform takes at most
     RADIX6_MAX_BYTES, merge the 2 and a 3 into a radix-6 pass at the end
     instead: [4, ..., 3 (one fewer), ..., 5, ..., 6]. In cache that was
     faster than the radix-8 merge (transform length N = 384, i.e.
     n = 96: f32 8-10% on two Xeons, 1% on an Apple M2; f64 3-9%), out of
     cache slower (N = 24576, f32 and f64: 1-5%). */
  if (nf >= 2 && ifac[2] == 2 && ifac[3] == 4) {
    int i3 = 0;
    if ((size_t)n * sizeof(v4sf) <= RADIX6_MAX_BYTES)
      for (k1 = 4; k1 <= nf+1; ++k1) if (ifac[k1] == 3) { i3 = k1; break; }
    if (i3) {
      for (k1 = 2; k1 < i3-1; ++k1) ifac[k1] = ifac[k1+1];
      for (k1 = i3-1; k1 < nf; ++k1) ifac[k1] = ifac[k1+2];
      ifac[nf] = 6;
    } else {
      for (k1 = 2; k1 < nf; ++k1) ifac[k1] = ifac[k1+2];
      ifac[nf] = 8;
    }
    ifac[1] = --nf;
  }

  float argh = (2*(float)M_PI) / n;
  int is = 0;
  int nfm1 = nf - 1;
  int l1 = 1;
  for (k1 = 1; k1 <= nfm1; k1++) {
    int ip = ifac[k1 + 1];
    int ld = 0;
    int l2 = l1*ip;
    int ido = n / l2;
    int ipm = ip - 1;
    for (j = 1; j <= ipm; ++j) {
      float argld;
      int i = is, fi=0;
      ld += l1;
      argld = ld*argh;
      for (ii = 3; ii <= ido; ii += 2) {
        i += 2;
        fi += 1;
        wa[i - 2] = FUNC_COS(fi*argld);
        wa[i - 1] = FUNC_SIN(fi*argld);
      }
      is += ido;
    }
    l1 = l2;
  }
} /* rffti1 */

static void cffti1_ps(int n, float *wa, int *ifac)
{
  static const int ntryh[] = { 5,3,4,2,0 };
  int k1, j, ii;

  int nf = decompose(n,ifac,ntryh);

  /* decompose() puts the lone factor 2, if there is one, first and the 4s
     last: [2, ..., 4]. Merge the 2 into the last 4 as one radix-8 pass:
     one pass fewer over the data, and the radix-2 pass that goes away is
     one with twiddles. The last pass has ido == 1 (idot == 2 in
     cfftf1_ps), hence no twiddles, so passf8_ps holds only its 8 inputs
     and temporaries. Factoring by {5,3,8,4,2} instead, which makes all
     but the last radix-8 pass one with twiddles (n = 16384: 8, 8, 8, 8),
     made n = 16384 15-27% slower (2ce8a89), likely because a radix-8
     butterfly with 7 twiddle multiplies holds about 40 vectors live.
     Where there is also a 3, the 2 and a 3 become a radix-6 pass at the
     end instead, as twiddle-free as the radix-8 one (passf6_ps):
     [2, 5.., 3, 3.., 4..] -> [5.., 3.., 4.., 6]. The radix-2 and a
     radix-3 pass with twiddles go away, and all the 4s stay. The other
     factors keep their order, so the 3s and 5s still run at an ido of
     at least 6 (idot > 2), as passf3_ps and passf5_ps assume. The loop
     below writes the 6 (ido == 1) only trivial twiddles, which passf6_ps
     does not read, as for the 8; the twiddles of all passes still take
     2*n floats, as sum (ip-1)*ido = n-1 for any order of the factors. */
  if (nf >= 2 && ifac[2] == 2) {
    int k3 = 0;
    for (k1 = 3; k1 <= nf+1; ++k1)
      if (ifac[k1] == 3) k3 = k1;
    if (k3) {
      int m = 2;
      for (k1 = 3; k1 <= nf+1; ++k1)
        if (k1 != k3) ifac[m++] = ifac[k1];
      ifac[m] = 6;  /* m == nf */
      ifac[1] = --nf;
    } else if (ifac[nf+1] == 4) {
      for (k1 = 2; k1 <= nf; ++k1) ifac[k1] = ifac[k1+1];
      ifac[nf] = 8;
      ifac[1] = --nf;
    }
  }

  float argh = (2*(float)M_PI) / n;
  int i = 1;
  int l1 = 1;
  for (k1=1; k1<=nf; k1++) {
    int ip = ifac[k1+1];
    int ld = 0;
    int l2 = l1*ip;
    int ido = n / l2;
    int idot = ido + ido + 2;
    int ipm = ip - 1;
    for (j=1; j<=ipm; j++) {
      float argld;
      int i1 = i, fi = 0;
      wa[i-1] = 1;
      wa[i] = 0;
      ld += l1;
      argld = ld*argh;
      for (ii = 4; ii <= idot; ii += 2) {
        i += 2;
        fi += 1;
        wa[i-1] = FUNC_COS(fi*argld);
        wa[i] = FUNC_SIN(fi*argld);
      }
      if (ip > 16) {
        wa[i1-1] = wa[i-1];
        wa[i1] = wa[i];
      }
    }
    l1 = l2;
  }
} /* cffti1 */


static v4sf *cfftf1_ps(int n, const v4sf *input_readonly, v4sf *work1, v4sf *work2, const float *wa, const int *ifac, int isign) {
  v4sf *in  = (v4sf*)input_readonly;
  v4sf *out = (in == work2 ? work1 : work2);
  int nf = ifac[1], k1;
  int l1 = 1;
  int iw = 0;
  const int pf = PF_USE_PREFETCH(2*n); /* the same for every pass */
  (void)pf;
  assert(in != out && work1 != work2);
  for (k1=2; k1<=nf+1; k1++) {
    int ip = ifac[k1];
    int l2 = ip*l1;
    int ido = n / l2;
    int idot = ido + ido;
    switch (ip) {
      case 8: {
        PF_CALL(pf, passf8, (idot, l1, in, out, isign));
      } break;
      case 6: {
        PF_CALL(pf, passf6, (idot, l1, in, out, isign));
      } break;
      case 5: {
        int ix2 = iw + idot;
        int ix3 = ix2 + idot;
        int ix4 = ix3 + idot;
        passf5_ps(idot, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4], isign);
      } break;
      case 4: {
        int ix2 = iw + idot;
        int ix3 = ix2 + idot;
        PF_CALL(pf, passf4, (idot, l1, in, out, &wa[iw], &wa[ix2], &wa[ix3], isign));
      } break;
      case 2: {
        passf2_ps(idot, l1, in, out, &wa[iw], isign);
      } break;
      case 3: {
        int ix2 = iw + idot;
        PF_CALL(pf, passf3, (idot, l1, in, out, &wa[iw], &wa[ix2], isign));
      } break;
      default:
        assert(0);
    }
    l1 = l2;
    iw += (ip - 1)*idot;
    if (out == work2) {
      out = work1; in = work2;
    } else {
      out = work2; in = work1;
    }
  }

  return in; /* this is in fact the output .. */
}

#undef PF_PREFETCH
#undef PF_USE_PREFETCH
#undef PF_CALL
#undef PF_ON


struct SETUP_STRUCT {
  int     N;
  int     Ncvec;  /* nb of complex simd vectors (N/4 if PFFFT_COMPLEX, N/8 if PFFFT_REAL) */
  /* hold the decomposition into small integers of N */
  int ifac[IFAC_MAX_SIZE]; /* N, number of factors, factors (admitted values: 2, 3, 4, 5, 6 or 8) */
  pffft_transform_t transform;
  v4sf *data;     /* allocated room for twiddle coefs */
  float *e;       /* points into 'data', N/4*3 elements */
  float *twiddle; /* points into 'data', N/4 elements */
};

SETUP_STRUCT *FUNC_NEW_SETUP(int N, pffft_transform_t transform) {
  SETUP_STRUCT *s = 0;
  int k, m;
  /* validate N for negative values or potential int overflow */
  if (N < 0) {
    return s;
  }
  if (N > (1<<26)) {
    /* higher values of N will make you enter in the integer overflow world... */
    assert(0);
    return s;
  }
  /* unfortunately, the fft size must be a multiple of 16 for complex FFTs
     and 32 for real FFTs -- a lot of stuff would need to be rewritten to
     handle other cases (or maybe just switch to a scalar fft, I don't know..) */
  if (transform == PFFFT_REAL)    { if ((N%(2*SIMD_SZ*SIMD_SZ)) || N<=0) return s; }
  if (transform == PFFFT_COMPLEX) { if ((N%(  SIMD_SZ*SIMD_SZ)) || N<=0) return s; }
  s = (SETUP_STRUCT*)malloc(sizeof(SETUP_STRUCT));
  /* assert((N % 32) == 0); */
  s->N = N;
  s->transform = transform;
  /* nb of complex simd vectors */
  s->Ncvec = (transform == PFFFT_REAL ? N/2 : N)/SIMD_SZ;
  s->data = (v4sf*)FUNC_ALIGNED_MALLOC(2*s->Ncvec * sizeof(v4sf));
  s->e = (float*)s->data;
  s->twiddle = (float*)(s->data + (2*s->Ncvec*(SIMD_SZ-1))/SIMD_SZ);

  for (k=0; k < s->Ncvec; ++k) {
    int i = k/SIMD_SZ;
    int j = k%SIMD_SZ;
    for (m=0; m < SIMD_SZ-1; ++m) {
      float A = -2*(float)M_PI*(m+1)*k / N;
      s->e[(2*(i*3 + m) + 0) * SIMD_SZ + j] = FUNC_COS(A);
      s->e[(2*(i*3 + m) + 1) * SIMD_SZ + j] = FUNC_SIN(A);
    }
  }

  if (transform == PFFFT_REAL) {
    rffti1_ps(N/SIMD_SZ, s->twiddle, s->ifac);
  } else {
    cffti1_ps(N/SIMD_SZ, s->twiddle, s->ifac);
  }

  /* check that N is decomposable with allowed prime factors */
  for (k=0, m=1; k < s->ifac[1]; ++k) { m *= s->ifac[2+k]; }
  if (m != N/SIMD_SZ) {
    FUNC_DESTROY(s); s = 0;
  }

  return s;
}


void FUNC_DESTROY(SETUP_STRUCT *s) {
  if (!s)
    return;
  FUNC_ALIGNED_FREE(s->data);
  free(s);
}

#if ( SIMD_SZ == 4 )    /* !defined(PFFFT_SIMD_DISABLE) */

/* [0 0 1 2 3 4 5 6 7 8] -> [0 8 7 6 5 4 3 2 1] */
static void reversed_copy(int N, const v4sf *in, int in_stride, v4sf *out) {
  v4sf g0, g1;
  int k;
  INTERLEAVE2(in[0], in[1], g0, g1); in += in_stride;

  *--out = VSWAPHL(g0, g1); /* [g0l, g0h], [g1l g1h] -> [g1l, g0h] */
  for (k=1; k < N; ++k) {
    v4sf h0, h1;
    INTERLEAVE2(in[0], in[1], h0, h1); in += in_stride;
    *--out = VSWAPHL(g1, h0);
    *--out = VSWAPHL(h0, h1);
    g1 = h1;
  }
  *--out = VSWAPHL(g1, g0);
}

static void unreversed_copy(int N, const v4sf *in, v4sf *out, int out_stride) {
  v4sf g0, g1, h0, h1;
  int k;
  g0 = g1 = in[0]; ++in;
  for (k=1; k < N; ++k) {
    h0 = *in++; h1 = *in++;
    g1 = VSWAPHL(g1, h0);
    h0 = VSWAPHL(h0, h1);
    UNINTERLEAVE2(h0, g1, out[0], out[1]); out += out_stride;
    g1 = h1;
  }
  h0 = *in++; h1 = g0;
  g1 = VSWAPHL(g1, h0);
  h0 = VSWAPHL(h0, h1);
  UNINTERLEAVE2(h0, g1, out[0], out[1]);
}

void FUNC_ZREORDER(const SETUP_STRUCT *setup, const float *in, float *out, pffft_direction_t direction) {
  int k, N = setup->N, Ncvec = setup->Ncvec;
  const v4sf *vin = (const v4sf*)in;
  v4sf *vout = (v4sf*)out;
  assert(in != out);
  if (setup->transform == PFFFT_REAL) {
    int dk = N/32;
    if (direction == PFFFT_FORWARD) {
      for (k=0; k < dk; ++k) {
        INTERLEAVE2(vin[k*8 + 0], vin[k*8 + 1], vout[2*(0*dk + k) + 0], vout[2*(0*dk + k) + 1]);
        INTERLEAVE2(vin[k*8 + 4], vin[k*8 + 5], vout[2*(2*dk + k) + 0], vout[2*(2*dk + k) + 1]);
      }
      reversed_copy(dk, vin+2, 8, (v4sf*)(out + N/2));
      reversed_copy(dk, vin+6, 8, (v4sf*)(out + N));
    } else {
      for (k=0; k < dk; ++k) {
        UNINTERLEAVE2(vin[2*(0*dk + k) + 0], vin[2*(0*dk + k) + 1], vout[k*8 + 0], vout[k*8 + 1]);
        UNINTERLEAVE2(vin[2*(2*dk + k) + 0], vin[2*(2*dk + k) + 1], vout[k*8 + 4], vout[k*8 + 5]);
      }
      unreversed_copy(dk, (v4sf*)(in + N/4), (v4sf*)(out + N - 6*SIMD_SZ), -8);
      unreversed_copy(dk, (v4sf*)(in + 3*N/4), (v4sf*)(out + N - 2*SIMD_SZ), -8);
    }
  } else {
    if (direction == PFFFT_FORWARD) {
      for (k=0; k < Ncvec; ++k) {
        int kk = (k/4) + (k%4)*(Ncvec/4);
        INTERLEAVE2(vin[k*2], vin[k*2+1], vout[kk*2], vout[kk*2+1]);
      }
    } else {
      for (k=0; k < Ncvec; ++k) {
        int kk = (k/4) + (k%4)*(Ncvec/4);
        UNINTERLEAVE2(vin[kk*2], vin[kk*2+1], vout[k*2], vout[k*2+1]);
      }
    }
  }
}

void FUNC_CPLX_FINALIZE(int Ncvec, const v4sf *in, v4sf *out, const v4sf *e) {
  int k, dk = Ncvec/SIMD_SZ; /* number of 4x4 matrix blocks */
  v4sf r0, i0, r1, i1, r2, i2, r3, i3;
  v4sf sr0, dr0, sr1, dr1, si0, di0, si1, di1;
  assert(in != out);
  for (k=0; k < dk; ++k) {
    r0 = in[8*k+0]; i0 = in[8*k+1];
    r1 = in[8*k+2]; i1 = in[8*k+3];
    r2 = in[8*k+4]; i2 = in[8*k+5];
    r3 = in[8*k+6]; i3 = in[8*k+7];
    VTRANSPOSE4(r0,r1,r2,r3);
    VTRANSPOSE4(i0,i1,i2,i3);
    VCPLXMUL(r1,i1,e[k*6+0],e[k*6+1]);
    VCPLXMUL(r2,i2,e[k*6+2],e[k*6+3]);
    VCPLXMUL(r3,i3,e[k*6+4],e[k*6+5]);

    sr0 = VADD(r0,r2); dr0 = VSUB(r0, r2);
    sr1 = VADD(r1,r3); dr1 = VSUB(r1, r3);
    si0 = VADD(i0,i2); di0 = VSUB(i0, i2);
    si1 = VADD(i1,i3); di1 = VSUB(i1, i3);

    /*
      transformation for each column is:

      [1   1   1   1   0   0   0   0]   [r0]
      [1   0  -1   0   0  -1   0   1]   [r1]
      [1  -1   1  -1   0   0   0   0]   [r2]
      [1   0  -1   0   0   1   0  -1]   [r3]
      [0   0   0   0   1   1   1   1] * [i0]
      [0   1   0  -1   1   0  -1   0]   [i1]
      [0   0   0   0   1  -1   1  -1]   [i2]
      [0  -1   0   1   1   0  -1   0]   [i3]
    */

    r0 = VADD(sr0, sr1); i0 = VADD(si0, si1);
    r1 = VADD(dr0, di1); i1 = VSUB(di0, dr1);
    r2 = VSUB(sr0, sr1); i2 = VSUB(si0, si1);
    r3 = VSUB(dr0, di1); i3 = VADD(di0, dr1);

    *out++ = r0; *out++ = i0; *out++ = r1; *out++ = i1;
    *out++ = r2; *out++ = i2; *out++ = r3; *out++ = i3;
  }
}

void FUNC_CPLX_PREPROCESS(int Ncvec, const v4sf *in, v4sf *out, const v4sf *e) {
  int k, dk = Ncvec/SIMD_SZ; /* number of 4x4 matrix blocks */
  v4sf r0, i0, r1, i1, r2, i2, r3, i3;
  v4sf sr0, dr0, sr1, dr1, si0, di0, si1, di1;
  assert(in != out);
  for (k=0; k < dk; ++k) {
    r0 = in[8*k+0]; i0 = in[8*k+1];
    r1 = in[8*k+2]; i1 = in[8*k+3];
    r2 = in[8*k+4]; i2 = in[8*k+5];
    r3 = in[8*k+6]; i3 = in[8*k+7];

    sr0 = VADD(r0,r2); dr0 = VSUB(r0, r2);
    sr1 = VADD(r1,r3); dr1 = VSUB(r1, r3);
    si0 = VADD(i0,i2); di0 = VSUB(i0, i2);
    si1 = VADD(i1,i3); di1 = VSUB(i1, i3);

    r0 = VADD(sr0, sr1); i0 = VADD(si0, si1);
    r1 = VSUB(dr0, di1); i1 = VADD(di0, dr1);
    r2 = VSUB(sr0, sr1); i2 = VSUB(si0, si1);
    r3 = VADD(dr0, di1); i3 = VSUB(di0, dr1);

    VCPLXMULCONJ(r1,i1,e[k*6+0],e[k*6+1]);
    VCPLXMULCONJ(r2,i2,e[k*6+2],e[k*6+3]);
    VCPLXMULCONJ(r3,i3,e[k*6+4],e[k*6+5]);

    VTRANSPOSE4(r0,r1,r2,r3);
    VTRANSPOSE4(i0,i1,i2,i3);

    *out++ = r0; *out++ = i0; *out++ = r1; *out++ = i1;
    *out++ = r2; *out++ = i2; *out++ = r3; *out++ = i3;
  }
}


static ALWAYS_INLINE(void) FUNC_REAL_FINALIZE_4X4(const v4sf *in0, const v4sf *in1, const v4sf *in,
                            const v4sf *e, v4sf *out) {
  v4sf r0, i0, r1, i1, r2, i2, r3, i3;
  v4sf sr0, dr0, sr1, dr1, si0, di0, si1, di1;
  r0 = *in0; i0 = *in1;
  r1 = *in++; i1 = *in++; r2 = *in++; i2 = *in++; r3 = *in++; i3 = *in++;
  VTRANSPOSE4(r0,r1,r2,r3);
  VTRANSPOSE4(i0,i1,i2,i3);

  /*
    transformation for each column is:

    [1   1   1   1   0   0   0   0]   [r0]
    [1   0  -1   0   0  -1   0   1]   [r1]
    [1   0  -1   0   0   1   0  -1]   [r2]
    [1  -1   1  -1   0   0   0   0]   [r3]
    [0   0   0   0   1   1   1   1] * [i0]
    [0  -1   0   1  -1   0   1   0]   [i1]
    [0  -1   0   1   1   0  -1   0]   [i2]
    [0   0   0   0  -1   1  -1   1]   [i3]
  */

  /* cerr << "matrix initial, before e , REAL:\n 1: " << r0 << "\n 1: " << r1 << "\n 1: " << r2 << "\n 1: " << r3 << "\n"; */
  /* cerr << "matrix initial, before e, IMAG :\n 1: " << i0 << "\n 1: " << i1 << "\n 1: " << i2 << "\n 1: " << i3 << "\n"; */

  VCPLXMUL(r1,i1,e[0],e[1]);
  VCPLXMUL(r2,i2,e[2],e[3]);
  VCPLXMUL(r3,i3,e[4],e[5]);

  /* cerr << "matrix initial, real part:\n 1: " << r0 << "\n 1: " << r1 << "\n 1: " << r2 << "\n 1: " << r3 << "\n"; */
  /* cerr << "matrix initial, imag part:\n 1: " << i0 << "\n 1: " << i1 << "\n 1: " << i2 << "\n 1: " << i3 << "\n"; */

  sr0 = VADD(r0,r2); dr0 = VSUB(r0,r2);
  sr1 = VADD(r1,r3); dr1 = VSUB(r3,r1);
  si0 = VADD(i0,i2); di0 = VSUB(i0,i2);
  si1 = VADD(i1,i3); di1 = VSUB(i3,i1);

  r0 = VADD(sr0, sr1);
  r3 = VSUB(sr0, sr1);
  i0 = VADD(si0, si1);
  i3 = VSUB(si1, si0);
  r1 = VADD(dr0, di1);
  r2 = VSUB(dr0, di1);
  i1 = VSUB(dr1, di0);
  i2 = VADD(dr1, di0);

  *out++ = r0;
  *out++ = i0;
  *out++ = r1;
  *out++ = i1;
  *out++ = r2;
  *out++ = i2;
  *out++ = r3;
  *out++ = i3;

}

static NEVER_INLINE(void) FUNC_REAL_FINALIZE(int Ncvec, const v4sf *in, v4sf *out, const v4sf *e) {
  int k, dk = Ncvec/SIMD_SZ; /* number of 4x4 matrix blocks */
  /* fftpack order is f0r f1r f1i f2r f2i ... f(n-1)r f(n-1)i f(n)r */

  v4sf_union cr, ci, *uout = (v4sf_union*)out;
  v4sf save = in[7], zero=VZERO();
  float xr0, xi0, xr1, xi1, xr2, xi2, xr3, xi3;
  static const float s = (float)(M_SQRT1_2);

  cr.v = in[0]; ci.v = in[Ncvec*2-1];
  assert(in != out);
  FUNC_REAL_FINALIZE_4X4(&zero, &zero, in+1, e, out);

  /*
    [cr0 cr1 cr2 cr3 ci0 ci1 ci2 ci3]

    [Xr(1)]  ] [1   1   1   1   0   0   0   0]
    [Xr(N/4) ] [0   0   0   0   1   s   0  -s]
    [Xr(N/2) ] [1   0  -1   0   0   0   0   0]
    [Xr(3N/4)] [0   0   0   0   1  -s   0   s]
    [Xi(1)   ] [1  -1   1  -1   0   0   0   0]
    [Xi(N/4) ] [0   0   0   0   0  -s  -1  -s]
    [Xi(N/2) ] [0  -1   0   1   0   0   0   0]
    [Xi(3N/4)] [0   0   0   0   0  -s   1  -s]
  */

  xr0=(cr.f[0]+cr.f[2]) + (cr.f[1]+cr.f[3]); uout[0].f[0] = xr0;
  xi0=(cr.f[0]+cr.f[2]) - (cr.f[1]+cr.f[3]); uout[1].f[0] = xi0;
  xr2=(cr.f[0]-cr.f[2]);                     uout[4].f[0] = xr2;
  xi2=(cr.f[3]-cr.f[1]);                     uout[5].f[0] = xi2;
  xr1= ci.f[0] + s*(ci.f[1]-ci.f[3]);        uout[2].f[0] = xr1;
  xi1=-ci.f[2] - s*(ci.f[1]+ci.f[3]);        uout[3].f[0] = xi1;
  xr3= ci.f[0] - s*(ci.f[1]-ci.f[3]);        uout[6].f[0] = xr3;
  xi3= ci.f[2] - s*(ci.f[1]+ci.f[3]);        uout[7].f[0] = xi3;

  for (k=1; k < dk; ++k) {
    v4sf save_next = in[8*k+7];
    FUNC_REAL_FINALIZE_4X4(&save, &in[8*k+0], in + 8*k+1,
                           e + k*6, out + k*8);
    save = save_next;
  }

}

static ALWAYS_INLINE(void) FUNC_REAL_PREPROCESS_4X4(const v4sf *in,
                                             const v4sf *e, v4sf *out, int first) {
  v4sf r0=in[0], i0=in[1], r1=in[2], i1=in[3], r2=in[4], i2=in[5], r3=in[6], i3=in[7];
  /*
    transformation for each column is:

    [1   1   1   1   0   0   0   0]   [r0]
    [1   0   0  -1   0  -1  -1   0]   [r1]
    [1  -1  -1   1   0   0   0   0]   [r2]
    [1   0   0  -1   0   1   1   0]   [r3]
    [0   0   0   0   1  -1   1  -1] * [i0]
    [0  -1   1   0   1   0   0   1]   [i1]
    [0   0   0   0   1   1  -1  -1]   [i2]
    [0   1  -1   0   1   0   0   1]   [i3]
  */

  v4sf sr0 = VADD(r0,r3), dr0 = VSUB(r0,r3);
  v4sf sr1 = VADD(r1,r2), dr1 = VSUB(r1,r2);
  v4sf si0 = VADD(i0,i3), di0 = VSUB(i0,i3);
  v4sf si1 = VADD(i1,i2), di1 = VSUB(i1,i2);

  r0 = VADD(sr0, sr1);
  r2 = VSUB(sr0, sr1);
  r1 = VSUB(dr0, si1);
  r3 = VADD(dr0, si1);
  i0 = VSUB(di0, di1);
  i2 = VADD(di0, di1);
  i1 = VSUB(si0, dr1);
  i3 = VADD(si0, dr1);

  VCPLXMULCONJ(r1,i1,e[0],e[1]);
  VCPLXMULCONJ(r2,i2,e[2],e[3]);
  VCPLXMULCONJ(r3,i3,e[4],e[5]);

  VTRANSPOSE4(r0,r1,r2,r3);
  VTRANSPOSE4(i0,i1,i2,i3);

  if (!first) {
    *out++ = r0;
    *out++ = i0;
  }
  *out++ = r1;
  *out++ = i1;
  *out++ = r2;
  *out++ = i2;
  *out++ = r3;
  *out++ = i3;
}

static NEVER_INLINE(void) FUNC_REAL_PREPROCESS(int Ncvec, const v4sf *in, v4sf *out, const v4sf *e) {
  int k, dk = Ncvec/SIMD_SZ; /* number of 4x4 matrix blocks */
  /* fftpack order is f0r f1r f1i f2r f2i ... f(n-1)r f(n-1)i f(n)r */

  v4sf_union Xr, Xi, *uout = (v4sf_union*)out;
  float cr0, ci0, cr1, ci1, cr2, ci2, cr3, ci3;
  static const float s = (float)M_SQRT2;
  assert(in != out);
  for (k=0; k < 4; ++k) {
    Xr.f[k] = ((float*)in)[8*k];
    Xi.f[k] = ((float*)in)[8*k+4];
  }

  FUNC_REAL_PREPROCESS_4X4(in, e, out+1, 1); /* will write only 6 values */

  /*
    [Xr0 Xr1 Xr2 Xr3 Xi0 Xi1 Xi2 Xi3]

    [cr0] [1   0   2   0   1   0   0   0]
    [cr1] [1   0   0   0  -1   0  -2   0]
    [cr2] [1   0  -2   0   1   0   0   0]
    [cr3] [1   0   0   0  -1   0   2   0]
    [ci0] [0   2   0   2   0   0   0   0]
    [ci1] [0   s   0  -s   0  -s   0  -s]
    [ci2] [0   0   0   0   0  -2   0   2]
    [ci3] [0  -s   0   s   0  -s   0  -s]
  */
  for (k=1; k < dk; ++k) {
    FUNC_REAL_PREPROCESS_4X4(in+8*k, e + k*6, out-1+k*8, 0);
  }

  cr0=(Xr.f[0]+Xi.f[0]) + 2*Xr.f[2]; uout[0].f[0] = cr0;
  cr1=(Xr.f[0]-Xi.f[0]) - 2*Xi.f[2]; uout[0].f[1] = cr1;
  cr2=(Xr.f[0]+Xi.f[0]) - 2*Xr.f[2]; uout[0].f[2] = cr2;
  cr3=(Xr.f[0]-Xi.f[0]) + 2*Xi.f[2]; uout[0].f[3] = cr3;
  ci0= 2*(Xr.f[1]+Xr.f[3]);                       uout[2*Ncvec-1].f[0] = ci0;
  ci1= s*(Xr.f[1]-Xr.f[3]) - s*(Xi.f[1]+Xi.f[3]); uout[2*Ncvec-1].f[1] = ci1;
  ci2= 2*(Xi.f[3]-Xi.f[1]);                       uout[2*Ncvec-1].f[2] = ci2;
  ci3=-s*(Xr.f[1]-Xr.f[3]) - s*(Xi.f[1]+Xi.f[3]); uout[2*Ncvec-1].f[3] = ci3;
}


void FUNC_TRANSFORM_INTERNAL(const SETUP_STRUCT *setup, const float *finput, float *foutput, v4sf *scratch,
                             pffft_direction_t direction, int ordered) {
  int k, Ncvec   = setup->Ncvec;
  int nf_odd = (setup->ifac[1] & 1);

  /* temporary buffer is allocated on the stack if the scratch pointer is NULL */
  int stack_allocate = (scratch == 0 ? Ncvec*2 : 1);
  VLA_ARRAY_ON_STACK(v4sf, scratch_on_stack, stack_allocate);

  const v4sf *vinput = (const v4sf*)finput;
  v4sf *voutput      = (v4sf*)foutput;
  v4sf *buff[2]      = { voutput, scratch ? scratch : scratch_on_stack };
  int ib = (nf_odd ^ ordered ? 1 : 0);

  assert(VALIGNED(finput) && VALIGNED(foutput));

  /* assert(finput != foutput); */
  if (direction == PFFFT_FORWARD) {
    ib = !ib;
    if (setup->transform == PFFFT_REAL) {
      ib = (rfftf1_ps(Ncvec*2, vinput, buff[ib], buff[!ib],
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);
      FUNC_REAL_FINALIZE(Ncvec, buff[ib], buff[!ib], (v4sf*)setup->e);
    } else {
      v4sf *tmp = buff[ib];
      for (k=0; k < Ncvec; ++k) {
        UNINTERLEAVE2(vinput[k*2], vinput[k*2+1], tmp[k*2], tmp[k*2+1]);
      }
      ib = (cfftf1_ps(Ncvec, buff[ib], buff[!ib], buff[ib],
                      setup->twiddle, &setup->ifac[0], -1) == buff[0] ? 0 : 1);
      FUNC_CPLX_FINALIZE(Ncvec, buff[ib], buff[!ib], (v4sf*)setup->e);
    }
    if (ordered) {
      FUNC_ZREORDER(setup, (float*)buff[!ib], (float*)buff[ib], PFFFT_FORWARD);
    } else ib = !ib;
  } else {
    if (vinput == buff[ib]) {
      ib = !ib; /* may happen when finput == foutput */
    }
    if (ordered) {
      FUNC_ZREORDER(setup, (float*)vinput, (float*)buff[ib], PFFFT_BACKWARD);
      vinput = buff[ib]; ib = !ib;
    }
    if (setup->transform == PFFFT_REAL) {
      FUNC_REAL_PREPROCESS(Ncvec, vinput, buff[ib], (v4sf*)setup->e);
      ib = (rfftb1_ps(Ncvec*2, buff[ib], buff[0], buff[1],
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);
    } else {
      FUNC_CPLX_PREPROCESS(Ncvec, vinput, buff[ib], (v4sf*)setup->e);
      ib = (cfftf1_ps(Ncvec, buff[ib], buff[0], buff[1],
                      setup->twiddle, &setup->ifac[0], +1) == buff[0] ? 0 : 1);
      for (k=0; k < Ncvec; ++k) {
        INTERLEAVE2(buff[ib][k*2], buff[ib][k*2+1], buff[ib][k*2], buff[ib][k*2+1]);
      }
    }
  }

  if (buff[ib] != voutput) {
    /* extra copy required -- this situation should only happen when finput == foutput */
    assert(finput==foutput);
    for (k=0; k < Ncvec; ++k) {
      v4sf a = buff[ib][2*k], b = buff[ib][2*k+1];
      voutput[2*k] = a; voutput[2*k+1] = b;
    }
    ib = !ib;
  }
  assert(buff[ib] == voutput);
}

void FUNC_ZCONVOLVE_ACCUMULATE(const SETUP_STRUCT *s, const float *a, const float *b, float *ab, float scaling) {
  int Ncvec = s->Ncvec;
  const v4sf *va = (const v4sf*)a;
  const v4sf *vb = (const v4sf*)b;
  v4sf *vab = (v4sf*)ab;

#ifdef __arm__
  __builtin_prefetch(va);
  __builtin_prefetch(vb);
  __builtin_prefetch(vab);
  __builtin_prefetch(va+2);
  __builtin_prefetch(vb+2);
  __builtin_prefetch(vab+2);
  __builtin_prefetch(va+4);
  __builtin_prefetch(vb+4);
  __builtin_prefetch(vab+4);
  __builtin_prefetch(va+6);
  __builtin_prefetch(vb+6);
  __builtin_prefetch(vab+6);
#endif

  float ar0, ai0, br0, bi0, abr0, abi0;
  v4sf vscal = LD_PS1(scaling);
  int i;

  assert(VALIGNED(a) && VALIGNED(b) && VALIGNED(ab));
  ar0 = ((v4sf_union*)va)[0].f[0];
  ai0 = ((v4sf_union*)va)[1].f[0];
  br0 = ((v4sf_union*)vb)[0].f[0];
  bi0 = ((v4sf_union*)vb)[1].f[0];
  abr0 = ((v4sf_union*)vab)[0].f[0];
  abi0 = ((v4sf_union*)vab)[1].f[0];

  for (i=0; i < Ncvec; i += 2) {
    v4sf ar, ai, br, bi;
    ar = va[2*i+0]; ai = va[2*i+1];
    br = vb[2*i+0]; bi = vb[2*i+1];
    VCPLXMUL(ar, ai, br, bi);
    vab[2*i+0] = VMADD(ar, vscal, vab[2*i+0]);
    vab[2*i+1] = VMADD(ai, vscal, vab[2*i+1]);
    ar = va[2*i+2]; ai = va[2*i+3];
    br = vb[2*i+2]; bi = vb[2*i+3];
    VCPLXMUL(ar, ai, br, bi);
    vab[2*i+2] = VMADD(ar, vscal, vab[2*i+2]);
    vab[2*i+3] = VMADD(ai, vscal, vab[2*i+3]);
  }
  if (s->transform == PFFFT_REAL) {
    ((v4sf_union*)vab)[0].f[0] = abr0 + ar0*br0*scaling;
    ((v4sf_union*)vab)[1].f[0] = abi0 + ai0*bi0*scaling;
  }
}

void FUNC_ZCONVOLVE_SCALE(const SETUP_STRUCT *s, const float *a, const float *b, float *ab, float scaling) {
  v4sf vscal = LD_PS1(scaling);
  const v4sf *va = (const v4sf*)a;
  const v4sf *vb = (const v4sf*)b;
  v4sf *vab = (v4sf*)ab;
  float sar, sai, sbr, sbi;
  const int NcvecMulTwo = 2*s->Ncvec;  /* int Ncvec = s->Ncvec; */
  int k; /* was i -- but always used "2*i" - except at for() */

#ifdef __arm__
  __builtin_prefetch(va);
  __builtin_prefetch(vb);
  __builtin_prefetch(vab);
  __builtin_prefetch(va+2);
  __builtin_prefetch(vb+2);
  __builtin_prefetch(vab+2);
  __builtin_prefetch(va+4);
  __builtin_prefetch(vb+4);
  __builtin_prefetch(vab+4);
  __builtin_prefetch(va+6);
  __builtin_prefetch(vb+6);
  __builtin_prefetch(vab+6);
#endif

  assert(VALIGNED(a) && VALIGNED(b) && VALIGNED(ab));
  sar = ((v4sf_union*)va)[0].f[0];
  sai = ((v4sf_union*)va)[1].f[0];
  sbr = ((v4sf_union*)vb)[0].f[0];
  sbi = ((v4sf_union*)vb)[1].f[0];

  /* default routine, works fine for non-arm cpus with current compilers */
  for (k=0; k < NcvecMulTwo; k += 4) {
    v4sf var, vai, vbr, vbi;
    var = va[k+0]; vai = va[k+1];
    vbr = vb[k+0]; vbi = vb[k+1];
    VCPLXMUL(var, vai, vbr, vbi);
    vab[k+0] = VMUL(var, vscal);
    vab[k+1] = VMUL(vai, vscal);
    var = va[k+2]; vai = va[k+3];
    vbr = vb[k+2]; vbi = vb[k+3];
    VCPLXMUL(var, vai, vbr, vbi);
    vab[k+2] = VMUL(var, vscal);
    vab[k+3] = VMUL(vai, vscal);
  }

  if (s->transform == PFFFT_REAL) {
    ((v4sf_union*)vab)[0].f[0] = sar*sbr*scaling;
    ((v4sf_union*)vab)[1].f[0] = sai*sbi*scaling;
  }
}


void FUNC_ZCONVOLVE(const SETUP_STRUCT *s, const float *a, const float *b, float *ab) {
  const v4sf *va = (const v4sf*)a;
  const v4sf *vb = (const v4sf*)b;
  v4sf *vab = (v4sf*)ab;
  float sar, sai, sbr, sbi;
  const int NcvecMulTwo = 2*s->Ncvec;  /* int Ncvec = s->Ncvec; */
  int k; /* was i -- but always used "2*i" - except at for() */

#ifdef __arm__
  __builtin_prefetch(va);
  __builtin_prefetch(vb);
  __builtin_prefetch(vab);
  __builtin_prefetch(va+2);
  __builtin_prefetch(vb+2);
  __builtin_prefetch(vab+2);
  __builtin_prefetch(va+4);
  __builtin_prefetch(vb+4);
  __builtin_prefetch(vab+4);
  __builtin_prefetch(va+6);
  __builtin_prefetch(vb+6);
  __builtin_prefetch(vab+6);
#endif

  assert(VALIGNED(a) && VALIGNED(b) && VALIGNED(ab));
  sar = ((v4sf_union*)va)[0].f[0];
  sai = ((v4sf_union*)va)[1].f[0];
  sbr = ((v4sf_union*)vb)[0].f[0];
  sbi = ((v4sf_union*)vb)[1].f[0];

  /* default routine, works fine for non-arm cpus with current compilers */
  for (k=0; k < NcvecMulTwo; k += 4) {
    v4sf var, vai, vbr, vbi;
    var = va[k+0]; vai = va[k+1];
    vbr = vb[k+0]; vbi = vb[k+1];
    VCPLXMUL(var, vai, vbr, vbi);
    vab[k+0] = var;
    vab[k+1] = vai;
    var = va[k+2]; vai = va[k+3];
    vbr = vb[k+2]; vbi = vb[k+3];
    VCPLXMUL(var, vai, vbr, vbi);
    vab[k+2] = var;
    vab[k+3] = vai;
  }

  if (s->transform == PFFFT_REAL) {
    ((v4sf_union*)vab)[0].f[0] = sar*sbr;
    ((v4sf_union*)vab)[1].f[0] = sai*sbi;
  }
}


int FUNC_ZCONVERT_ZP(const SETUP_STRUCT *s, const float *in, float *out, float scaling) {
  int j, nv = s->Ncvec * 2;  /* number of simd vectors: alternating real/imag vector pairs */
  /* no RESTRICT: in may alias out by contract */
  const v4sf *vin = (const v4sf*)in;
  v4sf *vout = (v4sf*)out;
  const float nscaling = scaling / (float)s->N;  /* unit-gain normalization */
  v4sf vs = LD_PS1(nscaling);
  if (s->transform != PFFFT_REAL) return -1;
  assert(VALIGNED(in) && VALIGNED(out));
  /*
    unordered forward REAL layout: even-indexed vectors hold the real
    parts of 4 coefficients each, odd-indexed vectors their imaginary
    parts -- except the very first odd lane, which packs the real
    Nyquist coefficient F(N/2) next to F(0) in vector 0 (same packing
    that pffft_zconvolve_accumulate() treats component-wise).
    Zero-phase form: keep/scale the real vectors, zero the imag
    vectors, keep/scale the Nyquist lane.

    The conversion applies a 1/N factor (times 'scaling'): with
    scaling == 1, transform(x, FORWARD) -> zconvolve_zp() ->
    transform(.., BACKWARD) yields the circular convolution of x
    with the filter at unit gain. */
  vout[0] = VMUL(vin[0], vs);
  {
    /* latch the scaled Nyquist coefficient before zeroing: with
       in == out the vector store below would destroy vin[1].f[0] */
    const float nyquist = ((const v4sf_union*)vin)[1].f[0] * nscaling;
    const v4sf vz = VZERO();
    vout[1] = vz;
    ((v4sf_union*)vout)[1].f[0] = nyquist;
    for (j = 2; j < nv; j += 2) {
      vout[j] = VMUL(vin[j], vs);
      vout[j+1] = vz;
    }
  }
  return 0;
}

int FUNC_ZCONVOLVE_ZP(const SETUP_STRUCT *s, const float *x, const float *hzp, float *ab) {
  int j, np = s->Ncvec;  /* number of real/imag vector pairs */
  /* no RESTRICT: any pair of pointers may alias by contract */
  const v4sf *vx = (const v4sf*)x;
  const v4sf *vh = (const v4sf*)hzp;
  v4sf *vab = (v4sf*)ab;
  if (s->transform != PFFFT_REAL) return -1;
  assert(VALIGNED(x) && VALIGNED(hzp) && VALIGNED(ab));
  /*
    with a zero-phase filter (imag vectors all zero):
    Y_re = X_re * H_re   and   Y_im = X_im * H_re.
    The first two vectors additionally carry the real DC and Nyquist
    coefficients in their f[0] lanes (see FUNC_ZCONVERT_ZP), so that
    pair is finished scalar-wise, component-wise, exactly like
    FUNC_ZCONVOLVE_ACCUMULATE() does. */
  {
    /* note: the Nyquist product reads hv1.f[0], not hv0.f[0] */
    v4sf_union xv1, hv0, hv1;
    xv1.v = vx[1];
    hv0.v = vh[0]; hv1.v = vh[1];
    vab[0] = VMUL(vx[0], vh[0]);
    ((v4sf_union*)vab)[1].f[0] = xv1.f[0] * hv1.f[0];
    ((v4sf_union*)vab)[1].f[1] = xv1.f[1] * hv0.f[1];
    ((v4sf_union*)vab)[1].f[2] = xv1.f[2] * hv0.f[2];
    ((v4sf_union*)vab)[1].f[3] = xv1.f[3] * hv0.f[3];
  }
  for (j = 1; j < np; ++j) {
    v4sf hr = vh[2*j];
    vab[2*j] = VMUL(vx[2*j], hr);
    vab[2*j+1] = VMUL(vx[2*j+1], hr);
  }
  return 0;
}

#else  /* #if ( SIMD_SZ == 4 )   * !defined(PFFFT_SIMD_DISABLE) */

/* standard routine using scalar floats, without SIMD stuff. */

#define pffft_zreorder_nosimd FUNC_ZREORDER
void pffft_zreorder_nosimd(const SETUP_STRUCT *setup, const float *in, float *out, pffft_direction_t direction) {
  int k, N = setup->N;
  if (setup->transform == PFFFT_COMPLEX) {
    for (k=0; k < 2*N; ++k) out[k] = in[k];
    return;
  }
  else if (direction == PFFFT_FORWARD) {
    float x_N = in[N-1];
    for (k=N-1; k > 1; --k) out[k] = in[k-1];
    out[0] = in[0];
    out[1] = x_N;
  } else {
    float x_N = in[1];
    for (k=1; k < N-1; ++k) out[k] = in[k+1];
    out[0] = in[0];
    out[N-1] = x_N;
  }
}

#define pffft_transform_internal_nosimd FUNC_TRANSFORM_INTERNAL
void pffft_transform_internal_nosimd(const SETUP_STRUCT *setup, const float *input, float *output, float *scratch,
                                    pffft_direction_t direction, int ordered) {
  int Ncvec   = setup->Ncvec;
  int nf_odd = (setup->ifac[1] & 1);

  /* temporary buffer is allocated on the stack if the scratch pointer is NULL */
  int stack_allocate = (scratch == 0 ? Ncvec*2 : 1);
  VLA_ARRAY_ON_STACK(v4sf, scratch_on_stack, stack_allocate);
  float *buff[2];
  int ib;
  if (scratch == 0) scratch = scratch_on_stack;
  buff[0] = output; buff[1] = scratch;

  if (setup->transform == PFFFT_COMPLEX) ordered = 0; /* it is always ordered. */
  ib = (nf_odd ^ ordered ? 1 : 0);

  if (direction == PFFFT_FORWARD) {
    if (setup->transform == PFFFT_REAL) {
      ib = (rfftf1_ps(Ncvec*2, input, buff[ib], buff[!ib],
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);
    } else {
      ib = (cfftf1_ps(Ncvec, input, buff[ib], buff[!ib],
                      setup->twiddle, &setup->ifac[0], -1) == buff[0] ? 0 : 1);
    }
    if (ordered) {
      FUNC_ZREORDER(setup, buff[ib], buff[!ib], PFFFT_FORWARD); ib = !ib;
    }
  } else {
    if (input == buff[ib]) {
      ib = !ib; /* may happen when finput == foutput */
    }
    if (ordered) {
      FUNC_ZREORDER(setup, input, buff[!ib], PFFFT_BACKWARD);
      input = buff[!ib];
    }
    if (setup->transform == PFFFT_REAL) {
      ib = (rfftb1_ps(Ncvec*2, input, buff[ib], buff[!ib],
                      setup->twiddle, &setup->ifac[0]) == buff[0] ? 0 : 1);
    } else {
      ib = (cfftf1_ps(Ncvec, input, buff[ib], buff[!ib],
                      setup->twiddle, &setup->ifac[0], +1) == buff[0] ? 0 : 1);
    }
  }
  if (buff[ib] != output) {
    int k;
    /* extra copy required -- this situation should happens only when finput == foutput */
    assert(input==output);
    for (k=0; k < Ncvec; ++k) {
      float a = buff[ib][2*k], b = buff[ib][2*k+1];
      output[2*k] = a; output[2*k+1] = b;
    }
    ib = !ib;
  }
  assert(buff[ib] == output);
}

#define pffft_zconvolve_accumulate_nosimd FUNC_ZCONVOLVE_ACCUMULATE
void pffft_zconvolve_accumulate_nosimd(const SETUP_STRUCT *s, const float *a, const float *b,
                                       float *ab, float scaling) {
  int NcvecMulTwo = 2*s->Ncvec;  /* int Ncvec = s->Ncvec; */
  int k; /* was i -- but always used "2*i" - except at for() */

  if (s->transform == PFFFT_REAL) {
    /* take care of the fftpack ordering */
    ab[0] += a[0]*b[0]*scaling;
    ab[NcvecMulTwo-1] += a[NcvecMulTwo-1]*b[NcvecMulTwo-1]*scaling;
    ++ab; ++a; ++b; NcvecMulTwo -= 2;
  }
  for (k=0; k < NcvecMulTwo; k += 2) {
    float ar, ai, br, bi;
    ar = a[k+0]; ai = a[k+1];
    br = b[k+0]; bi = b[k+1];
    VCPLXMUL(ar, ai, br, bi);
    ab[k+0] += ar*scaling;
    ab[k+1] += ai*scaling;
  }
}

#define pffft_zconvolve_scale_nosimd FUNC_ZCONVOLVE_SCALE
void pffft_zconvolve_scale_nosimd(const SETUP_STRUCT *s, const float *a, const float *b,
                                    float *ab, float scaling) {
  int NcvecMulTwo = 2*s->Ncvec;  /* int Ncvec = s->Ncvec; */
  int k; /* was i -- but always used "2*i" - except at for() */

  if (s->transform == PFFFT_REAL) {
    /* take care of the fftpack ordering */
    ab[0] = a[0]*b[0]*scaling;
    ab[NcvecMulTwo-1] = a[NcvecMulTwo-1]*b[NcvecMulTwo-1]*scaling;
    ++ab; ++a; ++b; NcvecMulTwo -= 2;
  }
  for (k=0; k < NcvecMulTwo; k += 2) {
    float ar, ai, br, bi;
    ar = a[k+0]; ai = a[k+1];
    br = b[k+0]; bi = b[k+1];
    VCPLXMUL(ar, ai, br, bi);
    ab[k+0] = ar*scaling;
    ab[k+1] = ai*scaling;
  }
}


#define pffft_zconvolve_nosimd FUNC_ZCONVOLVE
void pffft_zconvolve_nosimd(const SETUP_STRUCT *s, const float *a, const float *b,
                            float *ab) {
  int NcvecMulTwo = 2*s->Ncvec;  /* int Ncvec = s->Ncvec; */
  int k; /* was i -- but always used "2*i" - except at for() */

  if (s->transform == PFFFT_REAL) {
    /* take care of the fftpack ordering */
    ab[0] = a[0]*b[0];
    ab[NcvecMulTwo-1] = a[NcvecMulTwo-1]*b[NcvecMulTwo-1];
    ++ab; ++a; ++b; NcvecMulTwo -= 2;
  }
  for (k=0; k < NcvecMulTwo; k += 2) {
    float ar, ai, br, bi;
    ar = a[k+0]; ai = a[k+1];
    br = b[k+0]; bi = b[k+1];
    VCPLXMUL(ar, ai, br, bi);
    ab[k+0] = ar;
    ab[k+1] = ai;
  }
}


#define pffft_zconvert_zp_nosimd FUNC_ZCONVERT_ZP
int pffft_zconvert_zp_nosimd(const SETUP_STRUCT *s, const float *in, float *out, float scaling) {
  int k, N2 = s->Ncvec * 2;
  float nscaling = scaling / (float)s->N;  /* unit-gain normalization */
  if (s->transform != PFFFT_REAL) return -1;
  out[0] = in[0] * nscaling;
  out[N2-1] = in[N2-1] * nscaling;
  for (k = 1; k < N2-1; k += 2) {
    out[k] = in[k] * nscaling;
    out[k+1] = 0.f;
  }
  return 0;
}

#define pffft_zconvolve_zp_nosimd FUNC_ZCONVOLVE_ZP
int pffft_zconvolve_zp_nosimd(const SETUP_STRUCT *s, const float *x, const float *hzp, float *ab) {
  int k, N2 = s->Ncvec * 2;
  if (s->transform != PFFFT_REAL) return -1;
  ab[0] = x[0] * hzp[0];
  ab[N2-1] = x[N2-1] * hzp[N2-1];
  for (k = 1; k < N2-1; k += 2) {
    float hr = hzp[k];
    ab[k] = x[k] * hr;
    ab[k+1] = x[k+1] * hr;
  }
  return 0;
}

#endif /* #if ( SIMD_SZ == 4 )    * !defined(PFFFT_SIMD_DISABLE) */


void FUNC_TRANSFORM_UNORDRD(const SETUP_STRUCT *setup, const float *input, float *output, float *work, pffft_direction_t direction) {
  FUNC_TRANSFORM_INTERNAL(setup, input, output, (v4sf*)work, direction, 0);
}

void FUNC_TRANSFORM_ORDERED(const SETUP_STRUCT *setup, const float *input, float *output, float *work, pffft_direction_t direction) {
  FUNC_TRANSFORM_INTERNAL(setup, input, output, (v4sf*)work, direction, 1);
}


#if ( SIMD_SZ == 4 )

#define assertv4(v,f0,f1,f2,f3) assert(v.f[0] == (f0) && v.f[1] == (f1) && v.f[2] == (f2) && v.f[3] == (f3))

/* detect bugs with the vector support macros */
PFFFT_EXPORT void FUNC_VALIDATE_SIMD_A(void) {
  float f[16] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 };
  v4sf_union a0, a1, a2, a3, t, u;
  memcpy(a0.f, f, 4*sizeof(float));
  memcpy(a1.f, f+4, 4*sizeof(float));
  memcpy(a2.f, f+8, 4*sizeof(float));
  memcpy(a3.f, f+12, 4*sizeof(float));

  t = a0; u = a1; t.v = VZERO();
  printf("VZERO=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 0, 0, 0, 0);
  t.v = VADD(a1.v, a2.v);
  printf("VADD(4:7,8:11)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 12, 14, 16, 18);
  t.v = VMUL(a1.v, a2.v);
  printf("VMUL(4:7,8:11)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 32, 45, 60, 77);
  t.v = VMADD(a1.v, a2.v,a0.v);
  printf("VMADD(4:7,8:11,0:3)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]); assertv4(t, 32, 46, 62, 80);

  INTERLEAVE2(a1.v,a2.v,t.v,u.v);
  printf("INTERLEAVE2(4:7,8:11)=[%2g %2g %2g %2g] [%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3], u.f[0], u.f[1], u.f[2], u.f[3]);
  assertv4(t, 4, 8, 5, 9); assertv4(u, 6, 10, 7, 11);
  UNINTERLEAVE2(a1.v,a2.v,t.v,u.v);
  printf("UNINTERLEAVE2(4:7,8:11)=[%2g %2g %2g %2g] [%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3], u.f[0], u.f[1], u.f[2], u.f[3]);
  assertv4(t, 4, 6, 8, 10); assertv4(u, 5, 7, 9, 11);

  t.v=LD_PS1(f[15]);
  printf("LD_PS1(15)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]);
  assertv4(t, 15, 15, 15, 15);
  t.v = VSWAPHL(a1.v, a2.v);
  printf("VSWAPHL(4:7,8:11)=[%2g %2g %2g %2g]\n", t.f[0], t.f[1], t.f[2], t.f[3]);
  assertv4(t, 8, 9, 6, 7);
  VTRANSPOSE4(a0.v, a1.v, a2.v, a3.v);
  printf("VTRANSPOSE4(0:3,4:7,8:11,12:15)=[%2g %2g %2g %2g] [%2g %2g %2g %2g] [%2g %2g %2g %2g] [%2g %2g %2g %2g]\n",
         a0.f[0], a0.f[1], a0.f[2], a0.f[3], a1.f[0], a1.f[1], a1.f[2], a1.f[3],
         a2.f[0], a2.f[1], a2.f[2], a2.f[3], a3.f[0], a3.f[1], a3.f[2], a3.f[3]);
  assertv4(a0, 0, 4, 8, 12); assertv4(a1, 1, 5, 9, 13); assertv4(a2, 2, 6, 10, 14); assertv4(a3, 3, 7, 11, 15);
}


static void pffft_assert1( float result, float ref, const char * vartxt, const char * functxt, int * numErrs, const char * f, int lineNo )
{
  if ( !( fabs( result - ref ) < 0.01F ) )
  {
    fprintf(stderr, "%s: assert for %s at %s(%d)\n  expected %f  value %f\n", functxt, vartxt, f, lineNo, ref, result);
    ++(*numErrs);
  }
}

static void pffft_assert4(  vsfscalar v0, vsfscalar v1, vsfscalar v2, vsfscalar v3,
  float a, float b, float c, float d, const char * functxt, int * numErrs, const char * f, int lineNo )
{
  pffft_assert1( v0, a, "[0]", functxt, numErrs, f, lineNo );
  pffft_assert1( v1, b, "[1]", functxt, numErrs, f, lineNo );
  pffft_assert1( v2, c, "[2]", functxt, numErrs, f, lineNo );
  pffft_assert1( v3, d, "[3]", functxt, numErrs, f, lineNo );
}

#define PFFFT_ASSERT4( V, a, b, c, d, FUNCTXT )  pffft_assert4( (V).f[0], (V).f[1], (V).f[2], (V).f[3], a, b, c, d, FUNCTXT, &numErrs, __FILE__, __LINE__ )


PFFFT_EXPORT int FUNC_VALIDATE_SIMD_EX(FILE * DbgOut)
{
  int numErrs = 0;

  {
    v4sf_union C;
    int k;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: { }\n" );
    }
    C.v = VZERO();
    if (DbgOut) {
      fprintf(DbgOut, "VZERO(a) => C) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( C, 0.0F, 0.0F, 0.0F, 0.0F, "VZERO() Out C" );
  }

  {
    v4sf_union C;
    float a = 42.0F;
    int k;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: a = {\n" );
      fprintf(DbgOut, "  Inp a:  %f\n", a );
      fprintf(DbgOut, "}\n" );
    }
    C.v = LD_PS1(a);
    if (DbgOut) {
      fprintf(DbgOut, "LD_PS1(a) => C) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( C, 42.0F, 42.0F, 42.0F, 42.0F, "LD_PS1() Out C" );
  }

  {
    v4sf_union C;
    float a[16];
    int numAligned = 0, numUnaligned = 0;
    int k;
    const char * pUn;
    for ( k = 0; k < 16; ++k ) a[k] = k+1;

    for ( k = 0; k + 3 < 16; ++k )
    {
      const float * ptr = &a[k];
      if (DbgOut)
        fprintf(DbgOut, "\ninput: a = [ %f, %f, %f, %f ]\n", ptr[0], ptr[1], ptr[2], ptr[3] );
      if ( VALIGNED(ptr) )
      {
        C.v = VLOAD_ALIGNED( ptr );
        pUn = "";
        ++numAligned;
      }
      else
      {
        C.v = VLOAD_UNALIGNED( ptr );
        pUn = "UN";
        ++numUnaligned;
      }
      if (DbgOut) {
        fprintf(DbgOut, "C = VLOAD_%sALIGNED(&a[%d]) => {\n", pUn, k );
        fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
        fprintf(DbgOut, "}\n" );
      }
      //PFFFT_ASSERT4( C, 32.0F, 34.0F, 36.0F, 38.0F, "VADD(): Out C" );

      if ( numAligned >= 1 && numUnaligned >= 4 )
        break;
    }
    if ( numAligned < 1 ) {
      fprintf(stderr, "VALIGNED() should have found at least 1 occurence!");
      ++numErrs;
    }
    /* backends without alignment requirements (NEON/WASM float) may accept
       every float address -- then there is no unaligned path to exercise */
    if ( numUnaligned < 4 && !(VALIGNED(&a[1]) && VALIGNED(&a[2]) && VALIGNED(&a[3])) ) {
      fprintf(stderr, "!VALIGNED() should have found at least 4 occurences!");
      ++numErrs;
    }
  }

  {
    v4sf_union A, B, C;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = 20 + k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    C.v = VADD(A.v, B.v);
    if (DbgOut) {
      fprintf(DbgOut, "C = VADD(A,B) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "VADD(): Inp A" );
    PFFFT_ASSERT4( B, 21.0F, 22.0F, 23.0F, 24.0F, "VADD(): Inp B" );
    PFFFT_ASSERT4( C, 32.0F, 34.0F, 36.0F, 38.0F, "VADD(): Out C" );
  }

  {
    v4sf_union A, B, C;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 20 + 2*k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    C.v = VSUB(A.v, B.v);
    if (DbgOut) {
      fprintf(DbgOut, "C = VSUB(A,B) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 21.0F, 23.0F, 25.0F, 27.0F, "VSUB(): Inp A" );
    PFFFT_ASSERT4( B, 11.0F, 12.0F, 13.0F, 14.0F, "VSUB(): Inp B" );
    PFFFT_ASSERT4( C, 10.0F, 11.0F, 12.0F, 13.0F, "VSUB(): Out C" );
  }

  {
    v4sf_union A, B, C;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    C.v = VMUL(A.v, B.v);
    if (DbgOut) {
      fprintf(DbgOut, "C = VMUL(A,B) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "VMUL(): Inp A" );
    PFFFT_ASSERT4( B,  1.0F,  2.0F,  3.0F,  4.0F, "VMUL(): Inp B" );
    PFFFT_ASSERT4( C, 11.0F, 24.0F, 39.0F, 56.0F, "VMUL(): Out C" );
  }

  {
    v4sf_union A, B, C, D;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 10 + k;
    for ( k = 0; k < 4; ++k )  D.f[k] = 40 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B,C = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "  Inp C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    D.v = VMADD(A.v, B.v, C.v);
    if (DbgOut) {
      fprintf(DbgOut, "D = VMADD(A,B,C) => {\n" );
      fprintf(DbgOut, "  Out D:  %f, %f, %f, %f\n", D.f[0], D.f[1], D.f[2], D.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "VMADD(): Inp A" );
    PFFFT_ASSERT4( B,  1.0F,  2.0F,  3.0F,  4.0F, "VMADD(): Inp B" );
    PFFFT_ASSERT4( C, 10.0F, 11.0F, 12.0F, 13.0F, "VMADD(): Inp C" );
    PFFFT_ASSERT4( D, 21.0F, 35.0F, 51.0F, 69.0F, "VMADD(): Out D" );
  }

  {
    v4sf_union A, B, C, D;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = 20 + k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;
    for ( k = 0; k < 4; ++k )  D.f[k] = 40 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    INTERLEAVE2(A.v, B.v, C.v, D.v);
    if (DbgOut) {
      fprintf(DbgOut, "INTERLEAVE2(A,B, => C,D) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "  Out D:  %f, %f, %f, %f\n", D.f[0], D.f[1], D.f[2], D.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "INTERLEAVE2() Inp A" );
    PFFFT_ASSERT4( B, 21.0F, 22.0F, 23.0F, 24.0F, "INTERLEAVE2() Inp B" );
    PFFFT_ASSERT4( C, 11.0F, 21.0F, 12.0F, 22.0F, "INTERLEAVE2() Out C" );
    PFFFT_ASSERT4( D, 13.0F, 23.0F, 14.0F, 24.0F, "INTERLEAVE2() Out D" );
  }

  {
    v4sf_union A, B, C, D;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = 20 + k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;
    for ( k = 0; k < 4; ++k )  D.f[k] = 40 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    UNINTERLEAVE2(A.v, B.v, C.v, D.v);
    if (DbgOut) {
      fprintf(DbgOut, "UNINTERLEAVE2(A,B, => C,D) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "  Out D:  %f, %f, %f, %f\n", D.f[0], D.f[1], D.f[2], D.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "UNINTERLEAVE2() Inp A" );
    PFFFT_ASSERT4( B, 21.0F, 22.0F, 23.0F, 24.0F, "UNINTERLEAVE2() Inp B" );
    PFFFT_ASSERT4( C, 11.0F, 13.0F, 21.0F, 23.0F, "UNINTERLEAVE2() Out C" );
    PFFFT_ASSERT4( D, 12.0F, 14.0F, 22.0F, 24.0F, "UNINTERLEAVE2() Out D" );
  }

  {
    v4sf_union A, B, C, D;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = 20 + k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;
    for ( k = 0; k < 4; ++k )  D.f[k] = 40 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B,C,D = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "  Inp C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "  Inp D:  %f, %f, %f, %f\n", D.f[0], D.f[1], D.f[2], D.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    VTRANSPOSE4(A.v, B.v, C.v, D.v);
    if (DbgOut) {
      fprintf(DbgOut, "VTRANSPOSE4(A,B,C,D) => {\n" );
      fprintf(DbgOut, "  Out A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Out B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "  Out D:  %f, %f, %f, %f\n", D.f[0], D.f[1], D.f[2], D.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 21.0F, 31.0F, 41.0F, "VTRANSPOSE4(): Out A" );
    PFFFT_ASSERT4( B, 12.0F, 22.0F, 32.0F, 42.0F, "VTRANSPOSE4(): Out B" );
    PFFFT_ASSERT4( C, 13.0F, 23.0F, 33.0F, 43.0F, "VTRANSPOSE4(): Out C" );
    PFFFT_ASSERT4( D, 14.0F, 24.0F, 34.0F, 44.0F, "VTRANSPOSE4(): Out D" );
  }

  {
    v4sf_union A, B, C;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  B.f[k] = 20 + k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A,B = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "  Inp B:  %f, %f, %f, %f\n", B.f[0], B.f[1], B.f[2], B.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    C.v = VSWAPHL(A.v, B.v);
    if (DbgOut) {
      fprintf(DbgOut, "C = VSWAPHL(A,B) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "VSWAPHL(): Inp A" );
    PFFFT_ASSERT4( B, 21.0F, 22.0F, 23.0F, 24.0F, "VSWAPHL(): Inp B" );
    PFFFT_ASSERT4( C, 21.0F, 22.0F, 13.0F, 14.0F, "VSWAPHL(): Out C" );
  }

  {
    v4sf_union A, C;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;
    for ( k = 0; k < 4; ++k )  C.f[k] = 30 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    C.v = VREV_S(A.v);
    if (DbgOut) {
      fprintf(DbgOut, "C = VREV_S(A) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "VREV_S(): Inp A" );
    PFFFT_ASSERT4( C, 14.0F, 13.0F, 12.0F, 11.0F, "VREV_S(): Out C" );
  }

  {
    v4sf_union A, C;
    int k;
    for ( k = 0; k < 4; ++k )  A.f[k] = 10 + k+1;

    if (DbgOut) {
      fprintf(DbgOut, "\ninput: A = {\n" );
      fprintf(DbgOut, "  Inp A:  %f, %f, %f, %f\n", A.f[0], A.f[1], A.f[2], A.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    C.v = VREV_C(A.v);
    if (DbgOut) {
      fprintf(DbgOut, "C = VREV_C(A) => {\n" );
      fprintf(DbgOut, "  Out C:  %f, %f, %f, %f\n", C.f[0], C.f[1], C.f[2], C.f[3] );
      fprintf(DbgOut, "}\n" );
    }
    PFFFT_ASSERT4( A, 11.0F, 12.0F, 13.0F, 14.0F, "VREV_C(): Inp A" );
    PFFFT_ASSERT4( C, 13.0F, 14.0F, 11.0F, 12.0F, "VREV_C(): Out A" );
  }

  return numErrs;
}

#else  /* if ( SIMD_SZ == 4 ) */

PFFFT_EXPORT void FUNC_VALIDATE_SIMD_A()
{
}

PFFFT_EXPORT int FUNC_VALIDATE_SIMD_EX(FILE * DbgOut)
{
  return -1;
}

#endif  /* end if ( SIMD_SZ == 4 ) */

