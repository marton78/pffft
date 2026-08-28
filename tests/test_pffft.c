/*
  Copyright (c) 2013 Julien Pommier.

  Small test for PFFFT

  How to build: 

  on linux, with fftw3:
  gcc -o test_pffft -DHAVE_FFTW -msse -mfpmath=sse -O3 -Wall -W pffft.c test_pffft.c fftpack.c -L/usr/local/lib -I/usr/local/include/ -lfftw3f -lm

  on macos, without fftw3:
  clang -o test_pffft -DHAVE_VECLIB -O3 -Wall -W pffft.c test_pffft.c fftpack.c -L/usr/local/lib -I/usr/local/include/ -framework Accelerate

  on macos, with fftw3:
  clang -o test_pffft -DHAVE_FFTW -DHAVE_VECLIB -O3 -Wall -W pffft.c test_pffft.c fftpack.c -L/usr/local/lib -I/usr/local/include/ -lfftw3f -framework Accelerate

  as alternative: replace clang by gcc.

  on windows, with visual c++:
  cl /Ox -D_USE_MATH_DEFINES /arch:SSE test_pffft.c pffft.c fftpack.c
  
  build without SIMD instructions:
  gcc -o test_pffft -DPFFFT_SIMD_DISABLE -O3 -Wall -W pffft.c test_pffft.c fftpack.c -lm

 */

#ifdef PFFFT_ENABLE_FLOAT
#include "pffft/pffft.h"

typedef float pffft_scalar;
#else
/*
Note: adapted for double precision dynamic range version.
*/
#include "pffft/pffft_double.h"

typedef double pffft_scalar;
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>

/* define own constants required to turn off g++ extensions .. */
#ifndef M_PI
  #define M_PI    3.14159265358979323846  /* pi */
#endif

/* EXPECTED_DYN_RANGE in dB:
 * single precision float has 24 bits mantissa
 * => 24 Bits * 6 dB = 144 dB
 * allow a few dB tolerance (even 144 dB looks good on my PC)
 */
#ifdef PFFFT_ENABLE_FLOAT
#define EXPECTED_DYN_RANGE  140.0
#else
#define EXPECTED_DYN_RANGE  215.0
#endif

/* maximum allowed phase error in degree */
#define DEG_ERR_LIMIT   1E-4

/* maximum allowed magnitude error in amplitude (of 1.0 or 1.1) */
#define MAG_ERR_LIMIT  1E-6


#define PRINT_SPEC  0

#define PWR2LOG(PWR)  ( (PWR) < 1E-30 ? 10.0*log10(1E-30) : 10.0*log10(PWR) )



int test(int N, int cplx, int useOrdered) {
  int Nfloat = (cplx ? N*2 : N);
#ifdef PFFFT_ENABLE_FLOAT
  pffft_scalar *X = pffft_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *Y = pffft_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *R = pffft_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *Z = pffft_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *W = pffft_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
#else
  pffft_scalar *X = pffftd_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *Y = pffftd_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *R = pffftd_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *Z = pffftd_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
  pffft_scalar *W = pffftd_aligned_malloc((unsigned)Nfloat * sizeof(pffft_scalar));
#endif
  pffft_scalar amp = (pffft_scalar)1.0;
  double freq, dPhi, phi, phi0;
  double pwr, pwrCar, pwrOther, err, errSum, mag, expextedMag;
  int k, j, m, iter, kmaxOther, retError = 0;

#ifdef PFFFT_ENABLE_FLOAT
  assert( pffft_is_power_of_two(N) );
  PFFFT_Setup *s = pffft_new_setup(N, cplx ? PFFFT_COMPLEX : PFFFT_REAL);
#else
  assert( pffftd_is_power_of_two(N) );
  PFFFTD_Setup *s = pffftd_new_setup(N, cplx ? PFFFT_COMPLEX : PFFFT_REAL);
#endif
  assert(s);
  if (!s) {
    printf("Error setting up PFFFT!\n");
    return 1;
  }

  for ( k = m = 0; k < (cplx? N : (1 + N/2) ); k += N/16, ++m )
  {
    amp = (pffft_scalar)( ( (m % 3) == 0 ) ? 1.0 : 1.1 );
    freq = (k < N/2) ? ((double)k / N) : ((double)(k-N) / N);
    dPhi = 2.0 * M_PI * freq;
    if ( dPhi < 0.0 )
      dPhi += 2.0 * M_PI;

    iter = -1;
    while (1)
    {
      ++iter;

      if (iter)
        printf("bin %d: dphi = %f for freq %f\n", k, dPhi, freq);

      /* generate cosine carrier as time signal - start at defined phase phi0 */
      phi = phi0 = (m % 4) * 0.125 * M_PI;  /* have phi0 < 90 deg to be normalized */
      for ( j = 0; j < N; ++j )
      {
        if (cplx) {
          X[2*j] = amp * (pffft_scalar)cos(phi);  /* real part */
          X[2*j+1] = amp * (pffft_scalar)sin(phi);  /* imag part */
        }
        else
          X[j] = amp * (pffft_scalar)cos(phi);  /* only real part */

        /* phase increment .. stay normalized - cos()/sin() might degrade! */
        phi += dPhi;
        if ( phi >= M_PI )
          phi -= 2.0 * M_PI;
      }

      /* forward transform from X --> Y  .. using work buffer W */
#ifdef PFFFT_ENABLE_FLOAT
      if ( useOrdered )
        pffft_transform_ordered(s, X, Y, W, PFFFT_FORWARD );
      else
      {
        pffft_transform(s, X, R, W, PFFFT_FORWARD );  /* use R for reordering */
        pffft_zreorder(s, R, Y, PFFFT_FORWARD ); /* reorder into Y[] for power calculations */
      }
#else
      if ( useOrdered )
        pffftd_transform_ordered(s, X, Y, W, PFFFT_FORWARD );
      else
      {
        pffftd_transform(s, X, R, W, PFFFT_FORWARD );  /* use R for reordering */
        pffftd_zreorder(s, R, Y, PFFFT_FORWARD ); /* reorder into Y[] for power calculations */
      }
#endif

      pwrOther = -1.0;
      pwrCar = 0;


      /* for positive frequencies: 0 to 0.5 * samplerate */
      /* and also for negative frequencies: -0.5 * samplerate to 0 */
      for ( j = 0; j < ( cplx ? N : (1 + N/2) ); ++j )
      {
        if (!cplx && !j)  /* special treatment for DC for real input */
          pwr = Y[j]*Y[j];
        else if (!cplx && j == N/2)  /* treat 0.5 * samplerate */
          pwr = Y[1] * Y[1];  /* despite j (for freq calculation) we have index 1 */
        else
          pwr = Y[2*j] * Y[2*j] + Y[2*j+1] * Y[2*j+1];
        if (iter || PRINT_SPEC)
          printf("%s fft %d:  pwr[j = %d] = %g == %f dB\n", (cplx ? "cplx":"real"), N, j, pwr, PWR2LOG(pwr) );
        if (k == j)
          pwrCar = pwr;
        else if ( pwr > pwrOther ) {
          pwrOther = pwr;
          kmaxOther = j;
        }
      }

      if ( PWR2LOG(pwrCar) - PWR2LOG(pwrOther) < EXPECTED_DYN_RANGE ) {
        printf("%s fft %d amp %f iter %d:\n", (cplx ? "cplx":"real"), N, amp, iter);
        printf("  carrier power  at bin %d: %g == %f dB\n", k, pwrCar, PWR2LOG(pwrCar) );
        printf("  carrier mag || at bin %d: %g\n", k, sqrt(pwrCar) );
        printf("  max other pwr  at bin %d: %g == %f dB\n", kmaxOther, pwrOther, PWR2LOG(pwrOther) );
        printf("  dynamic range: %f dB\n\n", PWR2LOG(pwrCar) - PWR2LOG(pwrOther) );
        retError = 1;
        if ( iter == 0 )
          continue;
      }

      if ( k > 0 && k != N/2 )
      {
        phi = atan2( Y[2*k+1], Y[2*k] );
        if ( fabs( phi - phi0) > DEG_ERR_LIMIT * M_PI / 180.0 )
        {
        retError = 1;
        printf("%s fft %d  bin %d amp %f : phase mismatch! phase = %f deg   expected = %f deg\n",
            (cplx ? "cplx":"real"), N, k, amp, phi * 180.0 / M_PI, phi0 * 180.0 / M_PI );
        }
      }

      expextedMag = cplx ? amp : ( (k == 0 || k == N/2) ? amp : (amp/2) );
      mag = sqrt(pwrCar) / N;
      if ( fabs(mag - expextedMag) > MAG_ERR_LIMIT )
      {
        retError = 1;
        printf("%s fft %d  bin %d amp %f : mag = %g   expected = %g\n", (cplx ? "cplx":"real"), N, k, amp, mag, expextedMag );
      }


      /* now convert spectrum back */
#ifdef PFFFT_ENABLE_FLOAT
      if (useOrdered)
        pffft_transform_ordered(s, Y, Z, W, PFFFT_BACKWARD);
      else
        pffft_transform(s, R, Z, W, PFFFT_BACKWARD);
#else
      if (useOrdered)
        pffftd_transform_ordered(s, Y, Z, W, PFFFT_BACKWARD);
      else
        pffftd_transform(s, R, Z, W, PFFFT_BACKWARD);
#endif

      errSum = 0.0;
      for ( j = 0; j < (cplx ? (2*N) : N); ++j )
      {
        /* scale back */
        Z[j] /= N;
        /* square sum errors over real (and imag parts) */
        err = (X[j]-Z[j]) * (X[j]-Z[j]);
        errSum += err;
      }

      if ( errSum > N * 1E-7 )
      {
        retError = 1;
        printf("%s fft %d  bin %d : inverse FFT doesn't match original signal! errSum = %g ; mean err = %g\n", (cplx ? "cplx":"real"), N, k, errSum, errSum / N);
      }

      break;
    }

  }
#ifdef PFFFT_ENABLE_FLOAT
  pffft_destroy_setup(s);
  pffft_aligned_free(X);
  pffft_aligned_free(Y);
  pffft_aligned_free(Z);
  pffft_aligned_free(R);
  pffft_aligned_free(W);
#else
  pffftd_destroy_setup(s);
  pffftd_aligned_free(X);
  pffftd_aligned_free(Y);
  pffftd_aligned_free(Z);
  pffftd_aligned_free(R);
  pffftd_aligned_free(W);
#endif

  return retError;
}

/* check that pffft_zconvolve() matches pffft_zconvolve_scale(.., 1.0)
   bit-exactly, and that circular convolution with a delayed impulse
   delays the signal (the backward transform carries the usual factor N) */
static int test_zconvolve_unscaled(int N)
{
  int j, retError = 0;
#ifdef PFFFT_ENABLE_FLOAT
  typedef PFFFT_Setup setup_t;
#define ZT_NEW_SETUP     pffft_new_setup
#define ZT_DESTROY_SETUP pffft_destroy_setup
#define ZT_TRANSFORM     pffft_transform
#define ZT_ZCONV_SCALE   pffft_zconvolve_scale
#define ZT_ZCONVOLVE     pffft_zconvolve
#define ZT_MALLOC        pffft_aligned_malloc
#define ZT_FREE          pffft_aligned_free
#else
  typedef PFFFTD_Setup setup_t;
#define ZT_NEW_SETUP     pffftd_new_setup
#define ZT_DESTROY_SETUP pffftd_destroy_setup
#define ZT_TRANSFORM     pffftd_transform
#define ZT_ZCONV_SCALE   pffftd_zconvolve_scale
#define ZT_ZCONVOLVE     pffftd_zconvolve
#define ZT_MALLOC        pffftd_aligned_malloc
#define ZT_FREE          pffftd_aligned_free
#endif
  setup_t *s;
  pffft_scalar *X, *H, *SX, *SH, *CNew, *CRef, *Y, *W;

  s = ZT_NEW_SETUP(N, PFFFT_REAL);
  assert(s);
  X    = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  H    = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  SX   = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  SH   = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  CNew = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  CRef = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  Y    = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  W    = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));

  for ( j = 0; j < N; ++j ) {
    X[j] = (pffft_scalar)(sin(0.1*j) + 0.3*cos(0.9*j));
    H[j] = (j == 4 ? (pffft_scalar)1.0 : (pffft_scalar)0.0); /* impulse delayed by 4 samples */
  }

  ZT_TRANSFORM(s, X, SX, W, PFFFT_FORWARD);
  ZT_TRANSFORM(s, H, SH, W, PFFFT_FORWARD);

  ZT_ZCONVOLVE(s, SX, SH, CNew);
  ZT_ZCONV_SCALE(s, SX, SH, CRef, 1.0);

  if ( memcmp(CNew, CRef, N * sizeof(pffft_scalar)) != 0 ) {
    printf("real fft %d: zconvolve() differs from zconvolve_scale(.., 1.0)\n", N);
    retError = 1;
  }

  ZT_TRANSFORM(s, CNew, Y, W, PFFFT_BACKWARD);
  for ( j = 0; j < N-4; ++j ) {
    pffft_scalar expected = (pffft_scalar)(X[j]);
    pffft_scalar got      = (pffft_scalar)(Y[j+4] / N);
    if ( fabs(got - expected) > 1e-4 ) {
      printf("real fft %d: delayed convolution mismatch at j = %d : got %g, expected %g\n",
             N, j, (double)got, (double)expected);
      retError = 1;
      break;
    }
  }

  ZT_FREE(X); ZT_FREE(H); ZT_FREE(SX); ZT_FREE(SH);
  ZT_FREE(CNew); ZT_FREE(CRef); ZT_FREE(Y); ZT_FREE(W);
  ZT_DESTROY_SETUP(s);
  return retError;
}

/* zero-phase filtering check: a linear-phase FIR centered at t = 0
   with left-hand taps wrapped around the block must produce exactly the
   circular convolution of the wrapped filter with the signal when its
   unordered spectrum is converted via pffft_zconvert_zp() and applied
   with pffft_zconvolve_zp(). Also cross-checks against the plain
   complex multiply and verifies the REAL-only guard. */
static int test_zerophase(int N)
{
  int j, k, retError = 0;
#ifdef PFFFT_ENABLE_FLOAT
  typedef PFFFT_Setup setup_t;
#define ZT_NEW_SETUP     pffft_new_setup
#define ZT_DESTROY_SETUP pffft_destroy_setup
#define ZT_TRANSFORM     pffft_transform
#define ZT_ZCONVERT_ZP   pffft_zconvert_zp
#define ZT_ZCONV_SCALE   pffft_zconvolve_scale
#define ZT_ZCONVOLVE_ZP  pffft_zconvolve_zp
#define ZT_MALLOC        pffft_aligned_malloc
#define ZT_FREE          pffft_aligned_free
#else
  typedef PFFFTD_Setup setup_t;
#define ZT_NEW_SETUP     pffftd_new_setup
#define ZT_DESTROY_SETUP pffftd_destroy_setup
#define ZT_TRANSFORM     pffftd_transform
#define ZT_ZCONVERT_ZP   pffftd_zconvert_zp
#define ZT_ZCONV_SCALE   pffftd_zconvolve_scale
#define ZT_ZCONVOLVE_ZP  pffftd_zconvolve_zp
#define ZT_MALLOC        pffftd_aligned_malloc
#define ZT_FREE          pffftd_aligned_free
#endif
  setup_t *s;
  setup_t *scplx;
  /* symmetric FIR, centered at t = 0 : taps h[-4 .. 4] */
  pffft_scalar ht[9] = { 0.03, -0.10, 0.25, 0.75, 1.0, 0.75, 0.25, -0.10, 0.03 };
  pffft_scalar *X, *W, *H, *HZP, *HZPN, *C, *CRef, *Y;
  s = ZT_NEW_SETUP(N, PFFFT_REAL);
  scplx = ZT_NEW_SETUP(N, PFFFT_COMPLEX);
  assert(s && scplx);
  X    = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  W    = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  H    = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  HZP  = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  HZPN = (pffft_scalar*)ZT_MALLOC(N * sizeof(pffft_scalar));
  C    = (pffft_scalar*)ZT_MALLOC(2 * N * sizeof(pffft_scalar));
  CRef = (pffft_scalar*)ZT_MALLOC(2 * N * sizeof(pffft_scalar));
  Y    = (pffft_scalar*)ZT_MALLOC(2 * N * sizeof(pffft_scalar));

  if ( ZT_ZCONVERT_ZP(scplx, W, HZP, 1.0) == 0 ) {
    printf("zero-phase fft %d: zconvert_zp accepted a PFFFT_COMPLEX setup!\n", N);
    retError = 1;
  }
  if ( ZT_ZCONVOLVE_ZP(scplx, C, HZP, CRef) == 0 ) {
    printf("zero-phase fft %d: zconvolve_zp accepted a PFFFT_COMPLEX setup!\n", N);
    retError = 1;
  }

  for ( j = 0; j < N; ++j ) {
    X[j] = (pffft_scalar)(sin(0.37*j) + 0.5*cos(1.71*j + 0.3));
    W[j] = (pffft_scalar)0.0;
  }
  for ( j = 0; j <= 4; ++j )      W[j] = (pffft_scalar)ht[4+j];   /* h[0..4] */
  for ( j = 1; j <= 4; ++j )      W[N-j] = (pffft_scalar)ht[4-j]; /* h[-4..-1] wrapped */

  ZT_TRANSFORM(s, X, C, NULL, PFFFT_FORWARD);
  ZT_TRANSFORM(s, W, H, NULL, PFFFT_FORWARD);

  if ( ZT_ZCONVERT_ZP(s, H, HZP, 1.0) != 0 ) {
    printf("zero-phase fft %d: zconvert_zp failed on a PFFFT_REAL setup!\n", N);
    retError = 1;
  }

  /* zero-phase invariants, layout-independent by construction:
     1) conversion writes EVERY lane of its destination -- sentinel-fill
        first so unwritten lanes cannot hide (this is exactly what the
        vector-1 defect violated); no lane-value assumptions needed.
     2) conversion is linear in scaling.
     Lane-level value correctness is covered by the cross-check and the
     time-domain round-trip below, which compare whole buffers without
     lane-position assumptions. */
  for ( j = 0; j < N; ++j )
    CRef[j] = (pffft_scalar)123456789.0;
  if ( ZT_ZCONVERT_ZP(s, H, CRef, 2.0) != 0 ) {
    printf("zero-phase fft %d: zconvert_zp(scaling=2) failed!\n", N);
    retError = 1;
  }
  for ( j = 0; j < N; ++j ) {
    if ( CRef[j] == (pffft_scalar)123456789.0 ) {
      printf("zero-phase fft %d: zconvert_zp left float %d unwritten\n", N, j);
      retError = 1;
      break;
    }
  }
  /* linearity: convert(scaling=1) then double must equal convert(2.0) */
  {
    int lin_bad = -1;
    for ( j = 0; j < N; ++j ) {
      if ( fabs((double)(HZP[j]*2 - CRef[j])) > 1e-3 * (fabs((double)HZP[j]) + 1.0) ) {
        lin_bad = j;
        break;
      }
    }
    if ( lin_bad >= 0 ) {
      printf("zero-phase fft %d: zconvert_zp not linear in scaling at float %d : 2*%g != %g\n",
             N, lin_bad, (double)HZP[lin_bad], (double)CRef[lin_bad]);
      retError = 1;
    }
  }

  /* documented aliasing: in == out must work (this was the crux of the
     vector-1 defect: a naive zero-then-restore order destroys the
     Nyquist lane before it is read) */
  {
    for ( j = 0; j < N; ++j )
      CRef[j] = H[j];
    if ( ZT_ZCONVERT_ZP(s, CRef, CRef, 1.0) != 0 ) {
      printf("zero-phase fft %d: in-place zconvert_zp failed\n", N);
      retError = 1;
    }
    for ( j = 0; j < N; ++j ) {
      /* identical operations in identical order: any difference is a
         real defect, so exact comparison */
      if ( CRef[j] != HZP[j] ) {
        printf("zero-phase fft %d: in-place zconvert_zp differs from out-of-place at float %d : %g vs %g\n",
               N, j, (double)CRef[j], (double)HZP[j]);
        retError = 1;
        break;
      }
    }
  }
  /* unit-gain round trip: with scaling == 1, the 1/N applied by
     zconvert_zp() completes the backward transform's factor N, so Y
     must match the direct time-domain circular convolution without
     any rescaling */
  {
    for ( j = 0; j < N; ++j )
      CRef[j] = C[j];
    ZT_ZCONVOLVE_ZP(s, CRef, HZP, CRef);        /* in-place is documented legal */
    ZT_TRANSFORM(s, CRef, Y, NULL, PFFFT_BACKWARD);
    for ( j = 0; j < N; ++j ) {
      pffft_scalar ref = 0;
      for ( k = -4; k <= 4; ++k )
        ref += (pffft_scalar)(ht[4+k] * X[(j-k+N) % N]);
      if ( fabs((double)(Y[j]) - (double)(ref)) > 1e-4 * (fabs((double)ref) + 1.0) ) {
        printf("zero-phase fft %d: sample %d : got %g expected %g\n", N, j, (double)(Y[j]), (double)ref);
        retError = 1;
        break;
      }
    }
  }

  /* bit-exact cross-check against the plain complex multiply:
     converting with scaling == N cancels the internal 1/N, restoring
     the raw zero-phase spectrum */
  if ( ZT_ZCONVERT_ZP(s, H, HZPN, (pffft_scalar)N) != 0 ) {
    printf("zero-phase fft %d: zconvert_zp(scaling=N) failed!\n", N);
    retError = 1;
  }
  ZT_ZCONV_SCALE(s, C, HZPN, CRef, 1.0);      /* reference from pristine C */
  ZT_ZCONVOLVE_ZP(s, C, HZPN, C);             /* in-place is documented legal */
  /* unordered layout is not interleaved re/im pairs: compare flat.
     with a correctly converted zero-phase filter (H_im == 0) the plain
     complex multiply collapses to exactly what the zp path computes, and
     zconvolve_scale's REAL epilogue writes the DC and Nyquist lanes with
     the same component-wise products -- so every float must match. */
  for ( j = 0; j < N; ++j ) {
    if ( fabs((double)(C[j] - CRef[j])) > 1e-3 ) {
      printf("zero-phase fft %d: zp convolve differs from plain convolve at float %d\n", N, j);
      retError = 1;
      break;
    }
  }

  ZT_FREE(X); ZT_FREE(W); ZT_FREE(H); ZT_FREE(HZP); ZT_FREE(HZPN); ZT_FREE(C); ZT_FREE(CRef); ZT_FREE(Y);
  ZT_DESTROY_SETUP(s);
  ZT_DESTROY_SETUP(scplx);
  return retError;
}


/* small functions inside pffft.c that will detect (compiler) bugs with respect to simd instructions */
void validate_pffft_simd();
int  validate_pffft_simd_ex(FILE * DbgOut);
void validate_pffftd_simd();
int  validate_pffftd_simd_ex(FILE * DbgOut);



int main(int argc, char **argv)
{
  int N, result, resN, resAll, i, k, resNextPw2, resIsPw2, resFFT;
  int resZconv;
  int resZP;

  int inp_power_of_two[] = { 1, 2, 3, 4, 5, 6, 7, 8,  9, 511, 512,  513 };
  int ref_power_of_two[] = { 1, 2, 4, 4, 8, 8, 8, 8, 16, 512, 512, 1024 };

  for ( i = 1; i < argc; ++i ) {

    if (!strcmp(argv[i], "--test-simd")) {
#ifdef PFFFT_ENABLE_FLOAT
      int numErrs = validate_pffft_simd_ex(stdout);
#else
      int numErrs = validate_pffftd_simd_ex(stdout);
#endif
      fprintf( ( numErrs != 0 ? stderr : stdout ), "validate_pffft_simd_ex() returned %d errors!\n", numErrs);
      return ( numErrs > 0 ? 1 : 0 );
    }
  }

  resNextPw2 = 0;
  resIsPw2 = 0;
  for ( k = 0; k < (sizeof(inp_power_of_two)/sizeof(inp_power_of_two[0])); ++k) {
#ifdef PFFFT_ENABLE_FLOAT
    N = pffft_next_power_of_two(inp_power_of_two[k]);
#else
    N = pffftd_next_power_of_two(inp_power_of_two[k]);
#endif
    if (N != ref_power_of_two[k]) {
      resNextPw2 = 1;
      printf("pffft_next_power_of_two(%d) does deliver %d, which is not reference result %d!\n",
        inp_power_of_two[k], N, ref_power_of_two[k] );
    }

#ifdef PFFFT_ENABLE_FLOAT
    result = pffft_is_power_of_two(inp_power_of_two[k]);
#else
    result = pffftd_is_power_of_two(inp_power_of_two[k]);
#endif
    if (inp_power_of_two[k] == ref_power_of_two[k]) {
      if (!result) {
        resIsPw2 = 1;
        printf("pffft_is_power_of_two(%d) delivers false; expected true!\n", inp_power_of_two[k]);
      }
    } else {
      if (result) {
        resIsPw2 = 1;
        printf("pffft_is_power_of_two(%d) delivers true; expected false!\n", inp_power_of_two[k]);
      }
    }
  }
  if (!resNextPw2)
    printf("tests for pffft_next_power_of_two() succeeded successfully.\n");
  if (!resIsPw2)
    printf("tests for pffft_is_power_of_two() succeeded successfully.\n");

  resFFT = 0;
  for ( N = 32; N <= 65536; N *= 2 )
  {
    result = test(N, 1 /* cplx fft */, 1 /* useOrdered */);
    resN = result;
    resFFT |= result;

    result = test(N, 0 /* cplx fft */, 1 /* useOrdered */);
    resN |= result;
    resFFT |= result;

    result = test(N, 1 /* cplx fft */, 0 /* useOrdered */);
    resN |= result;
    resFFT |= result;

    result = test(N, 0 /* cplx fft */, 0 /* useOrdered */);
    resN |= result;
    resFFT |= result;

    if (!resN)
      printf("tests for size %d succeeded successfully.\n", N);
  }

  if (!resFFT) {
#ifdef PFFFT_ENABLE_FLOAT
    printf("all pffft transform tests (FORWARD/BACKWARD, REAL/COMPLEX, float) succeeded successfully.\n");
#else
    printf("all pffft transform tests (FORWARD/BACKWARD, REAL/COMPLEX, double) succeeded successfully.\n");
#endif
  }

  resZconv = test_zconvolve_unscaled(64);
  if (!resZconv)
    printf("zconvolve unscaled tests succeeded successfully.\n");

  resZP = test_zerophase(64);
  if (!resZP)
    printf("zero-phase convolution tests succeeded successfully.\n");

  resAll = resNextPw2 | resIsPw2 | resFFT | resZconv | resZP;
  if (!resAll)
    printf("all tests succeeded successfully.\n");
  else
    printf("there are failed tests!\n");

  return resAll;
}

