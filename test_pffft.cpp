/*
  Copyright (c) 2013 Julien Pommier.

  Small test & bench for PFFFT, comparing its performance with the scalar
  FFTPACK, FFTW, and Apple vDSP

  How to build:

  on linux, with fftw3:
  gcc -o test_pffft -DHAVE_FFTW -msse -mfpmath=sse -O3 -Wall -W pffft.c
  test_pffft.c fftpack.c -L/usr/local/lib -I/usr/local/include/ -lfftw3f -lm

  on macos, without fftw3:
  clang -o test_pffft -DHAVE_VECLIB -O3 -Wall -W pffft.c test_pffft.c fftpack.c
  -L/usr/local/lib -I/usr/local/include/ -framework Accelerate

  on macos, with fftw3:
  clang -o test_pffft -DHAVE_FFTW -DHAVE_VECLIB -O3 -Wall -W pffft.c
  test_pffft.c fftpack.c -L/usr/local/lib -I/usr/local/include/ -lfftw3f
  -framework Accelerate

  as alternative: replace clang by gcc.

  on windows, with visual c++:
  cl /Ox -D_USE_MATH_DEFINES /arch:SSE test_pffft.c pffft.c fftpack.c

  build without SIMD instructions:
  gcc -o test_pffft -DPFFFT_SIMD_DISABLE -O3 -Wall -W pffft.c test_pffft.c
  fftpack.c -lm

 */

#include "pffft.hpp"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* maximum allowed phase error in degree */
#define DEG_ERR_LIMIT 1E-4

/* maximum allowed magnitude error in amplitude (of 1.0 or 1.1) */
#define MAG_ERR_LIMIT 1E-6

#define PRINT_SPEC 0

#define PWR2LOG(PWR) ((PWR) < 1E-30 ? 10.0 * log10(1E-30) : 10.0 * log10(PWR))

template<typename T>
bool
Ttest(int N, bool useOrdered)
{
  using Fft = typename pffft::Fft<T>;
  using Scalar = typename Fft::Scalar;

  bool cplx = std::is_same<T, std::complex<float>>::value ||
              std::is_same<T, std::complex<double>>::value;

  double EXPECTED_DYN_RANGE =
    std::is_same<double, Scalar>::value ? 215.0 : 140.0;

  int Nsca = (cplx ? N * 2 : N);
  int Ncplx = (cplx ? N : N / 2);
  T* X = Fft::alignedAlloc<T>(Nsca);
  T* Z = Fft::alignedAlloc<T>(Nsca);
  Scalar* R = Fft::alignedAllocScalar(Nsca);
  std::complex<Scalar>* Y = Fft::alignedAllocComplex(Nsca);
  int k, j, m, iter, kmaxOther;
  bool retError = false;
  double freq, dPhi, phi, phi0;
  double pwr, pwrCar, pwrOther, err, errSum, mag, expextedMag;
  double amp = 1.0;

  assert(pffft::isPowerOfTwo(N));

  Fft fft = Fft(N);

  Scalar* Xs = reinterpret_cast<Scalar*>(X);
  Scalar* Ys = reinterpret_cast<Scalar*>(Y);
  Scalar* Zs = reinterpret_cast<Scalar*>(Z);

  for (k = m = 0; k < (cplx ? N : (1 + N / 2)); k += N / 16, ++m) {
    amp = ((m % 3) == 0) ? 1.0F : 1.1F;
    freq = (k < N / 2) ? ((double)k / N) : ((double)(k - N) / N);
    dPhi = 2.0 * M_PI * freq;
    if (dPhi < 0.0)
      dPhi += 2.0 * M_PI;

    iter = -1;
    while (1) {
      ++iter;

      if (iter)
        printf("bin %d: dphi = %f for freq %f\n", k, dPhi, freq);

      /* generate cosine carrier as time signal - start at defined phase phi0 */
      phi = phi0 =
        (m % 4) * 0.125 * M_PI; /* have phi0 < 90 deg to be normalized */
      for (j = 0; j < N; ++j) {
        if (cplx) {
          Xs[2 * j] = amp * cos(phi);     /* real part */
          Xs[2 * j + 1] = amp * sin(phi); /* imag part */
        } else
          Xs[j] = amp * cos(phi); /* only real part */

        /* phase increment .. stay normalized - cos()/sin() might degrade! */
        phi += dPhi;
        if (phi >= M_PI)
          phi -= 2.0 * M_PI;
      }

      /* forward transform from X --> Y  .. using work buffer W */
      if (useOrdered)
        fft.forward(X, Y);
      else {
        fft.forwardInternalLayout(X, R); /* temporarily use R for reordering */
        fft.reorderSpectrum(R, Y, PFFFT_FORWARD);
      }

      pwrOther = -1.0;
      pwrCar = 0;

      /* for positive frequencies: 0 to 0.5 * samplerate */
      /* and also for negative frequencies: -0.5 * samplerate to 0 */
      for (j = 0; j < (cplx ? N : (1 + N / 2)); ++j) {
        if (!cplx && !j) /* special treatment for DC for real input */
          pwr = Ys[j] * Ys[j];
        else if (!cplx && j == N / 2) /* treat 0.5 * samplerate */
          pwr = Ys[1] *
                Ys[1]; /* despite j (for freq calculation) we have index 1 */
        else
          pwr = Ys[2 * j] * Ys[2 * j] + Ys[2 * j + 1] * Ys[2 * j + 1];
        if (iter || PRINT_SPEC)
          printf("%s fft %d:  pwr[j = %d] = %g == %f dB\n",
                 (cplx ? "cplx" : "real"),
                 N,
                 j,
                 pwr,
                 PWR2LOG(pwr));
        if (k == j)
          pwrCar = pwr;
        else if (pwr > pwrOther) {
          pwrOther = pwr;
          kmaxOther = j;
        }
      }

      if (PWR2LOG(pwrCar) - PWR2LOG(pwrOther) < EXPECTED_DYN_RANGE) {
        printf("%s fft %d amp %f iter %d:\n",
               (cplx ? "cplx" : "real"),
               N,
               amp,
               iter);
        printf("  carrier power  at bin %d: %g == %f dB\n",
               k,
               pwrCar,
               PWR2LOG(pwrCar));
        printf("  carrier mag || at bin %d: %g\n", k, sqrt(pwrCar));
        printf("  max other pwr  at bin %d: %g == %f dB\n",
               kmaxOther,
               pwrOther,
               PWR2LOG(pwrOther));
        printf("  dynamic range: %f dB\n\n",
               PWR2LOG(pwrCar) - PWR2LOG(pwrOther));
        retError = true;
        if (iter == 0)
          continue;
      }

      if (k > 0 && k != N / 2) {
        phi = atan2(Ys[2 * k + 1], Ys[2 * k]);
        if (fabs(phi - phi0) > DEG_ERR_LIMIT * M_PI / 180.0) {
          retError = true;
          printf("%s fft %d  bin %d amp %f : phase mismatch! phase = %f deg   "
                 "expected = %f deg\n",
                 (cplx ? "cplx" : "real"),
                 N,
                 k,
                 amp,
                 phi * 180.0 / M_PI,
                 phi0 * 180.0 / M_PI);
        }
      }

      expextedMag = cplx ? amp : ((k == 0 || k == N / 2) ? amp : (amp / 2));
      mag = sqrt(pwrCar) / N;
      if (fabs(mag - expextedMag) > MAG_ERR_LIMIT) {
        retError = true;
        printf("%s fft %d  bin %d amp %f : mag = %g   expected = %g\n",
               (cplx ? "cplx" : "real"),
               N,
               k,
               amp,
               mag,
               expextedMag);
      }

      /* now convert spectrum back */
      fft.inverse(Y, Z);

      errSum = 0.0;
      for (j = 0; j < (cplx ? (2 * N) : N); ++j) {
        /* scale back */
        Z[j] /= N;
        /* square sum errors over real (and imag parts) */
        err = (Xs[j] - Zs[j]) * (Xs[j] - Zs[j]);
        errSum += err;
      }

      if (errSum > N * 1E-7) {
        retError = true;
        printf("%s fft %d  bin %d : inverse FFT doesn't match original signal! "
               "errSum = %g ; mean err = %g\n",
               (cplx ? "cplx" : "real"),
               N,
               k,
               errSum,
               errSum / N);
      }

      break;
    }
  }
  pffft::alignedFree(X);
  pffft::alignedFree(Y);
  pffft::alignedFree(Z);

  return retError;
}

bool
test(int N, bool useComplex, bool useOrdered)
{
  if (useComplex) {
    return Ttest<std::complex<float>>(N, useOrdered) &&
           Ttest<std::complex<double>>(N, useOrdered);
  } else {
    return Ttest<float>(N, useOrdered) && Ttest<double>(N, useOrdered);
  }
}

int
main(int argc, char** argv)
{
  int N, result, resN, resAll, k, resNextPw2, resIsPw2, resFFT;

  int inp_power_of_two[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 511, 512, 513 };
  int ref_power_of_two[] = { 1, 2, 4, 4, 8, 8, 8, 8, 16, 512, 512, 1024 };

  resNextPw2 = 0;
  resIsPw2 = 0;
  for (k = 0; k < (sizeof(inp_power_of_two) / sizeof(inp_power_of_two[0]));
       ++k) {
    N = pffft::nextPowerOfTwo(inp_power_of_two[k]);
    if (N != ref_power_of_two[k]) {
      resNextPw2 = 1;
      printf("pffft_next_power_of_two(%d) does deliver %d, which is not "
             "reference result %d!\n",
             inp_power_of_two[k],
             N,
             ref_power_of_two[k]);
    }

    result = pffft::isPowerOfTwo(inp_power_of_two[k]);
    if (inp_power_of_two[k] == ref_power_of_two[k]) {
      if (!result) {
        resIsPw2 = 1;
        printf("pffft_is_power_of_two(%d) delivers false; expected true!\n",
               inp_power_of_two[k]);
      }
    } else {
      if (result) {
        resIsPw2 = 1;
        printf("pffft_is_power_of_two(%d) delivers true; expected false!\n",
               inp_power_of_two[k]);
      }
    }
  }
  if (!resNextPw2)
    printf("tests for pffft_next_power_of_two() succeeded successfully.\n");
  if (!resIsPw2)
    printf("tests for pffft_is_power_of_two() succeeded successfully.\n");

  resFFT = 0;
  for (N = 32; N <= 65536; N *= 2) {
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

  if (!resFFT)
    printf("all pffft transform tests (FORWARD/BACKWARD, REAL/COMPLEX) "
           "succeeded successfully.\n");

  resAll = resNextPw2 | resIsPw2 | resFFT;
  if (!resAll)
    printf("all tests succeeded successfully.\n");
  else
    printf("there are failed tests!\n");

  return resAll;
}
