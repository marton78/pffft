/*
  Copyright (c) 2013 Julien Pommier.

  Small test & bench for PFFFT, comparing its performance with the scalar FFTPACK, FFTW, and Apple vDSP

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

#include "pffft.h"
#include "fftpack.h"

#ifdef HAVE_GREEN_FFTS
#include "fftext.h"
#endif

#ifdef HAVE_KISS_FFT
#include <kiss_fft.h>
#include <kiss_fftr.h>
#endif


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#ifdef HAVE_SYS_TIMES
#  include <sys/times.h>
#  include <unistd.h>
#endif

#ifdef HAVE_VECLIB
#  include <Accelerate/Accelerate.h>
#endif

#ifdef HAVE_FFTW
#  include <fftw3.h>
#endif


#define NUM_FFT_ALGOS  7
enum {
  ALGO_FFTPACK = 0,
  ALGO_VECLIB,
  ALGO_FFTW_ESTIM,
  ALGO_FFTW_AUTO,
  ALGO_GREEN,
  ALGO_KISS,
  ALGO_PFFFT /* = 6 */
};

#define NUM_TYPES      4
enum {
  TYPE_PREP = 0,         /* time for preparation in ms */
  TYPE_DUR_NS = 1,       /* time per fft in ns */
  TYPE_DUR_FASTEST = 2,  /* relative time to fastest */
  TYPE_REL_PFFFT = 3     /* relative time to ALGO_PFFFT */
};
// double tmeas[NUM_TYPES][NUM_FFT_ALGOS];

const char * algoName[NUM_FFT_ALGOS] = {
  "FFTPack    ",
  "vDSP (vec) ",
  "FFTW(estim)",
  "FFTW (auto)",
  "Green      ",
  "Kiss       ",
  "PFFFT(simd)"
};


int compiledInAlgo[NUM_FFT_ALGOS] = {
  1, /* "FFTPack    " */
#ifdef HAVE_VECLIB
  1, /* "vDSP (vec) " */
#else
  0,
#endif
#ifdef HAVE_FFTW
  1, /* "FFTW(estim)" */
  1, /* "FFTW (auto)" */
#else
  0, 0,
#endif
#ifdef HAVE_GREEN_FFTS
  1, /* "Green      " */
#else
  0,
#endif
#ifdef HAVE_KISS_FFT
  1, /* "Kiss       " */
#else
  0,
#endif
  1  /* "PFFFT      " */
};

const char * algoTableHeader[NUM_FFT_ALGOS][2] = {
{ "| real FFTPack ", "| cplx FFTPack " },
{ "|  real   vDSP ", "|  cplx   vDSP " },
{ "|real FFTWestim", "|cplx FFTWestim" },
{ "|real FFTWauto ", "|cplx FFTWauto " },
{ "|  real  Green ", "|  cplx  Green " },
{ "|  real   Kiss ", "|  cplx   Kiss " },
{ "|  real  PFFFT ", "|  cplx  PFFFT " } };

const char * typeText[NUM_TYPES] = {
  "preparation in ms",
  "time per fft in ns",
  "relative to fastest",
  "relative to pffft"
};


#define MAX(x,y) ((x)>(y)?(x):(y))

double frand() {
  return rand()/(double)RAND_MAX;
}

#if defined(HAVE_SYS_TIMES)
  inline double uclock_sec(void) {
    static double ttclk = 0.;
    if (ttclk == 0.) ttclk = sysconf(_SC_CLK_TCK);
    struct tms t; return ((double)times(&t)) / ttclk;
  }
# else
  double uclock_sec(void)
{ return (double)clock()/(double)CLOCKS_PER_SEC; }
#endif


/* compare results with the regular fftpack */
void pffft_validate_N(int N, int cplx) {
  int Nfloat = N*(cplx?2:1);
  int Nbytes = Nfloat * sizeof(float);
  float *ref, *in, *out, *tmp, *tmp2;
  PFFFT_Setup *s = pffft_new_setup(N, cplx ? PFFFT_COMPLEX : PFFFT_REAL);
  int pass;

  if (!s) { printf("Skipping N=%d, not supported\n", N); return; }
  ref = pffft_aligned_malloc(Nbytes);
  in = pffft_aligned_malloc(Nbytes);
  out = pffft_aligned_malloc(Nbytes);
  tmp = pffft_aligned_malloc(Nbytes);
  tmp2 = pffft_aligned_malloc(Nbytes);

  for (pass=0; pass < 2; ++pass) {
    float ref_max = 0;
    int k;
    //printf("N=%d pass=%d cplx=%d\n", N, pass, cplx);
    // compute reference solution with FFTPACK
    if (pass == 0) {
      float *wrk = malloc(2*Nbytes+15*sizeof(float));
      for (k=0; k < Nfloat; ++k) {
        ref[k] = in[k] = frand()*2-1; 
        out[k] = 1e30;
      }
      if (!cplx) {
        rffti(N, wrk);
        rfftf(N, ref, wrk);
        // use our ordering for real ffts instead of the one of fftpack
        {
          float refN=ref[N-1];
          for (k=N-2; k >= 1; --k) ref[k+1] = ref[k]; 
          ref[1] = refN;
        }
      } else {
        cffti(N, wrk);
        cfftf(N, ref, wrk);
      }
      free(wrk);
    }

    for (k = 0; k < Nfloat; ++k) ref_max = MAX(ref_max, fabs(ref[k]));

      
    // pass 0 : non canonical ordering of transform coefficients  
    if (pass == 0) {
      // test forward transform, with different input / output
      pffft_transform(s, in, tmp, 0, PFFFT_FORWARD);
      memcpy(tmp2, tmp, Nbytes);
      memcpy(tmp, in, Nbytes);
      pffft_transform(s, tmp, tmp, 0, PFFFT_FORWARD);
      for (k = 0; k < Nfloat; ++k) {
        assert(tmp2[k] == tmp[k]);
      }

      // test reordering
      pffft_zreorder(s, tmp, out, PFFFT_FORWARD);
      pffft_zreorder(s, out, tmp, PFFFT_BACKWARD);
      for (k = 0; k < Nfloat; ++k) {
        assert(tmp2[k] == tmp[k]);
      }
      pffft_zreorder(s, tmp, out, PFFFT_FORWARD);
    } else {
      // pass 1 : canonical ordering of transform coeffs.
      pffft_transform_ordered(s, in, tmp, 0, PFFFT_FORWARD);
      memcpy(tmp2, tmp, Nbytes);
      memcpy(tmp, in, Nbytes);
      pffft_transform_ordered(s, tmp, tmp, 0, PFFFT_FORWARD);
      for (k = 0; k < Nfloat; ++k) {
        assert(tmp2[k] == tmp[k]);
      }
      memcpy(out, tmp, Nbytes);
    }

    {
      for (k=0; k < Nfloat; ++k) {
        if (!(fabs(ref[k] - out[k]) < 1e-3*ref_max)) {
          printf("%s forward PFFFT mismatch found for N=%d\n", (cplx?"CPLX":"REAL"), N);
          exit(1);
        }
      }
        
      if (pass == 0) pffft_transform(s, tmp, out, 0, PFFFT_BACKWARD);
      else   pffft_transform_ordered(s, tmp, out, 0, PFFFT_BACKWARD);
      memcpy(tmp2, out, Nbytes);
      memcpy(out, tmp, Nbytes);
      if (pass == 0) pffft_transform(s, out, out, 0, PFFFT_BACKWARD);
      else   pffft_transform_ordered(s, out, out, 0, PFFFT_BACKWARD);
      for (k = 0; k < Nfloat; ++k) {
        assert(tmp2[k] == out[k]);
        out[k] *= 1.f/N;
      }
      for (k = 0; k < Nfloat; ++k) {
        if (fabs(in[k] - out[k]) > 1e-3 * ref_max) {
          printf("pass=%d, %s IFFFT does not match for N=%d\n", pass, (cplx?"CPLX":"REAL"), N); break;
          exit(1);
        }
      }
    }

    // quick test of the circular convolution in fft domain
    {
      float conv_err = 0, conv_max = 0;

      pffft_zreorder(s, ref, tmp, PFFFT_FORWARD);
      memset(out, 0, Nbytes);
      pffft_zconvolve_accumulate(s, ref, ref, out, 1.0);
      pffft_zreorder(s, out, tmp2, PFFFT_FORWARD);
      
      for (k=0; k < Nfloat; k += 2) {
        float ar = tmp[k], ai=tmp[k+1];
        if (cplx || k > 0) {
          tmp[k] = ar*ar - ai*ai;
          tmp[k+1] = 2*ar*ai;
        } else {
          tmp[0] = ar*ar;
          tmp[1] = ai*ai;
        }
      }
      
      for (k=0; k < Nfloat; ++k) {
        float d = fabs(tmp[k] - tmp2[k]), e = fabs(tmp[k]);
        if (d > conv_err) conv_err = d;
        if (e > conv_max) conv_max = e;
      }
      if (conv_err > 1e-5*conv_max) {
        printf("zconvolve error ? %g %g\n", conv_err, conv_max); exit(1);
      }
    }

  }

  printf("%s PFFFT is OK for N=%d\n", (cplx?"CPLX":"REAL"), N); fflush(stdout);
  
  pffft_destroy_setup(s);
  pffft_aligned_free(ref);
  pffft_aligned_free(in);
  pffft_aligned_free(out);
  pffft_aligned_free(tmp);
  pffft_aligned_free(tmp2);
}

void pffft_validate(int cplx) {
  static int Ntest[] = { 16, 32, 64, 96, 128, 160, 192, 256, 288, 384, 5*96, 512, 576, 5*128, 800, 864, 1024, 2048, 2592, 4000, 4096, 12000, 36864, 0};
  int k;
  for (k = 0; Ntest[k]; ++k) {
    int N = Ntest[k];
    if (N == 16 && !cplx) continue;
    pffft_validate_N(N, cplx);
  }
}

int array_output_format = 1;

double show_output(const char *name, int N, int cplx, float flops, float t0, float t1, int max_iter) {
  double T = (double)(t1-t0)/2/max_iter * 1e9;
  float mflops = flops/1e6/(t1 - t0 + 1e-16);
  if (array_output_format) {
    if (flops != -1) {
      printf("|%11.0f   ", mflops);
    } else printf("|      n/a     ");
  } else {
    if (flops != -1) {
      printf("N=%5d, %s %16s : %6.0f MFlops [t=%6.0f ns, %d runs]\n", N, (cplx?"CPLX":"REAL"), name, mflops, (t1-t0)/2/max_iter * 1e9, max_iter);
    }
  }
  fflush(stdout);
  return T;
}

void benchmark_ffts(int Nlog, int N, int cplx, int withFFTWmeas, double tmeas[NUM_TYPES][NUM_FFT_ALGOS], int haveAlgo[NUM_FFT_ALGOS] ) {
  int Nfloat = (cplx ? N*2 : N);
  int Nbytes = Nfloat * sizeof(float);
  float *X = pffft_aligned_malloc(Nbytes), *Y = pffft_aligned_malloc(Nbytes + 2*sizeof(float) ), *Z = pffft_aligned_malloc(Nbytes);
  double te, t0, t1, flops, Tfastest;

  int k;
  int max_iter = 5120000/N*4;
#ifdef __arm__
  max_iter /= 4;
#endif
  int iter;

  for (k = 0; k < Nfloat; ++k) {
    X[k] = sqrtf(k+1);
  }

  for ( k = 0; k < NUM_TYPES; ++k )
  {
    for ( iter = 0; iter < NUM_FFT_ALGOS; ++iter )
      tmeas[k][iter] = 0.0;
  }


  // FFTPack benchmark
  {
    float *wrk = malloc(2*Nbytes + 15*sizeof(float));
    int max_iter_ = max_iter/pffft_simd_size(); if (max_iter_ == 0) max_iter_ = 1;
    te = uclock_sec();  
    if (cplx) cffti(N, wrk);
    else      rffti(N, wrk);
    t0 = uclock_sec();  
    
    for (iter = 0; iter < max_iter_; ++iter) {
      if (cplx) {
        cfftf(N, X, wrk);
        cfftb(N, X, wrk);
      } else {
        rfftf(N, X, wrk);
        rfftb(N, X, wrk);
      }
    }
    t1 = uclock_sec();
    free(wrk);
    
    flops = (max_iter_*2) * ((cplx ? 5 : 2.5)*N*log((double)N)/M_LN2); // see http://www.fftw.org/speed/method.html
    tmeas[TYPE_DUR_NS][ALGO_FFTPACK] = show_output("FFTPack", N, cplx, flops, t0, t1, max_iter_);
    tmeas[TYPE_PREP][ALGO_FFTPACK] = (t0 - te) * 1e3;
    haveAlgo[ALGO_FFTPACK] = 1;
  }


#ifdef HAVE_VECLIB
  int log2N = (int)(log(N)/log(2) + 0.5f);
  te = uclock_sec();  
  if (N == (1<<log2N)) {
    FFTSetup setup;

    setup = vDSP_create_fftsetup(log2N, FFT_RADIX2);
    DSPSplitComplex zsamples;
    zsamples.realp = &X[0];
    zsamples.imagp = &X[Nfloat/2];
    t0 = uclock_sec();  
    for (iter = 0; iter < max_iter; ++iter) {
      if (cplx) {
        vDSP_fft_zip(setup, &zsamples, 1, log2N, kFFTDirection_Forward);
        vDSP_fft_zip(setup, &zsamples, 1, log2N, kFFTDirection_Inverse);
      } else {
        vDSP_fft_zrip(setup, &zsamples, 1, log2N, kFFTDirection_Forward); 
        vDSP_fft_zrip(setup, &zsamples, 1, log2N, kFFTDirection_Inverse);
      }
    }
    t1 = uclock_sec();
    vDSP_destroy_fftsetup(setup);

    flops = (max_iter*2) * ((cplx ? 5 : 2.5)*N*log((double)N)/M_LN2); // see http://www.fftw.org/speed/method.html
    tmeas[TYPE_DUR_NS][ALGO_VECLIB] = show_output("vDSP", N, cplx, flops, t0, t1, max_iter);
    tmeas[TYPE_PREP][ALGO_VECLIB] = (t0 - te) * 1e3;
    haveAlgo[ALGO_VECLIB] = 1;
  } else {
    show_output("vDSP", N, cplx, -1, -1, -1, -1);
  }
#endif

#ifdef HAVE_FFTW
  {
    /* int flags = (N <= (256*1024) ? FFTW_MEASURE : FFTW_ESTIMATE);  measure takes a lot of time on largest ffts */
    int flags = FFTW_ESTIMATE;
    te = uclock_sec();
    fftwf_plan planf, planb;
    fftw_complex *in = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * N);
    memset(in, 0, sizeof(fftw_complex) * N);
    if (cplx) {
      planf = fftwf_plan_dft_1d(N, (fftwf_complex*)in, (fftwf_complex*)out, FFTW_FORWARD, flags);
      planb = fftwf_plan_dft_1d(N, (fftwf_complex*)in, (fftwf_complex*)out, FFTW_BACKWARD, flags);
    } else {
      planf = fftwf_plan_dft_r2c_1d(N, (float*)in, (fftwf_complex*)out, flags);
      planb = fftwf_plan_dft_c2r_1d(N, (fftwf_complex*)in, (float*)out, flags);
    }

    t0 = uclock_sec();  
    for (iter = 0; iter < max_iter; ++iter) {
      fftwf_execute(planf);
      fftwf_execute(planb);
    }
    t1 = uclock_sec();

    fftwf_destroy_plan(planf);
    fftwf_destroy_plan(planb);
    fftwf_free(in); fftwf_free(out);

    flops = (max_iter*2) * ((cplx ? 5 : 2.5)*N*log((double)N)/M_LN2); // see http://www.fftw.org/speed/method.html
    tmeas[TYPE_DUR_NS][ALGO_FFTW_ESTIM] = show_output((flags == FFTW_MEASURE ? algoName[ALGO_FFTW_AUTO] : algoName[ALGO_FFTW_ESTIM]), N, cplx, flops, t0, t1, max_iter);
    tmeas[TYPE_PREP][ALGO_FFTW_ESTIM] = (t0 - te) * 1e3;
    haveAlgo[ALGO_FFTW_ESTIM] = 1;
  }
  do {
    /* int flags = (N <= (256*1024) ? FFTW_MEASURE : FFTW_ESTIMATE);  measure takes a lot of time on largest ffts */
    /* int flags = FFTW_MEASURE; */
    int flags = (N < 40000 ? FFTW_MEASURE : (withFFTWmeas ? FFTW_MEASURE : FFTW_ESTIMATE));
    if (flags == FFTW_ESTIMATE) {
      show_output((flags == FFTW_MEASURE ? algoName[ALGO_FFTW_AUTO] : algoName[ALGO_FFTW_ESTIM]), N, cplx, -1, -1, -1, -1);
      /* copy values from estimation */
      tmeas[TYPE_DUR_NS][ALGO_FFTW_AUTO] = tmeas[TYPE_DUR_NS][ALGO_FFTW_ESTIM];
      tmeas[TYPE_PREP][ALGO_FFTW_AUTO] = tmeas[TYPE_PREP][ALGO_FFTW_ESTIM];
      haveAlgo[ALGO_FFTW_AUTO] = 0;
    } else {
      te = uclock_sec();
      fftwf_plan planf, planb;
      fftw_complex *in = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * N);
      fftw_complex *out = (fftw_complex*) fftwf_malloc(sizeof(fftw_complex) * N);
      memset(in, 0, sizeof(fftw_complex) * N);
      if (cplx) {
        planf = fftwf_plan_dft_1d(N, (fftwf_complex*)in, (fftwf_complex*)out, FFTW_FORWARD, flags);
        planb = fftwf_plan_dft_1d(N, (fftwf_complex*)in, (fftwf_complex*)out, FFTW_BACKWARD, flags);
      } else {
        planf = fftwf_plan_dft_r2c_1d(N, (float*)in, (fftwf_complex*)out, flags);
        planb = fftwf_plan_dft_c2r_1d(N, (fftwf_complex*)in, (float*)out, flags);
      }

      t0 = uclock_sec();  
      for (iter = 0; iter < max_iter; ++iter) {
        fftwf_execute(planf);
        fftwf_execute(planb);
      }
      t1 = uclock_sec();

      fftwf_destroy_plan(planf);
      fftwf_destroy_plan(planb);
      fftwf_free(in); fftwf_free(out);

      flops = (max_iter*2) * ((cplx ? 5 : 2.5)*N*log((double)N)/M_LN2); // see http://www.fftw.org/speed/method.html
      tmeas[TYPE_DUR_NS][ALGO_FFTW_AUTO] = show_output((flags == FFTW_MEASURE ? algoName[ALGO_FFTW_AUTO] : algoName[ALGO_FFTW_ESTIM]), N, cplx, flops, t0, t1, max_iter);
      tmeas[TYPE_PREP][ALGO_FFTW_AUTO] = (t0 - te) * 1e3;
      haveAlgo[ALGO_FFTW_AUTO] = 1;
    }
  } while (0);
#else
  (void)withFFTWmeas;
#endif  

#ifdef HAVE_GREEN_FFTS
  {
    te = uclock_sec();
    fftInit(Nlog);
    
    {
      t0 = uclock_sec();  
      for (iter = 0; iter < max_iter; ++iter) {
        if (cplx) {
          ffts(X, Nlog, 1);
          iffts(X, Nlog, 1);
        } else {
          rffts(X, Nlog, 1);
          riffts(X, Nlog, 1);
        }
      }
      t1 = uclock_sec();
      fftFree();

      flops = (max_iter*2) * ((cplx ? 5 : 2.5)*N*log((double)N)/M_LN2); // see http://www.fftw.org/speed/method.html
      tmeas[TYPE_DUR_NS][ALGO_GREEN] = show_output("Green", N, cplx, flops, t0, t1, max_iter);
      tmeas[TYPE_PREP][ALGO_GREEN] = (t0 - te) * 1e3;
      haveAlgo[ALGO_GREEN] = 1;
    }
  }
#else
  (void)Nlog;
#endif

#ifdef HAVE_KISS_FFT
  {
    kiss_fft_cfg stf;
    kiss_fft_cfg sti;
    kiss_fftr_cfg stfr;
    kiss_fftr_cfg stir;

    te = uclock_sec();
    if (cplx) {
      stf = kiss_fft_alloc(N, 0, 0, 0);
      sti = kiss_fft_alloc(N, 1, 0, 0);
    } else {
      stfr = kiss_fftr_alloc(N, 0, 0, 0);
      stir = kiss_fftr_alloc(N, 1, 0, 0);
    }

    {
      t0 = uclock_sec();  
      for (iter = 0; iter < max_iter; ++iter) {
        if (cplx) {
          kiss_fft(stf, (const kiss_fft_cpx *)X, (kiss_fft_cpx *)Y);
          kiss_fft(sti, (const kiss_fft_cpx *)Y, (kiss_fft_cpx *)X);
        } else {
          kiss_fftr(stfr, X, (kiss_fft_cpx *)Y);
          kiss_fftri(stir, (const kiss_fft_cpx *)Y, X);
        }
      }
      t1 = uclock_sec();
      kiss_fft_cleanup();

      flops = (max_iter*2) * ((cplx ? 5 : 2.5)*N*log((double)N)/M_LN2); // see http://www.fftw.org/speed/method.html
      tmeas[TYPE_DUR_NS][ALGO_KISS] = show_output("Kiss", N, cplx, flops, t0, t1, max_iter);
      tmeas[TYPE_PREP][ALGO_KISS] = (t0 - te) * 1e3;
      haveAlgo[ALGO_KISS] = 1;
    }
  }
#endif



  // PFFFT benchmark
  {
    te = uclock_sec();
    PFFFT_Setup *s = pffft_new_setup(N, cplx ? PFFFT_COMPLEX : PFFFT_REAL);
    if (s) {
      t0 = uclock_sec();  
      for (iter = 0; iter < max_iter; ++iter) {
        pffft_transform(s, X, Z, Y, PFFFT_FORWARD);
        pffft_transform(s, X, Z, Y, PFFFT_BACKWARD);
      }
      t1 = uclock_sec();
      pffft_destroy_setup(s);
    
      flops = (max_iter*2) * ((cplx ? 5 : 2.5)*N*log((double)N)/M_LN2); // see http://www.fftw.org/speed/method.html
      tmeas[TYPE_DUR_NS][ALGO_PFFFT] = show_output("PFFFT", N, cplx, flops, t0, t1, max_iter);
      tmeas[TYPE_PREP][ALGO_PFFFT] = (t0 - te) * 1e3;
      haveAlgo[ALGO_PFFFT] = 1;
    }
  }


  if (!array_output_format)
  {
    printf("prepare/ms:     ");
    for ( iter = 0; iter < NUM_FFT_ALGOS; ++iter )
    {
      if ( haveAlgo[iter] && tmeas[TYPE_DUR_NS][iter] > 0.0 ) {
        printf("%s %.3f    ", algoName[iter], tmeas[TYPE_PREP][iter] );
      }
    }
    printf("\n");
  }
  Tfastest = 0.0;
  for ( iter = 0; iter < NUM_FFT_ALGOS; ++iter )
  {
    if ( Tfastest == 0.0 || ( tmeas[TYPE_DUR_NS][iter] != 0.0 && tmeas[TYPE_DUR_NS][iter] < Tfastest ) )
      Tfastest = tmeas[TYPE_DUR_NS][iter];
  }
  if ( Tfastest > 0.0 )
  {
    if (!array_output_format)
      printf("relative fast:  ");
    for ( iter = 0; iter < NUM_FFT_ALGOS; ++iter )
    {
      if ( haveAlgo[iter] && tmeas[TYPE_DUR_NS][iter] > 0.0 ) {
        tmeas[TYPE_DUR_FASTEST][iter] = tmeas[TYPE_DUR_NS][iter] / Tfastest;
        if (!array_output_format)
          printf("%s %.3f    ", algoName[iter], tmeas[TYPE_DUR_FASTEST][iter] );
      }
    }
    if (!array_output_format)
      printf("\n");
  }

  {
    if (!array_output_format)
      printf("relative pffft: ");
    for ( iter = 0; iter < NUM_FFT_ALGOS; ++iter )
    {
      if ( haveAlgo[iter] && tmeas[TYPE_DUR_NS][iter] > 0.0 ) {
        tmeas[TYPE_REL_PFFFT][iter] = tmeas[TYPE_DUR_NS][iter] / tmeas[TYPE_DUR_NS][ALGO_PFFFT];
        if (!array_output_format)
          printf("%s %.3f    ", algoName[iter], tmeas[TYPE_REL_PFFFT][iter] );
      }
    }
    if (!array_output_format)
      printf("\n");
  }

  if (!array_output_format) {
    printf("--\n");
  }

  pffft_aligned_free(X);
  pffft_aligned_free(Y);
  pffft_aligned_free(Z);
}

#ifndef PFFFT_SIMD_DISABLE
void validate_pffft_simd(); // a small function inside pffft.c that will detect compiler bugs with respect to simd instruction 
#endif



int main(int argc, char **argv) {
  /* unfortunately, the fft size must be a multiple of 16 for complex FFTs 
     and 32 for real FFTs -- a lot of stuff would need to be rewritten to
     handle other cases (or maybe just switch to a scalar fft, I don't know..) */

  /* int Nvalues[] = { 64, 96, 128, 160, 192, 256, 384, 5*96, 512, 5*128, 3*256, 800, 1024, 2048, 2400, 4096, 8192, 9*1024, 16384, 32768, 256*1024, 1024*1024, -1 };
   */
  int Nexp[] =    { 1, 2, 3,  4,  5,  6,   7,   8,   9,   10,
	               11, 12, 13, 14, 15, 16, 17, 18, 19,    20, -1 };
  int Nvalues[] = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,
                  2*1024, 4*1024, 8*1024, 16*1024, 32*1024, 64*1024, 128*1024, 256*1024, 512*1024,
                  1024*1024, -1 };
  int benchReal=1, benchCplx=1, withFFTWmeas=0;
  int realCplxIdx, typeIdx;
  int i, k;
  int smallestCplxN = pffft_simd_size()*pffft_simd_size();
  int smallestRealN = 2*smallestCplxN;

  double tmeas[2][sizeof(Nvalues)/sizeof(Nvalues[0])][NUM_TYPES][NUM_FFT_ALGOS];
  int haveAlgo[NUM_FFT_ALGOS];

  for ( i = 0; i < NUM_FFT_ALGOS; ++i )
    haveAlgo[i] = 0;

  for ( i = 1; i < argc; ++i ) {
    if (!strcmp(argv[i], "--array-format") || !strcmp(argv[i], "--table")) {
      array_output_format = 1;
    }
    if (!strcmp(argv[i], "--no-tab")) {
      array_output_format = 0;
    }
    if (!strcmp(argv[i], "--real")) {
      benchCplx = 0;
    }
    if (!strcmp(argv[i], "--cplx")) {
      benchReal = 0;
    }
    if (!strcmp(argv[i], "--fftw-measure")) {
      withFFTWmeas = 1;
    }
  }

#ifdef HAVE_FFTW
  if (withFFTWmeas)
  {
    compiledInAlgo[ALGO_FFTW_AUTO] = 1;
    algoName[ALGO_FFTW_AUTO] = "FFTW(meas.)"; /* "FFTW (auto)" */
    algoTableHeader[NUM_FFT_ALGOS][0] = "|real FFTWmeas "; /* "|real FFTWauto " */
    algoTableHeader[NUM_FFT_ALGOS][0] = "|cplx FFTWmeas "; /* "|cplx FFTWauto " */
  }
  /* else
    compiledInAlgo[ALGO_FFTW_AUTO] = 0; */
#endif

#ifdef PFFFT_SIMD_DISABLE
  algoName[ALGO_PFFFT] = "PFFFT(scal)";
#else
  validate_pffft_simd();
#endif
  pffft_validate(1);
  pffft_validate(0);
  if (!array_output_format) {
    if (benchReal) {
      for (i=0; Nvalues[i] > 0; ++i) {
        if ( Nvalues[i] >= smallestRealN )
          benchmark_ffts(Nexp[i], Nvalues[i], 0 /* real fft */, withFFTWmeas, tmeas[0][i], haveAlgo);
      }
    }
    if (benchCplx) {
      for (i=0; Nvalues[i] > 0; ++i) {
        if ( Nvalues[i] >= smallestCplxN )
          benchmark_ffts(Nexp[i], Nvalues[i], 1 /* cplx fft */, withFFTWmeas, tmeas[1][i], haveAlgo);
      }
    }
  } else {

    /* print table headers */
    {
      printf("| input len ");
      for (realCplxIdx = 0; realCplxIdx < 2; ++realCplxIdx)
      {
        if ( (realCplxIdx == 0 && !benchReal) || (realCplxIdx == 1 && !benchCplx) )
          continue;
        for (k=0; k < NUM_FFT_ALGOS; ++k)
        {
          if ( compiledInAlgo[k] )
            printf("%s", algoTableHeader[k][realCplxIdx]);
        }
      }
      printf("|\n");
    }
    /* print table value seperators */
    {
      printf("|----------");
      for (realCplxIdx = 0; realCplxIdx < 2; ++realCplxIdx)
      {
        if ( (realCplxIdx == 0 && !benchReal) || (realCplxIdx == 1 && !benchCplx) )
          continue;
        for (k=0; k < NUM_FFT_ALGOS; ++k)
        {
          if ( compiledInAlgo[k] )
            printf("%s", ":|-------------");
        }
      }
      printf(":|\n");
    }

    for (i=0; Nvalues[i] > 0; ++i) {
      if ( Nvalues[i] >= smallestRealN && Nvalues[i] >= smallestCplxN )
      {
        printf("|%9d  ", Nvalues[i]);
        if (benchReal) {
          benchmark_ffts(Nexp[i], Nvalues[i], 0, withFFTWmeas, tmeas[0][i], haveAlgo);
        }
        if (benchCplx) {
          benchmark_ffts(Nexp[i], Nvalues[i], 1, withFFTWmeas, tmeas[1][i], haveAlgo);
        }
        printf("|\n");
      }
    }
    printf(" (numbers are given in MFlops)\n");
  }


  printf("\n\n");
  printf("smallest cplx fft size: %d\n", smallestCplxN);
  printf("smallest real fft size: %d\n", smallestRealN);


  for (realCplxIdx = 0; realCplxIdx < 2; ++realCplxIdx)
  {
    if ( (benchReal && realCplxIdx == 0) || (benchCplx && realCplxIdx == 1) )
    {
      for (typeIdx = 0; typeIdx < NUM_TYPES; ++typeIdx)
      {
        printf("\n\n\n");
        printf("%s %s\n", (realCplxIdx == 0 ? "real" : "complex"), typeText[typeIdx]);
        {
          printf("size, log2, ");
          for (k=0; k < NUM_FFT_ALGOS; ++k)
            if ( haveAlgo[k] )
              printf("%s, ", algoName[k]);
          printf("\n");
        }
        for (i=0; Nvalues[i] > 0; ++i)
        {
          if ( (benchReal && Nvalues[i] >= smallestRealN) || (benchCplx && Nvalues[i] >= smallestCplxN) )
          {
            printf("%d, %.3f, ", Nvalues[i], log10((double)Nvalues[i])/log10(2.0) );
            for (k=0; k < NUM_FFT_ALGOS; ++k)
              if ( haveAlgo[k] )
                printf("%f, ", tmeas[realCplxIdx][i][typeIdx][k]);
            printf("\n");
          }
        }
      }
    }
  }

  return 0;
}

