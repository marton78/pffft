/*
This software is part of pffft/pfdsp, a set of simple DSP routines.

Copyright (c) 2021  Hayati Ayguen <h_ayguen@web.de>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ANDRAS RETZLER BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "pf_windowing.h"
#include "fmv.h"

#include <math.h>
#include <assert.h>
#include <stdio.h>


typedef float freal;
typedef double frealsum;

typedef freal freal2[2];

#define rcos(X)  cosf(X)
#define rsin(X)  sinf(X)
#define rsqrt(X) sqrtf(X)

#define TRIG_TAB_SIZE 4

//they dropped M_PI in C99, so we define it:
#define PI      ((freal)3.14159265358979323846)
#define TWO_PI  ((freal)(2.0 * 3.14159265358979323846))

#define c0          0.0F
#define c1          1.0F
#define c2          2.0F
#define c3          3.0F
#define c4          4.0F

#define c_0p5       0.5F


typedef struct pf_windowing_param_tag
{
  freal trig_tab[TRIG_TAB_SIZE +1][2];
  freal a[5];  // 0 .. 4
  freal N_div;
  int N_div_i;  // N_div as integer
  int num_a_coeff;
  int N;
  int w0;
} pf_windowing_param_t;


const char * pf_window_text(pf_windowing_T win)
{
    switch (win)
    {
    case pf_winRect:            return "Rectangle / Boxcar";
    case pf_winHann:            return "von Hann";
    case pf_winHamming:         return "Hamming";
    case pf_winBlackman:        return "Blackman";
    case pf_winBlackmanExact:   return "Blackman (Exact)";
    case pf_winNutall:          return "Nutall";
    case pf_winBlackmanHarris:  return "Blackman-Harris";
    case pf_winBlackmanNutall:  return "Blackman-Nutall";
    case pf_winFlatTop_MATLAB:  return "FlatTop (Matlab)";
    case pf_winFlatTop_SR785:   return "FlatTop (SR785)";
    case pf_winFlatTop_HFT70:   return "FlatTop (HFT70)";
    case pf_winFlatTop_HFT95:   return "FlatTop (HFT95)";
    case pf_winFlatTop_HFT90D:  return "FlatTop (HFT90D)";
    case pf_num_windows:
    default:  ;
    }
    return 0;
}


static
pf_windowing_param_t pf_window_param(int calc_cordic, pf_windowing_T win, int N, int w0, freal alpha)
{
  pf_windowing_param_t wp;
  //pf_windowing_param_t *p = &wp;
  // classical alpha for 3-term Blackman: 0.16
  const freal bl_alpha = (alpha >= 0) ? alpha : (freal)(0.16);
  wp.a[0] = (freal)(1.0);
  wp.a[1] = (freal)(0.0);
  wp.a[2] = (freal)(0.0);
  wp.a[3] = (freal)(0.0);
  wp.a[4] = (freal)(0.0);
  wp.N_div = (freal)(N -1);
  wp.N_div_i = N - 1;
  wp.num_a_coeff = 1;
  wp.N = N;
  wp.w0 = w0;

  switch( win )
  {
    case pf_winRect:
      break;
    default:    assert(0);
    case pf_winHann:
      wp.a[0] = (freal)(0.5);
      wp.a[1] = (freal)(-0.5);
      wp.num_a_coeff = 2;
      break;
    case pf_winHamming:
      wp.a[0] = (freal)( 25.0 / 46.0); // ~= 0.54
      wp.a[1] = (freal)(-21.0 / 46.0); // ~= 0.46
      wp.num_a_coeff = 2;
      break;
    case pf_winBlackman:
      wp.a[0] = c_0p5 * (c1 - bl_alpha);  // default: 0.42
      wp.a[1] = (freal)(-0.5);
      wp.a[2] = c_0p5 * bl_alpha;  // default: 0.08
      wp.num_a_coeff = 3;
      break;
    case pf_winBlackmanExact:
      wp.a[0] = (freal)( 7938.0/18608.0);  // = 0.42659
      wp.a[1] = (freal)(-9240.0/18608.0);  // = 0.49656
      wp.a[2] = (freal)( 1430.0/18608.0);  // = 0.076849
      wp.num_a_coeff = 3;
      break;
    case pf_winNutall:
      wp.a[0] = (freal)( 0.355768);
      wp.a[1] = (freal)(-0.487396);
      wp.a[2] = (freal)( 0.144232);
      wp.a[3] = (freal)(-0.012604);
      wp.num_a_coeff = 4;
      break;
    case pf_winBlackmanHarris:
      wp.a[0] = (freal)( 0.35875);
      wp.a[1] = (freal)(-0.48829);
      wp.a[2] = (freal)( 0.14128);
      wp.a[3] = (freal)(-0.01168);
      wp.num_a_coeff = 4;
      break;
    case pf_winBlackmanNutall:
      wp.a[0] = (freal)( 0.3635819);
      wp.a[1] = (freal)(-0.4891775);
      wp.a[2] = (freal)( 0.1365995);
      wp.a[3] = (freal)(-0.0106411);
      wp.num_a_coeff = 4;
      break;
    case pf_winFlatTop_MATLAB:
      wp.a[0] = (freal)( 0.21557895);
      wp.a[1] = (freal)(-0.41663158);
      wp.a[2] = (freal)( 0.277263158);
      wp.a[3] = (freal)(-0.083578947);
      wp.a[4] = (freal)( 0.006947368);
      //wp.N_div = (freal)(N -1);
      //wp.N_div_i = N -1;
      wp.num_a_coeff = 5;
      break;
    case pf_winFlatTop_SR785:
      wp.a[0] = (freal)( 1.0);
      wp.a[1] = (freal)(-1.93);
      wp.a[2] = (freal)( 1.29);
      wp.a[3] = (freal)(-0.388);
      wp.a[4] = (freal)( 0.028);
      wp.num_a_coeff = 5;
      break;
    case pf_winFlatTop_HFT70:
      wp.a[0] = (freal)( 1.0);
      wp.a[1] = (freal)(-1.90796);
      wp.a[2] = (freal)( 1.07349);
      wp.a[3] = (freal)(-0.18199);
      wp.num_a_coeff = 4;
      break;
    case pf_winFlatTop_HFT95:
      wp.a[0] = (freal)( 1.0);
      wp.a[1] = (freal)(-1.9383379);
      wp.a[2] = (freal)( 1.3045202);
      wp.a[3] = (freal)(-0.4028270);
      wp.a[4] = (freal)( 0.0350665);
      wp.num_a_coeff = 5;
      break;
    case pf_winFlatTop_HFT90D:
      wp.a[0] = (freal)( 1.0);
      wp.a[1] = (freal)(-1.942604);
      wp.a[2] = (freal)( 1.340318);
      wp.a[3] = (freal)(-0.440811);
      wp.a[4] = (freal)( 0.043097);
      wp.num_a_coeff = 5;
      break;
  }


  if (!calc_cordic)
      return wp;
  // pre-calculate small trig table - hope this allow SIMD parallelization
  freal2 *trig_tab = wp.trig_tab;
  const freal phi_inc = TWO_PI / wp.N_div;
  const freal c_inc_cos = rcos(phi_inc);
  const freal c_inc_sin = rsin(phi_inc);
  freal cos_phi = trig_tab[0][0] = c1;
  freal sin_phi = trig_tab[0][1] = c0;
#if 0
  printf("pf_window_param(N = %d, %s)\n", N, pf_window_text(win));
#endif
  for ( int k = 1; k <= TRIG_TAB_SIZE; ++k )
  {
    // rotate
    const freal cos_sav = cos_phi;
    cos_phi = c_inc_cos * cos_sav - c_inc_sin * sin_phi;
    sin_phi = c_inc_sin * cos_sav + c_inc_cos * sin_phi;
    // re-normalize cos_phi/sin_phi to unit-circle
    const freal mag = rsqrt(cos_phi * cos_phi + sin_phi * sin_phi);
    cos_phi /= mag;
    sin_phi /= mag;
    trig_tab[k][0] = cos_phi;
    trig_tab[k][1] = sin_phi;
#if 0
    printf("pf_window_param() trig_tab[%d] = %.3g + i * %.3g. %.3g + i %.3g. dphi = %g\n", k
        , trig_tab[k][0], trig_tab[k][1]
        , rcos(k*phi_inc), rsin(k*phi_inc)
        , phi_inc
        );
#endif
  }

  return wp;
}


const struct pf_windowing_param_tag * pf_window_alloc(pf_windowing_T win, int N, int w0, freal alpha)
{
  struct pf_windowing_param_tag * param = new pf_windowing_param_t();
  *param = pf_window_param(1, win, N, w0, alpha);
  return param;
}


void pf_get_window_param(const struct pf_windowing_param_tag * param, int *N_div, int *num_alpha_coeff)
{
  if (N_div)
      *N_div = param->N_div_i;
  if (num_alpha_coeff)
      *num_alpha_coeff = param->num_a_coeff;
}


void pf_window_free(const struct pf_windowing_param_tag * param)
{
  delete param;
}



PF_TARGET_CLONES
void pf_window_trig(freal * pafWin, pf_windowing_T win, int N, freal scale, int w0, freal alpha)
{
  const freal phi0 = w0 ? c0 : PI;
  const pf_windowing_param_t wp = pf_window_param(0, win, N, w0, alpha);

  switch( wp.num_a_coeff )
  {
    case 1:
      {
        if (scale == (freal)(1.0))
          break;
        const freal a0 = scale;
        for ( int n = 0; n < N; ++n )
          pafWin[n] *= a0;
      }
      break;
    default:    assert(0);
    case 2:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi);
          pafWin[n] *= scale * wn;
        }
      }
      break;
    case 3:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi)
                              + a2 * rcos(c2 * phi);
          pafWin[n] *= scale * wn;
        }
      }
      break;
    case 4:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi)
                              + a2 * rcos(c2 * phi)
                              + a3 * rcos(c3 * phi);
          pafWin[n] *= scale * wn;
        }
      }
      break;
    case 5:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        const freal a4 = wp.a[4];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi)
                              + a2 * rcos(c2 * phi)
                              + a3 * rcos(c3 * phi)
                              + a4 * rcos(c4 * phi);
          pafWin[n] *= scale * wn;
        }
      }
      break;
  }
}


PF_TARGET_CLONES
void pf_window_trig_cx(complexf * pacWin, pf_windowing_T win, int N, freal scale, int w0, freal alpha)
{
  const freal phi0 = w0 ? c0 : PI;
  const pf_windowing_param_t wp = pf_window_param(0, win, N, w0, alpha);

  switch( wp.num_a_coeff )
  {
    case 1:
      {
        if (scale == (freal)(1.0))
          break;
        const freal a0 = scale;
        for ( int n = 0; n < N; ++n )
        {
          pacWin[n].i *= a0;
          pacWin[n].q *= a0;
        }
      }
      break;
    default:    assert(0);
    case 2:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi);
          pacWin[n].i *= scale * wn;
          pacWin[n].q *= scale * wn;
        }
      }
      break;
    case 3:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi)
                              + a2 * rcos(c2 * phi);
          pacWin[n].i *= scale * wn;
          pacWin[n].q *= scale * wn;
        }
      }
      break;
    case 4:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi)
                              + a2 * rcos(c2 * phi)
                              + a3 * rcos(c3 * phi);
          pacWin[n].i *= scale * wn;
          pacWin[n].q *= scale * wn;
        }
      }
      break;
    case 5:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        const freal a4 = wp.a[4];
        for ( int n = 0; n < N; ++n )
        {
          const freal phi = phi0 + ( TWO_PI * n ) / wp.N_div;
          const freal wn = a0 + a1 * rcos(phi)
                              + a2 * rcos(c2 * phi)
                              + a3 * rcos(c3 * phi)
                              + a4 * rcos(c4 * phi);
          pacWin[n].i *= scale * wn;
          pacWin[n].q *= scale * wn;
        }
      }
      break;
  }
}


frealsum pf_window_sum(const struct pf_windowing_param_tag * param)
{
  const pf_windowing_param_t &wp = *param;
  const int N = wp.N;
  const freal2 *trig_tab = wp.trig_tab;
  // w0 != 0: start with phi = 0
  // w0 == 0: start with phi = ( TWO_PI * (-N/2) ) / N == -PI
  freal cos_phi0 = wp.w0 ? c1 : -c1;
  freal sin_phi0 = c0;
  frealsum wsum = (frealsum)(0.0);

  switch( wp.num_a_coeff )
  {
    case 1:
      wsum += N;
      break;
    default:    assert(0);
    case 2:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            //const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal wn = a0 + a1 * cos_phi;
            wsum += wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 3:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi;
            wsum += wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 4:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi + a3 * cos_3phi;
            wsum += wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 5:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        const freal a4 = wp.a[4];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal cos_4phi = cos_2phi * cos_2phi - sin_2phi * sin_2phi;
            const freal wn = a0 +a1 * cos_phi  +a2 * cos_2phi
                                +a3 * cos_3phi +a4 * cos_4phi;
            wsum += wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
  }
  return wsum;
}


PF_TARGET_CLONES
void pf_window_no_scale(freal * pafWin, const struct pf_windowing_param_tag * param)
{
  const pf_windowing_param_t &wp = *param;
  const int N = wp.N;
  const freal2 *trig_tab = wp.trig_tab;
  // w0 != 0: start with phi = 0
  // w0 == 0: start with phi = ( TWO_PI * (-N/2) ) / N == -PI
  freal cos_phi0 = wp.w0 ? c1 : -c1;
  freal sin_phi0 = c0;

  switch( wp.num_a_coeff )
  {
    case 1:
      break;
    default:    assert(0);
    case 2:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            //const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal wn = a0 + a1 * cos_phi;
            pafWin[n+k] *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 3:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi;
            pafWin[n+k] *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 4:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi + a3 * cos_3phi;
            pafWin[n+k] *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 5:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        const freal a4 = wp.a[4];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal cos_4phi = cos_2phi * cos_2phi - sin_2phi * sin_2phi;
            const freal wn = a0 +a1 * cos_phi  +a2 * cos_2phi
                                +a3 * cos_3phi +a4 * cos_4phi;
            pafWin[n+k] *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
  }
}


PF_TARGET_CLONES
void pf_window(freal * pafWin, const struct pf_windowing_param_tag * param, freal scale)
{
  if (scale == c1)
  {
    pf_window_no_scale(pafWin, param);
    return;
  }

  const pf_windowing_param_t &wp = *param;
  const int N = wp.N;
  const freal2 *trig_tab = wp.trig_tab;
  // w0 != 0: start with phi = 0
  // w0 == 0: start with phi = ( TWO_PI * (-N/2) ) / N == -PI
  freal cos_phi0 = wp.w0 ? c1 : -c1;
  freal sin_phi0 = c0;

  switch( wp.num_a_coeff )
  {
    case 1:
      {
        if (scale == (freal)(1.0))
          break;
        const freal a0 = scale;
        for ( int n =0; n < N; ++n )
          pafWin[n] *= a0;
      }
      break;
    default:    assert(0);
    case 2:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            //const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal wn = a0 + a1 * cos_phi;
            pafWin[n+k] *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 3:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi;
            pafWin[n+k] *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 4:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi + a3 * cos_3phi;
            pafWin[n+k] *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 5:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        const freal a4 = wp.a[4];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal cos_4phi = cos_2phi * cos_2phi - sin_2phi * sin_2phi;
            const freal wn = a0 +a1 * cos_phi  +a2 * cos_2phi
                                +a3 * cos_3phi +a4 * cos_4phi;
            pafWin[n+k] *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
  }
}


PF_TARGET_CLONES
void pf_window_cx_no_scale(complexf * pacWin, const struct pf_windowing_param_tag * param)
{
  const pf_windowing_param_t &wp = *param;
  const int N = wp.N;
  const freal2 *trig_tab = wp.trig_tab;
  // w0 != 0: start with phi = 0
  // w0 == 0: start with phi = ( TWO_PI * (-N/2) ) / N == -PI
  freal cos_phi0 = wp.w0 ? c1 : -c1;
  freal sin_phi0 = c0;

  switch( wp.num_a_coeff )
  {
    case 1:
      break;
    default:    assert(0);
    case 2:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            //const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal wn = a0 + a1 * cos_phi;
            pacWin[n+k].i *= wn;
            pacWin[n+k].q *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 3:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi;
            pacWin[n+k].i *= wn;
            pacWin[n+k].q *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 4:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi + a3 * cos_3phi;
            pacWin[n+k].i *= wn;
            pacWin[n+k].q *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 5:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        const freal a4 = wp.a[4];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal cos_4phi = cos_2phi * cos_2phi - sin_2phi * sin_2phi;
            const freal wn = a0 +a1 * cos_phi  +a2 * cos_2phi
                                +a3 * cos_3phi +a4 * cos_4phi;
            pacWin[n+k].i *= wn;
            pacWin[n+k].q *= wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
  }
}


PF_TARGET_CLONES
void pf_window_cx(complexf * pacWin, const struct pf_windowing_param_tag * param, freal scale)
{
  if (scale == c1)
  {
    pf_window_cx_no_scale(pacWin, param);
    return;
  }

  const pf_windowing_param_t &wp = *param;
  const int N = wp.N;
  const freal2 *trig_tab = wp.trig_tab;
  // w0 != 0: start with phi = 0
  // w0 == 0: start with phi = ( TWO_PI * (-N/2) ) / N == -PI
  freal cos_phi0 = wp.w0 ? c1 : -c1;
  freal sin_phi0 = c0;

  switch( wp.num_a_coeff )
  {
    case 1:
      {
        if (scale == (freal)(1.0))
          break;
        const freal a0 = scale;
        for ( int n =0; n < N; ++n )
        {
          pacWin[n].i *= a0;
          pacWin[n].q *= a0;
        }
      }
      break;
    default:    assert(0);
    case 2:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            //const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal wn = a0 + a1 * cos_phi;
            pacWin[n+k].i *= scale * wn;
            pacWin[n+k].q *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 3:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi;
            pacWin[n+k].i *= scale * wn;
            pacWin[n+k].q *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 4:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal wn = a0 + a1 * cos_phi + a2 * cos_2phi + a3 * cos_3phi;
            pacWin[n+k].i *= scale * wn;
            pacWin[n+k].q *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
    case 5:
      {
        const freal a0 = wp.a[0];
        const freal a1 = wp.a[1];
        const freal a2 = wp.a[2];
        const freal a3 = wp.a[3];
        const freal a4 = wp.a[4];
        for ( int n =0; n < N; )
        {
          const int M = (N - n >= TRIG_TAB_SIZE) ? TRIG_TAB_SIZE : (N - n);
          for ( int k = 0; k < M; ++k )
          {
            const freal cos_phi = trig_tab[k][0] * cos_phi0 - trig_tab[k][1] * sin_phi0;
            const freal sin_phi = trig_tab[k][1] * cos_phi0 + trig_tab[k][0] * sin_phi0;
            const freal cos_2phi = cos_phi * cos_phi - sin_phi * sin_phi;
            const freal sin_2phi = c2 * sin_phi * cos_phi;
            const freal cos_3phi = cos_2phi * cos_phi - sin_2phi * sin_phi;
            const freal cos_4phi = cos_2phi * cos_2phi - sin_2phi * sin_2phi;
            const freal wn = a0 +a1 * cos_phi  +a2 * cos_2phi
                                +a3 * cos_3phi +a4 * cos_4phi;
            pacWin[n+k].i *= scale * wn;
            pacWin[n+k].q *= scale * wn;
          }
          n += TRIG_TAB_SIZE;
          if ( n >= N )
            break;
          // rotate
          const freal cos_sav = cos_phi0;
          cos_phi0 = trig_tab[TRIG_TAB_SIZE][0] * cos_sav - trig_tab[TRIG_TAB_SIZE][1] * sin_phi0;
          sin_phi0 = trig_tab[TRIG_TAB_SIZE][1] * cos_sav + trig_tab[TRIG_TAB_SIZE][0] * sin_phi0;
          // re-normalize cos_phi/sin_phi to unit-circle
          const freal mag = rsqrt(cos_phi0 * cos_phi0 + sin_phi0 * sin_phi0);
          cos_phi0 /= mag;
          sin_phi0 /= mag;
        }
      }
      break;
  }
}

