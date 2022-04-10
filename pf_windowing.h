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

#pragma once

#include "pf_cplx.h"

#ifdef __cplusplus
extern "C" {
#endif

struct pf_windowing_param_tag;

/* https://en.wikipedia.org/wiki/Window_function */

typedef enum
{
    pf_winRect =0
  , pf_winHann
  , pf_winHamming
  , pf_winBlackman
  , pf_winBlackmanExact
  , pf_winNutall
  , pf_winBlackmanHarris
  , pf_winBlackmanNutall
  , pf_winFlatTop_MATLAB /* https://en.wikipedia.org/wiki/Window_function
                          * https://www.mathworks.com/help/signal/ref/flattopwin.html
                          */
  , pf_winFlatTop_SR785  /* Flat Top of spectrum-analyzer SR785 of Stanford Research Systems
                          * see https://de.wikipedia.org/wiki/Fensterfunktion
                          */
  , pf_winFlatTop_HFT70  /* "HFT70" from D.3.1 page 45 of
                          * https://pure.mpg.de/rest/items/item_152164_1/component/file_152163/content
                          * first zero at f = +/- 4.00 bins
                          * highest sidelobe: -70.4 dB at f = +/- 4.65 bins
                          * sidelobe drop at rate f^-1
                          * optimal overlap: 72.2%
                          * amplitude flatness: 0.964
                          * power flatness: 0.637
                          * overlap correlation: 0.041
                          */
  , pf_winFlatTop_HFT95  /* "HFT95" from D.3.2 page 46 of
                          * https://pure.mpg.de/rest/items/item_152164_1/component/file_152163/content
                          * first zero at f = +/- 5.00 bins
                          * highest sidelobe: -95.0 dB at f = +/- 7.49 bins
                          * sidelobe drop at rate f^-1
                          * optimal overlap: 75.6%
                          * amplitude flatness: 0.952
                          * power flatness: 0.647
                          * overlap correlation: 0.056
                          */
  , pf_winFlatTop_HFT90D /* "HFT90D" from D.3.3 page 46 of
                          * https://pure.mpg.de/rest/items/item_152164_1/component/file_152163/content
                          * first zero at f = +/- 5.00 bins
                          * highest sidelobe: -90.2 dB at f = +/- 5.58 bins
                          * sidelobe drop at rate f^-3
                          * optimal overlap: 76.0%
                          * amplitude flatness: 0.953
                          * power flatness: 0.646
                          * overlap correlation: 0.054
                          */
  , pf_num_windows
} pf_windowing_T;


const char * pf_window_text(pf_windowing_T win);

void pf_window_free(const struct pf_windowing_param_tag * param);

/*
 * scale: for additional scaling, e.g. normalization of window or FFT scale compensation
 *
 * w0: != 0 => top at center of coefficients; towards 0 at corners
 *     == 0 => top at corners; towards 0 at center
 *
 * alpha: for blackman window; negative value uses default of 0.16
 *
 */

#ifdef __cplusplus
// add default values for some parameters

// _trig functions utilize the (slow) trigonometric functions
void   pf_window_trig(      float * pafWin, pf_windowing_T win, int N, float scale =1.0F, int w0 =1, float alpha =-1.0F);
void   pf_window_trig_cx(complexf * pacWin, pf_windowing_T win, int N, float scale =1.0F, int w0 =1, float alpha =-1.0F);

const struct pf_windowing_param_tag * pf_window_alloc(pf_windowing_T win, int N, int w0 = 1, float alpha =-1.0F);
void   pf_get_window_param(const struct pf_windowing_param_tag * param, int *N_div =0, int *num_alpha_coeff =0);

double pf_window_sum(const struct pf_windowing_param_tag * param);
void   pf_window(    float * pafWin,    const struct pf_windowing_param_tag * param, float scale =1.0F);
void   pf_window_cx( complexf * pacWin, const struct pf_windowing_param_tag * param, float scale =1.0F);

#else
/* no overloads and no default parameters with pure C */

/* _trig functions utilize the (slow) trigonometric functions */
void   pf_window_trig(      float * pafWin, pf_windowing_T win, int N, float scale, int w0, float alpha);
void   pf_window_trig_cx(complexf * pacWin, pf_windowing_T win, int N, float scale, int w0, float alpha);

const struct pf_windowing_param_tag * pf_window_alloc(pf_windowing_T win, int N, int w0, float alpha);
void   pf_get_window_param(const struct pf_windowing_param_tag * param, int *N_div, int *num_alpha_coeff);

double pf_window_sum(const struct pf_windowing_param_tag * param);
void   pf_window(    float * pafWin,    const struct pf_windowing_param_tag * param, float scale);
void   pf_window_cx( complexf * pacWin, const struct pf_windowing_param_tag * param, float scale);
#endif



#ifdef __cplusplus
}


// additional C++ inline function pf_win_trig() with overloads:

inline void pf_win_trig(float * pafWin, pf_windowing_T win, int N, float scale =1.0F, int w0 =1, float alpha =-1.0F) {
  pf_window_trig(pafWin, win, N, scale, w0, alpha);
}

inline void pf_win_trig(complexf * pacWin, pf_windowing_T win, int N, float scale =1.0F, int w0 =1, float alpha =-1.0F) {
  pf_window_trig_cx(pacWin, win, N, scale, w0, alpha);
}


inline void pf_win(float * pafWin, const struct pf_windowing_param_tag * param, float scale = 1.0F ) {
  pf_window(pafWin, param, scale);
}

inline void pf_win(complexf * pacWin, const struct pf_windowing_param_tag * param, float scale = 1.0F ) {
  pf_window_cx(pacWin, param, scale);
}

#endif

