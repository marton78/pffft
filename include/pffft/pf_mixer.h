/*
This software is part of pffft/pfdsp, a set of simple DSP routines.

Copyright (c) 2014, Andras Retzler <randras@sdr.hu>
Copyright (c) 2020  Hayati Ayguen <h_ayguen@web.de>
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

#include <stdio.h>
#include <stdint.h>

#include "pf_cplx.h"

#ifndef PFDSP_EXPORT
  #ifdef PFDSP_STATIC_DEFINE
    #define PFDSP_EXPORT
  #elif defined(_WIN32) || defined(__CYGWIN__)
    #ifdef PFDSP_EXPORTS
      #define PFDSP_EXPORT __declspec(dllexport)
    #else
      #define PFDSP_EXPORT __declspec(dllimport)
    #endif
  #elif defined(PFDSP_EXPORTS) && (defined(__GNUC__) || defined(__clang__))
    #define PFDSP_EXPORT __attribute__((visibility("default")))
  #else
    #define PFDSP_EXPORT
  #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif


// =================================================================================

PFDSP_EXPORT int have_sse_shift_mixer_impl();


/*********************************************************************/

/**************/
/*** ALGO A ***/
/**************/

PFDSP_EXPORT float shift_math_cc(const complexf *input, complexf* output, int input_size, float rate, float starting_phase);


/*********************************************************************/

/**************/
/*** ALGO B ***/
/**************/

typedef struct shift_table_data_s
{
    float* table;
    int table_size;
} shift_table_data_t;

PFDSP_EXPORT void shift_table_deinit(shift_table_data_t table_data);
PFDSP_EXPORT shift_table_data_t shift_table_init(int table_size);
PFDSP_EXPORT float shift_table_cc(complexf* input, complexf* output, int input_size, float rate, shift_table_data_t table_data, float starting_phase);

/*********************************************************************/

/**************/
/*** ALGO C ***/
/**************/

typedef struct shift_addfast_data_s
{
    float dsin[4];
    float dcos[4];
    float phase_increment;
} shift_addfast_data_t;

PFDSP_EXPORT shift_addfast_data_t shift_addfast_init(float rate);
PFDSP_EXPORT float shift_addfast_cc(complexf *input, complexf* output, int input_size, shift_addfast_data_t* d, float starting_phase);
PFDSP_EXPORT float shift_addfast_inp_c(complexf *in_out, int N_cplx, shift_addfast_data_t* d, float starting_phase);


/*********************************************************************/

/**************/
/*** ALGO D ***/
/**************/

typedef struct shift_unroll_data_s
{
    float* dsin;
    float* dcos;
    float phase_increment;
    int size;
} shift_unroll_data_t;

PFDSP_EXPORT shift_unroll_data_t shift_unroll_init(float rate, int size);
PFDSP_EXPORT void shift_unroll_deinit(shift_unroll_data_t* d);
PFDSP_EXPORT float shift_unroll_cc(complexf *input, complexf* output, int size, shift_unroll_data_t* d, float starting_phase);
PFDSP_EXPORT float shift_unroll_inp_c(complexf* in_out, int size, shift_unroll_data_t* d, float starting_phase);


/*********************************************************************/

/**************/
/*** ALGO E ***/
/**************/

/* similar to shift_unroll_cc() - but, have fixed and limited precalc size
 * idea: smaller cache usage by table
 * size must be multiple of CSDR_SHIFT_LIMITED_SIMD (= 4)
 */
#define PF_SHIFT_LIMITED_UNROLL_SIZE  128
#define PF_SHIFT_LIMITED_SIMD_SZ  4

typedef struct shift_limited_unroll_data_s
{
    float dcos[PF_SHIFT_LIMITED_UNROLL_SIZE];
    float dsin[PF_SHIFT_LIMITED_UNROLL_SIZE];
    complexf complex_phase;
    float phase_increment;
} shift_limited_unroll_data_t;

PFDSP_EXPORT shift_limited_unroll_data_t shift_limited_unroll_init(float rate);
/* size must be multiple of PF_SHIFT_LIMITED_SIMD_SZ */
/* starting_phase for next call is kept internal in state */
PFDSP_EXPORT void shift_limited_unroll_cc(const complexf *input, complexf* output, int size, shift_limited_unroll_data_t* d);
PFDSP_EXPORT void shift_limited_unroll_inp_c(complexf* in_out, int size, shift_limited_unroll_data_t* d);


/*********************************************************************/

/**************/
/*** ALGO F ***/
/**************/

typedef struct shift_limited_unroll_A_sse_data_s
{
    /* small/limited trig table */
    float dcos[PF_SHIFT_LIMITED_UNROLL_SIZE+PF_SHIFT_LIMITED_SIMD_SZ];
    float dsin[PF_SHIFT_LIMITED_UNROLL_SIZE+PF_SHIFT_LIMITED_SIMD_SZ];
    /* 4 times complex phase */
    float phase_state_i[PF_SHIFT_LIMITED_SIMD_SZ];
    float phase_state_q[PF_SHIFT_LIMITED_SIMD_SZ];
    /* N_cplx_per_block times increment - for future parallel variants */
    float dcos_blk;
    float dsin_blk;
    /* */
    float phase_increment;
} shift_limited_unroll_A_sse_data_t;

PFDSP_EXPORT shift_limited_unroll_A_sse_data_t shift_limited_unroll_A_sse_init(float relative_freq, float phase_start_rad);
PFDSP_EXPORT void shift_limited_unroll_A_sse_inp_c(complexf* in_out, int N_cplx, shift_limited_unroll_A_sse_data_t* d);


/*********************************************************************/

/**************/
/*** ALGO G ***/
/**************/

typedef struct shift_limited_unroll_B_sse_data_s
{
    /* small/limited trig table */
    float dtrig[PF_SHIFT_LIMITED_UNROLL_SIZE+PF_SHIFT_LIMITED_SIMD_SZ];
    /* 4 times complex phase */
    float phase_state_i[PF_SHIFT_LIMITED_SIMD_SZ];
    float phase_state_q[PF_SHIFT_LIMITED_SIMD_SZ];
    /* N_cplx_per_block times increment - for future parallel variants */
    float dcos_blk;
    float dsin_blk;
    /* */
    float phase_increment;
} shift_limited_unroll_B_sse_data_t;

PFDSP_EXPORT shift_limited_unroll_B_sse_data_t shift_limited_unroll_B_sse_init(float relative_freq, float phase_start_rad);
PFDSP_EXPORT void shift_limited_unroll_B_sse_inp_c(complexf* in_out, int N_cplx, shift_limited_unroll_B_sse_data_t* d);

/*********************************************************************/

/**************/
/*** ALGO H ***/
/**************/

typedef struct shift_limited_unroll_C_sse_data_s
{
    /* small/limited trig table - interleaved: 4 cos, 4 sin, 4 cos, .. */
    float dinterl_trig[2*(PF_SHIFT_LIMITED_UNROLL_SIZE+PF_SHIFT_LIMITED_SIMD_SZ)];
    /* 4 times complex phase */
    float phase_state_i[PF_SHIFT_LIMITED_SIMD_SZ];
    float phase_state_q[PF_SHIFT_LIMITED_SIMD_SZ];
    /* N_cplx_per_block times increment - for future parallel variants */
    float dcos_blk;
    float dsin_blk;
    /* */
    float phase_increment;
} shift_limited_unroll_C_sse_data_t;

PFDSP_EXPORT shift_limited_unroll_C_sse_data_t shift_limited_unroll_C_sse_init(float relative_freq, float phase_start_rad);
PFDSP_EXPORT void shift_limited_unroll_C_sse_inp_c(complexf* in_out, int N_cplx, shift_limited_unroll_C_sse_data_t* d);



/*********************************************************************/

/**************/
/*** ALGO I ***/
/**************/

/* Recursive Quadrature Oscillator functions "recursive_osc"
 * see https://www.vicanek.de/articles/QuadOsc.pdf
 */
#define PF_SHIFT_RECURSIVE_SIMD_SZ  8
typedef struct shift_recursive_osc_s
{
    float u_cos[PF_SHIFT_RECURSIVE_SIMD_SZ];
    float v_sin[PF_SHIFT_RECURSIVE_SIMD_SZ];
} shift_recursive_osc_t;

typedef struct shift_recursive_osc_conf_s
{
    float k1;
    float k2;
} shift_recursive_osc_conf_t;

PFDSP_EXPORT void shift_recursive_osc_init(float rate, float starting_phase, shift_recursive_osc_conf_t *conf, shift_recursive_osc_t *state);
PFDSP_EXPORT void shift_recursive_osc_update_rate(float rate, shift_recursive_osc_conf_t *conf, shift_recursive_osc_t* state);

/* size must be multiple of PF_SHIFT_LIMITED_SIMD_SZ */
/* starting_phase for next call is kept internal in state */
PFDSP_EXPORT void shift_recursive_osc_cc(const complexf *input, complexf* output, int size, const shift_recursive_osc_conf_t *conf, shift_recursive_osc_t* state);
PFDSP_EXPORT void shift_recursive_osc_inp_c(complexf* output, int size, const shift_recursive_osc_conf_t *conf, shift_recursive_osc_t* state);
PFDSP_EXPORT void gen_recursive_osc_c(complexf* output, int size, const shift_recursive_osc_conf_t *conf, shift_recursive_osc_t* state);

/*********************************************************************/

/**************/
/*** ALGO J ***/
/**************/

#define PF_SHIFT_RECURSIVE_SIMD_SSE_SZ  4
typedef struct shift_recursive_osc_sse_s
{
    float u_cos[PF_SHIFT_RECURSIVE_SIMD_SSE_SZ];
    float v_sin[PF_SHIFT_RECURSIVE_SIMD_SSE_SZ];
} shift_recursive_osc_sse_t;

typedef struct shift_recursive_osc_sse_conf_s
{
    float k1;
    float k2;
} shift_recursive_osc_sse_conf_t;

PFDSP_EXPORT void shift_recursive_osc_sse_init(float rate, float starting_phase, shift_recursive_osc_sse_conf_t *conf, shift_recursive_osc_sse_t *state);
PFDSP_EXPORT void shift_recursive_osc_sse_update_rate(float rate, shift_recursive_osc_sse_conf_t *conf, shift_recursive_osc_sse_t* state);
PFDSP_EXPORT void shift_recursive_osc_sse_inp_c(complexf* in_out, int N_cplx, const shift_recursive_osc_sse_conf_t *conf, shift_recursive_osc_sse_t* state_ext);


#ifdef __cplusplus
}
#endif

