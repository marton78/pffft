
#include "pf_conv.h"

#include <string.h>
#include <assert.h>

#include <algorithm>
#include <type_traits>

#if 0
#include <stdio.h>

#define DPRINT(...) fprintf(stderr, __VA_ARGS__)

#else
#define DPRINT(...) do { } while (0)
#endif


#ifndef CONV_ARCH_POST
#error CONV_ARCH_POST not defined
#endif

/* pffft's SIMD macros, as pffft.c uses them. The 4xScalar backend
 * (PFFFT_SCALVEC_ENABLED) needs ALWAYS_INLINE, which pffft.c defines. */
#ifndef ALWAYS_INLINE
#  if defined(_MSC_VER)
#    define ALWAYS_INLINE(return_type) __forceinline return_type
#  elif defined(__GNUC__)
#    define ALWAYS_INLINE(return_type) inline return_type __attribute__ ((always_inline))
#  else
#    define ALWAYS_INLINE(return_type) inline return_type
#  endif
#endif
#include "simd/pf_float.h"

#define PP_STRINGIFY(X) #X
#define PP_TOSTRING(X)  PP_STRINGIFY(X)
#define PP_CONCAT_IMPL(x, y) x##y
#define PP_CONCAT(x, y) PP_CONCAT_IMPL( x, y )

#define ARCHFUNCNAME(X) PP_CONCAT(X##_,CONV_ARCH_POST)


const char * ARCHFUNCNAME(id)()
{
    return PP_TOSTRING(CONV_ARCH_POST);
}


int ARCHFUNCNAME(conv_float_simd_size)()
{
    return SIMD_SZ;
}


void ARCHFUNCNAME(conv_float_move_rest)(float * RESTRICT s, conv_buffer_state * RESTRICT state)
{
    int R = state->size - state->offset;    // this many samples from prev conv_float were not processed
    if (R > 0)
    {
        // memmove(s, &s[state->offset], R * sizeof(s[0]));   // move them to the begin
        ::std::copy(&s[state->offset], &s[state->size], s);
    }
    else
        R = 0;
    state->offset = 0;      // data - to be processed - is at begin
    state->size = R;        // this many unprocessed samples
}


void ARCHFUNCNAME(conv_cplx_move_rest)(complexf * RESTRICT s, conv_buffer_state * RESTRICT state)
{
    int R = state->size - state->offset;    // this many samples from prev conv_float were not processed
    if (R > 0)
    {
        // memmove(s, &s[state->offset], R * sizeof(s[0]));   // move them to the begin
        ::std::copy(&s[state->offset], &s[state->size], s);
    }
    else
        R = 0;
    state->offset = 0;      // data - to be processed - is at begin
    state->size = R;        // this many unprocessed samples
}


/* the sum of the SIMD_SZ lanes of a vector */
static inline float hsum(v4sf v)
{
    v4sf_union u;
    u.v = v;
    float sum = u.f[0];
    for (int k = 1; k < SIMD_SZ; ++k)
        sum += u.f[k];
    return sum;
}


/* The signal is loaded at every sample offset, so unaligned; the filter is
 * loaded from the start, so aligned. The filter length is a multiple of
 * SIMD_SZ, which conv_float_simd_size() returns.
 *
 * The kernels are limited by the loads, two per multiply-add if each output
 * is computed by itself. So they compute eight real outputs (four complex
 * ones) per pass, which share each load of the filter: 9 loads for 8
 * multiply-adds. Their accumulators also keep eight multiply-adds in flight.
 * The outputs left over at the end are computed one at a time, each with
 * four accumulators. */


/* out[j] = the sum of the lanes of a<j>, for j < 4 */
static inline void hsum4(v4sf a0, v4sf a1, v4sf a2, v4sf a3, float * RESTRICT out)
{
#if SIMD_SZ == 4
    VTRANSPOSE4(a0, a1, a2, a3);
    v4sf_union u;
    u.v = VADD(VADD(a0, a1), VADD(a2, a3));
    out[0] = u.f[0]; out[1] = u.f[1]; out[2] = u.f[2]; out[3] = u.f[3];
#else
    out[0] = hsum(a0); out[1] = hsum(a1); out[2] = hsum(a2); out[3] = hsum(a3);
#endif
}


/* out[j] = sum_k s[j+k] * filter[k], for k < sz_filter, for j < 8 */
static inline void dot8(const float * RESTRICT s, const float * RESTRICT filter, const int sz_filter, float * RESTRICT out)
{
    v4sf a0 = VZERO(), a1 = VZERO(), a2 = VZERO(), a3 = VZERO();
    v4sf a4 = VZERO(), a5 = VZERO(), a6 = VZERO(), a7 = VZERO();
    for (int k = 0; k < sz_filter; k += SIMD_SZ)
    {
        const v4sf h = VLOAD_ALIGNED(&filter[k]);
        const float * RESTRICT sk = &s[k];
        a0 = VMADD(VLOAD_UNALIGNED(sk    ), h, a0);
        a1 = VMADD(VLOAD_UNALIGNED(sk + 1), h, a1);
        a2 = VMADD(VLOAD_UNALIGNED(sk + 2), h, a2);
        a3 = VMADD(VLOAD_UNALIGNED(sk + 3), h, a3);
        a4 = VMADD(VLOAD_UNALIGNED(sk + 4), h, a4);
        a5 = VMADD(VLOAD_UNALIGNED(sk + 5), h, a5);
        a6 = VMADD(VLOAD_UNALIGNED(sk + 6), h, a6);
        a7 = VMADD(VLOAD_UNALIGNED(sk + 7), h, a7);
    }
    hsum4(a0, a1, a2, a3, out);
    hsum4(a4, a5, a6, a7, out + 4);
}


/* sum_k s[k] * filter[k], for k < sz_filter */
static inline float dot(const float * RESTRICT s, const float * RESTRICT filter, const int sz_filter)
{
    v4sf accu0 = VZERO(), accu1 = VZERO(), accu2 = VZERO(), accu3 = VZERO();
    int k = 0;
    for (; k + 4*SIMD_SZ <= sz_filter; k += 4*SIMD_SZ)
    {
        accu0 = VMADD(VLOAD_UNALIGNED(&s[k            ]), VLOAD_ALIGNED(&filter[k            ]), accu0);
        accu1 = VMADD(VLOAD_UNALIGNED(&s[k +   SIMD_SZ]), VLOAD_ALIGNED(&filter[k +   SIMD_SZ]), accu1);
        accu2 = VMADD(VLOAD_UNALIGNED(&s[k + 2*SIMD_SZ]), VLOAD_ALIGNED(&filter[k + 2*SIMD_SZ]), accu2);
        accu3 = VMADD(VLOAD_UNALIGNED(&s[k + 3*SIMD_SZ]), VLOAD_ALIGNED(&filter[k + 3*SIMD_SZ]), accu3);
    }
    for (; k < sz_filter; k += SIMD_SZ)
        accu0 = VMADD(VLOAD_UNALIGNED(&s[k]), VLOAD_ALIGNED(&filter[k]), accu0);
    return hsum(VADD(VADD(accu0, accu1), VADD(accu2, accu3)));
}

int ARCHFUNCNAME(conv_float_inplace)(
        float * RESTRICT s, conv_buffer_state * RESTRICT state,
        const float * RESTRICT filter, const int sz_filter
        )
{
    assert( (sz_filter % SIMD_SZ) == 0 );  // size of filter must be divisible by conv_float_simd_size()
    assert( VALIGNED(filter) );

    const int off0 = state->offset;
    const int sz_s = state->size;
    int offset;

    // s[offset] is read by no later offset: written after it is read
    for ( offset = off0; offset + 8 + sz_filter <= sz_s + 1; offset += 8)
    {
        float out[8];
        dot8(&s[offset], filter, sz_filter, out);
        for (int j = 0; j < 8; ++j)
            s[offset + j] = out[j];
    }
    for ( ; offset + sz_filter <= sz_s; ++offset)
        s[offset] = dot(&s[offset], filter, sz_filter);

    state->offset = offset;
    return offset - off0;
}


int ARCHFUNCNAME(conv_float_oop)(
        const float * RESTRICT s, conv_buffer_state * RESTRICT state,
        const float * RESTRICT filter, const int sz_filter,
        float * RESTRICT y
        )
{
    assert( (sz_filter % SIMD_SZ) == 0 );  // size of filter must be divisible by conv_float_simd_size()
    assert( VALIGNED(filter) );

    const int off0 = state->offset;
    const int sz_s = state->size;
    int offset;

    for ( offset = off0; offset + 8 + sz_filter <= sz_s + 1; offset += 8)
        dot8(&s[offset], filter, sz_filter, &y[offset]);
    for ( ; offset + sz_filter <= sz_s; ++offset)
        y[offset] = dot(&s[offset], filter, sz_filter);

    state->offset = offset;
    return offset - off0;
}


#if SIMD_SZ > 1
/* *y = the sums of the even lanes (real parts) and of the odd lanes
 * (imaginary parts) of v */
static inline void cplx_sum(v4sf v, complexf * RESTRICT y)
{
    v4sf_union u;
    u.v = v;
    float accu_re = u.f[0];
    float accu_im = u.f[1];
    for (int k = 2; k < SIMD_SZ; k += 2)
    {
        accu_re += u.f[k];
        accu_im += u.f[k+1];
    }
    y->i = accu_re;
    y->q = accu_im;
}
#endif


int ARCHFUNCNAME(conv_cplx_float_oop)(
        const complexf * RESTRICT s_cplx, conv_buffer_state * RESTRICT state,
        const float * RESTRICT filter, const int sz_filter,
        complexf * RESTRICT y_cplx
        )
{
    assert( (sz_filter % SIMD_SZ) == 0 );  // size of filter must be divisible by conv_float_simd_size()
    assert( VALIGNED(filter) );

    const int off0 = state->offset;
    const int sz_s = state->size;
    int offset;

#if SIMD_SZ == 1
    for ( offset = off0; offset + sz_filter <= sz_s; ++offset)
    {
        float accu_re = 0.0F;
        float accu_im = 0.0F;
        for (int k = 0; k < sz_filter; ++k)
        {
            accu_re += s_cplx[offset+k].i * filter[k];
            accu_im += s_cplx[offset+k].q * filter[k];
        }
        y_cplx[offset].i = accu_re;
        y_cplx[offset].q = accu_im;
    }
#else
    /* SIMD_SZ complex samples are two vectors of interleaved real and
     * imaginary parts, which are multiplied with SIMD_SZ filter taps, each
     * duplicated: the even lanes accumulate real parts, the odd lanes
     * imaginary parts. Four outputs per pass share the filter, eight
     * accumulators; one output at a time at the end, with four. */
    const float * RESTRICT s = &(s_cplx[0].i);

    for ( offset = off0; offset + 4 + sz_filter <= sz_s + 1; offset += 4)
    {
        const float * RESTRICT so = &s[2*offset];
        v4sf a0 = VZERO(), a1 = VZERO(), b0 = VZERO(), b1 = VZERO();
        v4sf c0 = VZERO(), c1 = VZERO(), d0 = VZERO(), d1 = VZERO();
        for (int k = 0; k < sz_filter; k += SIMD_SZ)
        {
            v4sf rH = VLOAD_ALIGNED(&filter[k]), rH0, rH1;
            INTERLEAVE2(rH, rH, rH0, rH1);
            const float * RESTRICT sk = &so[2*k];
            a0 = VMADD(VLOAD_UNALIGNED(sk        ), rH0, a0);
            a1 = VMADD(VLOAD_UNALIGNED(sk + SIMD_SZ    ), rH1, a1);
            b0 = VMADD(VLOAD_UNALIGNED(sk + 2      ), rH0, b0);
            b1 = VMADD(VLOAD_UNALIGNED(sk + SIMD_SZ + 2), rH1, b1);
            c0 = VMADD(VLOAD_UNALIGNED(sk + 4      ), rH0, c0);
            c1 = VMADD(VLOAD_UNALIGNED(sk + SIMD_SZ + 4), rH1, c1);
            d0 = VMADD(VLOAD_UNALIGNED(sk + 6      ), rH0, d0);
            d1 = VMADD(VLOAD_UNALIGNED(sk + SIMD_SZ + 6), rH1, d1);
        }
        cplx_sum(VADD(a0, a1), &y_cplx[offset    ]);
        cplx_sum(VADD(b0, b1), &y_cplx[offset + 1]);
        cplx_sum(VADD(c0, c1), &y_cplx[offset + 2]);
        cplx_sum(VADD(d0, d1), &y_cplx[offset + 3]);
    }

    for ( ; offset + sz_filter <= sz_s; ++offset)
    {
        const float * RESTRICT so = &s[2*offset];
        v4sf accu0 = VZERO(), accu1 = VZERO(), accu2 = VZERO(), accu3 = VZERO();
        int k = 0;
        for (; k + 2*SIMD_SZ <= sz_filter; k += 2*SIMD_SZ)
        {
            v4sf rH = VLOAD_ALIGNED(&filter[k]), rH0, rH1;
            v4sf rG = VLOAD_ALIGNED(&filter[k + SIMD_SZ]), rG0, rG1;
            INTERLEAVE2(rH, rH, rH0, rH1);
            INTERLEAVE2(rG, rG, rG0, rG1);
            accu0 = VMADD(VLOAD_UNALIGNED(&so[2*k            ]), rH0, accu0);
            accu1 = VMADD(VLOAD_UNALIGNED(&so[2*k +   SIMD_SZ]), rH1, accu1);
            accu2 = VMADD(VLOAD_UNALIGNED(&so[2*k + 2*SIMD_SZ]), rG0, accu2);
            accu3 = VMADD(VLOAD_UNALIGNED(&so[2*k + 3*SIMD_SZ]), rG1, accu3);
        }
        if (k < sz_filter)
        {
            v4sf rH = VLOAD_ALIGNED(&filter[k]), rH0, rH1;
            INTERLEAVE2(rH, rH, rH0, rH1);
            accu0 = VMADD(VLOAD_UNALIGNED(&so[2*k          ]), rH0, accu0);
            accu1 = VMADD(VLOAD_UNALIGNED(&so[2*k + SIMD_SZ]), rH1, accu1);
        }
        cplx_sum(VADD(VADD(accu0, accu1), VADD(accu2, accu3)), &y_cplx[offset]);
    }
#endif

    state->offset = offset;
    return offset - off0;
}


static const conv_f_ptrs conv_ptrs =
{
    PP_TOSTRING(CONV_ARCH_POST),
    VARCH,

    ARCHFUNCNAME(id),
    ARCHFUNCNAME(conv_float_simd_size),

    ARCHFUNCNAME(conv_float_move_rest),
    ARCHFUNCNAME(conv_float_inplace),
    ARCHFUNCNAME(conv_float_oop),

    ARCHFUNCNAME(conv_cplx_move_rest),
    ARCHFUNCNAME(conv_cplx_float_oop)
};


const conv_f_ptrs* ARCHFUNCNAME(conv_ptrs)()
{
    DPRINT("arch pointer for '%s': %s, float simd size %d\n", conv_ptrs.id, conv_ptrs.simd_arch, SIMD_SZ);
    return &conv_ptrs;
}

// conv_ptrs_<arch>() has the type that pf_conv.h declares
static_assert(std::is_same<decltype(&ARCHFUNCNAME(conv_ptrs)), f_conv_ptrs>::value,
              "conv_ptrs does not match f_conv_ptrs");
