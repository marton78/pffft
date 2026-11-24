
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
 * SIMD_SZ, which conv_float_simd_size() returns. */

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

    for ( offset = off0; offset + sz_filter <= sz_s; ++offset)
    {
        v4sf accu = VZERO();
        for (int k = 0; k < sz_filter; k += SIMD_SZ)
        {
            v4sf rS = VLOAD_UNALIGNED(&s[offset+k]);
            v4sf rH = VLOAD_ALIGNED(&filter[k]);
            accu = VMADD(rS, rH, accu);   // accu += rS * rH;
        }
        // s[offset] is read by no later offset
        s[offset] = hsum(accu);
    }

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

    for ( offset = off0; offset + sz_filter <= sz_s; ++offset)
    {
        v4sf accu = VZERO();
        for (int k = 0; k < sz_filter; k += SIMD_SZ)
        {
            v4sf rS = VLOAD_UNALIGNED(&s[offset+k]);
            v4sf rH = VLOAD_ALIGNED(&filter[k]);
            accu = VMADD(rS, rH, accu);   // accu += rS * rH;
        }
        y[offset] = hsum(accu);
    }

    state->offset = offset;
    return offset - off0;
}


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
     * imaginary parts. */
    const float * RESTRICT s = &(s_cplx[0].i);

    for ( offset = off0; offset + sz_filter <= sz_s; ++offset)
    {
        v4sf accu0 = VZERO();
        v4sf accu1 = VZERO();
        for (int k = 0; k < sz_filter; k += SIMD_SZ)
        {
            v4sf rS0 = VLOAD_UNALIGNED(&s[2*(offset+k)]);
            v4sf rS1 = VLOAD_UNALIGNED(&s[2*(offset+k) + SIMD_SZ]);
            v4sf rH = VLOAD_ALIGNED(&filter[k]);
            v4sf rH0, rH1;
            INTERLEAVE2(rH, rH, rH0, rH1);
            accu0 = VMADD(rS0, rH0, accu0);   // accu += rS * rH;
            accu1 = VMADD(rS1, rH1, accu1);
        }
        v4sf_union u;
        u.v = VADD(accu0, accu1);
        float accu_re = u.f[0];
        float accu_im = u.f[1];
        for (int k = 2; k < SIMD_SZ; k += 2)
        {
            accu_re += u.f[k];
            accu_im += u.f[k+1];
        }
        y_cplx[offset].i = accu_re;
        y_cplx[offset].q = accu_im;
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
