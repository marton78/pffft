/* Copyright (c) 2026  Márton Danóczy

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
*/

#ifndef PF_HELIUM_FLT_H
#define PF_HELIUM_FLT_H

/*
  Arm Helium (M-Profile Vector Extension, MVE) support macros, for
  Armv8.1-M cores such as Cortex-M55 and Cortex-M85.

  Experimental and opt-in via PFFFT_ENABLE_HELIUM: it is verified for
  correctness under QEMU, but not benchmarked on hardware, so it may be
  slower than the scalar code it replaces.

  The arithmetic maps 1:1 to NEON-style intrinsics. MVE has no zip, unzip,
  transpose or 64-bit half operations, though; the shuffles go through a
  small stack buffer using the interleaving loads/stores VST2/VLD2/VLD4, or
  move 64-bit halves as integer lanes. MVE has no double precision lanes,
  so this is float only.
*/
#if !defined(SIMD_SZ) && !defined(PFFFT_SIMD_DISABLE) && defined(PFFFT_ENABLE_HELIUM) && defined(__ARM_FEATURE_MVE) && (__ARM_FEATURE_MVE & 2)
#  include <arm_mve.h>
typedef float32x4_t v4sf;

#  define SIMD_SZ 4

typedef union v4sf_union {
  v4sf  v;
  float f[SIMD_SZ];
} v4sf_union;

#  define VARCH "Helium"
#  define VREQUIRES_ALIGN 0
#  define VZERO() vdupq_n_f32(0)
#  define VMUL(a,b) vmulq_f32(a,b)
#  define VADD(a,b) vaddq_f32(a,b)
#  define VMADD(a,b,c) vfmaq_f32(c,a,b)
#  define VMSUB(a,b,c) vfmsq_f32(c,a,b)
#  define VSUB(a,b) vsubq_f32(a,b)
#  define LD_PS1(p) vdupq_n_f32(p)
#  define VLOAD_UNALIGNED(ptr)  vld1q_f32((const float*)(ptr))
#  define VLOAD_ALIGNED(ptr)    vld1q_f32((const float*)(ptr))
#  define INTERLEAVE2(in1, in2, out1, out2) {                                \
    float tmp__[8] __attribute__((aligned(16)));                              \
    float32x4x2_t in__ = {{ in1, in2 }};                                      \
    vst2q_f32(tmp__, in__);                                                   \
    out1 = vld1q_f32(tmp__); out2 = vld1q_f32(tmp__ + 4);                     \
  }
#  define UNINTERLEAVE2(in1, in2, out1, out2) {                              \
    float tmp__[8] __attribute__((aligned(16)));                              \
    float32x4x2_t out__;                                                      \
    vst1q_f32(tmp__, in1); vst1q_f32(tmp__ + 4, in2);                         \
    out__ = vld2q_f32(tmp__);                                                 \
    out1 = out__.val[0]; out2 = out__.val[1];                                 \
  }
#  define VTRANSPOSE4(x0,x1,x2,x3) {                                          \
    float tmp__[16] __attribute__((aligned(16)));                             \
    float32x4x4_t t__;                                                        \
    vst1q_f32(tmp__, x0); vst1q_f32(tmp__ + 4, x1);                           \
    vst1q_f32(tmp__ + 8, x2); vst1q_f32(tmp__ + 12, x3);                      \
    t__ = vld4q_f32(tmp__);                                                   \
    x0 = t__.val[0]; x1 = t__.val[1]; x2 = t__.val[2]; x3 = t__.val[3];       \
  }
/* the 64-bit halves are moved as integer lanes */
static ALWAYS_INLINE(v4sf) VSWAPHL(v4sf a, v4sf b) {
  uint64x2_t r = vsetq_lane_u64(vgetq_lane_u64(vreinterpretq_u64_f32(b), 0), vreinterpretq_u64_f32(a), 0);
  return vreinterpretq_f32_u64(r);
}
/* reverse/flip complex floats: swap the 64-bit halves */
static ALWAYS_INLINE(v4sf) VREV_C(v4sf a) {
  uint64x2_t u = vreinterpretq_u64_f32(a);
  uint64x2_t r = vsetq_lane_u64(vgetq_lane_u64(u, 1), vsetq_lane_u64(vgetq_lane_u64(u, 0), u, 1), 0);
  return vreinterpretq_f32_u64(r);
}
/* reverse/flip all floats */
static ALWAYS_INLINE(v4sf) VREV_S(v4sf a) { return VREV_C(vrev64q_f32(a)); }

#  define VALIGNED(ptr) ((((uintptr_t)(ptr)) & 0x3) == 0)

#endif

#endif /* PF_HELIUM_FLT_H */
