/*
  test_pf_conv: the "slow" convolution kernels of benchmarks/pf_conv.cpp, of
  every architecture variant the dispatcher provides (those the CPU can run),
  against a plain convolution in double precision.

  The kernels compute y[n] = sum_k s[n+k] * h[k], for every n where the
  filter fits into the signal. Each is called on the whole signal at once,
  from a start offset, and in blocks the way bench_pf_conv_float streams the
  signal through a buffer with conv_*_move_rest().
*/

#include "pf_conv_dispatcher.h"
#include "pf_conv.h"

#include "pffft/pffft.hpp"  // pffft::AlignedVector

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <random>
#include <vector>


static const float GUARD = 1234.5F;

static int n_checks = 0;
static int n_failures = 0;

static void fail(const char * arch, const char * kernel, const char * what, int L, int detail)
{
    if (n_failures < 20)
        fprintf(stderr, "FAIL: arch '%s', %s, filter length %d: %s (%d)\n", arch, kernel, L, what, detail);
    ++n_failures;
}

/* y[n] == sum_k s[n+k] * h[k], within the rounding error of any order of
   summation in single precision */
static bool close_to_ref(float y, const float * s, int stride, const float * h, int L)
{
    double ref = 0.0, abs_sum = 0.0;
    for (int k = 0; k < L; ++k)
    {
        ref += double(s[k * stride]) * double(h[k]);
        abs_sum += fabs(double(s[k * stride]) * double(h[k]));
    }
    ++n_checks;
    return fabs(double(y) - ref) <= (L + 1) * double(FLT_EPSILON) * abs_sum + double(FLT_MIN);
}


static void test_arch(const conv_f_ptrs & c, int L, std::mt19937 & rng)
{
    std::uniform_real_distribution<float> dist(-1.0F, 1.0F);
    const int N = 4 * L + 61;   // real signal length
    const int M = 2 * L + 37;   // complex signal length

    pffft::AlignedVector<float> h(L);
    for (int k = 0; k < L; ++k)
        h[k] = dist(rng);

    // signals at every alignment: the kernels load them at every offset
    pffft::AlignedVector<float> s_buf(N + 8), y_buf(N + 8);
    pffft::AlignedVector<complexf> c_buf(M + 8), cy_buf(M + 8);
    for (int mis = 0; mis < 4; ++mis)
    {
        float * s = &s_buf[mis];
        float * y = &y_buf[mis];
        const int mc = mis & 1;
        complexf * cs = &c_buf[mc];
        complexf * cy = &cy_buf[mc];
        for (int n = 0; n < N; ++n)
            s[n] = dist(rng);
        s[N] = GUARD;
        for (int n = 0; n < M; ++n)
        {
            cs[n].i = dist(rng);
            cs[n].q = dist(rng);
        }
        cs[M].i = cs[M].q = GUARD;

        // 0 to 7: with the eight outputs per pass, every number of outputs
        // left over for the end, on every SIMD size
        for (int off0 = 0; off0 < 8; ++off0)
        {
            conv_buffer_state st;

            // out-of-place, real
            for (int n = 0; n <= N; ++n)
                y[n] = GUARD;
            st.offset = off0; st.size = N;
            int n_out = c.fp_conv_float_oop(s, &st, h.data(), L, y);
            if (n_out != N - L + 1 - off0 || st.offset != N - L + 1)
                fail(c.id, "conv_float_oop", "number of outputs", L, n_out);
            for (int n = 0; n <= N; ++n)
            {
                const bool computed = (n >= off0 && n <= N - L);
                if (computed ? !close_to_ref(y[n], &s[n], 1, h.data(), L) : y[n] != GUARD)
                    fail(c.id, "conv_float_oop", computed ? "wrong output" : "wrote outside the output", L, n);
            }
            if (s[N] != GUARD)
                fail(c.id, "conv_float_oop", "wrote to the signal", L, N);

            // in-place, real: the input is gone, so compare with a copy
            std::vector<float> s0(s, s + N);
            st.offset = off0; st.size = N;
            n_out = c.fp_conv_float_inplace(s, &st, h.data(), L);
            if (n_out != N - L + 1 - off0 || st.offset != N - L + 1)
                fail(c.id, "conv_float_inplace", "number of outputs", L, n_out);
            for (int n = 0; n <= N; ++n)
            {
                const bool computed = (n >= off0 && n <= N - L);
                const float expected = (n < N) ? s0[n] : GUARD;
                if (computed ? !close_to_ref(s[n], &s0[n], 1, h.data(), L) : s[n] != expected)
                    fail(c.id, "conv_float_inplace", computed ? "wrong output" : "wrote outside the output", L, n);
            }
            for (int n = 0; n < N; ++n)
                s[n] = s0[n];

            // out-of-place, complex signal and real filter
            for (int n = 0; n <= M; ++n)
                cy[n].i = cy[n].q = GUARD;
            st.offset = off0; st.size = M;
            n_out = c.fp_conv_cplx_float_oop(cs, &st, h.data(), L, cy);
            if (n_out != M - L + 1 - off0 || st.offset != M - L + 1)
                fail(c.id, "conv_cplx_float_oop", "number of outputs", L, n_out);
            for (int n = 0; n <= M; ++n)
            {
                const bool computed = (n >= off0 && n <= M - L);
                if (computed)
                {
                    if (!close_to_ref(cy[n].i, &cs[n].i, 2, h.data(), L))
                        fail(c.id, "conv_cplx_float_oop", "wrong real part", L, n);
                    if (!close_to_ref(cy[n].q, &cs[n].q, 2, h.data(), L))
                        fail(c.id, "conv_cplx_float_oop", "wrong imaginary part", L, n);
                }
                else if (cy[n].i != GUARD || cy[n].q != GUARD)
                    fail(c.id, "conv_cplx_float_oop", "wrote outside the output", L, n);
            }
            if (cs[M].i != GUARD || cs[M].q != GUARD)
                fail(c.id, "conv_cplx_float_oop", "wrote to the signal", L, M);
        }
    }

    // streamed in blocks through a buffer, as bench_pf_conv_float does
    const int block_lens[] = { 1, 3, L, 2 * L + 5 };
    for (int b = 0; b < 4; ++b)
    {
        const int B = block_lens[b];
        float * s = &s_buf[0];
        complexf * cs = &c_buf[0];

        pffft::AlignedVector<float> buf(B + L + 1);
        std::vector<float> y(N, GUARD);
        conv_buffer_state st;
        st.offset = st.size = 0;
        int n_out = 0;
        for (int off = 0; off + B <= N; off += B)
        {
            c.fp_conv_float_move_rest(buf.data(), &st);
            for (int k = 0; k < B; ++k)
                buf[st.size + k] = s[off + k];
            st.size += B;
            n_out += c.fp_conv_float_oop(buf.data(), &st, h.data(), L, &y[n_out]);
        }
        const int n_expected = (N / B) * B - L + 1;
        if (n_out != (n_expected > 0 ? n_expected : 0))
            fail(c.id, "conv_float_oop in blocks", "number of outputs", L, n_out);
        for (int n = 0; n < n_out && n < N; ++n)
            if (!close_to_ref(y[n], &s[n], 1, h.data(), L))
                fail(c.id, "conv_float_oop in blocks", "wrong output", L, n);

        pffft::AlignedVector<complexf> cbuf(B + L + 1);
        std::vector<complexf> cy(M);
        st.offset = st.size = 0;
        n_out = 0;
        for (int off = 0; off + B <= M; off += B)
        {
            c.fp_conv_cplx_move_rest(cbuf.data(), &st);
            for (int k = 0; k < B; ++k)
                cbuf[st.size + k] = cs[off + k];
            st.size += B;
            n_out += c.fp_conv_cplx_float_oop(cbuf.data(), &st, h.data(), L, &cy[n_out]);
        }
        const int m_expected = (M / B) * B - L + 1;
        if (n_out != (m_expected > 0 ? m_expected : 0))
            fail(c.id, "conv_cplx_float_oop in blocks", "number of outputs", L, n_out);
        for (int n = 0; n < n_out && n < M; ++n)
            if (!close_to_ref(cy[n].i, &cs[n].i, 2, h.data(), L) || !close_to_ref(cy[n].q, &cs[n].q, 2, h.data(), L))
                fail(c.id, "conv_cplx_float_oop in blocks", "wrong output", L, n);
    }
}


int main()
{
    int n_arch = 0;
    const ptr_to_conv_f_ptrs * arches = get_all_conv_arch_ptrs(&n_arch);
    std::mt19937 rng(42);
    int n_tested = 0;

    for (int a = 0; a < n_arch; ++a)
    {
        if (!arches[a])
        {
            // none and dflt are always there; the others where the CPU has
            // the instructions their -march allows
            if (a < 2)
            {
                fprintf(stderr, "FAIL: architecture %d is not available\n", a);
                ++n_failures;
            }
            else
                printf("arch %d: not supported by this CPU, skipped\n", a);
            continue;
        }
        const conv_f_ptrs & c = *arches[a];
        const int sz = c.fp_conv_float_simd_size();
        printf("arch '%s': SIMD backend %s, float SIMD size %d\n", c.id, c.simd_arch, sz);
        // filter lengths are multiples of the SIMD size; any number of
        // vectors, for the tails of the unrolled loops
        const int n_vecs[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33 };
        for (unsigned i = 0; i < sizeof(n_vecs) / sizeof(n_vecs[0]); ++i)
            test_arch(c, n_vecs[i] * sz, rng);
        ++n_tested;
    }

#ifdef PF_CONV_DFLT_IS_PFFFT_BACKEND
    // dflt is compiled like PFFFT: if pf_conv.cmake left out a definition
    // of PFFFT's, dflt would test another backend, e.g. 4xScalar for Helium
    if (n_arch > 1 && arches[1] && strcmp(arches[1]->simd_arch, pffft::Types<float>::simd_arch()) != 0)
    {
        fprintf(stderr, "FAIL: dflt uses the SIMD backend %s, PFFFT %s\n", arches[1]->simd_arch, pffft::Types<float>::simd_arch());
        ++n_failures;
    }
#endif

    printf("%d architectures, %d outputs checked, %d failures\n", n_tested, n_checks, n_failures);
    return (n_failures || !n_tested) ? 1 : 0;
}
