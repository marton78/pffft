/*
  Copyright (c) 2020  Hayati Ayguen ( h_ayguen@web.de )

  bench for mixer algorithm/implementations

 */

#include <pf_mixer.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <string.h>

#define HAVE_SYS_TIMES

#ifdef HAVE_SYS_TIMES
#  include <sys/times.h>
#  include <unistd.h>
#endif

#define BENCH_REF_TRIG_FUNC       1
#define BENCH_OUT_OF_PLACE_ALGOS  1
#define BENCH_INPLACE_ALGOS       1


#if defined(HAVE_SYS_TIMES)
    static double ttclk = 0.;

    static double uclock_sec(void)
    {
        struct tms t;
        if (ttclk == 0.)
            ttclk = sysconf(_SC_CLK_TCK);
        times(&t);
        /* use only the user time of this process - not realtime, which depends on OS-scheduler .. */
        return ((double)t.tms_utime) / ttclk;
    }
# else
    double uclock_sec(void)
    { return (double)clock()/(double)CLOCKS_PER_SEC; }
#endif


void save(complexf * d, int B, int N, const char * fn)
{
    if (!fn)
        fn = "/dev/shm/bench.bin";
    FILE * f = fopen(fn, "wb");
    if (!f) {
        fprintf(stderr, "error writing result to %s\n", fn);
        return;
    }
    for (int off = 0; off + B <= N; off += B)
    {
        fwrite(d+off, sizeof(complexf), B, f);
    }
    fclose(f);
}


double bench_shift_math_cc(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    complexf *output = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        phase = shift_math_cc(input+off, output+off, B, -0.0009F, phase);
        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(output, B, off, NULL);
    free(input);
    free(output);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_table_cc(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    int table_size=65536;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    complexf *output = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;

    shift_table_data_t table_data = shift_table_init(table_size);

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        phase = shift_table_cc(input+off, output+off, B, -0.0009F, table_data, phase);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(output, B, off, NULL);
    free(input);
    free(output);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_addfast(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    int table_size=65536;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    complexf *output = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;
    shift_addfast_data_t state = shift_addfast_init(-0.0009F);

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        phase = shift_addfast_cc(input+off, output+off, B, &state, phase);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(output, B, off, NULL);
    free(input);
    free(output);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_unroll_oop(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    complexf *output = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;
    shift_unroll_data_t state = shift_unroll_init(-0.0009F, B);

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        phase = shift_unroll_cc(input+off, output+off, B, &state, phase);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(output, B, off, NULL);
    free(input);
    free(output);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}

double bench_shift_unroll_inp(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;
    shift_unroll_data_t state = shift_unroll_init(-0.0009F, B);

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        phase = shift_unroll_inp_c(input+off, B, &state, phase);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(input, B, off, NULL);
    free(input);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}



double bench_shift_limited_unroll_oop(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    complexf *output = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;
    shift_limited_unroll_data_t state = shift_limited_unroll_init(-0.0009F);

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        shift_limited_unroll_cc(input+off, output+off, B, &state);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(output, B, off, NULL);
    free(input);
    free(output);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_limited_unroll_inp(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;
    shift_limited_unroll_data_t state = shift_limited_unroll_init(-0.0009F);

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        shift_limited_unroll_inp_c(input+off, B, &state);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(input, B, off, NULL);
    free(input);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_limited_unroll_sse_inp(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;
    shift_limited_unroll_sse_data_t *state = malloc(sizeof(shift_limited_unroll_sse_data_t));

    *state = shift_limited_unroll_sse_init(-0.0009F, 0.0F);

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        shift_limited_unroll_sse_inp_c(input+off, B, state);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(input, B, off, NULL);
    free(input);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_rec_osc_cc_oop(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    complexf *output = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state, shift_state;
    shift_recursive_osc_conf_t gen_conf, shift_conf;

    shift_recursive_osc_init(-0.0009F, 0.0F, &shift_conf, &shift_state);
    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        shift_recursive_osc_cc(input+off, output+off, B, &shift_conf, &shift_state);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    //save(input, B, off, "/dev/shm/bench_inp.bin");
    save(output, B, off, NULL);
    free(input);
    free(output);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_rec_osc_cc_inp(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state, shift_state;
    shift_recursive_osc_conf_t gen_conf, shift_conf;

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);
    shift_recursive_osc_init(-0.0009F, 0.0F, &shift_conf, &shift_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        shift_recursive_osc_inp_c(input+off, B, &shift_conf, &shift_state);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(input, B, off, NULL);
    free(input);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}


double bench_shift_rec_osc_sse_c_inp(int B, int N) {
    double t0, t1, tstop, T, nI;
    int iter, off;
    float phase = 0.0F;
    complexf *input = (complexf *)malloc(N * sizeof(complexf));
    shift_recursive_osc_t gen_state;
    shift_recursive_osc_conf_t gen_conf;

    shift_recursive_osc_sse_t *shift_state = malloc(sizeof(shift_recursive_osc_sse_t));
    shift_recursive_osc_sse_conf_t shift_conf;

    shift_recursive_osc_init(0.001F, 0.0F, &gen_conf, &gen_state);
    gen_recursive_osc_c(input, N, &gen_conf, &gen_state);

    shift_recursive_osc_sse_init(-0.0009F, 0.0F, &shift_conf, shift_state);

    iter = 0;
    off = 0;
    t0 = uclock_sec();
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */
    do {
        // work
        shift_recursive_osc_sse_inp_c(input+off, B, &shift_conf, shift_state);

        off += B;
        ++iter;
        t1 = uclock_sec();
    } while ( t1 < tstop && off + B < N );

    save(input, B, off, NULL);
    free(input);
    printf("processed %f Msamples\n", off * 1E-6);
    T = ( t1 - t0 );  /* duration per fft() */
    nI = ((double)iter) * B;  /* number of iterations "normalized" to O(N) = N */
    return (nI / T);    /* normalized iterations per second */
}



int main(int argc, char **argv)
{
    double rt;

    // process up to 64 MSample (512 MByte) in blocks of 8 kSamples (=64 kByte)
    int B = 8 * 1024;
    int N = 64 * 1024 * 1024;

#if BENCH_REF_TRIG_FUNC
    printf("\nstarting bench of shift_math_cc (out-of-place) with trig functions ..\n");
    rt = bench_shift_math_cc(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);
#endif

#if BENCH_OUT_OF_PLACE_ALGOS
    printf("starting bench of shift_table_cc (out-of-place) ..\n");
    rt = bench_shift_table_cc(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("starting bench of shift_addfast_cc (out-of-place) ..\n");
    rt = bench_shift_addfast(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("\nstarting bench of shift_unroll_cc (out-of-place) ..\n");
    rt = bench_shift_unroll_oop(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("\nstarting bench of shift_limited_unroll_cc (out-of-place) ..\n");
    rt = bench_shift_limited_unroll_oop(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("\nstarting bench of shift_recursive_osc_cc (out-of-place) ..\n");
    rt = bench_shift_rec_osc_cc_oop(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);
#endif

#if BENCH_INPLACE_ALGOS
    printf("starting bench of shift_unroll_cc in-place ..\n");
    rt = bench_shift_unroll_inp(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("starting bench of shift_limited_unroll_cc in-place ..\n");
    rt = bench_shift_limited_unroll_inp(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("starting bench of shift_limited_unroll_sse_inp_c in-place ..\n");
    rt = bench_shift_limited_unroll_sse_inp(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("starting bench of shift_recursive_osc_cc in-place ..\n");
    rt = bench_shift_rec_osc_cc_inp(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

    printf("starting bench of shift_recursive_osc_sse_c in-place ..\n");
    rt = bench_shift_rec_osc_sse_c_inp(B, N);
    printf("  %f MSamples/sec\n\n", rt * 1E-6);

#endif

    return 0;
}

