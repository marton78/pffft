/*
  Copyright (c) 2021  Hayati Ayguen ( h_ayguen@web.de )

  bench for windowing algorithm/implementations

 */

#include <pf_windowing.h>
#include <pf_helper_uclock.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


template <class TYPE>
double bench_windowing_trig_cached(int N, int normalize, pf_windowing_T w) {
    double t0, t1, tstop, T, nI, windSumI;
    int iter, off, k;
    int M = 512 * 1024 * 1024 / ( N * sizeof(TYPE) );
    float windScale;
    TYPE *data = (TYPE *)malloc(M * N * sizeof(TYPE));
    float *wind = (float *)malloc(N * sizeof(TYPE));

    for (k =0; k < M * N; ++k)
        data[k] = 1.0F;

    iter = 0;
    t0 = pf_uclock_sec(1);
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */

    for (k =0; k < N; ++k)
        wind[k] = 1.0F;

    /* calculate window */
    pf_win_trig(wind, w, N, 1.0F, 1, 0.16F);  /* pf_win_trig(.., int N, float scale, int w0, float alpha) */

    windScale = 1.0F;
    if (normalize)
    {
        /* window sum - for normalization */
        windSumI = 0.0;
        for (k =0; k < N; ++k)
            windSumI += wind[k];
        windScale = (float)N / windSumI;
        /* normalize window */
        for (k =0; k < N; ++k)
            wind[k] *= windScale;
    }

    int datapos = 0;
    do {
        TYPE *dataptr = data + datapos;
        /* work: multiply with cached and scaled window wind[] */
        for (k =0; k < N; ++k)
            dataptr[k] *= wind[k];
        datapos += N;
        if (datapos >= N*M)
            datapos = 0;
        ++iter;
        t1 = pf_uclock_sec(0);
    } while ( t1 < tstop );

    printf("windowing scale is %f\n", windScale);
    T = ( t1 - t0 );  /* duration per fft() */
    printf("processed %d windows in %f ms\n", iter, T*1E3);
    nI = iter / (T*1E3);
    printf("  %f windows per ms\n", nI);

    free(wind);
    free(data);
    return nI;
}


template <class TYPE, class SUMTYPE>
double bench_windowing_trig_uncached(int N, int normalize, pf_windowing_T w) {
    double t0, t1, tstop, T, nI;
    SUMTYPE windSumI;
    int iter, off, k;
    int M = 512 * 1024 * 1024 / ( N * sizeof(TYPE) );
    float windScale;
    TYPE *data = (TYPE *)malloc(M * N * sizeof(TYPE));
    for (k =0; k < M * N; ++k)
        data[k] = 1.0F;

    iter = 0;
    t0 = pf_uclock_sec(1);
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */

    /* calculate window */
    pf_win_trig(data, w, N, 1.0F, 1, 0.16F);  /* pf_win_trig(.., int N, float scale, int w0, float alpha) */
    /* window sum - for normalization */
    windSumI = 0.0;
    for (k =0; k < N; ++k)
        windSumI += data[k];
    windScale = normalize ? ((float)N / real(windSumI)) : 1.0F;

    int datapos = 0;
    do {
        TYPE *dataptr = data + datapos;
        /* work */
        pf_win_trig(dataptr, w, N, windScale, 1, 0.16F);
        datapos += N;
        if (datapos >= N*M)
            datapos = 0;
        ++iter;
        t1 = pf_uclock_sec(0);
    } while ( t1 < tstop );

    printf("windowing scale is %f\n", windScale);
    T = ( t1 - t0 );  /* duration per fft() */
    printf("processed %d windows in %f ms\n", iter, T*1E3);
    nI = iter / (T*1E3);
    printf("  %f windows per ms\n", nI);

    free(data);
    return nI;
}


template <class TYPE, class SUMTYPE>
double bench_windowing_uncached(int N, int normalize, pf_windowing_T w) {
    double t0, t1, tstop, T, nI;
    SUMTYPE windSumI;
    int iter, off, k;
    int M = 512 * 1024 * 1024 / ( N * sizeof(TYPE) );
    float windScale;
    TYPE *data = (TYPE *)malloc(M * N * sizeof(TYPE));
    for (k =0; k < M * N; ++k)
        data[k] = 1.0F;

    iter = 0;
    t0 = pf_uclock_sec(1);
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */

    /* calculate window */
    const struct pf_windowing_param_tag *wp = pf_window_alloc(w, N, 1, 0.16F);
    /* window sum - for normalization */
    windSumI = pf_window_sum(wp);
    windScale = normalize ? ((float)N / real(windSumI)) : 1.0F;

    int datapos = 0;
    do {
        TYPE *dataptr = data + datapos;
        /* work */
        pf_win(dataptr, wp, windScale);
        datapos += N;
        if (datapos >= N*M)
            datapos = 0;
        ++iter;
        t1 = pf_uclock_sec(0);
    } while ( t1 < tstop );

    printf("windowing scale is %f\n", windScale);
    T = ( t1 - t0 );  /* duration per fft() */
    printf("processed %d windows in %f ms\n", iter, T*1E3);
    nI = iter / (T*1E3);
    printf("  %f windows per ms\n", nI);

    pf_window_free(wp);
    free(data);
    return nI;
}


template <class TYPE, class SUMTYPE>
double bench_windowing_cached(int N, int normalize, pf_windowing_T w) {
    double t0, t1, tstop, T, nI;
    SUMTYPE windSumI;
    int iter, off, k;
    int M = 512 * 1024 * 1024 / ( N * sizeof(TYPE) );
    float windScale;
    TYPE *data = (TYPE *)malloc(M * N * sizeof(TYPE));
    float *wind = (float *)malloc(N * sizeof(TYPE));
    for (k =0; k < M * N; ++k)
        data[k] = 1.0F;

    iter = 0;
    t0 = pf_uclock_sec(1);
    tstop = t0 + 0.5;  /* benchmark duration: 500 ms */

    for (k =0; k < N; ++k)
        wind[k] = 1.0F;

    /* calculate window */
    const struct pf_windowing_param_tag *wp = pf_window_alloc(w, N, 1, 0.16F);
#if 0
    /* window sum - for normalization */
    windSumI = pf_window_sum(wp);
    windScale = normalize ? ((float)N / real(windSumI)) : 1.0F;
    pf_window(wind, wp, windScale);
#else
    pf_window(wind, wp);

    windScale = 1.0F;
    if (normalize)
    {
        /* window sum - for normalization */
        windSumI = 0.0;
        for (k =0; k < N; ++k)
            windSumI += wind[k];
        windScale = (float)N / windSumI;
        /* normalize window */
        for (k =0; k < N; ++k)
            wind[k] *= windScale;
    }
#endif

    int datapos = 0;
    do {
        TYPE *dataptr = data + datapos;
        /* work: multiply with cached and scaled window wind[] */
        for (k =0; k < N; ++k)
            dataptr[k] *= wind[k];
        datapos += N;
        if (datapos >= N*M)
            datapos = 0;
        ++iter;
        t1 = pf_uclock_sec(0);
    } while ( t1 < tstop );

    printf("windowing scale is %f\n", windScale);
    T = ( t1 - t0 );  /* duration per fft() */
    printf("processed %d windows in %f ms\n", iter, T*1E3);
    nI = iter / (T*1E3);
    printf("  %f windows per ms\n", nI);

    pf_window_free(wp);
    free(wind);
    free(data);
    return nI;
}


template <class TYPE, class SUMTYPE>
void bench(const char * type_str, const char * sum_str, const int N, const int normalize) {
    const pf_windowing_T aBenchWins[] = { pf_winHamming, pf_winBlackmanHarris, pf_winFlatTop_HFT90D };
    const char * acWinStrs[] = { "pf_winHamming", "pf_winBlackmanHarris", "pf_winFlatTop_HFT90D" };
    double rt, rtA, rtB, rtC, rtD;

    for ( int wno = 0; wno < 3; ++wno )
    {
        printf("\n==============================================================\n");
        printf("\nstarting bench of bench_windowing_trig_uncached<%s, %s>() with %s\n", type_str, "complexd", acWinStrs[wno]);
        rtA = bench_windowing_trig_uncached<TYPE, complexd>(N, normalize, aBenchWins[wno]);
        /* printf("  %f MSamples/sec\n\n", rt * 1E-6); */

        printf("\nstarting bench of bench_windowing_trig_cached<%s>() with %s\n", type_str, acWinStrs[wno]);
        rtB = bench_windowing_trig_cached<TYPE>(N, normalize, aBenchWins[wno]);

        printf("\nstarting bench of bench_windowing_uncached<%s, %s>() with %s\n", type_str, sum_str, acWinStrs[wno]);
        rtC = bench_windowing_uncached<TYPE, SUMTYPE>(N, normalize, aBenchWins[wno]);

        printf("\nstarting bench of bench_windowing_cached<%s, %s>() with %s\n", type_str, sum_str, acWinStrs[wno]);
        rtD = bench_windowing_cached<TYPE, SUMTYPE>(N, normalize, aBenchWins[wno]);

        printf("\nwindowing_trig_cached() / windowing()              = %f\n", rtB / rtC);
        printf("windowing() / windowing_trig_uncached()            = %f\n", rtC / rtA);
        printf("bench_windowing_cached() / windowing_trig_cached() = %f\n", rtD / rtB);
    }
    printf("\n==============================================================\n");
}


int main(int argc, char **argv)
{
    int N = 1024;
    int showUsage = 0;
    int normalize = 1;
    int cx = 1;

    if (argc == 1)
        showUsage = 1;
    if (1 < argc)
        N = atoi(argv[1]);
    if (2 < argc)
        normalize = atoi(argv[2]);
    if (3 < argc)
        cx = atoi(argv[3]);
    if ( !N || showUsage )
    {
        fprintf(stderr, "%s [<window_size> [<normalize: 0/1, default 1> [<complex: 0/1 default 1>] ] ]\n", argv[0]);
        return 0;
    }

    fprintf(stderr, "will benchmark window size N = %d with%s for %s data\n", N, normalize?"":"out", cx ? "complex":"real");
    pf_uclock_sec(0);

    if (!cx)
    {
        bench<float, double>("float", "double", N, normalize);
    }
    else
    {
        bench<complexf, double>("complexf", "double", N, normalize);
    }

    return 0;
}

