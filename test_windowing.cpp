
#include <pf_windowing.h>

#include <cmath>
#include <stdio.h>

#define WIN_SIZE 64
#define ERR_LIMIT 1E-5

int main(int argc, char** argv)
{
    float win_coeff_trig[WIN_SIZE], win_coeff_cordic[WIN_SIZE], delta[WIN_SIZE];
    int showUsage = 0;

    if (argc == 1)
        showUsage = 1;

    for (int winNo = 0; winNo < pf_num_windows; ++winNo)
    {
        for (int w0 = 0; w0 <= 1; ++w0 )
        {
            pf_windowing_T win = (pf_windowing_T)winNo;
            const int N = WIN_SIZE;

            for (int k = 0; k < N; ++k)
                win_coeff_trig[k] = 1.0F;
            //printf("calculate trig win ..\n");
            pf_window_trig(win_coeff_trig, win, N, 1.0F, w0);

            for (int k = 0; k < N; ++k)
                win_coeff_cordic[k] = 1.0F;
            //printf("calculate cordic win ..\n");
            const struct pf_windowing_param_tag *wp = pf_window_alloc(win, N, w0, 0.16F);
            pf_window(win_coeff_cordic, wp);

            double wind_sum = pf_window_sum(wp);
            int trig_max_idx = 0, cordic_max_idx = 0, delta_max_idx = 0;
            float trig_max = win_coeff_trig[0];
            float cordic_max = win_coeff_cordic[0];
            float delta_max = win_coeff_cordic[0] - win_coeff_trig[0];
            for (int k = 0; k < N; ++k)
            {
                delta[k] = std::abs(win_coeff_cordic[k] - win_coeff_trig[k]);
                if (delta_max < delta[k]) {
                    delta_max = delta[k];
                    delta_max_idx = k;
                }
                if (trig_max < win_coeff_trig[k]) {
                    trig_max = win_coeff_trig[k];
                    trig_max_idx = k;
                }
                if (cordic_max < win_coeff_cordic[k]) {
                    cordic_max = win_coeff_cordic[k];
                    cordic_max_idx = k;
                }
            }

            printf("\ntrig_max: %f @ %d\tcordic_max: %f @ %d\tdelta_max: %g @ %d\tfor %s  w0 = %d\n",
                trig_max, trig_max_idx, cordic_max, cordic_max_idx, delta_max, delta_max_idx,
                pf_window_text(win), w0 );
            printf("sum: %g\t0: %g / %g\tN-1: %g / %g\n", wind_sum
                 , win_coeff_trig[0], win_coeff_cordic[0]
                 , win_coeff_trig[N-1], win_coeff_cordic[N-1] );

            printf("trig:   ");
            for (int k=0; k < 8; ++k)
                printf("%d: %.2f\t", k * (N-1) / 7, win_coeff_trig[k * (N-1) / 7]);
            printf("\n");
            printf("cordic: ");
            for (int k=0; k < 8; ++k)
                printf("%d: %.2f\t", k * (N-1) / 7, win_coeff_cordic[k * (N-1) / 7]);
            printf("\n");

            if (delta_max > ERR_LIMIT)
            {
                printf("trig and cordic windows coefficient differ too much: %g @ idx %d -> test failed!\n", delta_max, delta_max_idx );
                return 1;
            }

            int N_div_i = 0;
            pf_get_window_param(wp, &N_div_i);

            int expected_peak_idx = (w0 == 0) ? 0 : (N / 2 -1);
            int expected_peak_idx2 = (w0 == 0) ? (N -1) : expected_peak_idx;
            //  && N_div_i == (N -1)
            if (trig_max_idx != expected_peak_idx && trig_max_idx != expected_peak_idx2 && win != pf_winRect)
            {
                printf("trig: expected maximum at index %d for w0 = %d: failed!\n", expected_peak_idx, w0);
                return 1;
            }
            if (cordic_max_idx != expected_peak_idx && cordic_max_idx != expected_peak_idx2 && win != pf_winRect)
            {
                printf("cordic: expected maximum at index %d for w0 = %d: failed!\n", expected_peak_idx, w0);
                return 1;
            }
        }
    }

    return 0;
}
