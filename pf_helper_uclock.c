/*
  Copyright (c) 2021  Hayati Ayguen ( h_ayguen@web.de )

  helper function(s)

 */

#include "pf_helper_uclock.h"

/* #include <math.h> */
#include <stdio.h>
#include <time.h>

#if defined(__linux__)
#define HAVE_SYS_TIMES
#endif

#ifdef HAVE_SYS_TIMES
#  include <sys/times.h>
#  include <unistd.h>
#endif

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>
#endif


#if defined(HAVE_SYS_TIMES)
    static double ttclk = 0.;

    double pf_uclock_sec(int find_start)
    {
        struct tms t0, t;
        if (ttclk == 0.)
        {
            ttclk = sysconf(_SC_CLK_TCK);
            fprintf(stderr, "sysconf(_SC_CLK_TCK) => %f\n", ttclk);
        }
        times(&t);
        if (find_start)
        {
            t0 = t;
            while (t0.tms_utime == t.tms_utime)
                times(&t);
        }
        /* use only the user time of this process - not realtime, which depends on OS-scheduler .. */
        return ((double)t.tms_utime) / ttclk;
    }

#elif defined(WIN32)
    // https://docs.microsoft.com/en-us/windows/win32/api/processthreadsapi/nf-processthreadsapi-getprocesstimes
    double pf_uclock_sec(int find_start)
    {
        FILETIME a, b, c, d;
        if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0)
        {
            //  Returns total user time.
            //  Can be tweaked to include kernel times as well.
            return
                (double)(d.dwLowDateTime |
                    ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
        }
        else {
            //  Handle error
            return 0;
        }
    }

#else
    double pf_uclock_sec(int find_start)
    { return (double)clock()/(double)CLOCKS_PER_SEC; }
#endif


