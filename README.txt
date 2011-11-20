PFFFT: a pretty fast FFT.

TL;DR
--

PFFFT does 1D Fast Fourier Transforms, of single precision real and
complex vectors. It tries do it fast, it tries to be correct, and it
tries to be small. Computations do take advantage of SSE1 instructions
on x86 cpus, Altivec on powerpc cpus, and NEON on ARM cpus. The
license is BSD-like.


Why does it exist:
--

I was in search of a good performing FFT library , preferably very
small and with a very liberal license.

When one says "fft library", FFTW ("Fastest Fourier Transform in the
West") is probably the first name that comes to mind -- I guess that
99% of open-source projects that need a FFT do use FFTW, and are happy
with it. However, it is quite a large library , which does everything
fft related (2d transforms, 3d transforms, other transformations such
as discrete cosine , or fast hartley). And it is licensed under the
GNU GPL , which means that it cannot be used in non open-source
products.

An alternative to FFTW that is really small, is the venerable FFTPACK
v4, which is available on NETLIB. A more recent version (v5) exists,
but it is larger as it deals with multi-dimensional transforms. This
is a library that is written in FORTRAN 77, a language that is now
considered as a bit antiquated by many. FFTPACKv4 was written in 1985,
by Dr Paul Swarztrauber of NCAR, more than 25 years ago ! And despite
its age, benchmarks show it that it still a very good performing FFT
library, see for example the 1d single precision benchmarks here:
http://www.fftw.org/speed/opteron-2.2GHz-32bit/ . It is however not
competitive with the fastest ones, such as FFTW, Intel MKL, AMD ACML,
Apple vDSP. The reason for that is that those libraries do take
advantage of the SSE SIMD instructions available on Intel CPUs,
available since the days of the Pentium III. These instructions deal
with small vectors of 4 floats at a time, instead of a single float
for a traditionnal FPU, so when using these instructions one may expect
a 4-fold performance improvement.

The idea was to take this fortran fftpack v4 code, translate to C,
modify it to deal with those SSE instructions, and check that the
final performance is not completely ridiculous when compared to other
SIMD FFT libraries. Translation to C was performed with f2c (
http://www.netlib.org/f2c/ ). The resulting file was a bit edited in
order to remove the thousands of gotos that were introduced by
f2c. You will find the fftpack.h and fftpack.c sources in the
repository, this a complete translation of
http://www.netlib.org/fftpack/ , with the discrete cosine transform
and the test program. There is no license information in the netlib
repository, but it was confirmed to me by the fftpack v5 curators that
the same terms do apply to fftpack v4:
http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html . This is a
"BSD-like" license, it is compatible with proprietary projects.

Adapting fftpack to deal with the SIMD 4-element vectors instead of
scalar single precision numbers was more complex than I originally
thought, especially with the real transforms, and I ended up writing
more code than I planned..


The code:
--

Only two files, in good old C, pffft.c and pffft.h . The API is very
very simple, just make sure that you read the comments in pffft.h.


Comparison with other FFTs:
--

The idea was not to break speed records, but to get a decently fast
fft that is at least 50% as fast as the fastest FFT -- especially on
slowest computers . I'm more focused on getting the best performance
on slow cpus (Atom, Intel Core 1, old Athlons, ARM Cortex-A9...), than
on getting top performance on today fastest cpus.

It can be used in a real-time context as the fft functions do not
perform any memory allocation -- that is why they accept a 'work'
array in their arguments.

It is also a bit focused on performing 1D convolutions, that is why it
provides "unordered" FFTs , and a fourier domain convolution
operation.


Benchmark results:
--

MacOS Lion, gcc 4.2, 64-bit, fftw 3.3 on a 3.4 GHz core i7 2600

Built with:

 gcc-4.2 -o test_pffft -arch x86_64 -O3 -Wall -W pffft.c test_pffft.c fftpack.c -L/usr/local/lib -I/usr/local/include/ -DHAVE_VECLIB -framework veclib -DHAVE_FFTW -lfftw3f 

| N (input length) | real FFTPack |   real vDSP  |   real FFTW  |  real PFFFT  | | cplx FFTPack |   cplx vDSP  |   cplx FFTW  |  cplx PFFFT  |
|------------------+--------------+--------------+--------------+--------------| |--------------+--------------+--------------+--------------|
|           64     |      2887    |      9124    |      7082    |      8140    | |      3416    |     15722    |     15705    |     11222    |
|           96     |      3356    |       n/a    |      8520    |      8158    | |      4038    |       n/a    |     16560    |     11141    |
|          128     |      3895    |     12135    |      9613    |     10321    | |      4292    |     17981    |     17054    |     12538    |
|          192     |      3941    |       n/a    |     10220    |     11167    | |      4388    |       n/a    |     16226    |     12964    |
|          256     |      4532    |     13669    |     11031    |     12779    | |      4628    |     19905    |     17259    |     14305    |
|          384     |      3512    |       n/a    |     11013    |     12278    | |      3645    |       n/a    |     16559    |     13370    |
|          512     |      3716    |     15236    |     11515    |     14376    | |      3737    |     20423    |     17050    |     14746    |
|          768     |      3756    |       n/a    |     11524    |     13659    | |      3748    |       n/a    |     16201    |     14891    |
|         1024     |      4060    |     15841    |     10393    |     15732    | |      3828    |     21555    |     15898    |     15883    |
|         2048     |      4646    |     16806    |     11888    |     15752    | |      4323    |     20802    |     15360    |     15219    |
|         4096     |      4794    |     17008    |     11866    |     15785    | |      4167    |     19842    |     14532    |     14723    |
|         8192     |      3887    |     16519    |     11290    |     12854    | |      3738    |     18923    |     12528    |     14164    |
|         9216     |      3924    |       n/a    |     10980    |     13211    | |      3684    |       n/a    |     12349    |     14474    |
|        16384     |      3907    |     16045    |     11146    |     13111    | |      3687    |     17628    |     12364    |     14176    |
|        32768     |      4279    |     15169    |     10946    |     11538    | |      3919    |     15179    |     11558    |     11911    |
|       262144     |      3423    |     11792    |      6753    |      9827    | |      2913    |     11989    |      8406    |     10960    |
|      1048576     |      3313    |     10613    |      5478    |      7142    | |      2683    |      8487    |      2826    |      5961    |
|------------------+--------------+--------------+--------------+--------------| |--------------+--------------+--------------+--------------|



Debian 6, gcc 4.4.5, 64-bit, fftw 3.3.1 on a 3.4 GHz core i7 2600

Built with:
gcc -o test_pffft -DHAVE_FFTW -msse2 -O3 -Wall -W pffft.c test_pffft.c fftpack.c -L$HOME/local/lib -I$HOME/local/include/ -lfftw3f -lm

| N (input length) | real FFTPack |   real FFTW  |  real PFFFT  | | cplx FFTPack |   cplx FFTW  |  cplx PFFFT  |
|------------------+--------------+--------------+--------------| |--------------+--------------+--------------|
|           64     |      3840    |      7680    |      8777    | |      4389    |     20480    |     11171    |
|           96     |      4214    |      9633    |      8429    | |      4816    |     22477    |     11238    |
|          128     |      3584    |     10240    |     10240    | |      5120    |     23893    |     11947    |
|          192     |      4854    |     11095    |     12945    | |      4854    |     22191    |     14121    |
|          256     |      4096    |     11703    |     16384    | |      5120    |     23406    |     13653    |
|          384     |      4395    |     14651    |     12558    | |      4884    |     19535    |     14651    |
|          512     |      5760    |     13166    |     15360    | |      4608    |     23040    |     15360    |
|          768     |      4907    |     14020    |     16357    | |      4461    |     19628    |     14020    |
|         1024     |      5120    |     14629    |     14629    | |      5120    |     20480    |     15754    |
|         2048     |      5632    |     14080    |     18773    | |      4693    |     12516    |     16091    |
|         4096     |      5120    |     13653    |     17554    | |      4726    |      7680    |     14456    |
|         8192     |      4160    |      7396    |     13312    | |      4437    |     14791    |     13312    |
|         9216     |      4210    |      6124    |     13473    | |      4491    |      7282    |     14970    |
|        16384     |      3976    |     11010    |     14313    | |      4210    |     11450    |     13631    |
|        32768     |      4260    |     10224    |     10954    | |      4260    |      6816    |     11797    |
|       262144     |      3736    |      6896    |      9961    | |      2359    |      8965    |      9437    |
|      1048576     |      2796    |      4534    |      6453    | |      1864    |      3078    |      5592    |
|------------------+--------------+--------------+--------------| |--------------+--------------+--------------|


Ubuntu 11.04, gcc 4.5, 32-bit, fftw 3.3 on a 1.6 GHz Atom

Built with:
gcc -o test_pffft -DHAVE_FFTW -msse -mfpmath=sse -O3 -Wall -W pffft.c test_pffft.c fftpack.c -L/usr/local/lib -I/usr/local/include/ -lfftw3f -lm

| N (input length) | real FFTPack |   real FFTW  |  real PFFFT  | | cplx FFTPack |   cplx FFTW  |  cplx PFFFT  |
|------------------+--------------+--------------+--------------| |--------------+--------------+--------------|
|           64     |       452    |      1041    |      1336    | |       549    |      2318    |      1781    |
|           96     |       444    |      1297    |      1297    | |       503    |      2408    |      1686    |
|          128     |       527    |      1525    |      1707    | |       543    |      2655    |      1886    |
|          192     |       498    |      1653    |      1849    | |       539    |      2678    |      1942    |
|          256     |       585    |      1862    |      2156    | |       594    |      2777    |      2244    |
|          384     |       499    |      1870    |      1998    | |       511    |      2586    |      1890    |
|          512     |       562    |      2095    |      2194    | |       542    |      2973    |      2194    |
|          768     |       545    |      2045    |      2133    | |       545    |      2365    |      2133    |
|         1024     |       595    |      2133    |      2438    | |       569    |      2695    |      2179    |
|         2048     |       587    |      2125    |      2347    | |       521    |      2230    |      1707    |
|         4096     |       495    |      1890    |      1834    | |       492    |      1876    |      1672    |
|         8192     |       469    |      1548    |      1729    | |       438    |      1740    |      1664    |
|         9216     |       468    |      1663    |      1663    | |       446    |      1585    |      1531    |
|        16384     |       453    |      1608    |      1767    | |       398    |      1476    |      1664    |
|        32768     |       456    |      1420    |      1503    | |       387    |      1388    |      1345    |
|       262144     |       309    |       385    |       726    | |       262    |       415    |       840    |
|      1048576     |       280    |       351    |       739    | |       261    |       313    |       797    |
|------------------+--------------+--------------+--------------| |--------------+--------------+--------------|

Windows 7, visual c++ 2010 on a 1.6 GHz Atom

Built with:
cl /Ox -D_USE_MATH_DEFINES /arch:SSE test_pffft.c pffft.c fftpack.c

(visual c++ is definitively not very good with SSE intrinsics...)

| N (input length) | real FFTPack |  real PFFFT  | | cplx FFTPack |  cplx PFFFT  |
|------------------+--------------+--------------| |--------------+--------------|
|           64     |       173    |      1009    | |       174    |      1159    |
|           96     |       169    |      1029    | |       188    |      1201    |
|          128     |       195    |      1242    | |       191    |      1275    |
|          192     |       178    |      1312    | |       184    |      1276    |
|          256     |       196    |      1591    | |       186    |      1281    |
|          384     |       172    |      1409    | |       181    |      1281    |
|          512     |       187    |      1640    | |       181    |      1313    |
|          768     |       171    |      1614    | |       176    |      1258    |
|         1024     |       186    |      1812    | |       178    |      1223    |
|         2048     |       190    |      1707    | |       186    |      1099    |
|         4096     |       182    |      1446    | |       177    |       975    |
|         8192     |       175    |      1345    | |       169    |      1034    |
|         9216     |       165    |      1271    | |       168    |      1023    |
|        16384     |       166    |      1396    | |       165    |       949    |
|        32768     |       172    |      1311    | |       161    |       881    |
|       262144     |       136    |       632    | |       134    |       629    |
|      1048576     |       134    |       698    | |       127    |       623    |
|------------------+--------------+--------------| |--------------+--------------|



Ubuntu 11.04 on Pandaboard, gcc-4.5.2, 32-bit, with fftw 3.3.1 beta (neon enabled), on a 1GHz ARM Cortex A9 (TI OMAP4430)

Built with:
gcc-4.5 -O3 -DHAVE_FFTW -march=armv7-a -mtune=cortex-a9 -mfloat-abi=softfp -mfpu=neon -ffast-math test_pffft.c pffft.c -o test_pffft_arm fftpack.c -lm -I/usr/local/include/ -L/usr/local/lib/ -lfftw3f

I must admit that the performance was a bit disappointing here...

| N (input length) | real FFTPack |   real FFTW  |  real PFFFT  | | cplx FFTPack |   cplx FFTW  |  cplx PFFFT  |
|------------------+--------------+--------------+--------------| |--------------+--------------+--------------|
|           64     |       384    |       614    |       591    | |       404    |      1024    |       549    |
|           96     |       324    |       702    |       562    | |       337    |       864    |       503    |
|          128     |       407    |       717    |       640    | |       407    |      1086    |       543    |
|          192     |       404    |       809    |       693    | |       388    |       903    |       547    |
|          256     |       465    |       788    |       788    | |       427    |       871    |       594    |
|          384     |       392    |       814    |       687    | |       343    |       862    |       543    |
|          512     |       411    |       768    |       794    | |       372    |       940    |       583    |
|          768     |       438    |       818    |       767    | |       383    |       846    |       584    |
|         1024     |       427    |       800    |       883    | |       400    |       883    |       602    |
|         2048     |       414    |       853    |       805    | |       343    |       828    |       477    |
|         4096     |       426    |       768    |       698    | |       341    |       808    |       469    |
|         8192     |       332    |       666    |       594    | |       297    |       765    |       438    |
|         9216     |       335    |       660    |       571    | |       294    |       687    |       432    |
|        16384     |       344    |       675    |       606    | |       314    |       709    |       456    |
|        32768     |       342    |       685    |       564    | |       295    |       634    |       399    |
|       262144     |       143    |       301    |       197    | |       160    |       321    |       251    |
|      1048576     |       138    |       238    |       174    | |       173    |       212    |       253    |
|------------------+--------------+--------------+--------------| |--------------+--------------+--------------|
