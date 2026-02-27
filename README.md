
---

# PFFFT: a pretty fast FFT and fast convolution with PFFASTCONV

---

<!-- toc -->

- [Brief Description](#brief-description)
- [Why does it exist?](#why-does-it-exist)
- [CMake](#cmake)
- [Using pffft in your CMake project](#using-pffft-in-your-cmake-project)
- [History / Origin / Changes](#history--origin--changes)
- [Comparison with other FFTs](#comparison-with-other-ffts)
- [Dependencies / Required Linux packages](#dependencies--required-linux-packages)
- [Benchmarks and results](#benchmarks-and-results)
- [Appendix: upstream commit map](#appendix-upstream-commit-map)

<!-- tocstop -->

---

## Brief description:

PFFFT does 1D Fast Fourier Transforms, of single precision real and
complex vectors. It tries do it fast, it tries to be correct, and it
tries to be small. Computations do take advantage of SSE1/AVX/AVX2 instructions
on x86 cpus, Altivec on powerpc cpus, and NEON on ARM cpus
(including Apple Silicon). The license is BSD-like.

PFFFT is a fork of [Julien Pommier's library on bitbucket](https://bitbucket.org/jpommier/pffft/)
with some changes and additions.


PFFASTCONV does fast convolution (FIR filtering), of single precision 
real vectors, utilizing the PFFFT library. The license is BSD-like.

PFDSP contains a few other signal processing functions.
Currently, mixing and carrier generation functions are contained.
It is work in progress - also the API!
The fast convolution from PFFASTCONV might get merged into PFDSP.


## Why does it exist:

I (Julien Pommier) was in search of a good performing FFT library ,
preferably very small and with a very liberal license.

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
library, see for example the 1d single precision benchmarks
[here](http://www.fftw.org/speed/opteron-2.2GHz-32bit/). It is however not
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
SIMD FFT libraries. Translation to C was performed with [f2c](
http://www.netlib.org/f2c/). The resulting file was a bit edited in
order to remove the thousands of gotos that were introduced by
f2c. You will find the fftpack.h and fftpack.c sources in the
repository, this a complete translation of [fftpack](
http://www.netlib.org/fftpack/), with the discrete cosine transform
and the test program. There is no license information in the netlib
repository, but it was confirmed to me by the fftpack v5 curators that
the [same terms do apply to fftpack v4]
(http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html). This is a
"BSD-like" license, it is compatible with proprietary projects.

Adapting fftpack to deal with the SIMD 4-element vectors instead of
scalar single precision numbers was more complex than I originally
thought, especially with the real transforms, and I ended up writing
more code than I planned..


## The code:

### Good old C:
The FFT API is very very simple, just make sure that you read the comments in `pffft.h`.

The Fast convolution's API is also very simple, just make sure that you read the comments 
in `pffastconv.h`.

### C++:
A simple C++ wrapper is available in `pffft.hpp`.

### Git:
This archive's source can be downloaded with git (without the submodules):
```
git clone https://github.com/marton78/pffft.git
```

### Only two files?:
_"Only two files, in good old C, pffft.c and pffft.h"_

This statement does **NO LONGER** hold!

With new functionality and support for AVX, there was need to restructure the sources.
But you can compile and link **pffft** as a static library.


## CMake:
There's now CMake support to build the static libraries `libPFFFT.a`
and `libPFFASTCONV.a` from the source files. The additional
`libFFTPACK.a` library is built only when tests or benchmarks are enabled.

There are several CMake options to modify library size and optimization.
You can explore all available options with `cmake-gui` or `ccmake`,
the console version - after having installed (on Debian/Ubuntu Linux) one of
```
sudo apt-get install cmake-qt-gui
sudo apt-get install cmake-curses-gui
```

Some of the options:
* `PFFFT_USE_TYPE_FLOAT` to activate single precision 'float' (default: ON)
* `PFFFT_USE_TYPE_DOUBLE` to activate 'double' precision float (default: ON)
* `PFFFT_USE_SIMD` to use SIMD (SSE/AVX/NEON/ALTIVEC) CPU features? (default: ON)
* `DISABLE_SIMD_AVX` to disable AVX CPU features (default: OFF)
* `PFFFT_USE_SIMD_NEON` to force using NEON on ARM (requires PFFFT_USE_SIMD) (default: OFF)
* `PFFFT_USE_SCALAR_VECT` to use 4-element vector scalar operations (if no other SIMD) (default: ON)

Options can be passed to `cmake` at command line, e.g.
```
cmake -DPFFFT_USE_TYPE_FLOAT=OFF -DPFFFT_USE_TYPE_DOUBLE=ON
```

My Linux distribution defaults to GCC. With installed CLANG and the bash shell, you can use it with
```
mkdir build
cd build
CC=/usr/bin/clang CXX=/usr/bin/clang++ cmake -DCMAKE_BUILD_TYPE=Debug ../
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=~ ../
ccmake .                          # or: cmake-gui .
cmake --build .                   # or simply: make
ctest                             # to execute some tests - including benchmarks
cmake --build . --target install  # or simply: [sudo] make install
```

With MSVC on Windows, you need some different options. Following ones to build a 64-bit Release with Visual Studio 2019:
```
mkdir build
cd build
cmake -G "Visual Studio 16 2019" -A x64 ..
cmake --build . --config Release
ctest -C Release
```

see [https://cmake.org/cmake/help/v3.15/manual/cmake-generators.7.html#visual-studio-generators](https://cmake.org/cmake/help/v3.15/manual/cmake-generators.7.html#visual-studio-generators)


## Using pffft in your CMake project

pffft can be included as a dependency in your CMake project via `FetchContent` or `add_subdirectory()`. When used this way, tests, benchmarks, and examples are automatically disabled, so pffft won't pollute your build or ctest runner with unnecessary targets.

### FetchContent (CMake 3.11+)

```cmake
include(FetchContent)

FetchContent_Declare(pffft
    GIT_REPOSITORY https://github.com/marton78/pffft.git
    GIT_TAG master
)
FetchContent_MakeAvailable(pffft)

target_link_libraries(my_app PRIVATE PFFFT::PFFFT)
```

### add_subdirectory

If you have pffft checked out as a subdirectory or git submodule:

```cmake
add_subdirectory(pffft)

target_link_libraries(my_app PRIVATE PFFFT::PFFFT)
```

### Available targets

| Target | Description |
|---|---|
| `PFFFT::PFFFT` (or `PFFFT`) | FFT library (float and/or double, depending on options) |
| `PFFASTCONV` | Fast convolution library |
| `PFDSP` | DSP utilities (mixer, carrier generation) |

### Build control options

When pffft is the top-level project, these default to `ON`. When included as a subdirectory, they default to `OFF`:

* `PFFFT_BUILD_TESTS` — Build test executables and register CTest tests
* `PFFFT_BUILD_BENCHMARKS` — Build benchmark executables
* `PFFFT_BUILD_EXAMPLES` — Build example programs


## History / Origin / Changes:
Origin for this code/fork is Julien Pommier's pffft on bitbucket:
[https://bitbucket.org/jpommier/pffft/](https://bitbucket.org/jpommier/pffft/)

Git history shows following first commits of the major contributors:
* Julien Pommier: November 19, 2011
* Marton Danoczy: September 30, 2015
* Hayati Ayguen: December 22, 2019
* Dario Mambro: March 24, 2020
* George Hilliard: March 15, 2025
* Gunj Joshi: April 1, 2025
* Aleksey Vaneev: December 20, 2025
* Sergi Granell (xerpi): February 26, 2026

There are a few other contributors not listed here.

The main changes on top of the original include:
* improved benchmarking, see [https://github.com/hayguen/pffft_benchmarks](https://github.com/hayguen/pffft_benchmarks)
* double support
* avx(2) support
* c++ headers (wrapper)
* additional API helper functions
* additional library for fast convolution
* cmake support
* ctest

The original repository has also continued to evolve independently, receiving
portability fixes and warning cleanups from various contributors. Those changes
have been cherry-picked into this fork where applicable. See the
[Appendix: upstream commit map](#appendix-upstream-commit-map) for a full
commit-by-commit comparison.


## Comparison with other FFTs:
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

Very interesting is [https://www.nayuki.io/page/free-small-fft-in-multiple-languages](https://www.nayuki.io/page/free-small-fft-in-multiple-languages).
It shows how small an FFT can be - including the Bluestein algorithm, but it's everything else than fast.
The whole C++ implementation file is 161 lines, including the Copyright header, see
[https://github.com/nayuki/Nayuki-web-published-code/blob/master/free-small-fft-in-multiple-languages/FftComplex.cpp](https://github.com/nayuki/Nayuki-web-published-code/blob/master/free-small-fft-in-multiple-languages/FftComplex.cpp)

## Dependencies / Required Linux packages

On Debian/Ubuntu Linux following packages should be installed:

```
sudo apt-get install build-essential gcc g++ cmake
```


## Benchmarks and results

#### Quicklink
Find results at [https://github.com/hayguen/pffft_benchmarks](https://github.com/hayguen/pffft_benchmarks).

#### General
My (Hayati Ayguen) first look at FFT-benchmarks was with [benchFFT](http://www.fftw.org/benchfft/)
and especially the results of the benchmarks [results](http://www.fftw.org/speed/),
which demonstrate the performance of the [FFTW](http://www.fftw.org/).
Looking at the benchmarked computer systems from todays view (2021), these are quite outdated.

Having a look into the [benchFFT source code](http://www.fftw.org/benchfft/benchfft-3.1.tar.gz),
the latest source changes, including competitive fft implementations, are dated November 2003.

In 2019, when pffft got my attention at [bitbucket](https://bitbucket.org/jpommier/pffft/src/master/),
there were also some benchmark results.
Unfortunately the results are tables with numbers - without graphical plots.
Without the plots, i could not get an impression. That was, why i started
[https://github.com/hayguen/pffft_benchmarks](https://github.com/hayguen/pffft_benchmarks),
which includes GnuPlot figures.

Today in June 2021, i realized the existence of [https://github.com/FFTW/benchfft](https://github.com/FFTW/benchfft).
This repository is much more up-to-date with a commit in December 2020.
Unfortunately, it looks not so simple to get it run - including the generation of plots.

Is there any website showing benchFFT results of more recent computer systems?

Of course, it's very important, that a benchmark can be compared with a bunch
of different FFT algorithms/implementations.
This requires to have these compiled/built and utilizable.


#### Git submodules for Green-, Kiss- and Pocket-FFT
Sources for [Green-](https://github.com/hayguen/greenffts),
[Kiss-](https://github.com/hayguen/kissfft)
and [Pocket-FFT](https://github.com/hayguen/pocketfft)
can be downloaded directly with the sources of this repository - using git submodules:
```
git clone --recursive https://github.com/marton78/pffft.git
```

Important is `--recursive`, that does also fetch the submodules directly.
But you might retrieve the submodules later, too:
```
git submodule update --init
```

#### Fastest Fourier Transform in the West: FFTW
To allow comparison with FFTW [http://www.fftw.org/](http://www.fftw.org/),
cmake option `-DPFFFT_USE_BENCH_FFTW=ON` has to be used with following commands.
The cmake option requires previous setup of following (debian/ubuntu) package:
```
sudo apt-get install libfftw3-dev
```

#### Intel Math Kernel Library: MKL
Intel's MKL [https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html)
currently looks even faster than FFTW.

On Ubuntu-Linux it's easy to setup with the package `intel-mkl`.
Similar on Debian: `intel-mkl-full`.

There are special repositories for following Linux distributions:
* Debian/apt: [https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html](https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html)
* RedHat/yum: [https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-yum-repo.html](https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-yum-repo.html)
* Gentoo/ebuild: [https://packages.gentoo.org/packages/sci-libs/mkl](https://packages.gentoo.org/packages/sci-libs/mkl)

#### Performing the benchmarks - with CMake
Benchmarks should be prepared by creating a special build folder
```
mkdir build_benches
cd build_benches
cmake ../bench
```

There are several CMake options to parametrize, which fft implementations should be benched.
You can explore all available options with `cmake-gui` or `ccmake`, see [CMake](#cmake).

Some of the options:
* `BENCH_ID`         name the benchmark - used in filename
* `BENCH_ARCH`       target architecture passed to compiler for code optimization
* `PFFFT_USE_BENCH_FFTW`   use (system-installed) FFTW3 in fft benchmark? (default: OFF)
* `PFFFT_USE_BENCH_GREEN`  use Green FFT in fft benchmark? (default: ON)
* `PFFFT_USE_BENCH_KISS`   use KissFFT in fft benchmark? (default: ON)
* `PFFFT_USE_BENCH_POCKET` use PocketFFT in fft benchmark? (default: ON)
* `PFFFT_USE_BENCH_MKL`    use Intel MKL in fft benchmark?  (default: OFF)

These options can be passed to `cmake` at command line, e.g.
```
cmake -DBENCH_ARCH=native -DPFFFT_USE_BENCH_FFTW=ON -DPFFFT_USE_BENCH_MKL=ON ../bench
```

The benchmarks are built and executed with
```
cmake --build .
```

You can also specify to use a different compiler/version with the cmake step, e.g.:

```
CC=/usr/bin/gcc-9 CXX=/usr/bin/g++-9 cmake -DBENCH_ID=gcc9 -DBENCH_ARCH=native -DPFFFT_USE_BENCH_FFTW=ON -DPFFFT_USE_BENCH_MKL=ON ../bench
```

```
CC=/usr/bin/clang-11 CXX=/usr/bin/clang++-11 cmake -DBENCH_ID=clang11 -DBENCH_ARCH=native -DPFFFT_USE_BENCH_FFTW=ON -DPFFFT_USE_BENCH_MKL=ON ../bench
```

For using MSVC/Windows, the cmake command requires/needs the generator and architecture options and to be called from the VS Developer prompt:
```
cmake -G "Visual Studio 16 2019" -A x64 ../bench/
```

see [https://cmake.org/cmake/help/v3.15/manual/cmake-generators.7.html#visual-studio-generators](https://cmake.org/cmake/help/v3.15/manual/cmake-generators.7.html#visual-studio-generators)



For running with different compiler version(s):
* copy the result file (.tgz), e.g. `cp *.tgz ../`
* delete the build directory: `rm -rf *`
* then continue with the cmake step


#### Benchmark results and contribution
You might contribute by providing us the results of your computer(s).

The benchmark results are stored in a separate git-repository:
See [https://github.com/hayguen/pffft_benchmarks](https://github.com/hayguen/pffft_benchmarks).

This is to keep this repositories' sources small.


## Appendix: upstream commit map

The table below lists every non-merge commit on
[jpommier/pffft master](https://bitbucket.org/jpommier/pffft/) starting from
the last common ancestor (`16490be`, _"added numbers for the tegra K1"_,
2014-09-03). For each commit the corresponding hash in this fork is given, or
an explanation of why the change was not ported.

| Date | Upstream | Author | Change | Ours | Notes |
|------|----------|--------|--------|------|-------|
| 2016-09-22 | `74d7261` | Julien Pommier | Add 64-bit ARM support | `3673ac0` | |
| 2020-10-15 | `e7de7c2` | Gregor Jasny | Fix 32-bit x86 iOS simulator define | `66ba262` | |
| 2021-05-07 | `82f5a59` | Julien Pommier | Add MKL to benchmark, update results | `e44faf5` | benchmark only |
| 2021-05-07 | `57b818a` | Julien Pommier | Remove old README.txt | -- | already absent |
| 2021-05-11 | `ed78751` | Julien Pommier | Fix compilation when SIMD not available | -- | handled by restructured code |
| 2021-12-01 | `ccd5628` | Mike Hommey | Add missing `#include <altivec.h>` | `b5f1aaa` | cherry-picked |
| 2021-12-02 | `233c266` | Dan Horák | Fix AltiVec vector initialization syntax | `ea9ccdf` | cherry-picked |
| 2021-12-02 | `d3fd82b` | Dan Horák | Fix ppc/ppc64 detection for GCC | `dbaaf40` | cherry-picked |
| 2022-02-15 | `7914df2` | Yair Chuchem | Fix Xcode 13.2.1 warnings | `617dbd9` | cherry-picked (comma-to-semicolons portion; `(void)` part was already done) |
| 2023-03-26 | `045c091` | Yair Chuchem | Fix `-Wshadow` in `pffft_zreorder` | `843d253` | cherry-picked |
| 2023-08-17 | `7ed8b4c` | Julien Pommier | Shorter test time on ARM | -- | test-only; our test infrastructure differs |
| 2023-08-17 | `1adf4b9` | Julien Pommier | Remove `f` suffixes for double-precision safety | `9c0c63c` | cherry-picked |
| 2023-10-11 | `f2455bf` | Nick Dowell | Fix clang `-Wstrict-prototypes` | `f41361d` | already done (void signatures) |
| 2024-03-19 | `17b6d50` | Yair Chuchem | Extract duplicated twiddle code | `2b5b658` | cherry-picked |
| 2024-04-08 | `180c2d6` | Yair Chuchem | Fix `-Wshadow` in `zconvolve_accumulate` | `c8259be` | cherry-picked |
| 2024-11-20 | `c43a62c` | Chris Robinson | Check AltiVec is enabled before using it | `b04ffaa` | cherry-picked |
| 2024-11-25 | `3e10d3b` | Fabian Greffrath | Better SIMD capability checks | `81e6fc3` | cherry-picked (SSE portion; NEON uses CMake option) |
| 2024-11-26 | `8873dd4` | Fabian Greffrath | Define `_USE_MATH_DEFINES` | `1e139d9` | mitigated differently (constants defined directly) |
| 2024-11-29 | `a9fdee9` | Julien Pommier | Add arm64ec support on Windows | `06c8002` | cherry-picked (adapted for split SIMD files) |
| 2024-11-29 | `c0ef461` | Julien Pommier | Replace `long long` with `size_t` | `20a6015` | already uses `uintptr_t` |
| 2024-11-29 | `a8a8013` | Julien Pommier | Enable WASM SIMD support | `ec031bd` | cherry-picked |
| 2024-11-29 | `de2c855` | Julien Pommier | Whitespace cleanup | -- | cosmetic only |
| 2024-11-29 | `6d6a123` | Julien Pommier | `ifac` overflow protection, N validation | `d16aab2` | cherry-picked (library portions) |
| 2024-11-29 | `02fe771` | Julien Pommier | Guard `_USE_MATH_DEFINES` redefinition | -- | not needed (constants defined directly) |
| 2025-02-12 | `0d7449a` | Dan Raviv | Fix MSVC `/fp:strict` C2099 errors | `4d1c78d` | cherry-picked |
| 2025-12-19 | `c306b13` | Julien Pommier | Fix implicit double-to-float conversions | `31be131` | cherry-picked (fftpack.c portion) |
| 2026-01-05 | `0979688` | Julien Pommier | Fix alignment for small `size_t` platforms | `a9786ad` | already uses `uintptr_t` |

