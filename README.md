
---

# PFFFT: A Pretty Fast FFT

---

<!-- toc -->

- [Brief description](#brief-description)
- [Origin and design](#origin-and-design)
- [The code](#the-code)
- [CMake](#cmake)
- [Using pffft in your CMake project](#using-pffft-in-your-cmake-project)
- [Performance and architecture notes](#performance-and-architecture-notes)
- [Benchmarks](#benchmarks)
- [History and changes](#history-and-changes)
- [Appendix: upstream commit map](#appendix-upstream-commit-map)

<!-- tocstop -->

---

## Brief description

PFFFT does 1D Fast Fourier Transforms, of single and double precision real
and complex vectors. It is fast, correct, and small. The license is BSD-like.

PFFFT was written by Julien Pommier; his original library is hosted on
[Bitbucket](https://bitbucket.org/jpommier/pffft/). This fork adds
double precision, further SIMD backends, a C++ wrapper, CMake support,
fast convolution, and a few DSP functions. Changes made to the original since
the fork have been cherry-picked where applicable, see the
[upstream commit map](#appendix-upstream-commit-map).

It supports these platforms and SIMD instruction sets:

- **x86 and x86-64** (PCs, Intel Macs, servers): SSE for float; SSE2, or AVX
  with `-mavx`, for double. Fused multiply-add with FMA3, e.g. via
  `-march=haswell` or `/arch:AVX2`.
- **64-bit ARM** (Apple Silicon Macs, iPhone and iPad, Android arm64-v8a,
  Raspberry Pi 3 and later with a 64-bit OS, ARM servers): NEON for float and
  double, with fused multiply-add.
- **32-bit ARM** (Android armeabi-v7a, Raspberry Pi 2 and later with a 32-bit
  OS): NEON for float, with fused multiply-add on VFPv4; double runs the scalar
  code.
- **WebAssembly** (browsers, Node.js, WASI): WASM SIMD for float and double,
  optionally with Relaxed SIMD multiply-add.
- **RISC-V** with the vector extension: RVV 1.0 for float and double.
- **PowerPC**: AltiVec for float; double runs the scalar code.
- **ARM Cortex-M microcontrollers**: scalar code; an experimental, opt-in
  Helium backend for float on Cortex-M55 and Cortex-M85.

On any other platform with a C99 compiler, e.g. ESP32 microcontrollers, it runs
portable scalar code. On ESP32, Espressif's
[ESP-DSP](https://github.com/espressif/esp-dsp) is likely faster, as it is
optimized for those chips.

CI runs the tests on Linux and Windows (MSVC and MinGW) on x86-64, on macOS on
64-bit ARM, on WebAssembly, and under emulation on RISC-V and Cortex-M55.

PFFASTCONV does fast convolution (FIR filtering), of single precision
real vectors, utilizing the PFFFT library.

PFDSP contains a few other signal processing functions: mixers, carrier
generation, and a CIC filter.


## Origin and design

PFFFT started from the venerable FFTPACK v4, which is available on NETLIB.
FFTPACK v4 was written in 1985, by Dr Paul Swarztrauber of NCAR. It performs
well, but it is not competitive with the fastest FFT libraries, which take
advantage of SIMD instructions.

The idea was to take this Fortran FFTPACK v4 code, translate it to C,
modify it to use SSE instructions, and check that the
final performance is not completely ridiculous when compared to other
SIMD FFT libraries. Translation to C was performed with [f2c](
http://www.netlib.org/f2c/). The resulting file was a bit edited in
order to remove the thousands of gotos that were introduced by
f2c. The result is in `src/fftpack.h` and `src/fftpack.c`: a complete
translation of [FFTPACK](
http://www.netlib.org/fftpack/), with the discrete cosine transform
and the test program. There is no license information in the netlib
repository, but the FFTPACK v5 curators confirmed to Julien Pommier that
the [same terms do apply to FFTPACK v4](
http://www.cisl.ucar.edu/css/software/fftpack5/ftpk.html). This is a
"BSD-like" license; it is compatible with proprietary projects.

The idea was not to break speed records, but to get a decently fast
FFT that is at least 50% as fast as the fastest FFT -- especially on
slow computers.

PFFFT can be used in a real-time context as the FFT functions do not
perform any memory allocation -- that is why they accept a 'work'
array in their arguments.

PFFFT also supports fast 1D convolution: it provides "unordered" FFTs, which
skip reordering the output (see `pffft_transform()` vs.
`pffft_transform_ordered()`), and a Fourier-domain convolution operation.

[README.pommier.md](README.pommier.md) tells the full story, in his words.


## The code

The FFT API is very very simple, just make sure that you read the comments in
`include/pffft/pffft.h` (float) and `include/pffft/pffft_double.h` (double).
A simple C++ wrapper is available in `include/pffft/pffft.hpp`; see
`examples/` for C and C++ examples.

Transform lengths must be of the form N = 2^a * 3^b * 5^c, and a multiple
of 32 for real or 16 for complex transforms (unless both `PFFFT_USE_SIMD`
and `PFFFT_USE_SCALAR_VECT` are off).
`pffft_nearest_transform_size()` finds the nearest valid length.

The PFFASTCONV API is also very simple, just make sure that you read the comments
in `include/pffft/pffastconv.h`.


## CMake

CMake builds the libraries `PFFFT`, `PFFASTCONV` and `PFDSP` (files
`libpffft`, `libpffastconv` and `libpfdsp`), static unless `BUILD_SHARED_LIBS`
is set. PFFASTCONV and PFDSP need `PFFFT_USE_TYPE_FLOAT`. The FFTPACK
libraries, used for validation, are built only when tests or benchmarks are
enabled and `PFFFT_USE_FFTPACK` is on.
On Debian/Ubuntu Linux, `sudo apt-get install build-essential cmake` is all
you need.

There are several CMake options to modify library size and optimization.
You can explore all available options with `cmake-gui` or `ccmake`.

Some of the options:
* `PFFFT_USE_TYPE_FLOAT` to activate single precision 'float' (default: ON)
* `PFFFT_USE_TYPE_DOUBLE` to activate 'double' precision float (default: ON)
* `PFFFT_USE_SIMD` to use SIMD (SSE/AVX/NEON/ALTIVEC/WASM SIMD/RVV) CPU features (default: ON)
* `PFFFT_USE_SCALAR_VECT` to emulate 4-wide vectors in scalar code when no SIMD is available (default: ON)
* `PFFFT_USE_WASM_RELAXED_SIMD` to use WebAssembly Relaxed SIMD multiply-add instructions, which may be fused when supported by the runtime (default: OFF)
* `PFFFT_USE_HELIUM` to use the experimental ARM Helium (MVE) backend (default: OFF)
* `INSTALL_PFFFT`, `INSTALL_PFFASTCONV`, `INSTALL_PFDSP` to install the
  respective library (default: ON, OFF, OFF)
* `TARGET_C_ARCH` / `TARGET_CXX_ARCH` to set the target architecture: a
  `-march` value for gcc/clang (e.g. `haswell`, `armv7-a`, `rv64gcv`, or
  `native`) or an `/arch` value for MSVC (e.g. `AVX2`)
* `TARGET_C_EXTRA` / `TARGET_CXX_EXTRA` to add a preset of gcc/clang options,
  e.g. `neon_vfpv4` or `neon_rpi4_a72` (32-bit ARM NEON), or `apple_m1`

The SIMD backend follows the compiler's target flags: NEON is used whenever
the compiler enables it (by default on AArch64; on 32-bit ARM e.g. with
`-mfpu=neon`, or one of the `TARGET_C_EXTRA=neon_*` presets), WASM SIMD with
`-msimd128`, and RVV with a `-march` that includes `v`.

Options can be passed to `cmake` at command line, e.g.
```
cmake -DPFFFT_USE_TYPE_FLOAT=OFF -DPFFFT_USE_TYPE_DOUBLE=ON
```

A typical build, here with clang instead of the default compiler:
```
mkdir build
cd build
CC=clang CXX=clang++ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=~ ..
cmake --build .                   # or simply: make
ctest                             # to execute some tests - including benchmarks
cmake --build . --target install  # or simply: [sudo] make install
```

With MSVC on Windows, the build type is chosen at build time:
```
cmake -B build
cmake --build build --config Release
ctest --test-dir build -C Release
```

With [Emscripten](https://emscripten.org/) for WebAssembly (requires the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html)):

```sh
mkdir build
cd build
emcmake cmake ..
cmake --build .
ctest
```

With [WASI-SDK](https://github.com/WebAssembly/wasi-sdk), or any other clang
toolchain targeting `wasm32`/`wasm64`:

```sh
cmake -S . -B build \
  -DCMAKE_TOOLCHAIN_FILE=$WASI_SDK_PATH/share/cmake/wasi-sdk.cmake \
  -DWASI_SDK_PREFIX=$WASI_SDK_PATH \
  -DPFFFT_BUILD_TESTS=OFF -DPFFFT_BUILD_BENCHMARKS=OFF -DPFFFT_BUILD_EXAMPLES=OFF
cmake --build build
```

The tests, benchmarks, and examples don't build out of the box under WASI.
They call `clock()`, which needs `-D_WASI_EMULATED_PROCESS_CLOCKS` in the C
flags and `-lwasi-emulated-process-clocks` in the linker flags. The C++
tests also use exceptions, which WASI-SDK doesn't support by default. With the
clock flags and `-DCMAKE_CROSSCOMPILING_EMULATOR=wasmtime`, the C tests build
and pass under [wasmtime](https://wasmtime.dev/).

WASM SIMD is enabled automatically for any WebAssembly target: Emscripten,
WASI-SDK, or plain clang with `CMAKE_SYSTEM_PROCESSOR` set to `wasm32` or
`wasm64`. See
[Performance and architecture notes](#performance-and-architecture-notes) for
the Relaxed SIMD option and the non-CMake build flags.

## Using pffft in your CMake project

pffft can be included as a dependency in your CMake project via `FetchContent` or `add_subdirectory()`. When used this way, tests, benchmarks, examples, and the standalone `uninstall` helper target are off by default, so pffft won't pollute your build or ctest runner with unnecessary targets or cause target-name collisions.

### FetchContent (CMake 3.14+)

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
| `PFDSP` | DSP utilities (mixer, carrier generation, CIC filter) |

### Build control options

When pffft is the top-level project, these default to `ON`. When included as a subdirectory, they default to `OFF`:

* `PFFFT_BUILD_TESTS` — Build test executables and register CTest tests. These
  include `test_simd_float`/`test_simd_double`, which run the SIMD backend's
  self-test (`test_pffft_float --test-simd`); they pass trivially when SIMD is
  compiled out
* `PFFFT_BUILD_BENCHMARKS` — Build benchmark executables
* `PFFFT_BUILD_EXAMPLES` — Build example programs


## Performance and architecture notes

Compiler- and architecture-specific tips for getting the best performance out
of pffft.

### x86: enable FMA

On x86 the fused multiply-add paths are **opt-in**, and the default build does
not get them. The complex multiply at the heart of `pffft_zconvolve*()` and of
every radix pass uses pffft's multiply-add macros `VMADD`/`VMSUB`. They compile
to a real FMA instruction only when the compiler guarantees FMA3:

| compiler | how to enable |
| --- | --- |
| gcc / clang | `-DTARGET_C_ARCH=haswell` (or a later `-march` value), likewise `TARGET_CXX_ARCH` |
| MSVC | `-DTARGET_C_ARCH=AVX2`, likewise `TARGET_CXX_ARCH` |

Without it, each fused operation costs a separate multiply plus add. Measured
on the `pffft_zconvolve_accumulate()` inner loop, identical for gcc 12 and
clang 14, x86-64 and i686, float (SSE) and double (AVX): **20 floating point
operations per iteration instead of 12**.

```
cmake -DTARGET_C_ARCH=haswell -DTARGET_CXX_ARCH=haswell ..
```

This is not the default on purpose: FMA3 requires an Intel Haswell (2013) or
AMD Piledriver CPU, so enabling it unconditionally would break older hardware.
`cmake` prints a reminder when the selected x86 architecture has no FMA.

Note that FMA changes results by up to one ULP, since it rounds once where the
separate multiply and add round twice.

### ARMv7 (32-bit ARM): prefer clang over gcc

On 32-bit ARM, **gcc emits very inefficient assembly**, with a lot more
operations than clang does. Counts of vector load/store instructions over the
whole float translation unit, both at `-O3 -march=armv7-a -mfpu=neon`:

| | 64-bit `vldr`/`vstr` | 128-bit `vld1`/`vst1` |
| --- | --- | --- |
| gcc 12.2 | 1600 | 344 |
| clang 14 | 181 | 797 |

GCC's bad code generation affects roughly 19 functions, including the hot ones
(`pffft_real_finalize`, `pffft_real_preprocess`, `pffft_cplx_finalize`,
`pffft_zreorder`, `passf4_ps`, `radf4_ps`, ...). This caveat applies to all versions
including GCC 16.1 (newest at time of writing), so this is not about old toolchains.

There is a GCC [ticket](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=43725)
reporting this since 2010 and it hasn't been fixed since. No flag changes it
(fourteen combinations were tried, including `-mtune=cortex-a9/a15`, `-mcpu=...`,
`-mfpu=neon-vfpv4`, `-mvectorize-with-neon-quad`, `-m[no-]unaligned-access`, `-O2`,
and the Raspberry Pi 3/4 presets). It doesn't seem like this problem is ever going
to be fixed, see [this comment](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=43725#c1)
by the ARM maintainer. So we recommend to **use Clang when targeting ARMv7**.

`cmake` prints a reminder when it sees gcc targeting ARMv7. Note that this applies
only to 32-bit ARM, on AArch64 both compilers are fine.

### WebAssembly: SIMD backend and Relaxed SIMD

pffft has a dedicated WASM SIMD backend using native `wasm_simd128.h`
intrinsics for both float and double precision.

> **Performance tip:** enable [Relaxed SIMD](https://github.com/WebAssembly/relaxed-simd) for a measured **6-11% throughput gain** on FFT sizes of 1024 and up, on [runtimes that support it](https://caniuse.com/mdn-webassembly_relaxed-simd) (e.g. V8/Chrome, Node.js):
>
> ```sh
> emcmake cmake -DPFFFT_USE_WASM_RELAXED_SIMD=ON ..
> ```
>
> This uses relaxed multiply-add instructions, which the runtime may fuse into a single-rounding FMA, or fall back to a two-rounding multiply followed by add/subtract. It is **off by default** because relaxed SIMD results are not guaranteed bit-identical across engines or hardware (unlike the default WASM SIMD path); enable it only if your application tolerates ULP-level differences between runtimes.

If you build for WebAssembly without the bundled CMake config (e.g. a raw
`emcc`/`clang`/Makefile/Bazel invocation), compile with `-msimd128` to select
the WASM SIMD backend.

### RISC-V: RVV vector backend

When the target is detected as `riscv64`, the build defaults to
`-march=rv64gcv` (the `v` is the vector extension, V), enabling the
RVV 1.0 vector backend in `src/simd/pf_rvv_{float,double}.h`. On
hardware with `VLEN >= 256` (e.g. Spacemit X100) the configure step
also probes the host VLEN and appends a `zvlNNNb` suffix so the
compiler can emit single 4-wide double-precision vector ops. If the
target turns out to lack V, the build falls back to `rv64gc` and the
scalar code. To opt out explicitly, pass `-DPFFFT_USE_SIMD=OFF` or set
both `-DTARGET_C_ARCH=rv64gc` and `-DTARGET_CXX_ARCH=rv64gc`.

When cross-compiling, the probe runs under `CMAKE_CROSSCOMPILING_EMULATOR`;
with QEMU, enable V and pick the VLEN explicitly, e.g.
`-DCMAKE_CROSSCOMPILING_EMULATOR="qemu-riscv64-static;-L;/usr/riscv64-linux-gnu;-cpu;rv64,v=true,vlen=256,vext_spec=v1.0"`.
Distributions that build their riscv64 libraries for the RVA23 profile (e.g.
Ubuntu 26.04) need an RVA23 CPU model to run the probe, e.g.
`-cpu;rva23u64,vlen=256`. Use `qemu-riscv64` from the `qemu-user` package
there, which no longer ships a `-static` binary. Without an emulator the build assumes V is present.

### PowerPC: AltiVec backend

The AltiVec backend (`src/simd/pf_altivec_float.h`) is **float only**: there
is no AltiVec double-precision backend, so double precision runs the scalar
code on PowerPC. It works on big- and little-endian targets and doesn't need
VSX. It is selected when the compiler defines `__ALTIVEC__`, which is the
default on little-endian `ppc64le`; elsewhere pass e.g.
`-DCMAKE_C_FLAGS=-maltivec`. GCC and clang warn that `vec_lvsl` is deprecated
on little-endian; this is harmless: the unaligned-load idiom it is used for
is correct there.

To check that the backend is active, look for `float: ALTIVEC` in CMake's
"Detected SIMD backend" line, and run `test_pffft_float --test-simd` (also
registered as the `test_simd_float` ctest test).

### Arm Cortex-M: experimental Helium backend

On Cortex-M cores pffft runs its portable scalar code by default. Armv8.1-M
cores with Helium (MVE), such as Cortex-M55 and Cortex-M85, can use the
**experimental** float backend in `src/simd/pf_helium_float.h` with
`-DPFFFT_USE_HELIUM=ON` (outside CMake: `-DPFFFT_ENABLE_HELIUM` with
`-mcpu=cortex-m55` or similar). It is tested for correctness under QEMU, but
has not been benchmarked on hardware, so it may well be slower than the scalar
code; measurements are welcome. Double precision stays scalar, as MVE has no
double-precision lanes.

The tests can be run on QEMU's Cortex-M55 board with `arm-none-eabi-gcc`
(with newlib) and `qemu-system-arm`:

```
cmake -B build_m55 -DCMAKE_TOOLCHAIN_FILE=tests/cortex-m55/toolchain.cmake -DPFFFT_USE_HELIUM=ON
cmake --build build_m55
ctest --test-dir build_m55 -E "^bench_"
```

The `bench_*` tests are excluded because they allocate about 1 GiB.

## Benchmarks

Benchmark results, with plots, are collected at [https://github.com/hayguen/pffft_benchmarks](https://github.com/hayguen/pffft_benchmarks).
They are stored in a separate git-repository to keep this repository's sources small.
You might contribute by providing us the results of your computer(s).

### Git submodules for Green-, Kiss-, and Pocket-FFT
The benchmarks compare against other FFT libraries. Sources for [Green-](https://github.com/hayguen/greenffts),
[Kiss-](https://github.com/hayguen/kissfft)
and [Pocket-FFT](https://github.com/hayguen/pocketfft)
can be downloaded directly with the sources of this repository - using git submodules:
```
git clone --recursive https://github.com/marton78/pffft.git
```

`--recursive` fetches the submodules too.
You can also retrieve them later:
```
git submodule update --init
```

### Fastest Fourier Transform in the West: FFTW
To compare with [FFTW](http://www.fftw.org/), install it
(on Debian/Ubuntu: `sudo apt-get install libfftw3-dev`) and use cmake option
`-DPFFFT_USE_BENCH_FFTW=ON` in the benchmark build below.

### Intel Math Kernel Library: MKL
To compare with Intel's [oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html),
install it (on Ubuntu-Linux the package `intel-mkl`, on Debian `intel-mkl-full`)
and use cmake option `-DPFFFT_USE_BENCH_MKL=ON`.

### Performing the benchmarks - with CMake
The main build's benchmark executables accept the same `PFFFT_USE_BENCH_*`
options. The separate project in `bench/` builds and runs the float and double
benchmarks and packs the results, with system information, into a `.tgz` for
the results repository. Use it in its own build folder:
```
mkdir build_benches
cd build_benches
cmake ../bench
```

Several CMake options select which FFT implementations are benchmarked.
You can explore all available options with `cmake-gui` or `ccmake`, see [CMake](#cmake).

Some of the options:
* `BENCH_ID`         name the benchmark - used in filename
* `BENCH_ARCH`       target architecture passed to compiler for code optimization
* `PFFFT_USE_BENCH_FFTW`   use (system-installed) FFTW3 in FFT benchmark (default: OFF)
* `PFFFT_USE_BENCH_GREEN`  use Green FFT in FFT benchmark (default: ON)
* `PFFFT_USE_BENCH_KISS`   use KissFFT in FFT benchmark (default: ON)
* `PFFFT_USE_BENCH_POCKET` use PocketFFT in FFT benchmark (default: ON)
* `PFFFT_USE_BENCH_MKL`    use Intel MKL in FFT benchmark (default: OFF)

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
CC=/usr/bin/clang-11 CXX=/usr/bin/clang++-11 cmake -DBENCH_ID=clang11 -DBENCH_ARCH=native -DPFFFT_USE_BENCH_FFTW=ON -DPFFFT_USE_BENCH_MKL=ON ../bench
```

On Windows, call cmake from the VS Developer prompt.

For running with different compiler version(s):
* copy the result file (.tgz), e.g. `cp *.tgz ../`
* empty the build directory: `rm -rf *`
* then continue with the cmake step

## History and changes
This is a fork of Julien Pommier's pffft on
[Bitbucket](https://bitbucket.org/jpommier/pffft/). His README, as of upstream
commit `0979688`, is kept in [README.pommier.md](README.pommier.md).

Git history shows the following first commits of the major contributors:
* Julien Pommier: November 19, 2011
* Marton Danoczy: September 30, 2015
* Hayati Ayguen: December 22, 2019
* Dario Mambro: March 24, 2020
* George Hilliard: March 15, 2025
* Gunj Joshi: April 1, 2025
* Aleksey Vaneev: December 20, 2025
* Sergi Granell: February 26, 2026
* Michael Neuling: May 26, 2026
* Piotr Kubaj: September 29, 2026

There are a few other contributors not listed here.

The main changes on top of the original include:
* improved benchmarking, see [https://github.com/hayguen/pffft_benchmarks](https://github.com/hayguen/pffft_benchmarks)
* double support
* AVX support for double
* NEON double, a dedicated WASM SIMD, RISC-V RVV, and ARM Helium backends
* C++ headers (wrapper)
* additional API helper functions
* additional library for fast convolution
* const-correct setup arguments for all transform functions
* zero-phase convolution support for linear-phase FIR filters
  (`pffft_zconvert_zp()` / `pffft_zconvolve_zp()`, see `include/pffft/pffft.h`);
  with default scaling the forward transform -> zconvolve_zp ->
  backward pipeline yields the circular convolution at unit gain
* unscaled frequency-domain multiply `pffft_zconvolve()` and `pffft_zconvolve_scale()`
* CMake support
* CTest

The original repository has also continued to evolve independently, receiving
portability fixes and warning cleanups from various contributors. Those changes
have been cherry-picked into this fork where applicable. See the
[Appendix: upstream commit map](#appendix-upstream-commit-map) for a full
commit-by-commit comparison.


## Appendix: upstream commit map

The table below lists every non-merge commit on
[jpommier/pffft master](https://bitbucket.org/jpommier/pffft/) after
the last common ancestor (`16490be`, _"added numbers for the tegra K1"_,
2014-09-03). The *Ours* column gives the commit in this fork that ports the
change or made it unnecessary; `--` means there is no such commit, and *Notes*
says why.

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
