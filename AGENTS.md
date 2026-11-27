# PFFFT Agent Guidelines

This guide provides build commands, testing procedures, and code style conventions for the PFFFT project.

## Build Commands

```bash
# Configure with CMake (Release build, default)
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..

# Configure with specific compiler (clang example)
CC=/usr/bin/clang CXX=/usr/bin/clang++ cmake -DCMAKE_BUILD_TYPE=Debug ..

# Configure with specific options
cmake -DPFFFT_USE_TYPE_FLOAT=ON -DPFFFT_USE_TYPE_DOUBLE=ON -DPFFFT_USE_SIMD=ON ..

# Build all targets
cmake --build . --parallel
# or
make -j$(nproc)

# Build specific target
cmake --build . --target PFFFT
make PFFFT
```

## Test Commands

```bash
# Run all tests via CTest
ctest

# Run all tests with verbose output
ctest -V

# Run a specific test by name
ctest -R test_fft_factors
ctest -R test_pffft_float
ctest -R test_pffft_double
ctest -R test_pfconv_lens_symetric
ctest -R test_pf_conv                  # the kernels of bench_pf_conv_float

# Run specific test executable directly (after build)
./test_fft_factors
./test_pffft_float
./test_pffft_double
./test_pffastconv --no-bench --quick

# Run the FFT benchmark (see benchmarks/README.md)
ctest -R bench_fft                                      # its smoke tests
./benchmarks/bench_fft --sizes 1024,1152 --prec f32     # pffft against the other libraries
./benchmarks/bench_pf_mixer_float
./benchmarks/bench_pf_conv_float

# Measure a change against master (from the repository root)
python3 benchmarks/compare.py master .

# Check its accuracy against master (see "Accuracy checks" below)
python3 benchmarks/accuracy.py master . --bits    # must give master's bits
python3 benchmarks/accuracy.py master .           # may differ, must not be less accurate

# Debug with AddressSanitizer (must rebuild)
cmake -DPFFFT_USE_DEBUG_ASAN=ON ..
cmake --build .
ctest
```

## Accuracy checks

`benchmarks/accuracy.py REV_A REV_B` builds both revisions and measures their
error against an exact (double-double) FFT on the same random inputs, forward
and backward, and checks bit for bit that split and unsplit builds, in-place
and out-of-place, ordered and transform + reorder agree; docs/accuracy.md
describes it. Accuracy depends on the instruction set and the compiler, not
on the CPU model, so unlike timing it needs neither the phones nor the GitHub
runners for most changes. Which check a change needs:

- **Every change:** `ctest` runs `test_accuracy` and `test_accuracy_split`
  (under a second): the bitwise equivalences and absolute error limits of the
  build at hand. CI runs them on every backend it builds, emulators included.
- **Changes that move data but compute nothing new** (reordering, loads and
  stores, prefetch, loop structure, unrolling, thresholds that choose between
  paths that give the same bits): `accuracy.py master . --bits`. B must give
  master's bits in every case. No exact transform is computed: seconds.
- **Changes to the arithmetic** (FMA placement, another radix or factor
  order, twiddles, scaling): while iterating, `accuracy.py master . --sizes
  LIST --max 65536` with the lengths the change touches (seconds); before the
  PR, the full default run (all lengths to 2^22, a few minutes on an M2),
  once per instruction set whose code changed: NEON natively on an Apple
  Mac, SSE/AVX/FMA on the same Mac under Rosetta (configure with
  `-DCMAKE_OSX_ARCHITECTURES=x86_64` and run with `ROSETTA_ADVERTISE_AVX=1`;
  accuracy needs no timing). Add one gcc run (an x86 or Arm Linux runner),
  as gcc contracts `a*b+c` into FMA by default where clang and MSVC do not.
  The phones run the same NEON code with the same compiler: not needed.
- **Twiddles or anything in the setup:** the full run in both precisions;
  they affect every length.

Judge the forward and backward errors only. The round trip is shown, not
judged: the backward pass uses the forward pass's twiddles conjugated, so
their errors largely cancel in it. Do not compare double precision against
another FFT library in double (PocketFFT, FFTW): their error is the size of
pffft's and correlated with it, which rewards resembling them rather than
being accurate.

## Thresholds that depend on the CPU

Some choices pay on one CPU and cost on another: from which array size a
transform is split into blocks that fit the cache, from which size the
passes prefetch what they write (x86), and up to which size a real
transform ends with a radix-6 pass. They are not compile-time constants:
`pffft_new_setup()` takes them from `pffft_default_tuning()`, which
identifies the CPU (x86: `cpuid`; Arm: the MIDR of the fastest core; Apple:
`sysctl`) and takes the first matching row of `src/pffft_tuning.h`
(`pffft_tuning_t` in `pffft_types.h`; `pffft_cpu_describe()` says which row). A
program can pass its own to `pffft_new_setup_tuned()`. Only one build option
is left: `PFFFT_DISABLE_PREFETCH`, which leaves the prefetching code out.

- **Measure thresholds in one build.** Every value of a threshold is a
  `--tune` of the same library (`bench_fft --tune LABEL=SPEC`,
  `compare_many.py base=HEAD x=HEAD --tune x=split=8m`), so that all values
  run the same code side by side; never one build per value, which mixes
  the threshold with code layout. `benchmarks/tune.py` sweeps all of them on
  one machine and prints the table row they make.
- **Values that were not measured:** a transform depends on the split
  threshold only through its split depth, so `benchmarks/split_sim.py`
  gives any split threshold from a sweep of others, exactly, without a new
  run (`--candidate NAME=SPEC`, `--cases REF` for the lengths that change).
- **Choosing a default:** every machine weighs the same; change a row only
  where a value beats the current one clearly (tune.py: by more than half a
  point of the geometric mean), and then compare it with master at the
  docs' lengths (the powers of two and mixed ones such as 96000, 192000,
  209952, 1500000) before deciding: a sweep's lengths (1, 3, 5, 9, 15 times
  2^k) once hid losses of 5 to 15% at mixed-radix lengths. Say in the table's
  comment which machines a row rests on, with the numbers.
- **Bits:** the split and prefetch thresholds choose between ways that give
  the same results (`accuracy.py master . --bits`); `radix6_max` changes the
  arithmetic, so a change to it needs the full accuracy run.
- **Code layout:** the same code at another address can move a short
  length by 5 to 16% on some machines (EPYC 7763 with either compiler,
  Neoverse N2, Core i7-8700B, EPYC 9V45; under 2% on the Xeons 6975P-C and
  8573C and on Apple). Call a loss layout only with the matched-layout
  control: `compare_many.py master . --pads spread` on that machine, next to
  an A/A run (`master m2=master --pads spread`) for the noise. A loss is
  real if it is `slower*` (every pad, beyond both bands), beyond the A/A
  noise of its case, and the kernels were sampled (the run's "layout:"
  lines); layout if the verdict is `layout`. On x86, judge with MSVC or
  Clang: GCC's own alignment leaves most kernels where they are at every
  pad. See docs/code-layout.md.

## Benchmarking on GitHub's runners

`.github/workflows/bench.yml` measures CPUs we do not own. It runs only when
started by hand, from the default branch's copy of the workflow, with
`--ref` picking the code (`benchmarks/README.md`, "Other CPUs"). Each label
lands on one of several CPUs:

| Label | CPUs (how often) |
|---|---|
| ubuntu-26.04 | Xeon 6975P-C (about 1 run in 30), Xeon Platinum 8573C, EPYC 9V45, 9V74, 7763 |
| windows-2025 | EPYC 7763, 9V74, Xeon Platinum 8573C, 8370C, 6973P-C (MSVC) |
| macos-26-intel | Core i7-8700B |
| ubuntu-26.04-arm | Neoverse N2 |
| windows-11-arm | Cobalt 100 (MSVC) |

Only the Linux runners count cycles (perf); the others are timed in ns, and
the VMs are shared, so expect noise and repeat a loss before working on it.

- **Use every landing once.** `benchmarks/ci_dispatch.py TABLE.json OUT`
  dispatches the workflow until every CPU of the table has done its entry:
  the sweep (`run.py --fftw`, for the docs), compares (`compare.py`, e.g.
  `origin/master HEAD` and `c14802c HEAD` for 1.1.0) and "rounds" (more
  `compare_many.py` runs). A run that lands on a CPU without an entry ends
  in seconds; a later run of an entry does only its missing parts. Plan
  everything a CPU must measure for the round before launching: a rare CPU
  may not come back for hours.
- **Deciding and measuring in one run:** if a run is to choose a default
  and also give the PR's numbers, put the candidates as `--tune` variants
  of one library next to `origin/master` in a round. Changing the table
  afterwards is then a change of data, not of the measured code. A default
  changed after the final compares otherwise needs those compares again.
- **On the runner,** `master` is `origin/master`.
- **Scratch branches:** push the code to measure as `scratch/<topic>`, with
  a commit on top for anything only the measurement needs, and delete the
  branch when the data is in. Keep the CSVs and the dispatch table
  (`build_bench_cache/` is ignored by git: back it up).
- **Local devices** (`--target android[:serial]`, `ios[:udid]`) run the same
  scripts; an iPhone must be unlocked to install. A 32-bit Raspberry Pi
  needs `--cmake="-DTARGET_C_ARCH=armv7-a -DTARGET_C_EXTRA=neon_vfpv4"` for
  NEON (tune.py) and takes hours per sweep.

## Code Style Guidelines

### Language Standards

- **C**: C99 standard (set in CMakeLists.txt), no C extensions
- **C++**: C++98 standard for most code, C++11 for some components (PF_DSP, convolution)
- Use `#pragma once` for C++ header guards
- Use standard include guards `#ifndef HEADER_H / #define HEADER_H / #endif` for C headers

### File Organization

- C headers: `.h` (e.g., `pffft.h`, `pffastconv.h`)
- C++ headers: `.hpp` (e.g., `pffft.hpp`)
- C sources: `.c` (e.g., `pffft.c`, `pffft_common.c`)
- C++ sources: `.cpp` (e.g., `pf_mixer.cpp`, `pf_conv.cpp`)
- SIMD headers: `simd/pf_*.h`
- Tests: `test_*.c` and `test_*.cpp`
- Benchmarks: `bench_*.c` and `bench_*.cpp`

### Naming Conventions

#### Types and Structs
```c
// C structs: typedef struct name_s { ... } name_t;
typedef struct PFFFT_Setup PFFFT_Setup;
typedef struct complexf_s { float i; float q; } complexf;

// Enums: snake_case
typedef enum { PFFFT_FORWARD, PFFFT_BACKWARD } pffft_direction_t;
typedef enum { PFFFT_REAL, PFFFT_COMPLEX } pffft_transform_t;
```

#### Functions
```c
// Public API: lowercase with underscores
pffft_new_setup()
pffft_transform()
pffft_destroy_setup()
pffastconv_apply()

// Internal functions: similar naming or with prefixes
```

#### Variables and Parameters
```c
// Parameters: camelCase for C++ functions, snake_case for C functions
float shift_math_cc(const complexf *input, complexf* output, int input_size);

// Local variables: camelCase or snake_case (be consistent within file)
float phaseIncrement;
int input_len;
```

#### Macros
```c
// UPPERCASE with underscores
#define PFFFT_FORWARD 0
#define M_PI 3.14159265358979323846
#define ALWAYS_INLINE(return_type) __attribute__((always_inline))
```

### Imports and Includes

```c
// Include own header first
#include "pffft.h"

// Then standard headers
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>

// Use <...> for system headers, "..." for project headers
```

### Formatting

- Indentation: 2 spaces (consistent with existing codebase)
- Max line length: ~80-100 characters (check existing patterns)
- Opening braces: Same line for functions, structures, etc.
- Space after keywords: `if (`, `while (`, `for (`, `switch (`
- No space before parentheses in function calls: `func(arg)`
- Spaces around operators: `a = b + c`
- Pointer declarations: `float *ptr` or `float* ptr` (be consistent per file)
- Alignment: 16-byte alignment for SIMD buffers (use `pffft_aligned_malloc()`)

### Memory Management

```c
// Use provided aligned allocation functions for SIMD buffers
void *pffft_aligned_malloc(size_t nb_bytes);
void pffft_aligned_free(void *ptr);

// For fast convolution
void *pffastconv_malloc(size_t nb_bytes);
void pffastconv_free(void *ptr);

// Always check for NULL after allocation
void *buf = pffft_aligned_malloc(size);
if (!buf) { /* handle error */ }
```

### Error Handling

- Return error codes or NULL for allocation failures
- Use assertions for internal consistency checks
- Document error conditions in function comments
- For user-facing APIs, return meaningful error codes

### Comments and Documentation

```c
/* Copyright header with license info at top of files */

/*
 * Multi-line comments for describing blocks of code
 * Use for explaining algorithms or complex logic
 */

// Single-line comments for brief explanations
// Use for describing variable purpose or short code blocks

// Function documentation: describe parameters, return values, behavior
/*
   Perform a Fourier transform.
   The z-domain data is stored in the most efficient order.
   Transforms are not scaled: PFFFT_BACKWARD(PFFFT_FORWARD(x)) = N*x.
   Input and output may alias.
*/
void pffft_transform(...);
```

### SIMD-Specific Guidelines

- Code is platform-independent via SIMD abstraction macros
- Use `V4SF`, `V2DF` or similar vector types from SIMD headers
- Never write intrinsics directly in main code - use macros
- SIMD size is retrieved via `pffft_simd_size()` (returns 4 or 1)
- Always align buffers to SIMD vector size requirements (16 bytes typically)

### Constants and Magic Numbers

```c
// Define constants rather than magic numbers
#define EXPECTED_DYN_RANGE  140.0
#define DEG_ERR_LIMIT   1E-4
#define MAG_ERR_LIMIT  1E-6

// For math constants
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif
```

### Compiler-Specific Code

```c
// Detect compiler with preprocessor
#if defined(_MSC_VER)
#  define COMPILER_MSVC
#elif defined(__GNUC__)
#  define COMPILER_GCC
#endif

// Use compiler attributes for optimization hints
#if defined(COMPILER_GCC)
#  define ALWAYS_INLINE(return_type) inline return_type __attribute__((always_inline))
#  define NEVER_INLINE(return_type) return_type __attribute__((noinline))
#elif defined(COMPILER_MSVC)
#  define ALWAYS_INLINE(return_type) __forceinline return_type
#  define NEVER_INLINE(return_type) __declspec(noinline)
#endif
```

### CMake Configuration

- Add new sources to appropriate lists in CMakeLists.txt
- Use `target_compile_definitions()` for preprocessor defines
- Use `target_link_libraries()` for dependencies
- Enable warnings with `target_activate_c_compiler_warnings()` or `target_activate_cxx_compiler_warnings()`
- Set C/C++ standard via `set_property(TARGET name PROPERTY C_STANDARD 99)`

### Testing Guidelines

- Write tests in `test_*.c` or `test_*.cpp` files
- Add tests via `add_test()` in CMakeLists.txt
- Use `assert()` for validation in test code
- Run full test suite before committing changes
- Test both float and double precision when applicable

## Common Build Configurations

```bash
# Default build with float and double
mkdir build && cd build && cmake .. && make

# Float-only build
cmake -DPFFFT_USE_TYPE_DOUBLE=OFF ..

# Double-only build
cmake -DPFFFT_USE_TYPE_FLOAT=OFF ..

# Build without SIMD
cmake -DPFFFT_USE_SIMD=OFF ..

# Debug build with AddressSanitizer
cmake -DCMAKE_BUILD_TYPE=Debug -DPFFFT_USE_DEBUG_ASAN=ON ..

# Build benchmarks with FFTW comparison (from the repository root; never a system FFTW)
python3 benchmarks/build_fftw.py --prefix build_bench_cache/fftw-local
cmake -S . -B build -DPFFFT_USE_BENCH_FFTW=ON -DFFTW3_ROOT=$PWD/build_bench_cache/fftw-local
```

## Platform-Specific Notes

- Windows MSVC: Use `/arch:SSE`, `/arch:AVX` etc. for SIMD
- Windows MinGW: Cross-compilation scripts available
- macOS: Works with native toolchain; Accelerate framework optional
- Linux: Works with GCC or Clang; standard `libm` linked by default
