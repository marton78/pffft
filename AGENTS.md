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

# Run specific test executable directly (after build)
./test_fft_factors
./test_pffft_float
./test_pffft_double
./test_pffastconv --no-bench --quick

# Run benchmark tests
./bench_pffft_float --max-len 128 --quick
./bench_pffft_float --non-pow2 --max-len 192 --quick
./bench_mixers
./bench_conv

# Debug with AddressSanitizer (must rebuild)
cmake -DPFFFT_USE_DEBUG_ASAN=ON ..
cmake --build .
ctest
```

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

# Build benchmarks with FFTW comparison
cmake -DPFFFT_USE_BENCH_FFTW=ON ..
```

## Platform-Specific Notes

- Windows MSVC: Use `/arch:SSE`, `/arch:AVX` etc. for SIMD
- Windows MinGW: Cross-compilation scripts available
- macOS: Works with native toolchain; Accelerate framework optional
- Linux: Works with GCC or Clang; standard `libm` linked by default
