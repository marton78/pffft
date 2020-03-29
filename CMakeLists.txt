cmake_minimum_required(VERSION 2.8)
project(PRETTY_FAST_FFT)

option(USE_SIMD  "use SIMD (SSE/NEON) CPU features?" ON)
option(USE_AVX   "use SIMD (AVX) CPU features? (requires SIMD)" OFF)
option(USE_NEON  "force using NEON on ARM? (requires SIMD)" OFF)
option(USE_FFTW  "use (system-installed) FFTW3 in fft benchmark?" OFF)
option(USE_GREEN "use Green FFT in fft benchmark? - if exists in subdir" ON)
option(USE_KISS  "use KissFFT in fft benchmark? - if exists in subdir" ON)
option(USE_ASAN  "use GCC's address sanitizer?" OFF)


# C90 requires the gcc extensions for function attributes like always_inline
# C99 provides the function attributes: no gcc extensions required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if (USE_ASAN)
  set(ASANLIB "asan")
else()
  set(ASANLIB "")
endif()


if (USE_GREEN)
  if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/greenffts/CMakeLists.txt")
    message(STATUS "found subdir greenffts")
    set(PATH_GREEN "${CMAKE_CURRENT_LIST_DIR}/greenffts")
    add_subdirectory( "${PATH_GREEN}" )
  else()
    message(WARNING "GreenFFT not found in subdir greenffts")
  endif()
endif()

if (USE_KISS)
  # git submodule add https://github.com/hayguen/kissfft.git
  if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/kissfft/CMakeLists.txt")
    message(STATUS "found subdir kissfft")
    set(PATH_KISS "${CMAKE_CURRENT_LIST_DIR}/kissfft")
    add_subdirectory( "${PATH_KISS}" )
  else()
    message(WARNING "KissFFT not found in subdir kissfft")
  endif()
endif()


########################################################################
# select the release build type by default to get optimization flags
########################################################################
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message(STATUS "INFO: detected MSVC: will not link math lib m")
  set(MATHLIB "")
else()
  message(STATUS "INFO: detected NO MSVC: ${CMAKE_C_COMPILER_ID}: will link math lib m")
  set(MATHLIB "m")
endif()

set( SIMD_FLOAT_HDRS pf_float.h pf_sse1_float.h pf_altivec_float.h pf_neon_float.h pf_scalar_float.h )
set( SIMD_DOUBLE_HDRS pf_double.h pf_avx_double.h pf_scalar_double.h )


add_library(PFFFT STATIC
  pffft.c pffft.h ${SIMD_FLOAT_HDRS}
  pffft_double.c pffft_double.h ${SIMD_DOUBLE_HDRS}
)
target_compile_definitions(PFFFT PRIVATE _USE_MATH_DEFINES)
if (USE_ASAN)
  target_compile_options(PFFFT PRIVATE "-fsanitize=address")
endif()
if (NOT USE_SIMD)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
if (USE_SIMD AND USE_NEON)
  target_compile_definitions(PFFFT PRIVATE PFFFT_ENABLE_NEON=1)
  target_compile_options(PFFFT PRIVATE "-mfpu=neon")
endif()
if (USE_SIMD AND USE_AVX)
  target_compile_definitions(PFFFT PRIVATE PFFFT_ENABLE_AVX=1)
  if(WIN32)
    #target_compile_options(PFFFT PRIVATE "/arch:AVX")
    set_property(SOURCE pffft_double.c PROPERTY COMPILE_FLAGS "/arch:AVX")
  else()
    #target_compile_options(PFFFT PRIVATE "-march=native")
    set_property(SOURCE pffft_double.c PROPERTY COMPILE_FLAGS "-march=native")
  endif()
endif()
target_link_libraries( PFFFT ${MATHLIB} )
set_property(TARGET PFFFT APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_library(FFTPACK STATIC fftpack.c fftpack.h)
target_compile_definitions(FFTPACK PRIVATE _USE_MATH_DEFINES)
target_link_libraries( FFTPACK ${MATHLIB} )
set_property(TARGET FFTPACK APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_library(PFFASTCONV STATIC pffastconv.c pffastconv.h pffft.h ${SIMD_FLOAT_HDRS} )
target_compile_definitions(PFFASTCONV PRIVATE _USE_MATH_DEFINES)
if (USE_ASAN)
  target_compile_options(PFFASTCONV PRIVATE "-fsanitize=address")
endif()
if (NOT USE_SIMD)
  target_compile_definitions(PFFASTCONV PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
if (USE_SIMD AND USE_NEON)
  target_compile_options(PFFASTCONV PRIVATE "-mfpu=neon")
endif()
target_link_libraries( PFFASTCONV PFFFT ${ASANLIB} ${MATHLIB} )
set_property(TARGET PFFASTCONV APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_executable( test_pffft  test_pffft.c )
target_compile_definitions(test_pffft PRIVATE _USE_MATH_DEFINES)
target_link_libraries( test_pffft  PFFFT ${ASANLIB} )

add_executable( test_pffft_double  test_pffft_double.c )
target_compile_definitions(test_pffft_double PRIVATE _USE_MATH_DEFINES)
target_link_libraries( test_pffft_double  PFFFT ${ASANLIB} )

add_executable( test_pffft_cpp test_pffft.cpp )
target_compile_definitions(test_pffft_cpp PRIVATE _USE_MATH_DEFINES)
target_link_libraries( test_pffft_cpp  PFFFT ${ASANLIB} )


add_executable(test_pffastconv   test_pffastconv.c )
target_compile_definitions(test_pffastconv PRIVATE _USE_MATH_DEFINES)
if (USE_ASAN)
  target_compile_options(test_pffastconv PRIVATE "-fsanitize=address")
endif()
if (NOT USE_SIMD)
  target_compile_definitions(test_pffastconv PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
if (USE_SIMD AND USE_NEON)
  target_compile_options(test_pffastconv PRIVATE "-mfpu=neon")
endif()
target_link_libraries( test_pffastconv  PFFASTCONV ${ASANLIB} ${MATHLIB} )


add_executable(bench_pffft   bench_pffft.c pffft.h fftpack.h)
target_compile_definitions(bench_pffft PRIVATE _USE_MATH_DEFINES)
if (NOT USE_SIMD)
  target_compile_definitions(bench_pffft PRIVATE PFFFT_SIMD_DISABLE=1)
endif()

target_link_libraries( bench_pffft  PFFFT FFTPACK ${ASANLIB} )

if (USE_FFTW)
  target_compile_definitions(bench_pffft PRIVATE HAVE_FFTW=1)
  target_link_libraries(bench_pffft  fftw3f)
endif()

if (PATH_GREEN AND USE_GREEN)
  target_compile_definitions(bench_pffft PRIVATE HAVE_GREEN_FFTS=1)
  target_link_libraries(bench_pffft  GreenFFT)
endif()

if (PATH_KISS AND USE_KISS)
  target_compile_definitions(bench_pffft PRIVATE HAVE_KISS_FFT=1)
  target_link_libraries(bench_pffft  KissFFT)
endif()


enable_testing()

add_test(NAME bench_pffft_pow2
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/bench_pffft"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME bench_pffft_non2
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/bench_pffft" "--non-pow2"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME bench_plots
  COMMAND bash "-c" "${CMAKE_CURRENT_SOURCE_DIR}/plots.sh"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_pfconv_lens_symetric
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_pffastconv" "--no-bench" "--quick" "--sym"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_pfconv_lens_non_sym
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_pffastconv" "--no-bench" "--quick"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME bench_pfconv_symetric
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_pffastconv" "--no-len" "--sym"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME bench_pfconv_non_sym
  COMMAND "${CMAKE_CURRENT_BINARY_DIR}/test_pffastconv" "--no-len"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

