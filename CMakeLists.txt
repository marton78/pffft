cmake_minimum_required(VERSION 3.0)
set(CMAKE_C_STANDARD 99)
project(PRETTY_FAST_FFT LANGUAGES C)

option(USE_SIMD "use SIMD (SSE/NEON) CPU features?" ON)
option(USE_FFTW "use fftw3 in benchmark?" OFF)

########################################################################
# select the release build type by default to get optimization flags
########################################################################
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)

add_library(PFFFT STATIC pffft.c)
if (NOT USE_SIMD)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
target_link_libraries(PFFFT m)
set_property(TARGET PFFFT APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

add_library(FFTPACK STATIC fftpack.c)
target_link_libraries(FFTPACK m)
set_property(TARGET FFTPACK APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_executable(bench_pffft   test_pffft.c )
if (NOT USE_SIMD)
  target_compile_definitions(bench_pffft PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
target_link_libraries(bench_pffft  PFFFT FFTPACK)

if (USE_FFTW)
  target_compile_definitions(bench_pffft PRIVATE HAVE_FFTW=1)
  target_link_libraries(bench_pffft  fftw3f)
endif()