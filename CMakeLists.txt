cmake_minimum_required(VERSION 2.8)
set(CMAKE_C_STANDARD 99)
project(PRETTY_FAST_FFT)

option(USE_SIMD  "use SIMD (SSE/NEON) CPU features?" ON)
option(USE_NEON  "force using NEON on ARM? (requires SIMD)" OFF)
option(USE_FFTW  "use (system-installed) FFTW3 in benchmark?" OFF)
option(USE_GREEN "use Green FFT in benchmark? - if exists in subdir" ON)
option(USE_KISS  "use KissFFT in benchmark? - if exists in subdir" ON)


if (USE_GREEN)
  if (IS_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/greenffts")
    message(STATUS "found subdir greenffts")
    set(PATH_GREEN "${CMAKE_CURRENT_LIST_DIR}/greenffts")
    add_subdirectory( "${PATH_GREEN}" )
  else()
    message(WARNING "GreenFFT not found in subdir greenffts")
  endif()
endif()

if (USE_KISS)
  # git submodule add https://github.com/hayguen/kissfft.git
  if (IS_DIRECTORY "${CMAKE_CURRENT_LIST_DIR}/kissfft")
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

add_library(PFFFT STATIC pffft.c)
if (NOT USE_SIMD)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
if (USE_SIMD AND USE_NEON)
  target_compile_options(PFFFT PRIVATE "-mfpu=neon")
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


add_executable(bench_pffft   bench_pffft.c )

if (NOT USE_SIMD)
  target_compile_definitions(bench_pffft PRIVATE PFFFT_SIMD_DISABLE=1)
endif()

target_link_libraries(bench_pffft  PFFFT FFTPACK)

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


