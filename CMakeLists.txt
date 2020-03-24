cmake_minimum_required(VERSION 2.8)
set(CMAKE_C_STANDARD 99)
project(PRETTY_FAST_FFT)

option(USE_SIMD  "use SIMD (SSE/NEON) CPU features?" ON)
option(USE_NEON  "force using NEON on ARM? (requires SIMD)" OFF)
option(USE_FFTW  "use (system-installed) FFTW3 in fft benchmark?" OFF)
option(USE_GREEN "use Green FFT in fft benchmark? - if exists in subdir" ON)
option(USE_KISS  "use KissFFT in fft benchmark? - if exists in subdir" ON)
option(USE_ASAN  "use GCC's address sanitizer?" OFF)

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


add_library(PFFFT STATIC pffft.c)
target_compile_definitions(PFFFT PRIVATE _USE_MATH_DEFINES)
if (USE_ASAN)
  target_compile_options(PFFFT PRIVATE "-fsanitize=address")
endif()
if (NOT USE_SIMD)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
if (USE_SIMD AND USE_NEON)
  target_compile_options(PFFFT PRIVATE "-mfpu=neon")
endif()
target_link_libraries( PFFFT ${MATHLIB} )
set_property(TARGET PFFFT APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

add_library(PFFFT_DOUBLE STATIC pffft_double.c)
target_compile_definitions(PFFFT_DOUBLE PRIVATE _USE_MATH_DEFINES)
if (USE_ASAN)
  target_compile_options(PFFFT_DOUBLE PRIVATE "-fsanitize=address")
endif()
if (NOT USE_SIMD)
  target_compile_definitions(PFFFT_DOUBLE PRIVATE PFFFTD_SIMD_DISABLE=1)
endif()
if (USE_SIMD AND USE_AVX)
if(WIN32)
  target_compile_options(PFFFT_DOUBLE PRIVATE "/arch=AVX")
  else(WIN32)
  target_compile_options(PFFFT_DOUBLE PRIVATE "-march=native")
  endif(WIN32)
endif()
target_link_libraries( PFFFT_DOUBLE ${MATHLIB} )
set_property(TARGET PFFFT_DOUBLE APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_library(FFTPACK STATIC fftpack.c)
target_compile_definitions(FFTPACK PRIVATE _USE_MATH_DEFINES)
target_link_libraries( FFTPACK ${MATHLIB} )
set_property(TARGET FFTPACK APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_library(PFFASTCONV STATIC pffastconv.c)
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
target_link_libraries( test_pffft_double  PFFFTD ${ASANLIB} )


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


add_executable(bench_pffft   bench_pffft.c )
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

