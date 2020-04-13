cmake_minimum_required(VERSION 2.8)
project(PRETTY_FAST_FFT)

# smaller library size?
option(USE_TYPE_FLOAT  "activate single precision 'float'?" ON)
option(USE_TYPE_DOUBLE "activate 'double' precision float?" ON)

# architecture/optimization options
option(USE_SIMD        "use SIMD (SSE/AVX/NEON/ALTIVEC) CPU features? - " ON)
option(USE_SIMD_NEON   "force using NEON on ARM? (requires USE_SIMD)" OFF)
option(USE_SCALAR_VECT "use 4-element vector scalar operations (if no other SIMD)" ON)

# test options
option(USE_BENCH_FFTW  "use (system-installed) FFTW3 in fft benchmark?" OFF)
option(USE_BENCH_GREEN "use Green FFT in fft benchmark? - if exists in subdir" ON)
option(USE_BENCH_KISS  "use KissFFT in fft benchmark? - if exists in subdir" ON)

option(USE_DEBUG_ASAN  "use GCC's address sanitizer?" OFF)


# C90 requires the gcc extensions for function attributes like always_inline
# C99 provides the function attributes: no gcc extensions required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if ( (NOT USE_TYPE_FLOAT) AND (NOT USE_TYPE_DOUBLE) )
  message(FATAL_ERROR "activate at least one of USE_TYPE_FLOAT or USE_TYPE_DOUBLE")
endif()


if (USE_DEBUG_ASAN)
  set(ASANLIB "asan")
else()
  set(ASANLIB "")
endif()


if (USE_BENCH_GREEN)
  if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/greenffts/CMakeLists.txt")
    message(STATUS "found subdir greenffts")
    set(PATH_GREEN "${CMAKE_CURRENT_LIST_DIR}/greenffts")
    add_subdirectory( "${PATH_GREEN}" )
  else()
    message(WARNING "GreenFFT not found in subdir greenffts")
  endif()
endif()

if (USE_BENCH_KISS)
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

set( SIMD_FLOAT_HDRS simd/pf_float.h simd/pf_sse1_float.h simd/pf_altivec_float.h simd/pf_neon_float.h simd/pf_scalar_float.h )
set( SIMD_DOUBLE_HDRS simd/pf_double.h simd/pf_avx_double.h simd/pf_scalar_double.h )

if (USE_TYPE_FLOAT)
  set( FLOAT_SOURCES pffft.c pffft.h ${SIMD_FLOAT_HDRS} )
else()
  set( FLOAT_SOURCES  )
endif()


if (USE_TYPE_DOUBLE)
  set( DOUBLE_SOURCES pffft_double.c pffft_double.h ${SIMD_DOUBLE_HDRS} )
else()
  set( DOUBLE_SOURCES )
endif()

######################################################

add_library(PFFFT STATIC ${FLOAT_SOURCES} ${DOUBLE_SOURCES} pffft_common.c pffft_priv_impl.h pffft.hpp )
target_compile_definitions(PFFFT PRIVATE _USE_MATH_DEFINES)
if (USE_SCALAR_VECT)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SCALVEC_ENABLED=1)
endif()
if (USE_DEBUG_ASAN)
  target_compile_options(PFFFT PRIVATE "-fsanitize=address")
endif()
if (NOT USE_SIMD)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
if (USE_SIMD AND USE_SIMD_NEON)
  target_compile_definitions(PFFFT PRIVATE PFFFT_ENABLE_NEON=1)
  target_compile_options(PFFFT PRIVATE "-mfpu=neon")
endif()
if (USE_SIMD AND USE_TYPE_DOUBLE)
  if(WIN32)
    set_property(SOURCE pffft_double.c PROPERTY COMPILE_FLAGS "/arch:AVX")
  else()
    set_property(SOURCE pffft_double.c PROPERTY COMPILE_FLAGS "-march=native")
  endif()
endif()
target_link_libraries( PFFFT ${MATHLIB} )
set_property(TARGET PFFFT APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

######################################################

add_library(FFTPACK STATIC fftpack.c fftpack.h)
target_compile_definitions(FFTPACK PRIVATE _USE_MATH_DEFINES)
target_link_libraries( FFTPACK ${MATHLIB} )
set_property(TARGET FFTPACK APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)

######################################################

if (USE_TYPE_FLOAT)
  # only 'float' supported in PFFASTCONV
  add_library(PFFASTCONV STATIC pffastconv.c pffastconv.h pffft.h )
  target_compile_definitions(PFFASTCONV PRIVATE _USE_MATH_DEFINES)
  if (USE_DEBUG_ASAN)
    target_compile_options(PFFASTCONV PRIVATE "-fsanitize=address")
  endif()
  target_link_libraries( PFFASTCONV PFFFT ${ASANLIB} ${MATHLIB} )
  set_property(TARGET PFFASTCONV APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )
endif()

#######################################################

if (USE_TYPE_FLOAT)
  add_executable( test_pffft_float  test_pffft.c )
  target_compile_definitions(test_pffft_float PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(test_pffft_float PRIVATE PFFFT_ENABLE_FLOAT)
  target_link_libraries( test_pffft_float  PFFFT ${ASANLIB} )
endif()

######################################################

if (USE_TYPE_DOUBLE)
  add_executable( test_pffft_double  test_pffft.c )
  target_compile_definitions(test_pffft_double PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(test_pffft_double PRIVATE PFFFT_ENABLE_DOUBLE)
  target_link_libraries( test_pffft_double  PFFFT ${ASANLIB} )
endif()

######################################################

add_executable( test_pffft_cpp test_pffft.cpp )
target_compile_definitions(test_pffft_cpp PRIVATE _USE_MATH_DEFINES)
if (USE_TYPE_FLOAT)
  target_compile_definitions(test_pffft_cpp PRIVATE PFFFT_ENABLE_FLOAT)
endif()
if (USE_TYPE_DOUBLE)
  target_compile_definitions(test_pffft_cpp PRIVATE PFFFT_ENABLE_DOUBLE)
endif()
target_link_libraries( test_pffft_cpp  PFFFT ${ASANLIB} )

######################################################

add_executable( test_pffft_cpp_11 test_pffft.cpp )
target_compile_definitions(test_pffft_cpp_11 PRIVATE _USE_MATH_DEFINES)
if (USE_TYPE_FLOAT)
  target_compile_definitions(test_pffft_cpp_11 PRIVATE PFFFT_ENABLE_FLOAT)
endif()
if (USE_TYPE_DOUBLE)
  target_compile_definitions(test_pffft_cpp_11 PRIVATE PFFFT_ENABLE_DOUBLE)
endif()
target_link_libraries( test_pffft_cpp_11  PFFFT ${ASANLIB} )

set_property(TARGET test_pffft_cpp_11 PROPERTY CXX_STANDARD 11)
set_property(TARGET test_pffft_cpp_11 PROPERTY CXX_STANDARD_REQUIRED ON)

######################################################

if (USE_TYPE_FLOAT)
  add_executable(test_pffastconv   test_pffastconv.c
    ${SIMD_FLOAT_HDRS} ${SIMD_DOUBLE_HDRS}
  )
  target_compile_definitions(test_pffastconv PRIVATE _USE_MATH_DEFINES)
  if (USE_DEBUG_ASAN)
    target_compile_options(test_pffastconv PRIVATE "-fsanitize=address")
  endif()
  if (NOT USE_SIMD)
    target_compile_definitions(test_pffastconv PRIVATE PFFFT_SIMD_DISABLE=1)
  endif()
  if (USE_SIMD AND USE_SIMD_NEON)
    target_compile_definitions(test_pffastconv PRIVATE PFFFT_ENABLE_NEON=1)
    target_compile_options(test_pffastconv PRIVATE "-mfpu=neon")
  endif()
  target_link_libraries( test_pffastconv  PFFASTCONV ${ASANLIB} ${MATHLIB} )
endif()

######################################################

if (USE_TYPE_FLOAT)
  add_executable(bench_pffft_float   bench_pffft.c pffft.h fftpack.h)
  target_compile_definitions(bench_pffft_float PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(bench_pffft_float PRIVATE PFFFT_ENABLE_FLOAT)

  target_link_libraries( bench_pffft_float  PFFFT FFTPACK ${ASANLIB} )

  if (USE_BENCH_FFTW)
    target_compile_definitions(bench_pffft_float PRIVATE HAVE_FFTW=1)
    target_link_libraries(bench_pffft_float  fftw3f)
  endif()

  if (PATH_GREEN AND USE_BENCH_GREEN)
    target_compile_definitions(bench_pffft_float PRIVATE HAVE_GREEN_FFTS=1)
    target_link_libraries(bench_pffft_float  GreenFFT)
  endif()

  if (PATH_KISS AND USE_BENCH_KISS)
    target_compile_definitions(bench_pffft_float PRIVATE HAVE_KISS_FFT=1)
    target_link_libraries(bench_pffft_float  KissFFT)
  endif()
endif()

if (USE_TYPE_DOUBLE)
  add_executable(bench_pffft_double   bench_pffft.c pffft.h fftpack.h)
  target_compile_definitions(bench_pffft_double PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(bench_pffft_double PRIVATE PFFFT_ENABLE_DOUBLE)
  target_link_libraries( bench_pffft_double  PFFFT ${ASANLIB} )

  if (USE_BENCH_FFTW)
    target_compile_definitions(bench_pffft_double PRIVATE HAVE_FFTW=1)
    target_link_libraries(bench_pffft_double  fftw3)
  endif()
endif()

######################################################

enable_testing()

if (USE_TYPE_FLOAT)

  add_test(NAME bench_pffft_pow2
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/bench_pffft_float"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(NAME bench_pffft_non2
    COMMAND "${CMAKE_CURRENT_BINARY_DIR}/bench_pffft_float" "--non-pow2"
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

endif()

