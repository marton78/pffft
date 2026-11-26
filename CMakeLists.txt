cmake_minimum_required(VERSION 3.10)
project(PRETTY_FAST_FFT)

# prefix this project's messages with [pffft] under cmake --log-context
list(APPEND CMAKE_MESSAGE_CONTEXT pffft)

include(GNUInstallDirs)

# detect if pffft is the top-level project or built as a dependency
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(PFFFT_IS_TOPLEVEL ON)
else()
  set(PFFFT_IS_TOPLEVEL OFF)
endif()

# smaller library size?
option(PFFFT_USE_TYPE_FLOAT  "activate single precision 'float'?" ON)
option(PFFFT_USE_TYPE_DOUBLE "activate 'double' precision float?" ON)

# architecture/optimization options
option(PFFFT_USE_SIMD        "use SIMD (SSE/AVX/NEON/ALTIVEC) CPU features? - " ON)
option(PFFFT_USE_SCALAR_VECT "use 4-element vector scalar operations (if no other SIMD)" ON)
option(PFFFT_USE_HELIUM "use the experimental, unbenchmarked Arm Helium (MVE) backend on Armv8.1-M?" OFF)
option(PFFFT_USE_WASM_RELAXED_SIMD "use WebAssembly Relaxed SIMD for higher performance? requires a browser/runtime with Relaxed SIMD support" OFF)
set(PFFFT_PREFETCH_DIST "" CACHE STRING "x86: prefetch distance in bytes for the FFT passes; 0 turns prefetching off (set it for AMD CPUs, see README); empty: the default")
set(PFFFT_PREFETCH_MIN "" CACHE STRING "x86: smallest pass in bytes that prefetches; empty: the default (24576)")
set(PFFFT_SPLIT_MAX "" CACHE STRING "largest array in bytes that the transforms do not split into blocks that fit the cache; empty: the default, by target (1 MiB on x86 but Zen and on Android, 8 MiB on Zen and Apple arm64, 16 MiB on other Arm)")

# what to install?
option(INSTALL_PFFFT      "install pffft to CMAKE_INSTALL_PREFIX?" ON)
option(INSTALL_PFDSP      "install pfdsp to CMAKE_INSTALL_PREFIX?" OFF)
option(INSTALL_PFFASTCONV "install pffastconv to CMAKE_INSTALL_PREFIX?" OFF)

# test options
option(PFFFT_USE_BENCH_FFTW   "use FFTW3 from FFTW3_ROOT (see benchmarks/build_fftw.py) in fft benchmark?" OFF)
option(PFFFT_USE_BENCH_GREEN  "use Green FFT in fft benchmark? - if exists in subdir" ON)
option(PFFFT_USE_BENCH_KISS   "use KissFFT in fft benchmark? - if exists in subdir" ON)
option(PFFFT_USE_BENCH_POCKET "use PocketFFT in fft benchmark? - if exists in subdir" ON)
option(PFFFT_USE_BENCH_OOURA  "use Ooura's FFT in fft benchmark? - if exists in subdir" ON)
option(PFFFT_USE_BENCH_MKL    "use Intel MKL in fft benchmark? needs to be installed" OFF)
option(PFFFT_USE_FFTPACK      "compile and use FFTPACK in fft benchmark & validation?" ON)
if (APPLE)
  option(PFFFT_USE_BENCH_VDSP "use vDSP of Apple's Accelerate framework in fft benchmark?" ON)
endif()

option(PFFFT_USE_DEBUG_ASAN  "use the address sanitizer (gcc, clang)?" OFF)

option(PFFFT_DISABLE_LINK_WITH_M "Disables linking with m library to build with clangCL from MSVC" OFF)

# build control options (default to ON for top-level builds, OFF for subdirectory/FetchContent)
option(PFFFT_BUILD_TESTS      "Build pffft tests?" ${PFFFT_IS_TOPLEVEL})
option(PFFFT_BUILD_BENCHMARKS "Build pffft benchmarks?" ${PFFFT_IS_TOPLEVEL})
option(PFFFT_BUILD_EXAMPLES   "Build pffft examples?" ${PFFFT_IS_TOPLEVEL})

# C90 requires the gcc extensions for function attributes like always_inline
# C99 provides the function attributes: no gcc extensions required
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_CXX_STANDARD 98)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# root source directory (for subdirectory CMakeLists.txt to reference source files)
set(PFFFT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# populate what to install
set(INSTALL_TARGETS "")
set(INSTALL_HEADERS "")


if ( (NOT PFFFT_USE_TYPE_FLOAT) AND (NOT PFFFT_USE_TYPE_DOUBLE) )
  message(FATAL_ERROR "activate at least one of PFFFT_USE_TYPE_FLOAT or PFFFT_USE_TYPE_DOUBLE")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(cmake/target_optimizations.cmake)
include(cmake/compiler_warnings.cmake)
include(cmake/simd_backend_probe.cmake)
find_package(PAPI)


if (PFFFT_USE_DEBUG_ASAN AND NOT MSVC)
  # the compiler driver links its own runtime: libasan for gcc, the
  # compiler-rt one for clang and AppleClang, which have no -lasan
  set(ASANLIB "-fsanitize=address")
else()
  # MSVC needs /fsanitize=address only when compiling
  set(ASANLIB "")
endif()

message(STATUS "INFO: CMAKE_C_COMPILER_ID is ${CMAKE_C_COMPILER_ID}")
message(STATUS "INFO: CMAKE_CXX_COMPILER_ID is ${CMAKE_CXX_COMPILER_ID}")
if (WIN32)
  message(STATUS "INFO: detected WIN32")
else()
  message(STATUS "INFO: NOT WIN32")
endif()
if (MINGW)
  message(STATUS "INFO: detected MINGW with compiler ${CMAKE_C_COMPILER_ID}")
else()
  message(STATUS "INFO: NOT MINGW")
endif()
if ( CMAKE_C_COMPILER_ID MATCHES "MSVC" )
  message(STATUS "INFO: detected MSVC with compiler ${CMAKE_C_COMPILER_ID}")
endif()


########################################################################
# select the release build type by default to get optimization flags
########################################################################
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)

if ( CMAKE_C_COMPILER_ID MATCHES "MSVC" )
  # using Visual Studio C++
  add_definitions("/D_CRT_SECURE_NO_WARNINGS")

  set(MSVC_DISABLED_WARNINGS_LIST
      "C4996"
  )
endif()

# The math library is only required on Unix
if(NOT WIN32 AND NOT APPLE)
  find_library(MATH_LIBRARY m)
  if(MATH_LIBRARY)
    if(PFFFT_DISABLE_LINK_WITH_M)
      set(MATHLIB "")
      message(STATUS "INFO: will not link math lib m as it is disabled")
    else()
      message(STATUS "INFO: detected Unix and found library: ${CMAKE_C_COMPILER_ID} will link math lib m")
      set(MATHLIB ${MATH_LIBRARY})
    endif()
  else()
    set(MATHLIB "")
    message(WARNING "Could not find math lib m")
  endif()
else()
  set(MATHLIB "")
  message(STATUS "INFO: will not link math lib m, no need for this target")
endif()

set(STDCXXLIB "")
if (MINGW)
  set(STDCXXLIB "stdc++")
endif()


set( SIMD_FLOAT_HDRS ${PFFFT_SOURCE_DIR}/src/simd/pf_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_sse1_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_altivec_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_neon_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_wasm_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_rvv_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_helium_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_scalar_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_fma_x86.h )
set( SIMD_DOUBLE_HDRS ${PFFFT_SOURCE_DIR}/src/simd/pf_double.h ${PFFFT_SOURCE_DIR}/src/simd/pf_avx_double.h ${PFFFT_SOURCE_DIR}/src/simd/pf_wasm_double.h ${PFFFT_SOURCE_DIR}/src/simd/pf_rvv_double.h ${PFFFT_SOURCE_DIR}/src/simd/pf_scalar_double.h ${PFFFT_SOURCE_DIR}/src/simd/pf_fma_x86.h )

if (PFFFT_USE_TYPE_FLOAT)
  set( FLOAT_SOURCES src/pffft.c include/pffft/pffft.h ${SIMD_FLOAT_HDRS} )
  if (INSTALL_PFFFT)
    set(INSTALL_HEADERS ${INSTALL_HEADERS} include/pffft/pffft.h)
  endif()
else()
  set( FLOAT_SOURCES  )
endif()


if (PFFFT_USE_TYPE_DOUBLE)
  set( DOUBLE_SOURCES src/pffft_double.c include/pffft/pffft_double.h ${SIMD_DOUBLE_HDRS} )
  if (INSTALL_PFFFT)
    set(INSTALL_HEADERS ${INSTALL_HEADERS} include/pffft/pffft_double.h)
  endif()
else()
  set( DOUBLE_SOURCES )
endif()

######################################################

add_library(PFFFT ${FLOAT_SOURCES} ${DOUBLE_SOURCES} src/pffft_common.c src/pffft_priv_impl.h include/pffft/pffft.hpp )
add_library(PFFFT::PFFFT ALIAS PFFFT)
set_target_properties(PFFFT PROPERTIES OUTPUT_NAME "pffft")
target_compile_definitions(PFFFT PRIVATE _USE_MATH_DEFINES)
target_activate_c_compiler_warnings(PFFFT)
if (PFFFT_USE_SCALAR_VECT)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SCALVEC_ENABLED=1)
endif()
if (PFFFT_USE_DEBUG_ASAN)
  target_compile_options(PFFFT PRIVATE "-fsanitize=address")
endif()
target_set_c_arch_flags(PFFFT)
# GCC's RISC-V vector code fuses multiplies and adds into FMAs even without
# floating-point contraction (the default with -std=c99), and not alike in
# every function: double-precision ordered transforms then differed in their
# last bits from pffft_transform followed by pffft_zreorder (gcc 14 and 15).
# With contraction, GCC fuses them in the source, the same way everywhere.
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_PROCESSOR MATCHES "riscv64")
  target_compile_options(PFFFT PRIVATE "-ffp-contract=fast")
endif()
if (NOT PFFFT_USE_SIMD)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
foreach(knob PFFFT_PREFETCH_DIST PFFFT_PREFETCH_MIN PFFFT_SPLIT_MAX)
  if (NOT ${knob} STREQUAL "")
    # decimal only: a leading 0 would make the preprocessor read octal
    if (NOT ${knob} MATCHES "^(0|[1-9][0-9]*)$")
      message(FATAL_ERROR "${knob} must be empty or a number of bytes, not '${${knob}}'")
    endif()
    target_compile_definitions(PFFFT PRIVATE "${knob}=${${knob}}")
  endif()
endforeach()
pffft_report_simd_backend(PFFFT)
target_link_libraries( PFFFT ${ASANLIB} ${MATHLIB} )
target_include_directories(PFFFT
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/include/pffft
)
set_target_properties(PFFFT PROPERTIES
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(PFFFT PUBLIC PFFFT_STATIC_DEFINE)
endif()
if (INSTALL_PFFFT)
  set(INSTALL_TARGETS ${INSTALL_TARGETS} PFFFT)
  set(INSTALL_HEADERS ${INSTALL_HEADERS} include/pffft/pffft.hpp)
endif()

######################################################

if (PFFFT_USE_TYPE_FLOAT)
  add_library(PFDSP src/pf_mixer.cpp include/pffft/pf_mixer.h include/pffft/pf_cplx.h src/pf_carrier.cpp include/pffft/pf_carrier.h src/pf_cic.cpp include/pffft/pf_cic.h src/fmv.h )
  set_property(TARGET PFDSP PROPERTY CXX_STANDARD 11)
  set_property(TARGET PFDSP PROPERTY CXX_STANDARD_REQUIRED ON)
  set_target_properties(PFDSP PROPERTIES OUTPUT_NAME "pfdsp")
  target_compile_definitions(PFDSP PRIVATE _USE_MATH_DEFINES)
  target_activate_cxx_compiler_warnings(PFDSP)
  if (PFFFT_USE_DEBUG_ASAN)
      target_compile_options(PFDSP PRIVATE "-fsanitize=address")
  endif()
  if (PFFFT_USE_SIMD)
      target_set_cxx_arch_flags(PFDSP)
  else()
      target_compile_definitions(PFDSP PRIVATE PFFFT_SIMD_DISABLE=1)
  endif()
  target_link_libraries( PFDSP ${MATHLIB} )
  target_include_directories(PFDSP
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/src
      ${CMAKE_CURRENT_SOURCE_DIR}/include/pffft
  )
  set_target_properties(PFDSP PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
  )
  if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(PFDSP PUBLIC PFDSP_STATIC_DEFINE)
  endif()
  if (INSTALL_PFDSP)
      set(INSTALL_TARGETS ${INSTALL_TARGETS} PFDSP)
      set(INSTALL_HEADERS ${INSTALL_HEADERS} include/pffft/pf_mixer.h include/pffft/pf_cplx.h include/pffft/pf_carrier.h include/pffft/pf_cic.h)
  endif()
endif()

######################################################

if (PFFFT_USE_TYPE_FLOAT)
  # only 'float' supported in PFFASTCONV
  add_library(PFFASTCONV src/pffastconv.c include/pffft/pffastconv.h include/pffft/pffft.h )
  set_target_properties(PFFASTCONV PROPERTIES OUTPUT_NAME "pffastconv")
  target_compile_definitions(PFFASTCONV PRIVATE _USE_MATH_DEFINES)
  target_activate_c_compiler_warnings(PFFASTCONV)
  if (PFFFT_USE_DEBUG_ASAN)
    target_compile_options(PFFASTCONV PRIVATE "-fsanitize=address")
  endif()
  target_link_libraries( PFFASTCONV PFFFT ${ASANLIB} ${MATHLIB} )
  target_include_directories(PFFASTCONV
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include/pffft
  )
  set_target_properties(PFFASTCONV PROPERTIES
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
  )
  if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(PFFASTCONV PUBLIC PFFASTCONV_STATIC_DEFINE)
  endif()
  if (INSTALL_PFFASTCONV)
    set(INSTALL_TARGETS ${INSTALL_TARGETS} PFFASTCONV)
    set(INSTALL_HEADERS ${INSTALL_HEADERS} include/pffft/pffastconv.h)
  endif()
endif()


######################################################

if (INSTALL_TARGETS)
  install( TARGETS ${INSTALL_TARGETS}
    EXPORT pffft-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    # a DLL is a RUNTIME artifact; CMake fills in a default destination for
    # it only when neither ARCHIVE nor RUNTIME is given
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
  install( FILES  ${INSTALL_HEADERS}  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pffft)

  # If ${INSTALL_TARGETS} is empty, then the export "pffft-targets" will
  # not exist and the following command will fail so we must make it conditional
  install(EXPORT pffft-targets
    FILE pffft-targets.cmake
    NAMESPACE PFFFT::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pffft
  )

  include(CMakePackageConfigHelpers)

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/pffft-config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/pffft-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pffft
  )

  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/pffft-config-version.cmake"
    VERSION 1.3.1
    COMPATIBILITY SameMajorVersion
  )

  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/pffft-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/pffft-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pffft
  )
  if (PFFFT_IS_TOPLEVEL AND NOT TARGET uninstall)
    add_custom_target(uninstall
        "${CMAKE_COMMAND}" -P "${CMAKE_SOURCE_DIR}/uninstall.cmake"
    )
  endif()
endif()

######################################################

if (PFFFT_IS_TOPLEVEL AND (PFFFT_BUILD_TESTS OR PFFFT_BUILD_BENCHMARKS))
  enable_testing()
endif()

if (PFFFT_BUILD_TESTS)
  add_subdirectory(tests)
endif()

if (PFFFT_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

if (PFFFT_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()
