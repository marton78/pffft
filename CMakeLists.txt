cmake_minimum_required(VERSION 3.10)
project(PRETTY_FAST_FFT)

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

# what to install?
option(INSTALL_PFFFT      "install pffft to CMAKE_INSTALL_PREFIX?" ON)
option(INSTALL_PFDSP      "install pfdsp to CMAKE_INSTALL_PREFIX?" OFF)
option(INSTALL_PFFASTCONV "install pffastconv to CMAKE_INSTALL_PREFIX?" OFF)

# test options
option(PFFFT_USE_BENCH_FFTW   "use (system-installed) FFTW3 in fft benchmark?" OFF)
option(PFFFT_USE_BENCH_GREEN  "use Green FFT in fft benchmark? - if exists in subdir" ON)
option(PFFFT_USE_BENCH_KISS   "use KissFFT in fft benchmark? - if exists in subdir" ON)
option(PFFFT_USE_BENCH_POCKET "use PocketFFT in fft benchmark? - if exists in subdir" ON)
option(PFFFT_USE_BENCH_MKL    "use Intel MKL in fft benchmark? needs to be installed" OFF)
option(PFFFT_USE_BENCH_FFTS      "use FFTS in fft benchmark? - if exists in subdir" ON)
option(PFFFT_FFTS_NO_CODEGEN     "build FFTS without dynamic code generation (x86_64 only)?" OFF)
option(PFFFT_USE_FFTPACK      "compile and use FFTPACK in fft benchmark & validation?" ON)

option(PFFFT_USE_DEBUG_ASAN  "use GCC's address sanitizer?" OFF)

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
find_package(PAPI)
find_package(MIPP)
if (MIPP_FOUND)
# if (TARGET MIPP)
    message(STATUS "found MIPP")
else()
    message(STATUS "NOT found MIPP")
endif()


if (PFFFT_USE_DEBUG_ASAN)
  set(ASANLIB "asan")
else()
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
  message(STATUS "INFO: detected MSVC: will not link math lib m")
  set(MATHLIB "")

  add_definitions("/D_CRT_SECURE_NO_WARNINGS")

  set(MSVC_DISABLED_WARNINGS_LIST
      "C4996"
  )

else()
  if(PFFFT_DISABLE_LINK_WITH_M)
  else()
    message(STATUS "INFO: detected NO MSVC: ${CMAKE_C_COMPILER_ID}: will link math lib m")
    set(MATHLIB "m")
  endif()
endif()

set(STDCXXLIB "")
if (MINGW)
  set(STDCXXLIB "stdc++")
endif()


set( SIMD_FLOAT_HDRS ${PFFFT_SOURCE_DIR}/src/simd/pf_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_sse1_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_altivec_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_neon_float.h ${PFFFT_SOURCE_DIR}/src/simd/pf_scalar_float.h )
set( SIMD_DOUBLE_HDRS ${PFFFT_SOURCE_DIR}/src/simd/pf_double.h ${PFFFT_SOURCE_DIR}/src/simd/pf_avx_double.h ${PFFFT_SOURCE_DIR}/src/simd/pf_scalar_double.h )

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
if (NOT PFFFT_USE_SIMD)
  target_compile_definitions(PFFFT PRIVATE PFFFT_SIMD_DISABLE=1)
endif()
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

if ((PFFFT_BUILD_TESTS OR PFFFT_BUILD_BENCHMARKS) AND PFFFT_USE_FFTPACK)

  # float / single precision
  add_library(FFTPACK_FLOAT STATIC src/fftpack.c src/fftpack.h)
  target_compile_definitions(FFTPACK_FLOAT PRIVATE _USE_MATH_DEFINES)
  target_activate_c_compiler_warnings(FFTPACK_FLOAT)
  target_link_libraries( FFTPACK_FLOAT ${MATHLIB} )
  target_include_directories(FFTPACK_FLOAT
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )

  # double precision
  add_library(FFTPACK_DOUBLE STATIC src/fftpack.c src/fftpack.h)
  target_compile_definitions(FFTPACK_DOUBLE PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(FFTPACK_DOUBLE PUBLIC FFTPACK_DOUBLE_PRECISION)
  target_activate_c_compiler_warnings(FFTPACK_DOUBLE)
  target_link_libraries( FFTPACK_DOUBLE ${MATHLIB} )
  target_include_directories(FFTPACK_DOUBLE
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )

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

install( TARGETS ${INSTALL_TARGETS}
  EXPORT pffft-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install( FILES  ${INSTALL_HEADERS}  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pffft)

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
  VERSION 1.0
  COMPATIBILITY SameMajorVersion
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/pffft-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/pffft-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pffft
)

add_custom_target(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_SOURCE_DIR}/uninstall.cmake"
)

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
