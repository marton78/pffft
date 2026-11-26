# bench_fft, the FFT benchmark: pffft, other FFT libraries and other builds of
# pffft, measured side by side in one process. See README.md in this directory.
#
# engines.c holds the adapters to all the libraries and is compiled once per
# precision. A HAVE_<LIBRARY> definition tells it that a library is there.

set(bench_precisions "")
if (PFFFT_USE_TYPE_FLOAT)
  list(APPEND bench_precisions F32)
endif()
if (PFFFT_USE_TYPE_DOUBLE)
  list(APPEND bench_precisions F64)
endif()

# filled in by the libraries below
set(bench_defs "")      # definitions for the engines of every precision
set(bench_defs_F32 "")  # definitions for the engines of one precision
set(bench_defs_F64 "")
set(bench_incs "")      # include directories for the engines
set(bench_libs PFFFT)   # libraries of bench_fft; those that are CMake targets
                        # also give the engines their include directories

######################################################

if (PFFFT_USE_FFTPACK)
  # not linked as a library: engines.c includes src/fftpack.c, to give the
  # functions of each precision their own names
  list(APPEND bench_defs HAVE_FFTPACK)
endif()

# Green, Kiss, Pocket and Ooura's FFT are git submodules: git submodule update --init
# Their CMakeLists.txt ask for CMake 2.8, which CMake 4 rejects, so they get
# the policies of this project instead.
if (NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
  set(CMAKE_POLICY_VERSION_MINIMUM 3.10)
endif()

# bench_submodule(<directory> <library target> <definition>)
macro(bench_submodule dir lib definition)
  if (EXISTS "${PFFFT_SOURCE_DIR}/${dir}/CMakeLists.txt")
    add_subdirectory("${PFFFT_SOURCE_DIR}/${dir}" "${CMAKE_CURRENT_BINARY_DIR}/${dir}")
    # the same -march / /arch as pffft (TARGET_C_ARCH), for a fair comparison
    target_set_c_arch_flags(${lib})
    if (MSVC)
      # Their CMakeLists.txt are written for GCC and Clang: MSVC rejects their
      # warning options (-Werror) and has no libm to link.
      set_property(TARGET ${lib} PROPERTY COMPILE_OPTIONS "")
      foreach (prop LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
        get_target_property(deps ${lib} ${prop})
        if (deps)
          list(REMOVE_ITEM deps m "$<LINK_ONLY:m>")
          set_property(TARGET ${lib} PROPERTY ${prop} "${deps}")
        endif()
      endforeach()
    endif()
    list(APPEND bench_defs ${definition})
    list(APPEND bench_libs ${lib})
  else()
    message(WARNING "${lib} not found in subdir ${dir}: fetch it with 'git submodule update --init'")
  endif()
endmacro()

if (PFFFT_USE_BENCH_GREEN AND PFFFT_USE_TYPE_FLOAT AND WIN32)
  # It crashes there (MSVC, x64 and ARM64), probably through its use of long,
  # which has 32 bits on Windows also on 64-bit systems.
  message(STATUS "Green FFT is left out of bench_fft on Windows, where it crashes")
elseif (PFFFT_USE_BENCH_GREEN AND PFFFT_USE_TYPE_FLOAT)  # single precision only
  # only its library: its test programs do not compile with current compilers
  # and its tests need Octave
  bench_submodule(greenffts/src GreenFFT HAVE_GREEN_FFTS)
endif()
if (PFFFT_USE_BENCH_KISS AND PFFFT_USE_TYPE_FLOAT)  # built for single precision
  bench_submodule(kissfft KissFFT HAVE_KISS_FFT)
endif()
if (PFFFT_USE_BENCH_POCKET)
  bench_submodule(pocketfft PocketFFT HAVE_POCKET_FFT)
endif()
if (PFFFT_USE_BENCH_OOURA)
  # not built with its own CMakeLists.txt, which makes many variants, tests and
  # flags of its own: engines.c includes ooura/fftsg.c, like FFTPACK
  if (EXISTS "${PFFFT_SOURCE_DIR}/ooura/fftsg.c")
    list(APPEND bench_defs HAVE_OOURA_FFT)
    list(APPEND bench_incs "${PFFFT_SOURCE_DIR}/ooura/include")
  else()
    message(WARNING "Ooura's FFT not found in subdir ooura: fetch it with 'git submodule update --init'")
  endif()
endif()

if (PFFFT_USE_BENCH_FFTW)
  # Never a system FFTW: how a package was built (with SIMD or without) would
  # decide how FFTW compares.
  set(FFTW3_ROOT "" CACHE PATH "FFTW3 install prefix with include/ and lib/, made by benchmarks/build_fftw.py")
  if (NOT FFTW3_ROOT)
    message(FATAL_ERROR "PFFFT_USE_BENCH_FFTW needs FFTW3_ROOT: python3 benchmarks/build_fftw.py --prefix <dir>, then -DFFTW3_ROOT=<dir>")
  endif()
  if (NOT "${FFTW3_ROOT}" STREQUAL "${BENCH_FFTW3_ROOT_SEARCHED}")
    # search again: the results in the cache belong to another FFTW3_ROOT
    unset(FFTW3_INCLUDE_DIR CACHE)
    unset(FFTW3_LIBRARY_F32 CACHE)
    unset(FFTW3_LIBRARY_F64 CACHE)
    set(BENCH_FFTW3_ROOT_SEARCHED "${FFTW3_ROOT}" CACHE INTERNAL "FFTW3_ROOT of the cached FFTW3 search results")
  endif()
  # a cross build has the prefix outside its sysroot
  set(bench_fftw_where PATHS "${FFTW3_ROOT}" NO_DEFAULT_PATH NO_CMAKE_FIND_ROOT_PATH)
  find_path(FFTW3_INCLUDE_DIR fftw3.h ${bench_fftw_where} PATH_SUFFIXES include)
  find_library(FFTW3_LIBRARY_F32 fftw3f ${bench_fftw_where} PATH_SUFFIXES lib)
  find_library(FFTW3_LIBRARY_F64 fftw3 ${bench_fftw_where} PATH_SUFFIXES lib)
  foreach (prec ${bench_precisions})
    if (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY_${prec})
      list(APPEND bench_defs_${prec} HAVE_FFTW)
      list(APPEND bench_libs "${FFTW3_LIBRARY_${prec}}")
      list(APPEND bench_incs "${FFTW3_INCLUDE_DIR}")
    else()
      message(WARNING "PFFFT_USE_BENCH_FFTW: FFTW3_ROOT=${FFTW3_ROOT} has no FFTW3 for ${prec} (fftw3.h, and libfftw3f for F32 or libfftw3 for F64), so bench_fft gets no fftw engines for ${prec}. python3 benchmarks/build_fftw.py --prefix <dir> builds both precisions.")
    endif()
  endforeach()
endif()

if (APPLE AND PFFFT_USE_BENCH_VDSP)
  list(APPEND bench_defs HAVE_VECLIB)
  list(APPEND bench_libs "-framework Accelerate")
endif()

if (PFFFT_USE_BENCH_MKL)
  if (NOT PFFFT_TARGET_IS_X86)
    message(WARNING "using Intel MKL on '${PFFFT_TARGET_PROCESSOR}' might fail.")
  endif()
  message(STATUS "In case compiling/linking with Intel MKL fails, check benchmarks/CMakeLists.txt or deactivate PFFFT_USE_BENCH_MKL")
  list(APPEND bench_defs HAVE_MKL)
  list(APPEND bench_libs mkl_intel_lp64 mkl_sequential -lmkl_core)
endif()

######################################################

# what the other C targets of the project get, too, and the same -march /
# /arch as pffft (TARGET_C_ARCH), which the CSV header records
function(bench_target_setup target)
  target_set_c_arch_flags(${target})
  target_compile_definitions(${target} PRIVATE "BENCH_TARGET_ARCH=\"${TARGET_C_ARCH}\"")
  target_compile_definitions(${target} PRIVATE _USE_MATH_DEFINES)
  target_activate_c_compiler_warnings(${target})
  if (PFFFT_USE_DEBUG_ASAN)
    target_compile_options(${target} PRIVATE "-fsanitize=address")
  endif()
endfunction()

set(bench_objects "")
foreach (prec ${bench_precisions})
  add_library(bench_engines_${prec} OBJECT engines.c)
  bench_target_setup(bench_engines_${prec})
  target_compile_definitions(bench_engines_${prec} PRIVATE BENCH_PREC_${prec} ${bench_defs} ${bench_defs_${prec}})
  # Before CMake 3.12 an object library cannot link to libraries, so it does
  # not get their include directories and definitions that way.
  foreach (lib ${bench_libs})
    if (TARGET "${lib}")
      target_include_directories(bench_engines_${prec} PRIVATE $<TARGET_PROPERTY:${lib},INTERFACE_INCLUDE_DIRECTORIES>)
      target_compile_definitions(bench_engines_${prec} PRIVATE $<TARGET_PROPERTY:${lib},INTERFACE_COMPILE_DEFINITIONS>)
    endif()
  endforeach()
  # last, so that nothing in a system include directory hides a header from above
  target_include_directories(bench_engines_${prec} PRIVATE ${bench_incs})
  list(APPEND bench_objects $<TARGET_OBJECTS:bench_engines_${prec}>)
endforeach()

add_executable(bench_fft bench_fft.c platform.c accuracy.c bench.h ${bench_objects})
bench_target_setup(bench_fft)
foreach (prec ${bench_precisions})
  target_compile_definitions(bench_fft PRIVATE BENCH_HAVE_${prec})
endforeach()
target_link_libraries(bench_fft ${bench_libs} ${ASANLIB} ${MATHLIB})
if (APPLE OR CMAKE_SYSTEM_NAME MATCHES "Linux|Android")
  # where platform.c uses dlopen; elsewhere -ldl may not exist (newlib)
  target_link_libraries(bench_fft ${CMAKE_DL_LIBS})
endif()

######################################################

# bench_pf_mixer and bench_pf_conv: the benchmarks of the mixers (PFDSP) and
# of the fast convolution (PFFASTCONV)

if (PFFFT_USE_TYPE_FLOAT)

    add_executable(bench_pf_mixer_float   bench_mixers.cpp papi_perf_counter.h)
    target_compile_definitions(bench_pf_mixer_float PRIVATE _USE_MATH_DEFINES)
    target_compile_definitions(bench_pf_mixer_float PRIVATE PFFFT_ENABLE_FLOAT)
    target_link_libraries( bench_pf_mixer_float  ${ASANLIB} )
    if (PFFFT_USE_DEBUG_ASAN)
      target_compile_options(bench_pf_mixer_float PRIVATE "-fsanitize=address")
    endif()
    if (PAPI_FOUND)
        target_compile_definitions(bench_pf_mixer_float PRIVATE HAVE_PAPI=1)
        target_link_libraries(bench_pf_mixer_float ${PAPI_LIBRARIES})
    endif()
    target_link_libraries( bench_pf_mixer_float  PFDSP $<$<CXX_COMPILER_ID:GNU>:stdc++> )


  ############################################################################

  include(${PFFFT_SOURCE_DIR}/benchmarks/pf_conv.cmake)

  ############################################################################

  add_executable(bench_pf_conv_float   bench_conv.cpp papi_perf_counter.h)
  set_property(TARGET bench_pf_conv_float PROPERTY CXX_STANDARD 11)
  set_property(TARGET bench_pf_conv_float PROPERTY CXX_STANDARD_REQUIRED ON)
  target_compile_definitions(bench_pf_conv_float PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(bench_pf_conv_float PRIVATE PFFFT_ENABLE_FLOAT)
  if (PFFFT_USE_DEBUG_ASAN)
      target_compile_options(bench_pf_conv_float PRIVATE "-fsanitize=address")
  endif()
  target_link_libraries( bench_pf_conv_float  ${ASANLIB} )
  if (PAPI_FOUND)
      target_compile_definitions(bench_pf_conv_float PRIVATE HAVE_PAPI=1)
      target_link_libraries(bench_pf_conv_float ${PAPI_LIBRARIES})
  endif()
  # PFFFT for pffft::AlignedVector
  target_link_libraries( bench_pf_conv_float  pf_conv_dispatcher PFFFT $<$<CXX_COMPILER_ID:GNU>:stdc++> )

endif()
