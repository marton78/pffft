# pf_conv: the "slow" convolution kernels of bench_pf_conv_float, compiled once
# per architecture option (pf_conv_arch_<opt>) and collected by
# pf_conv_dispatcher. Included by tests/ for test_pf_conv, so that every CI
# job that builds the tests checks the kernels on its SIMD backend, and by
# benchmarks/; the targets are created once.
include_guard(GLOBAL)

set(_pf_conv_dir ${PFFFT_SOURCE_DIR}/benchmarks)

# pf_conv.cpp is compiled once per architecture option, as
# pf_conv_arch_<opt>, with the SIMD macros of pffft (src/simd/pf_float.h)
# under the same definitions as PFFFT. pf_conv_arch_none is always scalar.
# The libraries are static also with BUILD_SHARED_LIBS: they export nothing,
# so MSVC would make no import library for them.
function(pf_conv_arch_target target arch_post)
  add_library(${target} STATIC ${_pf_conv_dir}/pf_conv.cpp ${_pf_conv_dir}/pf_conv.h ${PFFFT_SOURCE_DIR}/include/pffft/pf_cplx.h ${SIMD_FLOAT_HDRS})
  target_include_directories(${target} PRIVATE ${PFFFT_SOURCE_DIR}/include ${PFFFT_SOURCE_DIR}/src)
  target_compile_definitions(${target} PRIVATE CONV_ARCH_POST=${arch_post})
  if ("${arch_post}" STREQUAL "none" OR NOT PFFFT_USE_SIMD)
    target_compile_definitions(${target} PRIVATE PFFFT_SIMD_DISABLE=1)
  endif()
  if (PFFFT_USE_SCALAR_VECT AND NOT "${arch_post}" STREQUAL "none")
    target_compile_definitions(${target} PRIVATE PFFFT_SCALVEC_ENABLED=1)
  endif()
  if (PFFFT_USE_HELIUM)
    target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_HELIUM=1)
  endif()
  set_property(TARGET ${target} PROPERTY CXX_STANDARD 11)
  set_property(TARGET ${target} PROPERTY CXX_STANDARD_REQUIRED ON)
  target_activate_cxx_compiler_warnings(${target})
  if (PFFFT_USE_DEBUG_ASAN)
    target_compile_options(${target} PRIVATE "-fsanitize=address")
    target_link_libraries(${target} ${ASANLIB})
  endif()
endfunction()

pf_conv_arch_target(pf_conv_arch_none none)
add_library(pf_conv_dispatcher STATIC ${_pf_conv_dir}/pf_conv_dispatcher.cpp ${_pf_conv_dir}/pf_conv_dispatcher.h ${_pf_conv_dir}/pf_conv.h ${PFFFT_SOURCE_DIR}/include/pffft/pf_cplx.h)
target_include_directories(pf_conv_dispatcher PRIVATE ${PFFFT_SOURCE_DIR}/include)
set_property(TARGET pf_conv_dispatcher PROPERTY CXX_STANDARD 11)
set_property(TARGET pf_conv_dispatcher PROPERTY CXX_STANDARD_REQUIRED ON)
target_activate_cxx_compiler_warnings(pf_conv_dispatcher)

pf_conv_arch_target(pf_conv_arch_dflt dflt)
target_set_cxx_arch_flags(pf_conv_arch_dflt)

target_link_libraries(pf_conv_dispatcher pf_conv_arch_none pf_conv_arch_dflt)

if (PFFFT_TARGET_IS_X86 AND CMAKE_SIZEOF_VOID_P EQUAL 8)

  # The kernels are 4 wide (SSE) on every x86. Of the variants for core2,
  # nehalem and sandybridge, none was faster than avx2 (haswell: FMA) on the
  # x86 CPUs measured, and with gcc core2 and nehalem were slower than dflt,
  # the SSE2 baseline of x86-64. So avx2 is the only one besides dflt.
  if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
      set(PF_CONV_ARCHES "avx2")
      set(PF_CONV_OPT_avx2 "haswell")  # emulate a map
      target_compile_definitions(pf_conv_dispatcher PRIVATE CONV_ARCH_GCC_AMD64)
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
      set(PF_CONV_ARCHES "avx2")
      set(PF_CONV_OPT_avx2 "AVX2")  # emulate a map
      target_compile_definitions(pf_conv_dispatcher PRIVATE CONV_ARCH_MSVC_AMD64)
  else()
      set(PF_CONV_ARCHES "")
      message(WARNING "unknown compiler ${CMAKE_CXX_COMPILER_ID} on target processor ${PFFFT_TARGET_PROCESSOR}: can't do architecture specific compilation")
  endif()

elseif (PFFFT_TARGET_PROCESSOR MATCHES "^(aarch64|arm64)$"
        AND CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")

    set(PF_CONV_ARCHES "armv8a")
    set(PF_CONV_OPT_armv8a   "armv8-a")  # emulate a map for arch

    target_compile_definitions(pf_conv_dispatcher PRIVATE CONV_ARCH_GCC_AARCH64)

elseif (PFFFT_TARGET_PROCESSOR STREQUAL "armv7l")

  if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(PF_CONV_ARCHES "neon_vfpv4;neon_rpi3_a53;neon_rpi4_a72")
      set(PF_CONV_OPT_neon_vfpv4        "armv7-a")    # emulate a map for arch
      set(PF_CONV_EXTRA_neon_vfpv4      "neon_vfpv4") # emulate a map for additional options (EXTRA)
      set(PF_CONV_OPT_neon_rpi3_a53     "armv7-a")
      set(PF_CONV_EXTRA_neon_rpi3_a53   "neon_rpi3_a53")
      set(PF_CONV_OPT_neon_rpi4_a72     "armv7-a")
      set(PF_CONV_EXTRA_neon_rpi4_a72   "neon_rpi4_a72")

      target_compile_definitions(pf_conv_dispatcher PRIVATE CONV_ARCH_GCC_ARM32NEON)
  else()
      set(PF_CONV_ARCHES "")
      message(WARNING "unknown compiler ${CMAKE_CXX_COMPILER_ID} on target processor ${PFFFT_TARGET_PROCESSOR}: can't do architecture specific compilation")
  endif()

else()
    message(STATUS "pf_conv: no architecture-specific variants for target processor '${PFFFT_TARGET_PROCESSOR}', only none and dflt")
endif()

foreach (arch_opt ${PF_CONV_ARCHES})
    pf_conv_arch_target(pf_conv_arch_${arch_opt} ${arch_opt})
    target_set_cxx_arch_option(pf_conv_arch_${arch_opt} "${PF_CONV_OPT_${arch_opt}}" "${PF_CONV_EXTRA_${arch_opt}}"  "${PF_CONV_OPT_${arch_opt}}")
    target_link_libraries(pf_conv_dispatcher  pf_conv_arch_${arch_opt})
    message(STATUS "added library pf_conv_arch_${arch_opt}  with CONV_ARCH_POST=${arch_opt}")
endforeach()

if (PFFFT_USE_DEBUG_ASAN)
    target_compile_options(pf_conv_dispatcher  PRIVATE "-fsanitize=address")
    target_link_libraries(pf_conv_dispatcher ${ASANLIB})
endif()
