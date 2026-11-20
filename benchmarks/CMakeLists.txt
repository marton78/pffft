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

  add_library(pf_conv_arch_none pf_conv.cpp pf_conv.h ${PFFFT_SOURCE_DIR}/include/pffft/pf_cplx.h)
  target_include_directories(pf_conv_arch_none PRIVATE ${PFFFT_SOURCE_DIR}/include)
  target_compile_definitions(pf_conv_arch_none PRIVATE CONV_ARCH_POST=none MIPP_NO_INTRINSICS=1)
  set_property(TARGET pf_conv_arch_none PROPERTY CXX_STANDARD 11)
  set_property(TARGET pf_conv_arch_none PROPERTY CXX_STANDARD_REQUIRED ON)
  target_activate_cxx_compiler_warnings(pf_conv_arch_none)
  add_library(pf_conv_dispatcher  pf_conv_dispatcher.cpp pf_conv_dispatcher.h pf_conv.h ${PFFFT_SOURCE_DIR}/include/pffft/pf_cplx.h)
  target_include_directories(pf_conv_dispatcher PRIVATE ${PFFFT_SOURCE_DIR}/include)
  set_property(TARGET pf_conv_dispatcher PROPERTY CXX_STANDARD 11)
  set_property(TARGET pf_conv_dispatcher PROPERTY CXX_STANDARD_REQUIRED ON)
  target_activate_cxx_compiler_warnings(pf_conv_dispatcher)

  add_library(pf_conv_arch_dflt pf_conv.cpp pf_conv.h ${PFFFT_SOURCE_DIR}/include/pffft/pf_cplx.h)
  target_include_directories(pf_conv_arch_dflt PRIVATE ${PFFFT_SOURCE_DIR}/include)
  target_compile_definitions(pf_conv_arch_dflt PRIVATE CONV_ARCH_POST=dflt)
  set_property(TARGET pf_conv_arch_dflt PROPERTY CXX_STANDARD 11)
  set_property(TARGET pf_conv_arch_dflt PROPERTY CXX_STANDARD_REQUIRED ON)
  target_activate_cxx_compiler_warnings(pf_conv_arch_dflt)
  target_set_cxx_arch_flags(pf_conv_arch_dflt)

  target_link_libraries(pf_conv_dispatcher pf_conv_arch_none pf_conv_arch_dflt)

  if (PFFFT_TARGET_IS_X86 AND CMAKE_SIZEOF_VOID_P EQUAL 8)

    if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(PF_CONV_ARCHES "sse3;sse4;avx;avx2")
        set(PF_CONV_OPT_sse3 "core2")  # emulate a map
        set(PF_CONV_OPT_sse4 "nehalem")
        set(PF_CONV_OPT_avx  "sandybridge")
        set(PF_CONV_OPT_avx2 "haswell")
        target_compile_definitions(pf_conv_dispatcher PRIVATE CONV_ARCH_GCC_AMD64)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        set(PF_CONV_ARCHES "sse2;avx;avx2")
        set(PF_CONV_OPT_sse2 "SSE2")  # emulate a map
        set(PF_CONV_OPT_avx  "AVX")
        set(PF_CONV_OPT_avx2 "AVX2")
        target_compile_definitions(pf_conv_dispatcher PRIVATE CONV_ARCH_MSVC_AMD64)
    else()
        set(PF_CONV_ARCHES "")
        message(WARNING "unknown compiler ${CMAKE_CXX_COMPILER_ID} on target processor ${PFFFT_TARGET_PROCESSOR}: can't do architecture specific compilation")
    endif()

  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")

      if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
          set(PF_CONV_ARCHES "armv8a")
          set(PF_CONV_OPT_armv8a   "armv8-a")  # emulate a map for arch

          target_compile_definitions(pf_conv_dispatcher PRIVATE CONV_ARCH_GCC_AARCH64)
      else()
          set(PF_CONV_ARCHES "")
          message(WARNING "unknown compiler ${CMAKE_CXX_COMPILER_ID} on CMAKE_SYSTEM_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR}: can't do architecture specific compilation")
      endif()

  elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "armv7l")

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
        message(WARNING "unknown compiler ${CMAKE_CXX_COMPILER_ID} on CMAKE_SYSTEM_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR}: can't do architecture specific compilation")
    endif()

  else()
      message(STATUS "bench_pf_conv_float: no architecture-specific variants for target processor '${PFFFT_TARGET_PROCESSOR}', benchmarking the default build only")
  endif()

  foreach (arch_opt ${PF_CONV_ARCHES})
      add_library(pf_conv_arch_${arch_opt} pf_conv.cpp pf_conv.h ${PFFFT_SOURCE_DIR}/include/pffft/pf_cplx.h)
      target_include_directories(pf_conv_arch_${arch_opt} PRIVATE ${PFFFT_SOURCE_DIR}/include)
      set_property(TARGET pf_conv_arch_${arch_opt} PROPERTY CXX_STANDARD 11)
      set_property(TARGET pf_conv_arch_${arch_opt} PROPERTY CXX_STANDARD_REQUIRED ON)
      target_activate_cxx_compiler_warnings(pf_conv_arch_${arch_opt})
      target_compile_definitions(pf_conv_arch_${arch_opt} PRIVATE CONV_ARCH_POST=${arch_opt})

      target_set_cxx_arch_option(pf_conv_arch_${arch_opt} "${PF_CONV_OPT_${arch_opt}}" "${PF_CONV_EXTRA_${arch_opt}}"  "${PF_CONV_OPT_${arch_opt}}")
      target_link_libraries(pf_conv_dispatcher  pf_conv_arch_${arch_opt})
      message(STATUS "added library pf_conv_arch_${arch_opt}  with CONV_ARCH_POST=${arch_opt}")
  endforeach()

  if (PFFFT_USE_DEBUG_ASAN)
      foreach (arch_opt ${PF_CONV_ARCHES})
          target_compile_options(pf_conv_arch_${arch_opt} PRIVATE "-fsanitize=address")
          target_link_libraries( pf_conv_arch_${arch_opt} ${ASANLIB})
      endforeach()

      target_compile_options(pf_conv_arch_none  PRIVATE "-fsanitize=address")
      target_link_libraries( pf_conv_arch_none  ${ASANLIB})

      target_compile_options(pf_conv_dispatcher  PRIVATE "-fsanitize=address")
      target_link_libraries(pf_conv_dispatcher ${ASANLIB})
  endif()

  if(MIPP_FOUND)
      foreach (arch_opt ${PF_CONV_ARCHES})
          message(STATUS "link pf_conv_arch_${arch_opt} against MIPP")
          target_link_libraries(pf_conv_arch_${arch_opt} MIPP)
      endforeach()

      message(STATUS "link pf_conv_arch_none against MIPP")
      target_link_libraries(pf_conv_arch_none MIPP)
  endif()

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
  if(MIPP_FOUND)
      target_link_libraries(bench_pf_conv_float MIPP)
  endif()
  target_link_libraries( bench_pf_conv_float  pf_conv_dispatcher PFDSP $<$<CXX_COMPILER_ID:GNU>:stdc++> )

endif()
