if (PFFFT_USE_FFTPACK)
  add_executable(test_fftpack_float ${PFFFT_SOURCE_DIR}/src/fftpack.c ${PFFFT_SOURCE_DIR}/src/fftpack.h)
  target_compile_definitions(test_fftpack_float PRIVATE _USE_MATH_DEFINES TESTING_FFTPACK)
  target_link_libraries(test_fftpack_float ${MATHLIB})

  add_executable(test_fftpack_double ${PFFFT_SOURCE_DIR}/src/fftpack.c ${PFFFT_SOURCE_DIR}/src/fftpack.h)
  target_compile_definitions(test_fftpack_double PRIVATE _USE_MATH_DEFINES FFTPACK_DOUBLE_PRECISION TESTING_FFTPACK)
  target_link_libraries(test_fftpack_double ${MATHLIB})
endif()

######################################################

if (PFFFT_USE_TYPE_FLOAT)
  add_executable( test_pffft_float  test_pffft.c )
  target_compile_definitions(test_pffft_float PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(test_pffft_float PRIVATE PFFFT_ENABLE_FLOAT)
  target_link_libraries( test_pffft_float  PFFFT ${ASANLIB} )
endif()

######################################################

if (PFFFT_USE_TYPE_DOUBLE)
  add_executable( test_pffft_double  test_pffft.c )
  target_compile_definitions(test_pffft_double PRIVATE _USE_MATH_DEFINES)
  target_compile_definitions(test_pffft_double PRIVATE PFFFT_ENABLE_DOUBLE)
  target_link_libraries( test_pffft_double  PFFFT ${ASANLIB} )
endif()

######################################################

add_executable( test_fft_factors  test_fft_factors.c )
if (PFFFT_USE_TYPE_FLOAT)
  target_compile_definitions(test_fft_factors PRIVATE PFFFT_ENABLE_FLOAT)
endif()
if (PFFFT_USE_TYPE_DOUBLE)
  target_compile_definitions(test_fft_factors PRIVATE PFFFT_ENABLE_DOUBLE)
endif()
target_link_libraries(test_fft_factors PFFFT ${ASANLIB} ${MATHLIB})

######################################################

add_executable( test_pffft_cpp test_pffft.cpp )
target_compile_definitions(test_pffft_cpp PRIVATE _USE_MATH_DEFINES)
if (PFFFT_USE_TYPE_FLOAT)
  target_compile_definitions(test_pffft_cpp PRIVATE PFFFT_ENABLE_FLOAT)
endif()
if (PFFFT_USE_TYPE_DOUBLE)
  target_compile_definitions(test_pffft_cpp PRIVATE PFFFT_ENABLE_DOUBLE)
endif()
target_link_libraries( test_pffft_cpp  PFFFT ${STDCXXLIB} ${ASANLIB} )
if (EMSCRIPTEN)
  # the large-size tests overflow Emscripten's default 64 KB stack
  target_link_options(test_pffft_cpp PRIVATE "-sSTACK_SIZE=1MB")
endif()

######################################################

add_executable( test_pffft_cpp_11 test_pffft.cpp )
target_compile_definitions(test_pffft_cpp_11 PRIVATE _USE_MATH_DEFINES)
if (PFFFT_USE_TYPE_FLOAT)
  target_compile_definitions(test_pffft_cpp_11 PRIVATE PFFFT_ENABLE_FLOAT)
endif()
if (PFFFT_USE_TYPE_DOUBLE)
  target_compile_definitions(test_pffft_cpp_11 PRIVATE PFFFT_ENABLE_DOUBLE)
endif()
target_link_libraries( test_pffft_cpp_11  PFFFT ${STDCXXLIB} ${ASANLIB} )
if (EMSCRIPTEN)
  # the large-size tests overflow Emscripten's default 64 KB stack
  target_link_options(test_pffft_cpp_11 PRIVATE "-sSTACK_SIZE=1MB")
endif()

set_property(TARGET test_pffft_cpp_11 PROPERTY CXX_STANDARD 11)
set_property(TARGET test_pffft_cpp_11 PROPERTY CXX_STANDARD_REQUIRED ON)

######################################################

if (PFFFT_USE_TYPE_FLOAT)
  add_executable(test_pffastconv   test_pffastconv.c
    ${SIMD_FLOAT_HDRS} ${SIMD_DOUBLE_HDRS}
  )
  target_compile_definitions(test_pffastconv PRIVATE _USE_MATH_DEFINES)
  if (PFFFT_USE_DEBUG_ASAN)
    target_compile_options(test_pffastconv PRIVATE "-fsanitize=address")
  endif()
  target_set_c_arch_flags(test_pffastconv)
  if (NOT PFFFT_USE_SIMD)
    target_compile_definitions(test_pffastconv PRIVATE PFFFT_SIMD_DISABLE=1)
  endif()
  target_include_directories(test_pffastconv PRIVATE ${PFFFT_SOURCE_DIR}/src)
  target_link_libraries( test_pffastconv  PFFASTCONV ${ASANLIB} ${MATHLIB} )
  if (EMSCRIPTEN)
    target_link_options(test_pffastconv PRIVATE "-sALLOW_MEMORY_GROWTH=1")
  endif()
endif()

######################################################

if (PFFFT_USE_TYPE_FLOAT)
  add_test(NAME test_pffft_float
    COMMAND test_pffft_float
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  # self-test of the SIMD macros; passes trivially when SIMD is compiled out
  add_test(NAME test_simd_float
    COMMAND test_pffft_float --test-simd
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

if (PFFFT_USE_TYPE_DOUBLE)
  add_test(NAME test_pffft_double
    COMMAND test_pffft_double
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_test(NAME test_simd_double
    COMMAND test_pffft_double --test-simd
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

add_test(NAME test_pffft_cpp
  COMMAND test_pffft_cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_pffft_cpp_11
  COMMAND test_pffft_cpp_11
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_test(NAME test_fft_factors
  COMMAND test_fft_factors
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

if (PFFFT_USE_FFTPACK)
  add_test(NAME test_fftpack_float
    COMMAND test_fftpack_float
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(NAME test_fftpack_double
    COMMAND test_fftpack_double
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()


if (PFFFT_USE_TYPE_FLOAT)

  add_test(NAME test_pfconv_lens_symetric
    COMMAND test_pffastconv "--no-bench" "--quick" "--sym"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(NAME test_pfconv_lens_non_sym
    COMMAND test_pffastconv "--no-bench" "--quick"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(NAME bench_pfconv_symetric
    COMMAND test_pffastconv "--no-len" "--quick" "--sym"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(NAME bench_pfconv_non_sym
    COMMAND test_pffastconv "--no-len" "--quick"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

endif()

if (PFFFT_USE_TYPE_FLOAT)
  # the "slow" convolution kernels of bench_pf_conv_float, every architecture
  # variant of benchmarks/pf_conv.cpp against a plain convolution; built here
  # also without the benchmarks, for the SIMD backends of every CI job
  include(${PFFFT_SOURCE_DIR}/benchmarks/pf_conv.cmake)
  add_executable(test_pf_conv test_pf_conv.cpp)
  set_property(TARGET test_pf_conv PROPERTY CXX_STANDARD 11)
  set_property(TARGET test_pf_conv PROPERTY CXX_STANDARD_REQUIRED ON)
  target_compile_definitions(test_pf_conv PRIVATE PFFFT_ENABLE_FLOAT)
  target_include_directories(test_pf_conv PRIVATE ${PFFFT_SOURCE_DIR}/benchmarks)
  target_activate_cxx_compiler_warnings(test_pf_conv)
  # Is pf_conv_arch_dflt (C++) compiled for the same target as PFFFT (C)?
  # Then it must get the same SIMD backend. That needs the same arch options
  # and the same flags, of which MSVC gives C++ /GR and /EHsc on top.
  set(pf_conv_same_target TRUE)
  if (NOT (("${TARGET_C_ARCH}" STREQUAL "${TARGET_CXX_ARCH}") AND ("${TARGET_C_EXTRA}" STREQUAL "${TARGET_CXX_EXTRA}")))
    set(pf_conv_same_target FALSE)
  endif()
  foreach (cfg "" _DEBUG _RELEASE _RELWITHDEBINFO _MINSIZEREL)
    # as lists of options, which drops the spaces between them
    string(REGEX REPLACE "[ \t]+" ";" c_flags "${CMAKE_C_FLAGS${cfg}}")
    string(REGEX REPLACE "[ \t]+" ";" cxx_flags "${CMAKE_CXX_FLAGS${cfg}}")
    list(REMOVE_ITEM c_flags "")
    list(REMOVE_ITEM cxx_flags "")
    if (MSVC AND cxx_flags)
      list(REMOVE_ITEM cxx_flags /GR -GR /EHsc -EHsc)
    endif()
    if (NOT "${c_flags}" STREQUAL "${cxx_flags}")
      set(pf_conv_same_target FALSE)
    endif()
  endforeach()
  if (pf_conv_same_target)
    target_compile_definitions(test_pf_conv PRIVATE PF_CONV_DFLT_IS_PFFFT_BACKEND=1)
  endif()
  if (PFFFT_USE_DEBUG_ASAN)
    target_compile_options(test_pf_conv PRIVATE "-fsanitize=address")
  endif()
  target_link_libraries(test_pf_conv pf_conv_dispatcher PFFFT ${STDCXXLIB} ${ASANLIB} ${MATHLIB})

  add_test(NAME test_pf_conv
    COMMAND test_pf_conv
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

if (PFFFT_BUILD_BENCHMARKS)
  # Smoke tests of the FFT benchmark: every engine that is built in checks its
  # result and is timed briefly, at powers of two and at another size. One
  # test prints the table, the other the CSV. FFTW_MEASURE is limited to one
  # size: an FFTW without a cycle counter (as Debian's on arm64) takes seconds
  # for every plan.
  add_test(NAME bench_fft_smoke
    COMMAND bench_fft --sizes 64,96,128 --rounds 10 --fftw-measure-max 64
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_test(NAME bench_fft_csv
    COMMAND bench_fft --sizes 64,96 --rounds 10 --fftw-measure-max 0 --csv -
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()
