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

  if (PFFFT_BUILD_BENCHMARKS)
    add_test(NAME bench_pffft_pow2
      COMMAND bench_pffft_float "--max-len" "128" "--quick"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_test(NAME bench_pffft_non2
      COMMAND bench_pffft_float "--non-pow2" "--max-len" "192" "--quick"
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()

endif()
