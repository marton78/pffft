cmake_minimum_required(VERSION 3.1)
project(examples)

if ( CMAKE_C_COMPILER_ID MATCHES "MSVC" )
  # using Visual Studio C++
  message(STATUS "INFO: detected MSVC: will not link math lib m")
  set(MATHLIB "")
  add_definitions("/D_CRT_SECURE_NO_WARNINGS")
  set(MSVC_DISABLED_WARNINGS_LIST "C4996")
else()
  if(DISABLE_LINK_WITH_M)
  else()
    message(STATUS "INFO: detected NO MSVC: ${CMAKE_C_COMPILER_ID}: will link math lib m")
    set(MATHLIB "m")
  endif()
endif()


set(CMAKE_CXX_EXTENSIONS OFF)


add_executable(example_cpp11_real_dbl_fwd example_cpp11_real_dbl_fwd.cpp)
target_compile_definitions(example_cpp11_real_dbl_fwd PRIVATE PFFFT_ENABLE_DOUBLE PFFFT_ENABLE_FLOAT)
target_link_libraries(example_cpp11_real_dbl_fwd PFFFT ${MATHLIB})
set_property(TARGET example_cpp11_real_dbl_fwd PROPERTY CXX_STANDARD 11)
set_property(TARGET example_cpp11_real_dbl_fwd PROPERTY CXX_STANDARD_REQUIRED ON)


add_executable(example_cpp11_cplx_dbl_fwd example_cpp11_cplx_dbl_fwd.cpp)
target_compile_definitions(example_cpp11_cplx_dbl_fwd PRIVATE PFFFT_ENABLE_DOUBLE PFFFT_ENABLE_FLOAT)
target_link_libraries(example_cpp11_cplx_dbl_fwd PFFFT ${MATHLIB})
set_property(TARGET example_cpp11_cplx_dbl_fwd PROPERTY CXX_STANDARD 11)
set_property(TARGET example_cpp11_cplx_dbl_fwd PROPERTY CXX_STANDARD_REQUIRED ON)


add_executable(example_cpp98_real_flt_fwd example_cpp98_real_flt_fwd.cpp)
target_link_libraries(example_cpp98_real_flt_fwd PFFFT ${MATHLIB})
# float is default - when not setting any of PFFFT_ENABLE_DOUBLE or PFFFT_ENABLE_FLOAT
# target_compile_definitions(example_cpp98_cplx_flt_fwd PRIVATE PFFFT_ENABLE_FLOAT)
set_property(TARGET example_cpp98_real_flt_fwd PROPERTY CXX_STANDARD 98)
set_property(TARGET example_cpp98_real_flt_fwd PROPERTY CXX_STANDARD_REQUIRED ON)


add_executable(example_cpp98_cplx_flt_fwd example_cpp98_cplx_flt_fwd.cpp)
target_link_libraries(example_cpp98_cplx_flt_fwd PFFFT ${MATHLIB})
set_property(TARGET example_cpp98_cplx_flt_fwd PROPERTY CXX_STANDARD 98)
set_property(TARGET example_cpp98_cplx_flt_fwd PROPERTY CXX_STANDARD_REQUIRED ON)
