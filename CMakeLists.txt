cmake_minimum_required(VERSION 3.0)
set(CMAKE_C_STANDARD 99)
project(PRETTY_FAST_FFT LANGUAGES C)

########################################################################
# select the release build type by default to get optimization flags
########################################################################
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
   message(STATUS "Build type not specified: defaulting to release.")
endif(NOT CMAKE_BUILD_TYPE)


add_library(PFFFT STATIC pffft.c)
target_link_libraries(PFFFT m)
set_property(TARGET PFFFT APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_library(FFTPACK STATIC fftpack.c)
target_link_libraries(FFTPACK m)
set_property(TARGET FFTPACK APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
)


add_executable(bench_pffft   test_pffft.c )
target_link_libraries(bench_pffft  PFFFT FFTPACK)
