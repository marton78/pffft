# Build FFmpeg's libavutil (which contains the AVTx FFT implementation)
# as a static library for benchmarking against pffft.
#
# Usage (standalone):
#   cmake -S benchmarks/ffmpeg_tx -B build/ffmpeg_tx \
#         -DFFMPEG_SOURCE_DIR=/path/to/FFmpeg
#   cmake --build build/ffmpeg_tx
#
# Usage (from parent CMakeLists.txt):
#   set(FFMPEG_SOURCE_DIR "/path/to/FFmpeg")
#   add_subdirectory(benchmarks/ffmpeg_tx)
#   # Then use the imported target:  ffmpeg_tx::avutil
#
# The resulting target exposes:
#   - libavutil.a  with full SIMD support (x86 SSE/AVX/AVX2/FMA3, AArch64 NEON)
#   - Include paths for <libavutil/tx.h> and generated headers (config.h, avconfig.h)
#
# Public API header of interest:
#   #include <libavutil/tx.h>    // av_tx_init(), av_tx_fn, AV_TX_FLOAT_FFT, …
#
# Requirements:
#   - FFmpeg source tree (git clone https://git.ffmpeg.org/ffmpeg.git)
#   - C compiler, GNU make
#   - nasm >= 2.13  (for x86/x86-64 SSE/AVX assembly; auto-detected by FFmpeg configure)
#   - On macOS: Xcode command line tools
#
# To force a reconfigure (e.g. after changing FFMPEG_SOURCE_DIR or the compiler):
#   rm -rf <build_dir>/ffmpeg_tx/ffmpeg_build

cmake_minimum_required(VERSION 3.16)

# Allow being used as a subdirectory without clobbering the parent project.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(ffmpeg_tx C)
endif()

# ──────────────────────────────────────────────────────────────────────────────
# Input validation
# ──────────────────────────────────────────────────────────────────────────────

if(NOT DEFINED FFMPEG_SOURCE_DIR)
  message(FATAL_ERROR
    "FFMPEG_SOURCE_DIR is not set.\n"
    "Pass it on the cmake command line:\n"
    "  cmake -DFFMPEG_SOURCE_DIR=/path/to/FFmpeg ...")
endif()

get_filename_component(FFMPEG_SOURCE_DIR "${FFMPEG_SOURCE_DIR}" ABSOLUTE)

if(NOT EXISTS "${FFMPEG_SOURCE_DIR}/configure")
  message(FATAL_ERROR
    "FFmpeg configure script not found at: ${FFMPEG_SOURCE_DIR}/configure\n"
    "Make sure FFMPEG_SOURCE_DIR points to the FFmpeg source tree root.")
endif()

# ──────────────────────────────────────────────────────────────────────────────
# Parallel build jobs
# ──────────────────────────────────────────────────────────────────────────────

include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
  set(NPROC 4)
endif()

# ──────────────────────────────────────────────────────────────────────────────
# Build paths
# ──────────────────────────────────────────────────────────────────────────────

set(FFMPEG_BUILD_DIR   "${CMAKE_CURRENT_BINARY_DIR}/ffmpeg_build")
set(FFMPEG_AVUTIL_A    "${FFMPEG_BUILD_DIR}/libavutil/libavutil.a")

# Headers live in both the source tree (public API) and the build dir (generated
# headers: config.h, libavutil/avconfig.h, libavutil/ffversion.h).
set(FFMPEG_INCLUDE_DIRS
  "${FFMPEG_SOURCE_DIR}"
  "${FFMPEG_BUILD_DIR}"
)

# ──────────────────────────────────────────────────────────────────────────────
# FFmpeg configure flags
#
# Goals:
#  - Fast configure (~10-15 s) and build (~5-10 s on modern hardware)
#  - No external library dependencies beyond a C compiler and make
#  - Full platform SIMD support:
#      x86-64: SSE2, SSE3, AVX, AVX2, FMA3  (requires nasm >= 2.13)
#      AArch64: NEON
#    All detected and enabled automatically by FFmpeg's configure script.
# ──────────────────────────────────────────────────────────────────────────────

set(FFMPEG_CONFIGURE_FLAGS
  --prefix=${CMAKE_CURRENT_BINARY_DIR}/ffmpeg_install

  # Disable all optional components (encoders, decoders, muxers, filters, …)
  --disable-everything

  # Disable libraries beyond libavutil
  --disable-avdevice
  --disable-avcodec
  --disable-avformat
  --disable-swscale
  --disable-swresample

  # Disable binary programs (ffmpeg, ffprobe, ffplay)
  --disable-programs

  # Disable documentation
  --disable-doc

  # Disable network support
  --disable-network

  # Static library only
  --enable-static
  --disable-shared

  # Optimised build (no debug symbols)
  --disable-debug

  # Optimise for the host CPU (e.g. -mcpu=apple-m1 on Apple Silicon)
  --extra-cflags=-mcpu=native
)

# ──────────────────────────────────────────────────────────────────────────────
# ExternalProject: run FFmpeg configure + make
#
# The configure script auto-detects the source path from its own location, so
# no --srcdir flag is required.  It writes config.h and a Makefile into
# FFMPEG_BUILD_DIR, then `make` compiles libavutil (the only enabled library).
# ──────────────────────────────────────────────────────────────────────────────

include(ExternalProject)

ExternalProject_Add(ffmpeg_avutil_build
  SOURCE_DIR   "${FFMPEG_SOURCE_DIR}"
  BINARY_DIR   "${FFMPEG_BUILD_DIR}"

  # FFmpeg uses its own autoconf-style configure, not CMake.
  CONFIGURE_COMMAND
    "${FFMPEG_SOURCE_DIR}/configure" ${FFMPEG_CONFIGURE_FLAGS}

  # Build just libavutil.  Because avcodec/avformat/… are disabled, `make`
  # only ever builds libavutil anyway; the explicit target makes intent clear.
  BUILD_COMMAND
    ${CMAKE_MAKE_PROGRAM} -j${NPROC} libavutil/libavutil.a

  # No install step – we reference the build tree directly via IMPORTED_LOCATION.
  INSTALL_COMMAND ""

  # Declare the byproduct so that CMake/Ninja can track the dependency correctly.
  BUILD_BYPRODUCTS "${FFMPEG_AVUTIL_A}"
)

# ──────────────────────────────────────────────────────────────────────────────
# Exported imported target: ffmpeg_tx::avutil
# ──────────────────────────────────────────────────────────────────────────────

add_library(ffmpeg_tx_avutil STATIC IMPORTED GLOBAL)
add_dependencies(ffmpeg_tx_avutil ffmpeg_avutil_build)

set_target_properties(ffmpeg_tx_avutil PROPERTIES
  IMPORTED_LOCATION             "${FFMPEG_AVUTIL_A}"
  INTERFACE_INCLUDE_DIRECTORIES "${FFMPEG_INCLUDE_DIRS}"
)

# Namespace alias: target_link_libraries(my_exe ffmpeg_tx::avutil)
add_library(ffmpeg_tx::avutil ALIAS ffmpeg_tx_avutil)

# ──────────────────────────────────────────────────────────────────────────────
# Status summary
# ──────────────────────────────────────────────────────────────────────────────

message(STATUS "FFmpeg TX: source dir  = ${FFMPEG_SOURCE_DIR}")
message(STATUS "FFmpeg TX: build dir   = ${FFMPEG_BUILD_DIR}")
message(STATUS "FFmpeg TX: output lib  = ${FFMPEG_AVUTIL_A}")
message(STATUS "FFmpeg TX: parallelism = -j${NPROC}")
message(STATUS "FFmpeg TX: link target = ffmpeg_tx::avutil")
message(STATUS "FFmpeg TX: public API  = #include <libavutil/tx.h>")
