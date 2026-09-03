
# cmake options: TARGET_C_ARCH / TARGET_CPP_ARCH:
#   and optionally:  TARGET_C_EXTRA TARGET_CXX_EXTRA
#
# provided:
#   - function: target_set_c_arch_flags(<target>)    # uses options TARGET_C_ARCH and TARGET_C_EXTRA
#   - function: target_set_cxx_arch_flags(<target>)  # uses options TARGET_CXX_ARCH and TARGET_CXX_EXTRA
#   - macro:    target_set_cxx_arch_option(<target> <gcc/clang_march> <gcc/clang_extra> <msvc_arch>)
#
# see https://en.wikichip.org/wiki/x86/extensions
# and https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
#   for gcc specific architecture options
# and https://docs.microsoft.com/en-us/cpp/build/reference/arch-x64
# or  https://docs.microsoft.com/en-us/cpp/build/reference/arch-x86
#   for msvc specific architecture options

# https://en.wikichip.org/wiki/arm/versions
# https://en.wikipedia.org/wiki/Raspberry_Pi
# https://gcc.gnu.org/onlinedocs/gcc/ARM-Options.html#ARM-Options
# https://en.wikipedia.org/wiki/Comparison_of_ARMv7-A_cores
# https://en.wikipedia.org/wiki/Comparison_of_ARMv8-A_cores

# arm32_rpi1 untested
#   -mcpu=arm1176jzf-s -mfloat-abi=hard -mfpu=vfp         -mtune=arm1176jzf-s
# arm32_rpi2 untested
#   "-march=armv7-a"   "-mfloat-abi=hard" "-mfpu=neon-vfpv4"
#   "-march=armv8-a"   "-mfloat-abi=hard" "-mfpu=neon-vfpv4"
# arm32_rpi3 with "armv7-a" tested on Raspbian GNU/Linux 10 (buster), 32-bit  => MIPP test reports: NEONv1, 128 bits
#   "-march=armv7-a"   "-mfloat-abi=hard" "-mfpu=neon-vfpv4"
# arm32_rpi3 with "armv8-a" tested on Raspbian GNU/Linux 10 (buster), 32-bit  => MIPP test reports: NEONv1, 128 bits
#   "-march=armv8-a"   "-mfloat-abi=hard" "-mfpu=neon-vfpv4"
# arm32_rpi3 with "armv8-a" tested on Raspbian GNU/Linux 10 (buster), 32-bit  => MIPP test reports: NEONv1, 128 bits
#   "-march=armv8-a"   "-mfloat-abi=hard" "-mfpu=neon-vfpv4" "-mtune=cortex-a53"
# arm32_rpi4 untested
#   RPi 4 Model B:    Cortex-A72  =>  "-mtune=cortex-a72"  ?
#   "-mcpu=cortex-a72 -mfloat-abi=hard -mfpu=neon-fp-armv8 -mneon-for-64bits  -mtune=cortex-a72"

set(MSVC_EXTRA_OPT_none "")
set(GCC_EXTRA_OPT_none "")
set(GCC_EXTRA_OPT_neon_vfpv4    "-mfloat-abi=hard" "-mfpu=neon-vfpv4")
set(GCC_EXTRA_OPT_neon_rpi3_a53 "-mfloat-abi=hard" "-mfpu=neon-vfpv4" "-mtune=cortex-a53")
set(GCC_EXTRA_OPT_neon_rpi4_a72 "-mfloat-abi=hard" "-mfpu=neon-fp-armv8" "-mtune=cortex-a72")
set(GCC_EXTRA_OPT_apple_m1      "-mcpu=apple-m1")
set(GCC_EXTRA_OPT_ios_arm64     "-arch arm64")

string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" CMAKE_SYSTEM_PROCESSOR_LOWER)
if ( (CMAKE_SYSTEM_PROCESSOR_LOWER STREQUAL "i686") OR (CMAKE_SYSTEM_PROCESSOR_LOWER STREQUAL "x86_64")
    # On Windows CMake emits "AMD64" (64-bit) or "x86" (32-bit) rather than the above
    OR (CMAKE_SYSTEM_PROCESSOR_LOWER STREQUAL "amd64") OR (CMAKE_SYSTEM_PROCESSOR_LOWER STREQUAL "x86") )
    set(GCC_MARCH_DESC "native/SSE2:pentium4/SSE3:core2/SSE4:nehalem/AVX:sandybridge/AVX2:haswell")
    set(GCC_MARCH_VALUES "none;native;pentium4;core2;nehalem;sandybridge;haswell" CACHE INTERNAL "List of possible architectures")
    set(GCC_EXTRA_VALUES "" CACHE INTERNAL "List of possible EXTRA options")
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")
    set(GCC_MARCH_DESC "native/ARMwNEON:armv8-a")
    set(GCC_MARCH_VALUES "none;native;armv8-a" CACHE INTERNAL "List of possible architectures")
    if (APPLE)
        set(GCC_EXTRA_VALUES "none;apple_m1;ios_arm64" CACHE INTERNAL "List of possible additional options")
    else()
        set(GCC_EXTRA_VALUES "" CACHE INTERNAL "List of possible additional options")
    endif()
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "armv7l")
    set(GCC_MARCH_DESC "native/ARMwNEON:armv7-a")
    set(GCC_MARCH_VALUES "none;native;armv7-a" CACHE INTERNAL "List of possible architectures")
    set(GCC_EXTRA_VALUES "none;neon_vfpv4;neon_rpi3_a53;neon_rpi4_a72" CACHE INTERNAL "List of possible additional options")
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "riscv64")
    # RISC-V with the V (vector) extension; rv64gcv enables the standard
    # G+C profile plus V (RVV 1.0). pf_rvv_*.h is gated on __riscv_vector.
    set(GCC_MARCH_DESC "native/RISCVwRVV:rv64gcv")
    set(GCC_MARCH_VALUES "none;native;rv64gc;rv64gcv" CACHE INTERNAL "List of possible architectures")
    set(GCC_EXTRA_VALUES "" CACHE INTERNAL "List of possible additional options")
elseif (EMSCRIPTEN)
    # Emscripten WASM SIMD is handled automatically in target_set_c/cxx_arch_flags
    set(GCC_MARCH_DESC "wasm-simd")
    set(GCC_MARCH_VALUES "none" CACHE INTERNAL "List of possible architectures")
    set(GCC_EXTRA_VALUES "" CACHE INTERNAL "List of possible additional options")
else()
    message(WARNING "unsupported CMAKE_SYSTEM_PROCESSOR '${CMAKE_SYSTEM_PROCESSOR}'")
    # other PROCESSORs could be "ppc", "ppc64",  "arm" - or something else?!
    set(GCC_MARCH_DESC "native")
    set(GCC_MARCH_VALUES "none;native" CACHE INTERNAL "List of possible architectures")
    set(GCC_EXTRA_VALUES "" CACHE INTERNAL "List of possible additional options")
endif()

# cmake options - depending on C/C++ compiler
# how are chances, that C and C++ compilers are from different vendors?
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(TARGET_C_ARCH "none" CACHE STRING "gcc target C architecture (-march): ${GCC_MARCH_DESC}")
    set_property(CACHE TARGET_C_ARCH PROPERTY STRINGS ${GCC_MARCH_VALUES})
    if ( NOT (GCC_EXTRA_VALUES STREQUAL "") )
        set(TARGET_C_EXTRA "none" CACHE STRING "gcc additional options for C")
        set_property(CACHE TARGET_C_EXTRA PROPERTY STRINGS ${GCC_EXTRA_VALUES})
    endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang" OR CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    set(TARGET_C_ARCH "none" CACHE STRING "clang target C architecture (-march): ${GCC_MARCH_DESC}")
    set_property(CACHE TARGET_C_ARCH PROPERTY STRINGS ${GCC_MARCH_VALUES})
    if ( NOT (GCC_EXTRA_VALUES STREQUAL "") )
        set(TARGET_C_EXTRA "none" CACHE STRING "gcc additional options for C")
        set_property(CACHE TARGET_C_EXTRA PROPERTY STRINGS ${GCC_EXTRA_VALUES})
    endif()
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(TARGET_C_ARCH "none" CACHE STRING "msvc target C architecture (/arch): SSE2/AVX/AVX2/AVX512")
    set(TARGET_C_EXTRA "none" CACHE STRING "msvc additional options")
else()
    message(WARNING "unsupported C compiler '${CMAKE_C_COMPILER_ID}', see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(TARGET_CXX_ARCH "none" CACHE STRING "gcc target C++ architecture (-march): ${GCC_MARCH_DESC}")
    set_property(CACHE TARGET_CXX_ARCH PROPERTY STRINGS ${GCC_MARCH_VALUES})
    if ( NOT (GCC_EXTRA_VALUES STREQUAL "") )
        set(TARGET_CXX_EXTRA "none" CACHE STRING "gcc additional options for C++")
        set_property(CACHE TARGET_CXX_EXTRA PROPERTY STRINGS ${GCC_EXTRA_VALUES})
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(TARGET_CXX_ARCH "none" CACHE STRING "clang target C++ architecture (-march): ${GCC_MARCH_DESC}")
    set_property(CACHE TARGET_CXX_ARCH PROPERTY STRINGS ${GCC_MARCH_VALUES})
    if ( NOT (GCC_EXTRA_VALUES STREQUAL "") )
        set(TARGET_CXX_EXTRA "none" CACHE STRING "clang additional options for C++")
        set_property(CACHE TARGET_CXX_EXTRA PROPERTY STRINGS ${GCC_EXTRA_VALUES})
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(TARGET_CXX_ARCH "none" CACHE STRING "msvc target C++ architecture (/arch): SSE2/AVX/AVX2/AVX512")
    set(TARGET_CXX_EXTRA "none" CACHE STRING "msvc additional options")
else()
    message(WARNING "unsupported C++ compiler '${CMAKE_CXX_COMPILER_ID}', see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
endif()

######################################################

# Pick an -march string for RISC-V and report whether the V extension
# is available. We compile and run a tiny probe against -march=rv64gcv:
#
#   - if compile + run succeed, the target supports V; the probe prints
#     the host VLEN, and if VLEN >= 256 we append the matching zvlNNNb
#     so the compiler can emit single-instruction 4-wide double ops
#     instead of splitting them into VLEN_min=128 pieces;
#   - if compile or run fails, the target lacks V; we fall back to
#     plain rv64gc and the caller is expected to skip PFFFT_ENABLE_RVV
#     so the headers stay on the scalar path.
#
# The probe runs natively or under CMAKE_CROSSCOMPILING_EMULATOR (e.g.
# qemu-riscv64-static). When neither applies (bare cross-compile) we
# optimistically default to rv64gcv with has_v=TRUE; users targeting
# non-V cores should set TARGET_C_ARCH=rv64gc explicitly.
function(_pffft_riscv_march out_march out_has_v)
    if (NOT DEFINED PFFFT_RVV_HAS_V
        AND (NOT CMAKE_CROSSCOMPILING OR CMAKE_CROSSCOMPILING_EMULATOR))
        set(_probe "${CMAKE_CURRENT_BINARY_DIR}/pffft_rvv_vlen_probe.c")
        file(WRITE "${_probe}"
            "#include <stdio.h>\n"
            "#include <riscv_vector.h>\n"
            "int main(void){ printf(\"%zu\\n\", __riscv_vsetvlmax_e32m1() * 32); return 0; }\n")
        try_run(_run_res _compile_ok
            "${CMAKE_CURRENT_BINARY_DIR}/rvv_vlen_probe"
            "${_probe}"
            COMPILE_DEFINITIONS "-march=rv64gcv"
            RUN_OUTPUT_VARIABLE _vlen)
        if (_compile_ok AND _run_res EQUAL 0)
            string(STRIP "${_vlen}" _vlen)
            set(PFFFT_RVV_HAS_V TRUE CACHE INTERNAL "Target supports the RISC-V V extension")
            set(PFFFT_RVV_VLEN "${_vlen}" CACHE INTERNAL "Detected RISC-V VLEN in bits")
        else()
            set(PFFFT_RVV_HAS_V FALSE CACHE INTERNAL "Target lacks the RISC-V V extension")
        endif()
    endif()

    if (DEFINED PFFFT_RVV_HAS_V AND NOT PFFFT_RVV_HAS_V)
        set(${out_march} "rv64gc" PARENT_SCOPE)
        set(${out_has_v} FALSE PARENT_SCOPE)
    elseif (PFFFT_RVV_VLEN AND PFFFT_RVV_VLEN GREATER_EQUAL 256)
        set(${out_march} "rv64gcv_zvl${PFFFT_RVV_VLEN}b" PARENT_SCOPE)
        set(${out_has_v} TRUE PARENT_SCOPE)
    else()
        set(${out_march} "rv64gcv" PARENT_SCOPE)
        set(${out_has_v} TRUE PARENT_SCOPE)
    endif()
endfunction()

# Matches RISC-V march strings that include the V extension. Examples:
#   rv64gcv, rv64gcv_zvl256b, rv64imafdcv, rv32gc_v, rv64gc_v_zvl128b.
# Deliberately does NOT match plain rv64gc / rv64imafdc — the leading
# "rv" alone is not enough to claim the V extension.
set(PFFFT_RVV_MARCH_REGEX "rv[0-9]+[a-z_]*v($|_)")

function(target_set_c_arch_flags target)
    # Emscripten WASM SIMD via NEON emulation
    if (EMSCRIPTEN)
        message(STATUS "Emscripten detected: enabling WASM SIMD with NEON emulation for C target ${target}")
        target_compile_options(${target} PRIVATE "-msimd128")
        target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_NEON=1)
        return()
    endif()
    # On RISC-V the V extension is opt-in via -march; auto-enable RVV
    # when we can verify the target supports V, and skip it otherwise
    # so a non-V target still gets a usable rv64gc build.
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "riscv64")
        if ( ("${TARGET_C_ARCH}" STREQUAL "") OR ("${TARGET_C_ARCH}" STREQUAL "none") )
            _pffft_riscv_march(_rvv_march _rvv_has_v)
            target_compile_options(${target} PRIVATE "-march=${_rvv_march}")
            if (_rvv_has_v)
                target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_RVV=1)
                message(STATUS "RISC-V detected: defaulting C target ${target} to -march=${_rvv_march} with PFFFT_ENABLE_RVV")
            else()
                message(STATUS "RISC-V detected (no V extension): defaulting C target ${target} to -march=${_rvv_march}")
            endif()
            return()
        elseif ("${TARGET_C_ARCH}" MATCHES "${PFFFT_RVV_MARCH_REGEX}")
            target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_RVV=1)
            message(STATUS "RISC-V detected: enabling PFFFT_ENABLE_RVV for C target ${target}")
        endif()
    endif()
    if ( ("${TARGET_C_ARCH}" STREQUAL "") OR ("${TARGET_C_ARCH}" STREQUAL "none") )
        message(STATUS "C ARCH for target ${target} is not set!")
    else()
        if ( (CMAKE_C_COMPILER_ID STREQUAL "GNU") OR (CMAKE_C_COMPILER_ID STREQUAL "Clang") OR (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") )
            target_compile_options(${target} PRIVATE "-march=${TARGET_C_ARCH}")
            message(STATUS "C ARCH for target ${target} set: ${TARGET_C_ARCH}")
        elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
            target_compile_options(${target} PRIVATE "/arch:${TARGET_C_ARCH}")
            message(STATUS "C ARCH for target ${target} set: ${TARGET_C_ARCH}")
        else()
            message(WARNING "unsupported C compiler '${CMAKE_C_COMPILER_ID}' for target_set_c_arch_flags(), see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
        endif()
    endif()

    if ( ("${TARGET_C_EXTRA}" STREQUAL "") OR ("${TARGET_C_EXTRA}" STREQUAL "none") )
        message(STATUS "C additional options for target ${target} is not set!")
    else()
        if ( (CMAKE_C_COMPILER_ID STREQUAL "GNU") OR (CMAKE_C_COMPILER_ID STREQUAL "Clang") OR (CMAKE_C_COMPILER_ID STREQUAL "AppleClang"))
            target_compile_options(${target} PRIVATE "${GCC_EXTRA_OPT_${TARGET_C_EXTRA}}")
            message(STATUS "C additional options for target ${target} set: ${GCC_EXTRA_OPT_${TARGET_C_EXTRA}}")
        elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
            # target_compile_options(${target} PRIVATE "${MSVC_EXTRA_OPT_${TARGET_C_EXTRA}}")
            message(STATUS "C additional options for target ${target} not usable with MSVC")
        else()
            message(WARNING "unsupported C compiler '${CMAKE_C_COMPILER_ID}' for target_set_c_arch_flags(), see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
        endif()
    endif()

    # Enable NEON on ARM: either explicitly requested via TARGET_C_EXTRA, or auto-detected on aarch64/arm64
    if ( ("${TARGET_C_EXTRA}" MATCHES "^neon_.*") OR (CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64") OR (CMAKE_SYSTEM_PROCESSOR MATCHES "arm64"))
        message(STATUS "setting PFFFT_ENABLE_NEON for C target ${target}")
        target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_NEON=1)
    endif()
endfunction()

function(target_set_cxx_arch_flags target)
    # Emscripten WASM SIMD via NEON emulation
    if (EMSCRIPTEN)
        message(STATUS "Emscripten detected: enabling WASM SIMD with NEON emulation for C++ target ${target}")
        target_compile_options(${target} PRIVATE "-msimd128")
        target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_NEON=1)
        return()
    endif()
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "riscv64")
        if ( ("${TARGET_CXX_ARCH}" STREQUAL "") OR ("${TARGET_CXX_ARCH}" STREQUAL "none") )
            _pffft_riscv_march(_rvv_march _rvv_has_v)
            target_compile_options(${target} PRIVATE "-march=${_rvv_march}")
            if (_rvv_has_v)
                target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_RVV=1)
                message(STATUS "RISC-V detected: defaulting C++ target ${target} to -march=${_rvv_march} with PFFFT_ENABLE_RVV")
            else()
                message(STATUS "RISC-V detected (no V extension): defaulting C++ target ${target} to -march=${_rvv_march}")
            endif()
            return()
        elseif ("${TARGET_CXX_ARCH}" MATCHES "${PFFFT_RVV_MARCH_REGEX}")
            target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_RVV=1)
            message(STATUS "RISC-V detected: enabling PFFFT_ENABLE_RVV for C++ target ${target}")
        endif()
    endif()
    if ( ("${TARGET_CXX_ARCH}" STREQUAL "") OR ("${TARGET_CXX_ARCH}" STREQUAL "none") )
        message(STATUS "C++ ARCH for target ${target} is not set!")
    else()
        if ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") )
            target_compile_options(${target} PRIVATE "-march=${TARGET_CXX_ARCH}")
            message(STATUS "C++ ARCH for target ${target} set: ${TARGET_CXX_ARCH}")
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            target_compile_options(${target} PRIVATE "/arch:${TARGET_CXX_ARCH}")
            message(STATUS "C++ ARCH for target ${target} set: ${TARGET_CXX_ARCH}")
        else()
            message(WARNING "unsupported C++ compiler '${CMAKE_CXX_COMPILER_ID}' for target_set_cxx_arch_flags(), see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
        endif()
    endif()
    if ( ("${TARGET_CXX_EXTRA}" STREQUAL "") OR ("${TARGET_CXX_EXTRA}" STREQUAL "none") )
        message(STATUS "C++ additional options for target ${target} is not set!")
    else()
        if ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") )
            target_compile_options(${target} PRIVATE "${GCC_EXTRA_OPT_${TARGET_CXX_EXTRA}}")
            message(STATUS "C++ additional options for target ${target} set: ${GCC_EXTRA_OPT_${TARGET_CXX_EXTRA}}")
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            # target_compile_options(${target} PRIVATE "${MSVC_EXTRA_OPT_${TARGET_CXX_EXTRA}}")
            message(STATUS "C++ additional options for target ${target} not usable with MSVC")
        else()
            message(WARNING "unsupported C compiler '${CMAKE_C_COMPILER_ID}' for target_set_c_arch_flags(), see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
        endif()
    endif()
    # Enable NEON on ARM: either explicitly requested via TARGET_CXX_EXTRA, or auto-detected on aarch64/arm64
    if ( ("${TARGET_CXX_EXTRA}" MATCHES "^neon_.*") OR (CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64") OR (CMAKE_SYSTEM_PROCESSOR MATCHES "arm64"))
        message(STATUS "setting PFFFT_ENABLE_NEON for C++ target ${target}")
        target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_NEON=1)
    endif()
endfunction()


macro(target_set_cxx_arch_option target gcc_clang_arch gcc_clang_extra msvc_arch )
    if ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") )

        if ( NOT (("${gcc_clang_arch}" STREQUAL "") OR ("${gcc_clang_arch}" STREQUAL "none") ) )
            target_compile_options(${target} PRIVATE "-march=${gcc_clang_arch}")
            message(STATUS "C++ ARCH for target ${target}: ${gcc_clang_arch}")
        endif()
        if (NOT ( ("${gcc_clang_extra}" STREQUAL "") OR ("${gcc_clang_extra}" STREQUAL "none") ) )
            target_compile_options(${target} PRIVATE "${GCC_EXTRA_OPT_${gcc_clang_extra}}")
            message(STATUS "C++ additional options for target ${target}: ${GCC_EXTRA_OPT_${gcc_clang_extra}}")
        endif()
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        if ( NOT (("${msvc_arch}" STREQUAL "") OR ("${msvc_arch}" STREQUAL "none") ) )
            target_compile_options(${target} PRIVATE "/arch:${msvc_arch}")
            message(STATUS "C++ ARCH for target ${target} set: ${msvc_arch}")
        endif()
    else()
        message(WARNING "unsupported C++ compiler '${CMAKE_CXX_COMPILER_ID}' for target_set_cxx_arch_option(), see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
    endif()
endmacro()
