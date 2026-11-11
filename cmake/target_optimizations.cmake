
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

# pffft is compiled with the toolchain's default target flags, plus
# TARGET_C_ARCH / TARGET_C_EXTRA if set. On WebAssembly, and on RISC-V while
# TARGET_C_ARCH is unset, target_set_c_arch_flags() adds its own flags
# instead (see the list below).
# pf_float.h / pf_double.h select the
# SIMD backend from the compiler's predefined macros (__SSE__, __ARM_NEON,
# __ALTIVEC__, ...), and pffft_report_simd_backend() prints which one. E.g.
# the Android NDK enables NEON for armeabi-v7a by default, so pffft uses NEON
# there; Debian's armhf gcc needs -mfpu=neon or a neon_* preset for it.
#
# These targets get extra handling, recognised by CMAKE_SYSTEM_PROCESSOR:
#   - WebAssembly: target_set_c_arch_flags() adds -msimd128 for WASM SIMD,
#     and -mrelaxed-simd with PFFFT_USE_WASM_RELAXED_SIMD, in place of
#     TARGET_C_ARCH / TARGET_C_EXTRA. Checked first, as Emscripten's
#     toolchain file reports "x86".
#   - x86: the FMA reminder below, and the /arch choices offered for MSVC.
#   - 32-bit ARM A-profile, i.e. Cortex-A cores: the note on gcc's NEON
#     loads and stores below.
# RISC-V gets its -march from a probe in target_set_c_arch_flags(), unless
# TARGET_C_ARCH is set.
#
# CMAKE_SYSTEM_PROCESSOR is whatever the platform or toolchain file calls the
# target: uname -m on Linux ("x86_64", "i686", "armv7l"), uname -p on FreeBSD
# ("amd64", "i386"), "AMD64" / "x86" on Windows, "armv7-a" / "i686" for the
# Android NDK, "armv8l" for 32-bit ARM code on a 64-bit Linux kernel.
# Caveat: Visual Studio builds for -A ARM64 / -A ARM64EC without a
# toolchain file keep the host's "AMD64" here, so they count as x86.
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" CMAKE_SYSTEM_PROCESSOR_LOWER)
if (EMSCRIPTEN OR (CMAKE_SYSTEM_PROCESSOR_LOWER MATCHES "^wasm"))
    set(PFFFT_TARGET_IS_WASM TRUE)
elseif (CMAKE_SYSTEM_PROCESSOR_LOWER MATCHES "^(x86|x86_64|x64|amd64|i[3-6]86)$")
    set(PFFFT_TARGET_IS_X86 TRUE)
elseif (CMAKE_SYSTEM_PROCESSOR_LOWER MATCHES "^armv7(-?a|h?l)?$|^armv8l$")
    set(PFFFT_TARGET_IS_ARMV7 TRUE)
endif()

# cmake options - depending on C/C++ compiler
# TARGET_C_ARCH takes any -march value of the compiler as free text;
# TARGET_C_EXTRA offers the presets defined above, each for its own target.
# HELPSTRING and STRINGS are set explicitly, as set(... CACHE ...) keeps
# those of an existing build directory.
set(GCC_EXTRA_VALUES "none;neon_vfpv4;neon_rpi3_a53;neon_rpi4_a72;apple_m1")
set(MSVC_ARCH_VALUES "none;SSE2;AVX;AVX2;AVX512")
set(_pffft_march_desc "target architecture (-march), e.g. native, haswell, armv7-a, rv64gcv; none: the compiler's default")
set(_pffft_extra_desc "preset of additional options: neon_vfpv4, neon_rpi3_a53, neon_rpi4_a72 for 32-bit ARM hard-float; apple_m1 for Apple arm64; none")

# _pffft_gcc_clang_options(<C|CXX> <display name>)
# WebAssembly gets neither option, as target_set_c_arch_flags() applies its
# own flags there.
macro(_pffft_gcc_clang_options lang langname)
    if (NOT PFFFT_TARGET_IS_WASM)
        # A build directory configured by an older pffft holds a list of
        # choices for TARGET_<lang>_ARCH. CMake has no way to remove a single
        # cache property, so the entry is recreated with its current value.
        get_property(_pffft_has_list CACHE TARGET_${lang}_ARCH PROPERTY STRINGS SET)
        if (_pffft_has_list)
            get_property(_pffft_value CACHE TARGET_${lang}_ARCH PROPERTY VALUE)
            unset(TARGET_${lang}_ARCH CACHE)
            set(TARGET_${lang}_ARCH "${_pffft_value}" CACHE STRING "")
        endif()
        set(TARGET_${lang}_ARCH "none" CACHE STRING "")
        set_property(CACHE TARGET_${lang}_ARCH PROPERTY HELPSTRING "gcc/clang ${langname} ${_pffft_march_desc}")
        set(TARGET_${lang}_EXTRA "none" CACHE STRING "")
        set_property(CACHE TARGET_${lang}_EXTRA PROPERTY HELPSTRING "gcc/clang ${langname} ${_pffft_extra_desc}")
        set_property(CACHE TARGET_${lang}_EXTRA PROPERTY STRINGS ${GCC_EXTRA_VALUES})
    endif()
endmacro()

if (CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    _pffft_gcc_clang_options(C "C")
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(TARGET_C_ARCH "none" CACHE STRING "msvc target C architecture (/arch): SSE2/AVX/AVX2/AVX512")
    if (PFFFT_TARGET_IS_X86)
        set_property(CACHE TARGET_C_ARCH PROPERTY STRINGS ${MSVC_ARCH_VALUES})
    endif()
    set(TARGET_C_EXTRA "none" CACHE STRING "msvc additional options")
else()
    message(WARNING "unsupported C compiler '${CMAKE_C_COMPILER_ID}', see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    _pffft_gcc_clang_options(CXX "C++")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(TARGET_CXX_ARCH "none" CACHE STRING "msvc target C++ architecture (/arch): SSE2/AVX/AVX2/AVX512")
    if (PFFFT_TARGET_IS_X86)
        set_property(CACHE TARGET_CXX_ARCH PROPERTY STRINGS ${MSVC_ARCH_VALUES})
    endif()
    set(TARGET_CXX_EXTRA "none" CACHE STRING "msvc additional options")
else()
    message(WARNING "unsupported C++ compiler '${CMAKE_CXX_COMPILER_ID}', see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
endif()

######################################################
# x86: FMA is opt-in and easy to miss.
#
# The SSE-float and AVX-double VMADD/VMSUB macros lower to true FMA only when
# the compiler guarantees FMA3 -- gcc/clang via -march=haswell or later, MSVC
# via /arch:AVX2. Without it every fused op costs a separate multiply plus add,
# which in the frequency-domain multiply loops means 20 floating point
# operations per iteration instead of 12 (measured on gcc 12 and clang 14,
# x86-64 and i686, float and double alike).
#
# The default is deliberately left alone: FMA3 needs Haswell (2013) or
# Piledriver, and raising the baseline would break older CPUs. So just say so.
if (PFFFT_TARGET_IS_X86)
    set(_pffft_have_fma FALSE)
    if (CMAKE_C_COMPILER_ID MATCHES "MSVC")
        set(_pffft_fma_arch "AVX2")
        if ("${TARGET_C_ARCH}" MATCHES "AVX2|AVX512|AVX10")
            set(_pffft_have_fma TRUE)
        endif()
    else()
        set(_pffft_fma_arch "haswell")
        # "native" is resolved by the compiler on the build machine, so whether
        # it implies FMA is unknown here -- stay quiet rather than guess wrong
        if ( ("${TARGET_C_ARCH}" STREQUAL "haswell") OR ("${TARGET_C_ARCH}" STREQUAL "native") )
            set(_pffft_have_fma TRUE)
        endif()
    endif()
    if (NOT _pffft_have_fma)
        message(STATUS "TARGET_C_ARCH='${TARGET_C_ARCH}' has no FMA;"
                       " use -DTARGET_C_ARCH=${_pffft_fma_arch} for ~40% fewer FP ops"
                       " in the convolution loops (requires Haswell/Piledriver or later)")
    endif()
endif()

######################################################
# ARMv7: gcc moves 128-bit NEON vectors in 64-bit pieces.
#
# Not a gcc bug. float32x4_t is specified to be laid out as a NEON register
# stored with VSTM, which on big-endian differs from what VST1 produces, so for
# a plain v4sf array access gcc restricts itself to the layout-compatible forms:
# vld1.64/vst1.64, or a vldr/vstr pair whenever a nonzero immediate offset is
# needed, since vld1.64 has no offset addressing mode. Any loop touching a
# pointer more than once per iteration hits the latter. clang takes a different
# position and uses vld1.32 for the same code. See GCC PR43725 comment 1.
#
# Measured with gcc 12.2 over the whole float translation unit: 1600 vldr/vstr
# and 344 vld1/vst1, against clang 14's 181 and 797 -- about twice the
# load/store operations for the same bytes, spread over ~19 functions. No flag
# changes it (14 combinations tried), and gcc 16.1.0 behaves like 12.2. The
# alternative to switching compilers would be to load/store through
# vld1q_f32()/vst1q_f32() intrinsics, which this library does not do.
if (PFFFT_TARGET_IS_ARMV7 AND CMAKE_C_COMPILER_ID STREQUAL "GNU")
    message(STATUS "gcc moves 128-bit NEON vectors in 64-bit pieces on ARMv7"
                   " (see GCC PR43725), roughly doubling load/store work across the"
                   " transform; clang is recommended here (CC=clang CXX=clang++)")
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
#     plain rv64gc, so the headers stay on the scalar path.
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
        try_run(PFFFT_RVV_PROBE_RUN PFFFT_RVV_PROBE_COMPILED
            "${CMAKE_CURRENT_BINARY_DIR}/rvv_vlen_probe"
            "${_probe}"
            COMPILE_DEFINITIONS "-march=rv64gcv"
            RUN_OUTPUT_VARIABLE _vlen)
        # the output may also hold emulator chatter (qemu warns on stderr
        # when no vext_spec is given); the VLEN is the last number printed
        string(REGEX MATCH "[0-9]+[ \t\r\n]*$" _vlen "${_vlen}")
        string(STRIP "${_vlen}" _vlen)
        if (PFFFT_RVV_PROBE_COMPILED AND PFFFT_RVV_PROBE_RUN EQUAL 0 AND _vlen)
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

function(target_set_c_arch_flags target)
    if (PFFFT_USE_HELIUM)
        target_compile_definitions(${target} PRIVATE PFFFT_ENABLE_HELIUM=1)
    endif()
    # WebAssembly SIMD (Emscripten or WASI-SDK / plain clang)
    if (PFFFT_TARGET_IS_WASM)
        target_compile_options(${target} PRIVATE "-msimd128")
        if (PFFFT_USE_WASM_RELAXED_SIMD)
            message(STATUS "WebAssembly target: enabling WASM SIMD + Relaxed SIMD for C target ${target}")
            target_compile_options(${target} PRIVATE "-mrelaxed-simd")
        else()
            message(STATUS "WebAssembly target: enabling WASM SIMD for C target ${target}")
        endif()
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
                message(STATUS "RISC-V detected: defaulting C target ${target} to -march=${_rvv_march}")
            else()
                message(STATUS "RISC-V detected (no V extension): defaulting C target ${target} to -march=${_rvv_march}")
            endif()
            return()
        endif()
    endif()
    if ( ("${TARGET_C_ARCH}" STREQUAL "") OR ("${TARGET_C_ARCH}" STREQUAL "none") )
        message(STATUS "TARGET_C_ARCH not set for target ${target}, using default")
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
        message(STATUS "TARGET_C_EXTRA not set for target ${target}, using default")
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

endfunction()

function(target_set_cxx_arch_flags target)
    # WebAssembly SIMD (Emscripten or WASI-SDK / plain clang)
    if (PFFFT_TARGET_IS_WASM)
        target_compile_options(${target} PRIVATE "-msimd128")
        if (PFFFT_USE_WASM_RELAXED_SIMD)
            message(STATUS "WebAssembly target: enabling WASM SIMD + Relaxed SIMD for C++ target ${target}")
            target_compile_options(${target} PRIVATE "-mrelaxed-simd")
        else()
            message(STATUS "WebAssembly target: enabling WASM SIMD for C++ target ${target}")
        endif()
        return()
    endif()
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "riscv64")
        if ( ("${TARGET_CXX_ARCH}" STREQUAL "") OR ("${TARGET_CXX_ARCH}" STREQUAL "none") )
            _pffft_riscv_march(_rvv_march _rvv_has_v)
            target_compile_options(${target} PRIVATE "-march=${_rvv_march}")
            if (_rvv_has_v)
                message(STATUS "RISC-V detected: defaulting C++ target ${target} to -march=${_rvv_march}")
            else()
                message(STATUS "RISC-V detected (no V extension): defaulting C++ target ${target} to -march=${_rvv_march}")
            endif()
            return()
        endif()
    endif()
    if ( ("${TARGET_CXX_ARCH}" STREQUAL "") OR ("${TARGET_CXX_ARCH}" STREQUAL "none") )
        message(STATUS "TARGET_CXX_ARCH not set for target ${target}, using default")
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
        message(STATUS "TARGET_CXX_EXTRA not set for target ${target}, using default")
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
