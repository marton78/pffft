
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
# arm32_rpi3 with "armv7-a" tested on Raspbian GNU/Linux 10 (buster), 32-bit: NEON, 128 bits
#   "-march=armv7-a"   "-mfloat-abi=hard" "-mfpu=neon-vfpv4"
# arm32_rpi3 with "armv8-a" tested on Raspbian GNU/Linux 10 (buster), 32-bit: NEON, 128 bits
#   "-march=armv8-a"   "-mfloat-abi=hard" "-mfpu=neon-vfpv4"
# arm32_rpi3 with "armv8-a" tested on Raspbian GNU/Linux 10 (buster), 32-bit: NEON, 128 bits
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
# TARGET_C_ARCH / TARGET_C_EXTRA if set. On WebAssembly,
# target_set_c_arch_flags() adds its own flags instead (see the list below).
# pf_float.h / pf_double.h select the
# SIMD backend from the compiler's predefined macros (__SSE__, __ARM_NEON,
# __ALTIVEC__, ...), and pffft_report_simd_backend() prints which one. E.g.
# the Android NDK enables NEON for armeabi-v7a by default, so pffft uses NEON
# there; Debian's armhf gcc needs -mfpu=neon or a neon_* preset for it.
#
# These targets get extra handling, recognised by PFFFT_TARGET_PROCESSOR:
#   - WebAssembly: target_set_c_arch_flags() adds -msimd128 for WASM SIMD,
#     and -mrelaxed-simd with PFFFT_USE_WASM_RELAXED_SIMD, in place of
#     TARGET_C_ARCH / TARGET_C_EXTRA. Checked first, as Emscripten's
#     toolchain file reports "x86".
#   - x86: the FMA reminder below, and the /arch choices offered for MSVC.
#   - 32-bit ARM A-profile, i.e. Cortex-A cores: the note on gcc's NEON
#     loads and stores below.
#   - 64-bit RISC-V: TARGET_C_ARCH=auto, which enables the vector extension
#     if the build machine has it (see _pffft_riscv_auto_march() below).
#
# PFFFT_TARGET_PROCESSOR is the lower-case name of the target processor:
#   - compilers for the MSVC ABI (cl, clang-cl, clang --target=*-windows-msvc):
#     the compiler's own target, CMAKE_C_COMPILER_ARCHITECTURE_ID: "x64",
#     "x86", "arm64", "arm64ec", "armv7". Unless a toolchain file sets it,
#     CMAKE_SYSTEM_PROCESSOR is the host's there, e.g. "AMD64" in a Visual
#     Studio -A ARM64 build on x64. CMake 3.20 and later report ARM64EC;
#     older ones report it as "x64".
#   - Apple platforms with CMAKE_OSX_ARCHITECTURES: its one architecture,
#     "x86_64" or "arm64". CMAKE_SYSTEM_PROCESSOR is the host's there unless
#     cross-compiling, e.g. "arm64" for an x86_64 build on Apple Silicon.
#     With several architectures (a universal binary), "universal", which
#     gets no extra handling, as no single target is built.
#   - otherwise CMAKE_SYSTEM_PROCESSOR, whatever the platform or toolchain file
#     calls the target: uname -m on Linux ("x86_64", "i686", "armv7l"),
#     uname -p on FreeBSD ("amd64", "i386"), "armv7-a" / "i686" for the
#     Android NDK, "armv8l" for 32-bit ARM code on a 64-bit Linux kernel.
list(LENGTH CMAKE_OSX_ARCHITECTURES _pffft_n_osx_archs)
if ((CMAKE_C_COMPILER_ID STREQUAL "MSVC" OR CMAKE_C_SIMULATE_ID STREQUAL "MSVC")
    AND CMAKE_C_COMPILER_ARCHITECTURE_ID)
    string(TOLOWER "${CMAKE_C_COMPILER_ARCHITECTURE_ID}" PFFFT_TARGET_PROCESSOR)
elseif (APPLE AND _pffft_n_osx_archs EQUAL 1)
    string(TOLOWER "${CMAKE_OSX_ARCHITECTURES}" PFFFT_TARGET_PROCESSOR)
elseif (APPLE AND _pffft_n_osx_archs GREATER 1)
    set(PFFFT_TARGET_PROCESSOR "universal")
else()
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" PFFFT_TARGET_PROCESSOR)
endif()
if (EMSCRIPTEN OR (PFFFT_TARGET_PROCESSOR MATCHES "^wasm"))
    set(PFFFT_TARGET_IS_WASM TRUE)
elseif (PFFFT_TARGET_PROCESSOR MATCHES "^(x86|x86_64|x64|amd64|i[3-6]86)$")
    set(PFFFT_TARGET_IS_X86 TRUE)
elseif (PFFFT_TARGET_PROCESSOR MATCHES "^armv7(-?a|h?l)?$|^armv8l$")
    set(PFFFT_TARGET_IS_ARMV7 TRUE)
elseif (PFFFT_TARGET_PROCESSOR STREQUAL "riscv64")
    set(PFFFT_TARGET_IS_RISCV64 TRUE)
endif()

# cmake options - depending on C/C++ compiler
# TARGET_C_ARCH takes any -march value of the compiler as free text;
# TARGET_C_EXTRA offers the presets defined above, each for its own target.
# HELPSTRING and STRINGS are set explicitly, as set(... CACHE ...) keeps
# those of an existing build directory.
set(GCC_EXTRA_VALUES "none;neon_vfpv4;neon_rpi3_a53;neon_rpi4_a72;apple_m1")
set(MSVC_ARCH_VALUES "none;SSE2;AVX;AVX2;AVX512")
set(_pffft_march_desc "target architecture (-march), e.g. native, haswell, armv7-a, rv64gcv; none: the compiler's default; auto (64-bit RISC-V only): the vector extension if the build machine has it")
set(_pffft_extra_desc "preset of additional options: neon_vfpv4, neon_rpi3_a53, neon_rpi4_a72 for 32-bit ARM hard-float; apple_m1 for Apple arm64; none")

# _pffft_drop_choice_list(<cache entry>)
# A build directory configured by an older pffft, or for another target,
# may hold a list of choices for the entry. CMake has no way to remove a
# single cache property, so the entry is recreated with its current value.
macro(_pffft_drop_choice_list entry)
    get_property(_pffft_has_list CACHE ${entry} PROPERTY STRINGS SET)
    if (_pffft_has_list)
        get_property(_pffft_value CACHE ${entry} PROPERTY VALUE)
        unset(${entry} CACHE)
        set(${entry} "${_pffft_value}" CACHE STRING "")
    endif()
endmacro()

# _pffft_gcc_clang_options(<C|CXX> <display name>)
# WebAssembly gets neither option, as target_set_c_arch_flags() applies its
# own flags there.
macro(_pffft_gcc_clang_options lang langname)
    if (NOT PFFFT_TARGET_IS_WASM)
        _pffft_drop_choice_list(TARGET_${lang}_ARCH)
        set(TARGET_${lang}_ARCH "none" CACHE STRING "")
        set_property(CACHE TARGET_${lang}_ARCH PROPERTY HELPSTRING "gcc/clang ${langname} ${_pffft_march_desc}")
        set(TARGET_${lang}_EXTRA "none" CACHE STRING "")
        set_property(CACHE TARGET_${lang}_EXTRA PROPERTY HELPSTRING "gcc/clang ${langname} ${_pffft_extra_desc}")
        set_property(CACHE TARGET_${lang}_EXTRA PROPERTY STRINGS ${GCC_EXTRA_VALUES})
    endif()
endmacro()

# _pffft_msvc_options(<C|CXX> <display name>)
# The /arch choices are offered for x86 only.
macro(_pffft_msvc_options lang langname)
    if (NOT PFFFT_TARGET_IS_X86)
        _pffft_drop_choice_list(TARGET_${lang}_ARCH)
    endif()
    set(TARGET_${lang}_ARCH "none" CACHE STRING "")
    set_property(CACHE TARGET_${lang}_ARCH PROPERTY HELPSTRING "msvc target ${langname} architecture (/arch), e.g. AVX2 on x86, armv8.2 on ARM64")
    if (PFFFT_TARGET_IS_X86)
        set_property(CACHE TARGET_${lang}_ARCH PROPERTY STRINGS ${MSVC_ARCH_VALUES})
    endif()
    set(TARGET_${lang}_EXTRA "none" CACHE STRING "msvc additional options")
endmacro()

if (CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    _pffft_gcc_clang_options(C "C")
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    _pffft_msvc_options(C "C")
else()
    message(WARNING "unsupported C compiler '${CMAKE_C_COMPILER_ID}', see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    _pffft_gcc_clang_options(CXX "C++")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    _pffft_msvc_options(CXX "C++")
else()
    message(WARNING "unsupported C++ compiler '${CMAKE_CXX_COMPILER_ID}', see https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html")
endif()

# "auto" is resolved for 64-bit RISC-V with gcc/clang only; elsewhere
# target_set_c_arch_flags() treats it as "none"
foreach(_pffft_lang C CXX)
    if (("${TARGET_${_pffft_lang}_ARCH}" STREQUAL "auto")
        AND NOT (PFFFT_TARGET_IS_RISCV64 AND CMAKE_${_pffft_lang}_COMPILER_ID MATCHES "^(GNU|Clang)$"))
        message(WARNING "TARGET_${_pffft_lang}_ARCH=auto is supported for 64-bit RISC-V with gcc or clang only;"
                        " using the compiler's default")
    endif()
endforeach()

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

# _pffft_riscv_auto_march(<C|CXX> <out_march> <out_reason>)
# Resolves TARGET_<lang>_ARCH=auto on 64-bit RISC-V: an -march with the
# vector extension, V, if the CPU the build runs on has it, or empty to keep
# the compiler's default. Like -march=native, the result is for CPUs like
# the build machine, so it is opt-in.
#
#   - If the compiler already enables V, with CMAKE_<lang>_FLAGS and those of
#     the build type (e.g. a toolchain for the RVA23 profile, or a -march
#     there), nothing is added.
#   - Otherwise a probe built with -march=rv64gcv runs natively, or under
#     CMAKE_CROSSCOMPILING_EMULATOR (e.g. qemu-riscv64). If it runs, the
#     result is rv64gcv. In a native build with VLEN >= 256 it also gets the
#     matching zvlNNNb, so the compiler can emit 4-wide double ops as one
#     instruction instead of splitting them into VLEN_min=128 pieces; an
#     emulator's VLEN is only what it was configured with, so it is not used.
#   - If the probe fails, or a cross-compile has no emulator to run it, the
#     compiler's default is kept.
function(_pffft_riscv_auto_march lang out_march out_reason)
    set(${out_march} "" PARENT_SCOPE)

    # a cross toolchain may not link at configure time
    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    set(CMAKE_REQUIRED_QUIET ON)
    # try_compile() uses the flags of CMAKE_TRY_COMPILE_CONFIGURATION, else
    # of Debug; follow the build type of single-config generators, as
    # _pffft_probe_simd_arch() does
    get_property(_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if (NOT CMAKE_TRY_COMPILE_CONFIGURATION AND NOT _multi_config AND CMAKE_BUILD_TYPE)
        set(CMAKE_TRY_COMPILE_CONFIGURATION "${CMAKE_BUILD_TYPE}")
    endif()
    if (CMAKE_TRY_COMPILE_CONFIGURATION)
        string(TOUPPER "${CMAKE_TRY_COMPILE_CONFIGURATION}" _config)
    else()
        set(_config "DEBUG")
    endif()
    # the result is cached; check again when the compiler or its flags change
    set(_has_v PFFFT_RISCV_${lang}_COMPILER_HAS_V)
    set(_sig "${CMAKE_${lang}_COMPILER};${CMAKE_${lang}_COMPILER_TARGET};${CMAKE_SYSROOT};${CMAKE_${lang}_FLAGS};${_config};${CMAKE_${lang}_FLAGS_${_config}}")
    if (NOT "${${_has_v}_SIG}" STREQUAL "${_sig}")
        unset(${_has_v} CACHE)
    endif()
    set(_src "#ifndef __riscv_vector\n#error no V\n#endif\nint pffft_rvv_check;\n")
    if (lang STREQUAL "C")
        include(CheckCSourceCompiles)
        check_c_source_compiles("${_src}" ${_has_v})
    else()
        include(CheckCXXSourceCompiles)
        check_cxx_source_compiles("${_src}" ${_has_v})
    endif()
    set(${_has_v}_SIG "${_sig}" CACHE INTERNAL "inputs of ${_has_v}")
    if (${_has_v})
        set(${out_reason} "the compiler already enables the vector extension" PARENT_SCOPE)
        return()
    endif()

    if (CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
        set(${out_reason} "cross-compiling without CMAKE_CROSSCOMPILING_EMULATOR, so the CPU cannot be probed" PARENT_SCOPE)
        return()
    endif()

    # what the probe found about the CPU is cached; run it again for another
    # compiler, C flags, or emulator, including another qemu CPU model
    set(_sig "${CMAKE_C_COMPILER};${CMAKE_C_FLAGS};${_config};${CMAKE_C_FLAGS_${_config}};${CMAKE_CROSSCOMPILING_EMULATOR};$ENV{QEMU_CPU}")
    if (NOT "${PFFFT_RVV_HAS_V_SIG}" STREQUAL "${_sig}")
        unset(PFFFT_RVV_HAS_V CACHE)
        unset(PFFFT_RVV_VLEN CACHE)
    endif()
    if (NOT DEFINED PFFFT_RVV_HAS_V)
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
        if (NOT PFFFT_RVV_PROBE_COMPILED)
            # no answer about the CPU, e.g. flags that break the link: cache
            # nothing, so that the next configure probes again
            set(${out_reason} "the vector extension probe could not be built (see the CMake configure log)" PARENT_SCOPE)
            return()
        endif()
        # the output may also hold emulator chatter (qemu warns on stderr
        # when no vext_spec is given); the VLEN is the last number printed
        string(REGEX MATCH "[0-9]+[ \t\r\n]*$" _vlen "${_vlen}")
        string(STRIP "${_vlen}" _vlen)
        if (PFFFT_RVV_PROBE_RUN EQUAL 0 AND _vlen)
            set(PFFFT_RVV_HAS_V TRUE CACHE INTERNAL "The build machine supports the RISC-V V extension")
            set(PFFFT_RVV_VLEN "${_vlen}" CACHE INTERNAL "VLEN of the build machine in bits")
        else()
            set(PFFFT_RVV_HAS_V FALSE CACHE INTERNAL "The build machine lacks the RISC-V V extension")
        endif()
        set(PFFFT_RVV_HAS_V_SIG "${_sig}" CACHE INTERNAL "inputs of PFFFT_RVV_HAS_V")
    endif()

    if (NOT PFFFT_RVV_HAS_V)
        set(${out_reason} "the build machine has no vector extension" PARENT_SCOPE)
    elseif (NOT CMAKE_CROSSCOMPILING AND PFFFT_RVV_VLEN GREATER_EQUAL 256)
        set(${out_march} "rv64gcv_zvl${PFFFT_RVV_VLEN}b" PARENT_SCOPE)
        set(${out_reason} "the build machine has the vector extension, VLEN ${PFFFT_RVV_VLEN}" PARENT_SCOPE)
    else()
        set(${out_march} "rv64gcv" PARENT_SCOPE)
        set(${out_reason} "the build machine has the vector extension" PARENT_SCOPE)
    endif()
endfunction()

# _pffft_active_configs(<out_var>)
# The configurations being built, in upper case: CMAKE_BUILD_TYPE for
# single-config generators, CMAKE_CONFIGURATION_TYPES for multi-config ones.
function(_pffft_active_configs out_var)
    get_property(_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if (_multi_config)
        set(_configs ${CMAKE_CONFIGURATION_TYPES})
    else()
        set(_configs ${CMAKE_BUILD_TYPE})
    endif()
    string(TOUPPER "${_configs}" _configs)
    set(${out_var} "${_configs}" PARENT_SCOPE)
endfunction()

# _pffft_eval_compile_options(<C|CXX> <configs> <options> <out_applied> <out_unevaluated>)
# Of a COMPILE_OPTIONS list, the options the <lang> compiler gets, as far as
# that is known at configure time: SHELL: groups are split as CMake does, and
# $<$<COMPILE_LANGUAGE:...>:...> and $<$<CONFIG:...>:...> around a single
# option are evaluated. A CONFIG option applies if it names one of <configs>,
# a list of configurations in upper case. What other
# generator expressions hold depends on what CMake knows only when it
# generates the build; their list items go to <out_unevaluated> unchanged.
function(_pffft_eval_compile_options lang configs options out_applied out_unevaluated)
    set(_applied "")
    set(_unevaluated "")
    # one generator expression can span several list items, so track the
    # nesting depth, as _pffft_strip_genex() does
    set(_depth 0)
    foreach(_opt IN LISTS options)
        if (_depth EQUAL 0 AND _opt MATCHES "^\\$<\\$<(COMPILE_LANGUAGE|CONFIG):([A-Za-z0-9_,]+)>:([^$<>]*)>$")
            # CONFIG compares case-insensitively; the option applies if one
            # of the configurations being built matches
            if (CMAKE_MATCH_1 STREQUAL "CONFIG")
                string(TOUPPER "${CMAKE_MATCH_2}" _values)
                set(_wanted ${configs})
            else()
                set(_values "${CMAKE_MATCH_2}")
                set(_wanted ${lang})
            endif()
            set(_opt "${CMAKE_MATCH_3}")
            string(REPLACE "," ";" _values "${_values}")
            set(_applies FALSE)
            foreach(_value IN LISTS _values)
                if (_value IN_LIST _wanted)
                    set(_applies TRUE)
                endif()
            endforeach()
            if (NOT _applies)
                continue()
            endif()
        elseif (_depth GREATER 0 OR _opt MATCHES "\\$<")
            string(REGEX MATCHALL "\\$<" _open "${_opt}")
            string(REGEX MATCHALL ">" _close "${_opt}")
            list(LENGTH _open _n_open)
            list(LENGTH _close _n_close)
            math(EXPR _depth "${_depth} + ${_n_open} - ${_n_close}")
            list(APPEND _unevaluated "${_opt}")
            continue()
        endif()
        if (_opt MATCHES "^SHELL:(.*)")
            separate_arguments(_args UNIX_COMMAND "${CMAKE_MATCH_1}")
            list(APPEND _applied ${_args})
        elseif (NOT _opt STREQUAL "")
            list(APPEND _applied "${_opt}")
        endif()
    endforeach()
    set(${out_applied} "${_applied}" PARENT_SCOPE)
    set(${out_unevaluated} "${_unevaluated}" PARENT_SCOPE)
endfunction()

# _pffft_find_march(<C|CXX> <target> <out_flag> <out_where>)
# The first -march or -mcpu (which implies an -march on RISC-V) that the
# <lang> compiler gets for <target> from outside pffft's arch options: the
# flags, CMAKE_<lang>_FLAGS and those of the configurations being built, split
# as a command line is, and the compile options of <target>, which include
# those of its directory, as _pffft_eval_compile_options() sees them. Empty if
# none. An -march in a generator expression that it cannot evaluate still
# wins, as auto puts its -march first among the target's compile options.
function(_pffft_find_march lang target out_flag out_where)
    set(${out_flag} "" PARENT_SCOPE)
    _pffft_active_configs(_configs)
    set(_vars CMAKE_${lang}_FLAGS)
    foreach(_config IN LISTS _configs)
        list(APPEND _vars CMAKE_${lang}_FLAGS_${_config})
    endforeach()
    foreach(_var IN LISTS _vars)
        # split as the compiler's command line is, which removes quotes
        separate_arguments(_args NATIVE_COMMAND "${${_var}}")
        foreach(_arg IN LISTS _args)
            if (_arg MATCHES "^-m(arch|cpu)=")
                set(${out_flag} "${_arg}" PARENT_SCOPE)
                set(${out_where} "${_var}" PARENT_SCOPE)
                return()
            endif()
        endforeach()
    endforeach()
    # a target starts with the compile options of its directory, so these
    # include e.g. add_compile_options() of a parent project
    get_target_property(_opts ${target} COMPILE_OPTIONS)
    if (NOT _opts)
        set(_opts "")
    endif()
    # any configuration being built would get the -march
    _pffft_eval_compile_options(${lang} "${_configs}" "${_opts}" _applied _unevaluated)
    foreach(_opt IN LISTS _applied)
        if (_opt MATCHES "^-m(arch|cpu)=")
            set(${out_flag} "${_opt}" PARENT_SCOPE)
            set(${out_where} "the compile options of target ${target}" PARENT_SCOPE)
            return()
        endif()
    endforeach()
endfunction()

# _pffft_prepend_compile_option(<target> <option>)
# Puts <option> first among the compile options of <target>, which start with
# those of its directory. target_compile_options(BEFORE PRIVATE) would do the
# same, but ignores BEFORE under the OLD behaviour of policy CMP0101, which
# CMake 3.17 introduced and cmake_minimum_required(VERSION 3.10) keeps.
function(_pffft_prepend_compile_option target option)
    get_target_property(_opts ${target} COMPILE_OPTIONS)
    if (NOT _opts)
        set(_opts "")
    endif()
    set_property(TARGET ${target} PROPERTY COMPILE_OPTIONS "${option}" ${_opts})
endfunction()

# _pffft_resolve_arch(<C|CXX> <target> <out_arch>)
# TARGET_<lang>_ARCH, with "auto" resolved and "none" as empty.
function(_pffft_resolve_arch lang target out_arch)
    set(_arch "${TARGET_${lang}_ARCH}")
    if (_arch STREQUAL "none")
        set(_arch "")
    elseif (_arch STREQUAL "auto")
        set(_arch "")
        # elsewhere "auto" acts as "none", with a warning at configure time
        if (PFFFT_TARGET_IS_RISCV64 AND CMAKE_${lang}_COMPILER_ID MATCHES "^(GNU|Clang)$")
            # auto chooses the -march itself: one from elsewhere contradicts it
            _pffft_find_march(${lang} ${target} _flag _where)
            if (_flag)
                message(FATAL_ERROR "TARGET_${lang}_ARCH=auto conflicts with ${_flag} in ${_where}:"
                                    " auto chooses the -march itself. Set TARGET_${lang}_ARCH=none to"
                                    " build with ${_flag}, or remove ${_flag} to let auto choose.")
            endif()
            _pffft_riscv_auto_march(${lang} _arch _reason)
            if (_arch)
                message(STATUS "TARGET_${lang}_ARCH=auto for target ${target}: -march=${_arch}, as ${_reason}")
            else()
                message(STATUS "TARGET_${lang}_ARCH=auto for target ${target}: the compiler's default, as ${_reason}")
            endif()
        endif()
    endif()
    set(${out_arch} "${_arch}" PARENT_SCOPE)
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
    _pffft_resolve_arch(C ${target} _arch)
    if ("${_arch}" STREQUAL "")
        message(STATUS "TARGET_C_ARCH not set for target ${target}, using default")
    else()
        if ( (CMAKE_C_COMPILER_ID STREQUAL "GNU") OR (CMAKE_C_COMPILER_ID STREQUAL "Clang") OR (CMAKE_C_COMPILER_ID STREQUAL "AppleClang") )
            # auto's -march goes first, so that one from the target's compile
            # options, which auto could not see, comes later and wins
            if ("${TARGET_C_ARCH}" STREQUAL "auto")
                _pffft_prepend_compile_option(${target} "-march=${_arch}")
            else()
                target_compile_options(${target} PRIVATE "-march=${_arch}")
            endif()
            message(STATUS "C ARCH for target ${target} set: ${_arch}")
        elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
            target_compile_options(${target} PRIVATE "/arch:${_arch}")
            message(STATUS "C ARCH for target ${target} set: ${_arch}")
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
    _pffft_resolve_arch(CXX ${target} _arch)
    if ("${_arch}" STREQUAL "")
        message(STATUS "TARGET_CXX_ARCH not set for target ${target}, using default")
    else()
        if ( (CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang") )
            # auto's -march goes first, so that one from the target's compile
            # options, which auto could not see, comes later and wins
            if ("${TARGET_CXX_ARCH}" STREQUAL "auto")
                _pffft_prepend_compile_option(${target} "-march=${_arch}")
            else()
                target_compile_options(${target} PRIVATE "-march=${_arch}")
            endif()
            message(STATUS "C++ ARCH for target ${target} set: ${_arch}")
        elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            target_compile_options(${target} PRIVATE "/arch:${_arch}")
            message(STATUS "C++ ARCH for target ${target} set: ${_arch}")
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
