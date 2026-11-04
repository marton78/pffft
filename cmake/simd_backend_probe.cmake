# provided:
#   - function: pffft_report_simd_backend(<target>)
#
# The SIMD backend is chosen by the preprocessor in pf_float.h / pf_double.h.
# To report it without guessing, compile a probe with the target's flags:
#   const char pffft_simd_probe[] = "PFFFT_SIMD_ARCH[" VARCH "]";
# The literals concatenate to e.g. "PFFFT_SIMD_ARCH[NEON]", stored verbatim in
# the object's read-only data. The probe is built as a static library, which is
# only the uncompressed object files strung together, so file(STRINGS) can scan
# its raw bytes for printable runs, like the Unix strings tool, and the marker
# regex picks out ours. No object format is parsed and nothing is run, so this
# works for any cross target. CMake finds e.g. the pointer size of the target
# the same way, see Modules/CMakeCompilerABI.h.
#
# Single-config generators probe CMAKE_BUILD_TYPE. Multi-config generators
# probe one configuration, CMAKE_TRY_COMPILE_CONFIGURATION or Debug, and the
# report names it, since CMAKE_C_FLAGS_<CONFIG> may differ between them.

# One generator expression can span several list items, e.g.
#   $<$<C_COMPILER_ID:GNU>:-Wall;-Wextra;-pedantic>
# so track the nesting depth. A leftover "-pedantic>" reaches the shell in
# Makefile builds, which takes ">" as a redirect and eats the next flag.
function(_pffft_strip_genex in_list out_var)
    set(_out "")
    set(_depth 0)
    foreach(_item IN LISTS in_list)
        if (_depth GREATER 0 OR _item MATCHES "\\$<")
            string(REGEX MATCHALL "\\$<" _open "${_item}")
            string(REGEX MATCHALL ">" _close "${_item}")
            list(LENGTH _open _n_open)
            list(LENGTH _close _n_close)
            math(EXPR _depth "${_depth} + ${_n_open} - ${_n_close}")
        elseif (NOT _item STREQUAL "")
            list(APPEND _out "${_item}")
        endif()
    endforeach()
    set(${out_var} "${_out}" PARENT_SCOPE)
endfunction()

function(_pffft_probe_simd_arch target header out_var)
    get_target_property(_opts ${target} COMPILE_OPTIONS)
    get_target_property(_defs ${target} COMPILE_DEFINITIONS)
    if (NOT _opts)
        set(_opts "")
    endif()
    if (NOT _defs)
        set(_defs "")
    endif()
    # add_compile_definitions() of this or a parent directory, e.g. a project
    # using pffft via add_subdirectory(); the build applies them on top
    get_directory_property(_dir_defs COMPILE_DEFINITIONS)
    list(APPEND _defs ${_dir_defs})
    _pffft_strip_genex("${_opts}" _flags)
    _pffft_strip_genex("${_defs}" _defs)
    foreach(_d IN LISTS _defs)
        list(APPEND _flags "-D${_d}")
    endforeach()

    # try_compile() builds CMAKE_TRY_COMPILE_CONFIGURATION, Debug if unset,
    # not CMAKE_BUILD_TYPE; follow the build type of single-config generators
    get_property(_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
    if (NOT CMAKE_TRY_COMPILE_CONFIGURATION AND NOT _multi_config AND CMAKE_BUILD_TYPE)
        set(CMAKE_TRY_COMPILE_CONFIGURATION "${CMAKE_BUILD_TYPE}")
    endif()
    if (CMAKE_TRY_COMPILE_CONFIGURATION)
        set(_config "${CMAKE_TRY_COMPILE_CONFIGURATION}")
    else()
        set(_config "Debug")
    endif()
    string(TOUPPER "${_config}" _config_upper)

    # besides flags, try_compile() forwards these, which can change the target
    set(_sig "${CMAKE_C_COMPILER};${CMAKE_C_COMPILER_TARGET};${CMAKE_OSX_ARCHITECTURES};${CMAKE_OSX_SYSROOT};${CMAKE_SYSROOT}")
    foreach(_v IN LISTS CMAKE_TRY_COMPILE_PLATFORM_VARIABLES)
        string(APPEND _sig ";${_v}=${${_v}}")
    endforeach()
    string(APPEND _sig ";${CMAKE_C_FLAGS};${_config};${CMAKE_C_FLAGS_${_config_upper}};${_flags};${header}")
    # and the headers themselves, e.g. after updating pffft in a build tree
    file(GLOB _simd_headers "${CMAKE_CURRENT_SOURCE_DIR}/src/simd/*.h")
    foreach(_h IN LISTS _simd_headers)
        file(SHA1 "${_h}" _h_sha1)
        string(APPEND _sig ";${_h_sha1}")
    endforeach()
    set(_cache "PFFFT_SIMD_PROBE_${out_var}")
    set(${out_var}_CONFIG "${_config}" PARENT_SCOPE)
    if (DEFINED ${_cache} AND "${${_cache}_SIG}" STREQUAL "${_sig}")
        set(${out_var} "${${_cache}}" PARENT_SCOPE)
        return()
    endif()

    set(_dir "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/pffft_simd_probe")
    set(_src "${_dir}/probe_${out_var}.c")
    # the headers expect ALWAYS_INLINE from the including .c file
    file(WRITE "${_src}"
        "#if defined(_MSC_VER)\n"
        "#  define ALWAYS_INLINE(return_type) __forceinline return_type\n"
        "#else\n"
        "#  define ALWAYS_INLINE(return_type) inline return_type\n"
        "#endif\n"
        "#include \"simd/${header}\"\n"
        "const char pffft_simd_probe[] = \"PFFFT_SIMD_ARCH[\" VARCH \"]\";\n")

    # no link step: a cross toolchain may not link at configure time
    set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
    # LTO stores bitcode instead of plain data bytes, which file(STRINGS)
    # cannot see into; LTO does not change the predefined SIMD macros
    set(_probe_flags ${_flags})
    if (CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        list(APPEND _probe_flags "-fno-lto")
    endif()
    # a universal Apple build has one string per slice, in no reliable order;
    # probe each architecture on its own to report which slice gets what
    set(_slices "-")
    list(LENGTH CMAKE_OSX_ARCHITECTURES _n_slices)
    if (APPLE AND _n_slices GREATER 1)
        set(_slices ${CMAKE_OSX_ARCHITECTURES})
    endif()
    set(_results "")
    foreach(_slice IN LISTS _slices)
        set(_out "${_dir}/probe_${out_var}")
        if (NOT _slice STREQUAL "-")
            set(CMAKE_OSX_ARCHITECTURES "${_slice}")
            string(APPEND _out "_${_slice}")
        endif()
        try_compile(_ok "${_dir}/build_${out_var}" "${_src}"
            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${CMAKE_CURRENT_SOURCE_DIR}/src"
            COMPILE_DEFINITIONS ${_probe_flags}
            OUTPUT_VARIABLE _log
            COPY_FILE "${_out}.bin")

        if (NOT _ok)
            file(WRITE "${_out}.log" "${_log}")
            set(_arch "unknown (probe did not compile, see ${_out}.log)")
        else()
            set(_arch "unknown")
            file(STRINGS "${_out}.bin" _hits REGEX "PFFFT_SIMD_ARCH\\[[^]]*\\]")
            if (_hits MATCHES "PFFFT_SIMD_ARCH\\[([^]]*)\\]")
                set(_arch "${CMAKE_MATCH_1}")
            endif()
        endif()
        if (NOT _slice STREQUAL "-")
            string(APPEND _arch " (${_slice})")
        endif()
        list(APPEND _results "${_arch}")
    endforeach()
    string(REPLACE ";" ", " _arch "${_results}")

    set(${_cache} "${_arch}" CACHE INTERNAL "pffft SIMD backend from ${header}")
    set(${_cache}_SIG "${_sig}" CACHE INTERNAL "inputs of the pffft SIMD probe")
    set(${out_var} "${_arch}" PARENT_SCOPE)
endfunction()

function(pffft_report_simd_backend target)
    set(_parts "")
    if (PFFFT_USE_TYPE_FLOAT)
        _pffft_probe_simd_arch(${target} pf_float.h FLOAT)
        list(APPEND _parts "float: ${FLOAT}")
        set(_config "${FLOAT_CONFIG}")
    endif()
    if (PFFFT_USE_TYPE_DOUBLE)
        _pffft_probe_simd_arch(${target} pf_double.h DOUBLE)
        list(APPEND _parts "double: ${DOUBLE}")
        set(_config "${DOUBLE_CONFIG}")
    endif()
    if (_parts)
        string(REPLACE ";" ", " _parts "${_parts}")
        get_property(_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
        if (_multi_config)
            set(_parts "${_parts} (${_config} configuration)")
        endif()
        message(STATUS "Detected SIMD backend for ${target}: ${_parts}")
    endif()
endfunction()
