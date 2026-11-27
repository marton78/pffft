# Moves the code of the pffft library by a pad, for comparing revisions at
# matched code layouts (compare_many.py --pads). It works on any revision, as
# it changes no source of pffft: compare_many.py passes it to pffft's CMake as
#
#   -DCMAKE_PROJECT_PRETTY_FAST_FFT_INCLUDE=<this file> -DPFFFT_LAYOUT_PAD=<bytes>
#
# and CMake runs it right after pffft's project(). With a pad of 0 it changes
# nothing in the library: the build is the plain one.
#
# GCC and Clang: the first C file of the library gets a header forced in
# (-include) whose top-level asm puts PFFFT_LAYOUT_PAD bytes at the start of
# its object's code section. The functions after it, in this object and the
# ones linked after it, move by the pad, rounded up only to their own
# alignment. A pad in front of the objects, as a
# separate object, would not do: the linker starts each object's code
# section at that section's alignment (32 or 64 bytes with GCC for some
# CPUs), which would round every pad to a multiple of it, so that loops
# never moved relative to the 32- and 64-byte blocks that decoders see.
#
# MSVC: one more object, first in the library, with a function of that many
# bytes (no-operation instructions) in section .text$a. The linker sorts the
# parts of a section by the name after the $, so it comes before the code of
# every object (.text$mn), whose functions are aligned to 16 bytes each.
# /include keeps the function, which nothing calls. (A pad of data would go to
# a section of its own: the linker merges only parts of the same kind.)
#
# The objects' code is the same at every pad. The library's symbol addresses
# show where the functions went: nm reads them (layout.py), and with MSVC the
# linker writes them to pffft.map in the build directory (/MAP, which changes
# nothing in the library).
#
# -DPFFFT_LAYOUT_ALIGN=16 (GCC only) builds the library with functions and
# loops aligned to 16 bytes exactly: -falign-functions=16 -falign-loops=16
# -falign-jumps=1 -falign-labels=1. GCC aligns them by default to 32 or 64
# for many CPUs, and only where that skips few bytes (".p2align 5,,15"), so
# that the alignment takes up much of a pad: the hot functions then keep
# their offsets modulo 32 or 64 at every pad, and the pads vary only the
# address bits above. Give it to every label of a comparison, pad 0 too: it
# is another build, with its own layouts.
#
# With CMake 3.19 or later it acts once pffft's CMakeLists.txt is done
# (cmake_language(DEFER)); older ones (from 3.13) get add_library wrapped,
# and so does -DPFFFT_LAYOUT_PAD_WRAP=ON, which tests that way.

if (NOT DEFINED PFFFT_LAYOUT_PAD)
  set(PFFFT_LAYOUT_PAD 0)
endif()
if (NOT PFFFT_LAYOUT_PAD MATCHES "^[0-9]+$")
  message(FATAL_ERROR "PFFFT_LAYOUT_PAD: expected a number of bytes, got '${PFFFT_LAYOUT_PAD}'")
endif()

# write a file only when its text changes, so that it causes no rebuild
function(pffft_layout_pad_write path text)
  set(old "")
  if (EXISTS "${path}")
    file(READ "${path}" old)
  endif()
  if (NOT old STREQUAL text)
    file(WRITE "${path}" "${text}")
  endif()
endfunction()

set(PFFFT_LAYOUT_PAD_FILE "")
if (PFFFT_LAYOUT_PAD GREATER 0)
  set(head "/* made by benchmarks/layout_pad.cmake: ${PFFFT_LAYOUT_PAD} bytes before the code of pffft */\n")
  if (MSVC)
    if (CMAKE_C_COMPILER_ARCHITECTURE_ID MATCHES "ARM" OR CMAKE_SYSTEM_PROCESSOR MATCHES "^(ARM|arm|aarch)")
      math(EXPR nops "${PFFFT_LAYOUT_PAD} / 4 - 1")  # 4 bytes each, and the return
    else()
      # 1 byte each, and the return; MSVC begins an x64 function with an
      # instruction of 2 bytes at least (here 66 90), so that it can be patched
      math(EXPR nops "${PFFFT_LAYOUT_PAD} - 3")
    endif()
    if (nops LESS 0)
      set(nops 0)
    endif()
    string(REPEAT "  __nop();\n" ${nops} body)
    set(PFFFT_LAYOUT_PAD_FILE "${CMAKE_BINARY_DIR}/pffft_layout_pad.c")
    pffft_layout_pad_write("${PFFFT_LAYOUT_PAD_FILE}" "${head}#include <intrin.h>
#pragma code_seg(\".text$a\")
__declspec(noinline) void pffft_layout_pad(void) {
${body}}
#pragma code_seg()
#pragma comment(linker, \"/include:pffft_layout_pad\")
")
  else()
    # A label of its own makes it an atom for ld64 (Mach-O), which keeps it.
    # With ELF and -ffunction-sections (the Android NDK's flags), the object's
    # .text holds only the pad, and --gc-sections drops it, as nothing
    # references it: there it gets a section of its own, the first of its
    # object, with the flag R (SHF_GNU_RETAIN), where the assembler has it.
    # Without function sections the pad stays in .text: GNU ld would put a
    # section of its own after the object's .text, so that the first
    # object's code did not move. %progbits: on 32-bit Arm @ begins a comment.
    set(pad_section ".text")
    if (NOT APPLE AND "${CMAKE_C_FLAGS}" MATCHES "-ffunction-sections" AND NOT DEFINED PFFFT_LAYOUT_PAD_RETAIN)
      set(probe "${CMAKE_BINARY_DIR}/pffft_layout_pad_retain.c")
      file(WRITE "${probe}" [[__asm__(".section .text.pffft_layout_pad_probe,\"axR\",%progbits\n.text\n");
int main(void) { return 0; }
]])
      try_compile(retain "${CMAKE_BINARY_DIR}/pffft_layout_pad_retain" SOURCES "${probe}")
      set(PFFFT_LAYOUT_PAD_RETAIN ${retain} CACHE INTERNAL "the assembler takes the section flag R")
    endif()
    if (PFFFT_LAYOUT_PAD_RETAIN AND "${CMAKE_C_FLAGS}" MATCHES "-ffunction-sections")
      set(pad_section ".section .text.pffft_layout_pad,\\\"axR\\\",%progbits")
    endif()
    set(PFFFT_LAYOUT_PAD_FILE "${CMAKE_BINARY_DIR}/pffft_layout_pad.h")
    pffft_layout_pad_write("${PFFFT_LAYOUT_PAD_FILE}" "${head}__asm__(\"${pad_section}\\n\"
        \"pffft_layout_pad:\\n\"
        \".skip ${PFFFT_LAYOUT_PAD}\\n\"
        \".text\\n\");
")
  endif()
endif()

function(pffft_layout_pad_apply)
  if (NOT TARGET PFFFT)
    message(FATAL_ERROR "layout_pad.cmake: this revision has no library target PFFFT")
  endif()
  if (PFFFT_LAYOUT_PAD_FILE AND MSVC)
    get_target_property(srcs PFFFT SOURCES)
    set_property(TARGET PFFFT PROPERTY SOURCES "${PFFFT_LAYOUT_PAD_FILE}" ${srcs})
  elseif (PFFFT_LAYOUT_PAD_FILE)
    # the first object only: everything after it moves by the pad
    get_target_property(srcs PFFFT SOURCES)
    set(first "")
    foreach (src ${srcs})
      if (src MATCHES "\\.c$" AND NOT first)
        set(first "${src}")
      endif()
    endforeach()
    set_property(SOURCE "${first}" APPEND PROPERTY COMPILE_OPTIONS -include "${PFFFT_LAYOUT_PAD_FILE}")
  endif()
  if (MSVC)
    target_link_options(PFFFT PRIVATE "/MAP:${CMAKE_BINARY_DIR}/pffft.map")
  endif()
  if (PFFFT_LAYOUT_ALIGN)
    if (NOT CMAKE_C_COMPILER_ID STREQUAL "GNU")
      message(FATAL_ERROR "PFFFT_LAYOUT_ALIGN: for GCC only")
    endif()
    target_compile_options(PFFFT PRIVATE -falign-functions=${PFFFT_LAYOUT_ALIGN} -falign-loops=${PFFFT_LAYOUT_ALIGN}
                           -falign-jumps=1 -falign-labels=1)
    # x86: newer GCC aligns small loops to a cache line regardless of -falign-loops
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-mno-align-tight-loops PFFFT_LAYOUT_HAS_TIGHT_LOOPS)
    if (PFFFT_LAYOUT_HAS_TIGHT_LOOPS)
      target_compile_options(PFFFT PRIVATE -mno-align-tight-loops)
    endif()
  endif()
endfunction()

if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.19 AND NOT PFFFT_LAYOUT_PAD_WRAP)
  # after the CMakeLists.txt of pffft has made its targets
  cmake_language(DEFER DIRECTORY "${CMAKE_SOURCE_DIR}" CALL pffft_layout_pad_apply)
else()
  # older CMake (3.13 on): right after add_library(PFFFT ...) itself
  macro(add_library name)
    _add_library(${ARGV})
    if ("${name}" STREQUAL "PFFFT")
      pffft_layout_pad_apply()
    endif()
  endmacro()
endif()
