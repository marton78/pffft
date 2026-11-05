# Cross-compile for Cortex-M55 (Armv8.1-M with Helium) with arm-none-eabi-gcc
# and newlib, and run the tests on QEMU's mps3-an547 board via semihosting:
#
#   cmake -B build -DCMAKE_TOOLCHAIN_FILE=tests/cortex-m55/toolchain.cmake
#   cmake --build build && ctest --test-dir build -E "^bench_"
#
# The bench_* tests allocate about 1 GiB and do not fit.
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)
set(CMAKE_C_COMPILER arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(_cpu_flags -mcpu=cortex-m55 -mfloat-abi=hard -mthumb)
string(REPLACE ";" " " _cpu_flags_str "${_cpu_flags}")
set(CMAKE_C_FLAGS_INIT "${_cpu_flags_str}")
set(CMAKE_CXX_FLAGS_INIT "${_cpu_flags_str}")
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

# the startup code is linked into every executable
set(_startup "${CMAKE_BINARY_DIR}/cortex_m55_startup.o")
execute_process(COMMAND ${CMAKE_C_COMPILER} ${_cpu_flags} -O2 -c "${CMAKE_CURRENT_LIST_DIR}/startup.c" -o "${_startup}"
                RESULT_VARIABLE _rc)
if (NOT _rc EQUAL 0)
    message(FATAL_ERROR "could not compile ${CMAKE_CURRENT_LIST_DIR}/startup.c")
endif()
set(CMAKE_EXE_LINKER_FLAGS_INIT "--specs=rdimon.specs -T ${CMAKE_CURRENT_LIST_DIR}/mps3-an547.ld ${_startup} -Wl,--no-warn-rwx-segments")
set(CMAKE_C_STANDARD_LIBRARIES "-lm")
set(CMAKE_CXX_STANDARD_LIBRARIES "-lm")
set(CMAKE_CROSSCOMPILING_EMULATOR "${CMAKE_CURRENT_LIST_DIR}/qemu-run.sh")
