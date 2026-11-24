
#include "pf_conv_dispatcher.h"

#if 0
#include <stdio.h>

#define DPRINT(...) fprintf(stderr, __VA_ARGS__)

#else
#define DPRINT(...) do { } while (0)
#endif


/* Which of the architecture variants can this CPU run? Each is compiled
 * with its own -march / /arch (benchmarks/pf_conv.cmake), so the compiler
 * may use every instruction that implies. A variant the CPU lacks is left
 * as nullptr: running it would end in an illegal instruction. */

#if defined(CONV_ARCH_GCC_AMD64) || defined(CONV_ARCH_MSVC_AMD64)
#  if defined(_MSC_VER)
#    include <intrin.h>
#    include <immintrin.h>
#  endif

struct x86_features { bool ssse3, sse42, popcnt, avx, avx2, fma, bmi1, bmi2; };

static x86_features x86_cpu_features()
{
    x86_features f;
#  if defined(_MSC_VER)
    int r[4];
    __cpuid(r, 0);
    const int max_leaf = r[0];
    __cpuid(r, 1);
    const int ecx1 = r[2];
    // AVX also needs the OS to save the YMM registers: XCR0 bits 1 and 2
    const bool os_avx = (ecx1 & (1 << 27)) && ((_xgetbv(0) & 6) == 6);
    int ebx7 = 0;
    if (max_leaf >= 7)
    {
        __cpuidex(r, 7, 0);
        ebx7 = r[1];
    }
    f.ssse3  = (ecx1 & (1 << 9)) != 0;
    f.sse42  = (ecx1 & (1 << 20)) != 0;
    f.popcnt = (ecx1 & (1 << 23)) != 0;
    f.avx    = os_avx && (ecx1 & (1 << 28));
    f.fma    = os_avx && (ecx1 & (1 << 12));
    f.avx2   = os_avx && (ebx7 & (1 << 5));
    f.bmi1   = (ebx7 & (1 << 3)) != 0;
    f.bmi2   = (ebx7 & (1 << 8)) != 0;
#  else
    // gcc's and clang's runtime checks AVX's OS support as well
    __builtin_cpu_init();
    f.ssse3  = __builtin_cpu_supports("ssse3");
    f.sse42  = __builtin_cpu_supports("sse4.2");
    f.popcnt = __builtin_cpu_supports("popcnt");
    f.avx    = __builtin_cpu_supports("avx");
    f.avx2   = __builtin_cpu_supports("avx2");
    f.fma    = __builtin_cpu_supports("fma");
    f.bmi1   = __builtin_cpu_supports("bmi");
    f.bmi2   = __builtin_cpu_supports("bmi2");
#  endif
    return f;
}
#endif

#if defined(CONV_ARCH_GCC_ARM32NEON) && defined(__linux__)
#  include <sys/auxv.h>
#  include <string.h>

static bool arm32_has_vfpv4()
{
    const unsigned long HWCAP_NEON_ = 1UL << 12, HWCAP_VFPv4_ = 1UL << 16;
    const unsigned long hwcap = getauxval(AT_HWCAP);
    return (hwcap & HWCAP_NEON_) && (hwcap & HWCAP_VFPv4_);
}

/* -mfpu=neon-fp-armv8 needs an ARMv8 core. A 32-bit kernel reports those as
 * "v7l" (AT_PLATFORM), but sets the HWCAP2 bits of the ARMv8 extensions,
 * AES, PMULL, SHA1, SHA2 and CRC32, where the core has one: CRC32 is in
 * every Cortex-A53 and A72 (Raspberry Pi 3 and 4). A 64-bit kernel running
 * 32-bit code reports "v8l". An ARMv8.0 core with none of the extensions
 * under a 32-bit kernel is taken for ARMv7, which is safe. */
static bool arm32_is_armv8()
{
    const char * platform = reinterpret_cast<const char *>(getauxval(AT_PLATFORM));
    bool armv8 = platform && !strncmp(platform, "v8", 2);
#  ifdef AT_HWCAP2
    const unsigned long HWCAP2_ARMV8_ = 0x1F;  // AES, PMULL, SHA1, SHA2, CRC32
    armv8 = armv8 || (getauxval(AT_HWCAP2) & HWCAP2_ARMV8_) != 0;
#  endif
    return arm32_has_vfpv4() && armv8;
}
#elif defined(CONV_ARCH_GCC_ARM32NEON)
// no way to ask: assume the CPU is what the build was configured for
static bool arm32_has_vfpv4() { return true; }
static bool arm32_is_armv8() { return true; }
#endif


#define N_DEFAULT_ARCHES  2
// 0 is "none"
// 1 "dflt"
// the others may be nullptr, where the CPU cannot run them

ptr_to_conv_f_ptrs * get_all_conv_arch_ptrs(int * p_num_arch)
{
    static ptr_to_conv_f_ptrs * all_arches = nullptr;
    static int n_arch = 0;
    if (!all_arches)
    {
        n_arch = N_DEFAULT_ARCHES;
#if defined(CONV_ARCH_GCC_AMD64)
        static const conv_f_ptrs *conv_arch_ptrs[N_DEFAULT_ARCHES+4] = {0};
        DPRINT("CONV_ARCH_GCC_AMD64: sse3, sse4, avx, avx2\n");
        const x86_features f = x86_cpu_features();
        // -march=core2, nehalem, sandybridge, haswell
        conv_arch_ptrs[n_arch++] = f.ssse3 ? CONV_FN_ARCH(conv_ptrs, sse3)() : nullptr;
        conv_arch_ptrs[n_arch++] = (f.sse42 && f.popcnt) ? CONV_FN_ARCH(conv_ptrs, sse4)() : nullptr;
        conv_arch_ptrs[n_arch++] = (f.sse42 && f.popcnt && f.avx) ? CONV_FN_ARCH(conv_ptrs, avx) () : nullptr;
        // haswell also has F16C, LZCNT and MOVBE, which every CPU with AVX2,
        // FMA and BMI2 has too, and which not every compiler can ask for
        conv_arch_ptrs[n_arch++] = (f.sse42 && f.popcnt && f.avx2 && f.fma && f.bmi1 && f.bmi2) ? CONV_FN_ARCH(conv_ptrs, avx2)() : nullptr;
#elif defined(CONV_ARCH_MSVC_AMD64)
        static const conv_f_ptrs *conv_arch_ptrs[N_DEFAULT_ARCHES+3] = {0};
        DPRINT("CONV_ARCH_MSVC_AMD64: sse2, avx, avx2\n");
        const x86_features f = x86_cpu_features();
        // /arch:SSE2, AVX, AVX2 (which includes FMA and BMI)
        conv_arch_ptrs[n_arch++] = CONV_FN_ARCH(conv_ptrs, sse2)();
        conv_arch_ptrs[n_arch++] = f.avx ? CONV_FN_ARCH(conv_ptrs, avx) () : nullptr;
        conv_arch_ptrs[n_arch++] = (f.avx2 && f.fma && f.bmi1 && f.bmi2) ? CONV_FN_ARCH(conv_ptrs, avx2)() : nullptr;
#elif defined(CONV_ARCH_GCC_ARM32NEON)
        static const conv_f_ptrs *conv_arch_ptrs[N_DEFAULT_ARCHES+3] = {0};
        DPRINT("CONV_ARCH_GCC_ARM32NEON: neon_vfpv4, neon_rpi3_a53, neon_rpi4_a72\n");
        // -mfpu=neon-vfpv4, neon-vfpv4, neon-fp-armv8
        conv_arch_ptrs[n_arch++] = arm32_has_vfpv4() ? CONV_FN_ARCH(conv_ptrs, neon_vfpv4)() : nullptr;
        conv_arch_ptrs[n_arch++] = arm32_has_vfpv4() ? CONV_FN_ARCH(conv_ptrs, neon_rpi3_a53)() : nullptr;
        conv_arch_ptrs[n_arch++] = arm32_is_armv8() ? CONV_FN_ARCH(conv_ptrs, neon_rpi4_a72)() : nullptr;
#elif defined(CONV_ARCH_GCC_AARCH64)
        static const conv_f_ptrs *conv_arch_ptrs[N_DEFAULT_ARCHES+1] = {0};
        DPRINT("CONV_ARCH_GCC_AARCH64: -\n");
        conv_arch_ptrs[n_arch++] = CONV_FN_ARCH(conv_ptrs, armv8a)();
#else
        static const conv_f_ptrs *conv_arch_ptrs[N_DEFAULT_ARCHES] = {0};
        DPRINT("unknown CONV_ARCH: -\n");
#endif
        conv_arch_ptrs[0] = CONV_FN_ARCH(conv_ptrs, none)();
        conv_arch_ptrs[1] = CONV_FN_ARCH(conv_ptrs, dflt)();
        all_arches = conv_arch_ptrs;
    }
    if (p_num_arch)
        *p_num_arch = n_arch;
    return all_arches;
}

