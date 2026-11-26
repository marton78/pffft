/*
  The thresholds that pffft_new_setup() takes on each kind of CPU
  (pffft_tuning_t in pffft_types.h), as pffft_default_tuning() picks them.

  pffft_common.c includes this file, identifies the CPU at run time and
  takes the first row of the table below that matches it:
  - x86: the vendor of cpuid leaf 0 and family << 8 | model of leaf 1,
    e.g. 0x1901 for an AMD Zen 3 EPYC (family 0x19, model 0x01). The cache
    sizes come from cpuid (leaf 4, or 0x8000001D on AMD and Hygon).
  - Arm: the implementer and the part number of the MIDR register of the
    fastest core (Linux and Android: /proc/cpuinfo, Windows: the registry),
    e.g. implementer 0x41 (Arm) and part 0xd49 (Neoverse N2). The cache
    sizes come from /sys/devices/system/cpu on Linux, which Android often
    hides from apps. Apple's chips match implementer 0x61 (sysctl gives
    their caches, no part number).
  pffft_cpu_describe() prints what it found, and which row matched.

  Each threshold is either a number of bytes, BYTES(b), or a fraction of
  a cache, CACHE(level, num, den, b): num/den of the L1 data cache
  (level 1), the L2 (2) or the L3 (3) of a core (of its cluster or CCD
  where they share one), and b where that size is unknown. NEVER is
  (size_t)-1: never split, never prefetch. KIB (1024 bytes) and MIB
  (1048576 bytes) make sizes readable: BYTES(8 * MIB), CACHE(1, 1, 2, 24 * KIB).

  To add a CPU: measure it with benchmarks/tune.py, which sweeps every
  threshold in one build and prints the row to add, and put that row
  above the rows that would match it now. Before a pull request, please
  run tune.py twice and give both outputs: a row should rest on more
  than one run, and the rows below say which machines they rest on.
*/

#ifndef PFFFT_TUNING_H
#define PFFFT_TUNING_H

#include <stddef.h> /* for size_t */

/* the kind of CPU (cpu_info.arch, tuning_row.arch) */
#define PFFFT_CPU_OTHER 0
#define PFFFT_CPU_X86   1
#define PFFFT_CPU_ARM   2

/* x86 vendors; Arm's are the implementer codes of the MIDR */
#define X86_INTEL 1
#define X86_AMD   2
#define X86_HYGON 3

/* A threshold: bytes, or num/den of the cache of a level (1 to 3), bytes
   where that size is unknown */
typedef struct tuning_rule {
  unsigned char level, num, den;
  size_t bytes;
} tuning_rule;

#define BYTES(b)          {0, 0, 1, (size_t)(b)}
#define CACHE(l, n, d, b) {(l), (n), (d), (size_t)(b)}
#define NEVER             BYTES((size_t)-1)

/* sizes for BYTES and CACHE: BYTES(8 * MIB) */
#define KIB ((size_t)1 << 10)
#define MIB ((size_t)1 << 20)

typedef struct tuning_row {
  const char *name;
  int arch;
  unsigned vendor, id_min, id_max;
  tuning_rule split_max_real, split_max_complex, prefetch_min, radix6_max;
} tuning_row;

/* The columns (tuning_row):
   name, arch, vendor (0: any), first and last id, then the thresholds
   (bytes per array of the transform):
   split_max_real, split_max_complex, prefetch_min, radix6_max.
   radix6_max is 32 KiB everywhere: 0, 8, 16 or 64 KiB were not faster on any
   machine measured (Xeon 6975P-C, Core i7-8700B, EPYC 7763, 9V45, Neoverse
   N2, Apple M2: within -1.4 and +0.4% at the real lengths with a factor 3
   from 1536 to 98304). */
static const tuning_row tuning_rows[] = {
  /* ---- x86 ----
     AMD Zen 1 to 5 (families 0x17 to 0x1A) and Hygon (Zen 1): splitting from
     8 MiB, as their L3 is fast. At n = 65536 to 4194304 against never splitting,
     the best of 256 KiB to 16 MiB (EPYC 7763 Zen 3 with MSVC, 9V74 Zen 4, 9V45
     Zen 5: +7.0, +10.3, +7.9%; 1 MiB +2.8, +8.4, +3.0%). No prefetches: they cost
     the 7763 4 to 9% from 4K to 64K points (#129, 2026 sweep). */
  {"AMD Zen", PFFFT_CPU_X86, X86_AMD, 0x1700, 0x1FFF,
   BYTES(8 * MIB), BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  {"Hygon Dhyana", PFFFT_CPU_X86, X86_HYGON, 0x1800, 0x18FF,
   BYTES(8 * MIB), BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  /* Intel, and every other x86: splitting from 1 MiB (against never splitting,
     n = 65536 to 4194304: Xeon 6975P-C +31.1%, Platinum 8573C +19.5%, 8370C
     +17.3%, Core i7-8700B +14.0%). 512 KiB did as well on average at those
     lengths, but lost 5 to 15% at real arrays of 512 KiB to 1 MiB, mixed radix
     included (f64 real 73728 to 98304, 96000; f32 real 147456 to 192000). The
     prefetches pay from half the L1 data cache on (#129): 24 KiB on the Xeons
     (48 KiB L1d), 16 KiB on the Core i7 (32 KiB, +0.6% from 1024 to 8192 points
     against 24 KiB, f64 real +1.6%). */
  {"x86", PFFFT_CPU_X86, 0, 0, 0xFFFFFFFF,
   BYTES(1 * MIB), BYTES(1 * MIB), CACHE(1, 1, 2, 24 * KIB), BYTES(32 * KIB)},

  /* ---- Arm ----
     Apple: from 8 MiB (M2: real +10.7%, complex +17.8% against never splitting;
     A13 alike). Complex transforms split from 4 MiB gained 1% at 65536 to 10^6
     points and lost 1 to 5% above, where they split once more (M2 and A13). */
  {"Apple", PFFFT_CPU_ARM, 0x61, 0, 0xFFFFFFFF,
   BYTES(8 * MIB), BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  /* Arm Neoverse servers (N1 0xd0c, V1 0xd40, N2 0xd49, V2 0xd4f, V3AE and V3
     0xd83-0xd84, N3 0xd8e; the Cobalt 100 is an N2 by Microsoft): a large, fast
     system cache. Real transforms never split (Neoverse N2, Cobalt 100: no
     threshold of 256 KiB to 16 MiB gained), complex ones from 8 MiB (+2.2 and
     +4.3% against never splitting; 2026: +2.7 and +2.2%). */
  {"Arm Neoverse N1", PFFFT_CPU_ARM, 0x41, 0xd0c, 0xd0c,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  {"Arm Neoverse", PFFFT_CPU_ARM, 0x41, 0xd40, 0xd40,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  {"Arm Neoverse", PFFFT_CPU_ARM, 0x41, 0xd49, 0xd49,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  {"Arm Neoverse", PFFFT_CPU_ARM, 0x41, 0xd4f, 0xd4f,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  {"Arm Neoverse", PFFFT_CPU_ARM, 0x41, 0xd83, 0xd84,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  {"Arm Neoverse", PFFFT_CPU_ARM, 0x41, 0xd8e, 0xd8e,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  {"Microsoft Cobalt 100", PFFFT_CPU_ARM, 0x6d, 0xd49, 0xd49,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
  /* Arm's Cortex-A, Cortex-X and C1 cores (part numbers from Linux's cputype.h),
     as in phones: from 256 KiB (MediaTek MT6878, Cortex-A78: real +22.6%,
     complex +21.1% against never splitting, the best of 128 KiB to 16 MiB; at the
     docs' lengths +0.7% against 1 MiB, mixed-radix ones between 256 KiB and 1 MiB
     -0.2%). A Raspberry Pi 2 (Cortex-A7, 32 bits) gained 17.7% (real) and 7.7%
     (complex) at 256 KiB, and 1.3% more at 128 KiB: small in-order cores may
     deserve a row of their own once more of them are measured. */
  {"Arm Cortex", PFFFT_CPU_ARM, 0x41, 0xc00, 0xd0b,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
  {"Arm Cortex", PFFFT_CPU_ARM, 0x41, 0xd0d, 0xd3f,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
  {"Arm Cortex", PFFFT_CPU_ARM, 0x41, 0xd41, 0xd48,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
  {"Arm Cortex", PFFFT_CPU_ARM, 0x41, 0xd4a, 0xd4e,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
  {"Arm Cortex", PFFFT_CPU_ARM, 0x41, 0xd80, 0xd82,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
  {"Arm Cortex", PFFFT_CPU_ARM, 0x41, 0xd85, 0xd8d,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
  {"Arm C1", PFFFT_CPU_ARM, 0x41, 0xd90, 0xd90,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
  /* other Arm CPUs: on Android as the phones, elsewhere as the servers */
#if defined(__ANDROID__)
  {"Arm, other (Android)", PFFFT_CPU_ARM, 0, 0, 0xFFFFFFFF,
   BYTES(256 * KIB), BYTES(256 * KIB), NEVER, BYTES(32 * KIB)},
#else
  {"Arm, other", PFFFT_CPU_ARM, 0, 0, 0xFFFFFFFF,
   NEVER, BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
#endif

  /* ---- anything else ---- */
  {"other", PFFFT_CPU_OTHER, 0, 0, 0xFFFFFFFF,
   BYTES(8 * MIB), BYTES(8 * MIB), NEVER, BYTES(32 * KIB)},
};

#endif /* PFFFT_TUNING_H */
