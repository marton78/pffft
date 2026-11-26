
#include "pffft.h"

#include <stdint.h>
#include <stdlib.h>

/* SSE and co like 16-bytes aligned pointers
 * with a 64-byte alignment, we are even aligned on L2 cache lines... */
#define MALLOC_V4SF_ALIGNMENT 64
#define EXTRA_BYTES_FOR_ALIGNMENT_AND_A_POINTER (MALLOC_V4SF_ALIGNMENT - 1 + sizeof(void*))

static void * Valigned_malloc(size_t nb_bytes) {
  void *p, *p0 = malloc(nb_bytes + EXTRA_BYTES_FOR_ALIGNMENT_AND_A_POINTER);
  if (!p0) return (void *) 0;
  p = (void *) (((uintptr_t) p0 + EXTRA_BYTES_FOR_ALIGNMENT_AND_A_POINTER) & (~((uintptr_t) (MALLOC_V4SF_ALIGNMENT-1))));
  *((void **) p - 1) = p0;
  return p;
}

static void Valigned_free(void *p) {
  if (p) free(*((void **) p - 1));
}


static int next_power_of_two(int N) {
  /* https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2 */
  /* compute the next highest power of 2 of 32-bit v */
  unsigned v = N;
  v--;
  v |= v >> 1;
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v++;
  return v;
}

static int is_power_of_two(int N) {
  /* https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2 */
  int f = N && !(N & (N - 1));
  return f;
}



PFFFT_EXPORT void *pffft_aligned_malloc(size_t nb_bytes) { return Valigned_malloc(nb_bytes); }
PFFFT_EXPORT void pffft_aligned_free(void *p) { Valigned_free(p); }
PFFFT_EXPORT int pffft_next_power_of_two(int N) { return next_power_of_two(N); }
PFFFT_EXPORT int pffft_is_power_of_two(int N) { return is_power_of_two(N); }

PFFFT_EXPORT void *pffftd_aligned_malloc(size_t nb_bytes) { return Valigned_malloc(nb_bytes); }
PFFFT_EXPORT void pffftd_aligned_free(void *p) { Valigned_free(p); }
PFFFT_EXPORT int pffftd_next_power_of_two(int N) { return next_power_of_two(N); }
PFFFT_EXPORT int pffftd_is_power_of_two(int N) { return is_power_of_two(N); }

/* ------------------------------------------------------------------------
   The default tuning (pffft_default_tuning(), pffftd_default_tuning()): the
   row of pffft_tuning.h that matches the CPU this runs on.
   ------------------------------------------------------------------------ */

#include <stdio.h>
#include <string.h>

#include "pffft_tuning.h"

typedef struct cpu_info {
  int arch;          /* PFFFT_CPU_* */
  unsigned vendor;   /* x86: X86_*, 0 for another; Arm: the implementer, 0 if unknown */
  unsigned id;       /* x86: family << 8 | model; Arm: the part number; 0 if unknown */
  size_t cache[3];   /* bytes of the L1 data cache, the L2 and the L3 of the
                        core (of its cluster or CCD where shared); 0 if unknown */
  char vendor_text[13]; /* x86: cpuid's vendor string */
  int hybrid;        /* x86: cores of two kinds, whose caches differ (cache[] is 0) */
} cpu_info;

/* ---- x86: cpuid ----
   not ARM64EC, which defines _M_X64 too, and whose cpuid would describe the
   emulated x86 CPU */
#if (defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)) && !defined(_M_ARM64EC)
#  if defined(_MSC_VER)
#    include <intrin.h>
static void cpuid(unsigned leaf, unsigned sub, unsigned r[4]) {
  int v[4];
  __cpuidex(v, (int)leaf, (int)sub);
  r[0] = (unsigned)v[0]; r[1] = (unsigned)v[1]; r[2] = (unsigned)v[2]; r[3] = (unsigned)v[3];
}
#  else
#    include <cpuid.h>
static void cpuid(unsigned leaf, unsigned sub, unsigned r[4]) {
  __cpuid_count(leaf, sub, r[0], r[1], r[2], r[3]);
}
#  endif

/* the cache sizes from the deterministic cache parameters, leaf 4 or
   0x8000001D: 1 if any were found */
static int x86_caches(unsigned leaf, cpu_info *c) {
  unsigned i, r[4];
  int found = 0;
  for (i = 0; i < 16; i++) {
    unsigned type, level;
    cpuid(leaf, i, r);
    type = r[0] & 0x1f;
    level = (r[0] >> 5) & 7;
    if (!type) break;
    if (type == 2 || level < 1 || level > 3) continue; /* an instruction cache */
    c->cache[level - 1] = (size_t)((r[1] >> 22) + 1) * (((r[1] >> 12) & 0x3ff) + 1) *
                          ((r[1] & 0xfff) + 1) * ((size_t)r[2] + 1);
    found = 1;
  }
  return found;
}

static void cpu_detect(cpu_info *c) {
  unsigned r[4], max, ext;
  cpuid(0, 0, r);
  max = r[0];
  memcpy(c->vendor_text, &r[1], 4);
  memcpy(c->vendor_text + 4, &r[3], 4);
  memcpy(c->vendor_text + 8, &r[2], 4);
  c->vendor_text[12] = 0;
  c->arch = PFFFT_CPU_X86;
  c->vendor = !strcmp(c->vendor_text, "GenuineIntel") ? X86_INTEL
            : !strcmp(c->vendor_text, "AuthenticAMD") ? X86_AMD
            : !strcmp(c->vendor_text, "HygonGenuine") ? X86_HYGON : 0;
  if (max >= 1) {
    unsigned family, model;
    cpuid(1, 0, r);
    family = (r[0] >> 8) & 0xf;
    model = (r[0] >> 4) & 0xf;
    if (family == 0xf) family += (r[0] >> 20) & 0xff;
    if (family == 6 || family >= 0xf) model |= ((r[0] >> 16) & 0xf) << 4;
    c->id = family << 8 | model;
  }
  cpuid(0x80000000, 0, r);
  ext = r[0];
  if (c->vendor == X86_AMD || c->vendor == X86_HYGON) {
    /* 0x8000001D with the topology extensions; the L3 is that of a CCD */
    int topoext = 0;
    if (ext >= 0x80000001) {
      cpuid(0x80000001, 0, r);
      topoext = (r[2] >> 22) & 1;
    }
    if ((!topoext || ext < 0x8000001D || !x86_caches(0x8000001D, c)) && ext >= 0x80000006) {
      cpuid(0x80000005, 0, r);
      c->cache[0] = (size_t)(r[2] >> 24) << 10;
      cpuid(0x80000006, 0, r);
      c->cache[1] = (size_t)(r[2] >> 16) << 10;
      c->cache[2] = (size_t)(r[3] >> 18) << 19;
    }
  } else if (max >= 4) {
    x86_caches(4, c);
  }
  /* Intel's hybrid parts (P-cores and E-cores): leaf 4 describes the core
     this thread runs on, which may not be the one that runs the transforms,
     and the L1 data caches of the two kinds differ (48 and 32 KiB). The
     rows' fixed values serve them better than a size from either. */
  if (max >= 7) {
    cpuid(7, 0, r);
    if ((r[3] >> 15) & 1) {
      c->hybrid = 1;
      c->cache[0] = c->cache[1] = c->cache[2] = 0;
    }
  }
}

/* ---- Arm, Linux and Android: /proc/cpuinfo and sysfs ---- */
#elif (defined(__aarch64__) || defined(__arm__)) && defined(__linux__)

/* a number from /sys/devices/system/cpu/cpu<cpu>/<leaf>, with a suffix K
   or M multiplied in; -1 if there is none */
static long sysfs_number(int cpu, const char *leaf, char *text, size_t n) {
  char path[96];
  FILE *f;
  long v;
  char *end;
  text[0] = 0;
  snprintf(path, sizeof path, "/sys/devices/system/cpu/cpu%d/%s", cpu, leaf);
  f = fopen(path, "r");
  if (!f) return -1;
  if (!fgets(text, (int)n, f)) text[0] = 0;
  fclose(f);
  v = strtol(text, &end, 0);
  if (end == text) return -1;
  if (*end == 'K') v <<= 10;
  if (*end == 'M') v <<= 20;
  return v;
}

static void cpu_detect(cpu_info *c) {
  char line[256];
  FILE *f = fopen("/proc/cpuinfo", "r");
  int proc = -1, best = -1, i;
  long best_key = -2;
  unsigned imp = 0;
  c->arch = PFFFT_CPU_ARM;
  if (!f) return;
  /* the fastest core: of the largest cpu_capacity, else of the highest
     maximum frequency, else the first */
  while (fgets(line, sizeof line, f)) {
    const char *colon = strchr(line, ':');
    if (!colon) continue;
    if (!strncmp(line, "processor", 9)) {
      proc = atoi(colon + 1);
      imp = 0;
    } else if (!strncmp(line, "CPU implementer", 15)) {
      imp = (unsigned)strtoul(colon + 1, NULL, 0);
    } else if (!strncmp(line, "CPU part", 8) && proc >= 0) {
      char text[32];
      long key = sysfs_number(proc, "cpu_capacity", text, sizeof text);
      if (key < 0) key = sysfs_number(proc, "cpufreq/cpuinfo_max_freq", text, sizeof text);
      if (key > best_key) {
        best_key = key;
        best = proc;
        c->vendor = imp;
        c->id = (unsigned)strtoul(colon + 1, NULL, 0);
      }
    }
  }
  fclose(f);
  if (best < 0) return;
  for (i = 0; i < 10; i++) {
    char leaf[48], text[32];
    long level, bytes;
    snprintf(leaf, sizeof leaf, "cache/index%d/level", i);
    level = sysfs_number(best, leaf, text, sizeof text);
    if (level < 1 || level > 3) continue;
    snprintf(leaf, sizeof leaf, "cache/index%d/type", i);
    if (sysfs_number(best, leaf, text, sizeof text) == -1 && !strncmp(text, "Instruction", 11)) continue;
    snprintf(leaf, sizeof leaf, "cache/index%d/size", i);
    bytes = sysfs_number(best, leaf, text, sizeof text);
    if (bytes > 0) c->cache[level - 1] = (size_t)bytes;
  }
}

/* ---- Apple's arm64 chips: sysctl ---- */
#elif defined(__APPLE__) && (defined(__aarch64__) || defined(__arm64__))
#  include <sys/sysctl.h>

/* the cache of the performance cores, if there are two kinds */
static size_t sysctl_cache(const char *leaf) {
  static const char *const where[2] = {"hw.perflevel0.", "hw."};
  int i;
  for (i = 0; i < 2; i++) {
    char name[64];
    int64_t v = 0; /* 4 or 8 bytes wide; little-endian */
    size_t len = sizeof v;
    snprintf(name, sizeof name, "%s%s", where[i], leaf);
    if (sysctlbyname(name, &v, &len, NULL, 0) == 0 && v > 0) return (size_t)v;
  }
  return 0;
}

static void cpu_detect(cpu_info *c) {
  c->arch = PFFFT_CPU_ARM;
  c->vendor = 0x61; /* Apple */
  c->cache[0] = sysctl_cache("l1dcachesize");
  c->cache[1] = sysctl_cache("l2cachesize");
  c->cache[2] = sysctl_cache("l3cachesize");
}

/* ---- Arm, Windows: the MIDR of each core in the registry ---- */
#elif defined(_WIN32) && (defined(_M_ARM64) || defined(_M_ARM64EC) || defined(_M_ARM) || defined(__aarch64__) || \
                          defined(__arm__))
#  include <windows.h>
#  if defined(_MSC_VER)
#    pragma comment(lib, "advapi32.lib")
#  endif

static void cpu_detect(cpu_info *c) {
  DWORD best_mhz = 0;
  int i;
  c->arch = PFFFT_CPU_ARM;
  /* the core of the highest frequency, the first of them */
  for (i = 0; i < 1024; i++) {
    char key[96];
    unsigned long long midr = 0;
    DWORD mhz = 0, len = sizeof midr, mlen = sizeof mhz;
    snprintf(key, sizeof key, "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\%d", i);
    if (RegGetValueA(HKEY_LOCAL_MACHINE, key, "CP 4000", RRF_RT_REG_QWORD, NULL, &midr, &len) != ERROR_SUCCESS)
      break;
    RegGetValueA(HKEY_LOCAL_MACHINE, key, "~MHz", RRF_RT_REG_DWORD, NULL, &mhz, &mlen);
    if (i == 0 || mhz > best_mhz) {
      best_mhz = mhz;
      c->vendor = (unsigned)(midr >> 24) & 0xff;
      c->id = (unsigned)(midr >> 4) & 0xfff;
    }
  }
}

#else
static void cpu_detect(cpu_info *c) {
#  if defined(__aarch64__) || defined(__arm__) || defined(_M_ARM64) || defined(_M_ARM64EC) || defined(_M_ARM)
  c->arch = PFFFT_CPU_ARM;
#  else
  c->arch = PFFFT_CPU_OTHER;
#  endif
}
#endif

/* The CPU, found once per process where atomics allow it: a thread that
   finds it first publishes it, the others use it or find it themselves. */
#if defined(__GNUC__) || defined(__clang__)
static cpu_info g_cpu;
static int g_cpu_state; /* 0: not found, 1: being published, 2: published */
static void cpu_get(cpu_info *c) {
  int expected = 0;
  if (__atomic_load_n(&g_cpu_state, __ATOMIC_ACQUIRE) == 2) {
    *c = g_cpu;
    return;
  }
  memset(c, 0, sizeof *c);
  cpu_detect(c);
  if (__atomic_compare_exchange_n(&g_cpu_state, &expected, 1, 0, __ATOMIC_ACQUIRE, __ATOMIC_RELAXED)) {
    g_cpu = *c;
    __atomic_store_n(&g_cpu_state, 2, __ATOMIC_RELEASE);
  }
}
#elif defined(_MSC_VER)
#  include <intrin.h>
static cpu_info g_cpu;
static volatile long g_cpu_state;
static void cpu_get(cpu_info *c) {
  if (_InterlockedCompareExchange(&g_cpu_state, 2, 2) == 2) {
    *c = g_cpu;
    return;
  }
  memset(c, 0, sizeof *c);
  cpu_detect(c);
  if (_InterlockedCompareExchange(&g_cpu_state, 1, 0) == 0) {
    g_cpu = *c;
    _InterlockedExchange(&g_cpu_state, 2);
  }
}
#else
static void cpu_get(cpu_info *c) {
  memset(c, 0, sizeof *c);
  cpu_detect(c);
}
#endif

static const tuning_row *cpu_row(const cpu_info *c) {
  size_t i, n = sizeof tuning_rows / sizeof tuning_rows[0];
  for (i = 0; i < n; i++) {
    const tuning_row *r = &tuning_rows[i];
    if (r->arch == c->arch && (!r->vendor || r->vendor == c->vendor) && c->id >= r->id_min &&
        c->id <= r->id_max)
      return r;
  }
  return &tuning_rows[n - 1];
}

static size_t rule_value(const tuning_rule *r, const cpu_info *c) {
  if (r->level >= 1 && r->level <= 3 && c->cache[r->level - 1])
    return c->cache[r->level - 1] / r->den * r->num;
  return r->bytes;
}

/* for pffft_priv_impl.h: the default thresholds on this CPU (not exported) */
void pffft_priv_cpu_tuning(pffft_tuning_t *t) {
  cpu_info c;
  const tuning_row *r;
  cpu_get(&c);
  r = cpu_row(&c);
  t->size = sizeof *t;
  t->split_max_real = rule_value(&r->split_max_real, &c);
  t->split_max_complex = rule_value(&r->split_max_complex, &c);
  t->prefetch_min = rule_value(&r->prefetch_min, &c);
  t->radix6_max = rule_value(&r->radix6_max, &c);
}

static void cache_text(char *buf, size_t n, const char *name, size_t bytes) {
  if (!bytes) snprintf(buf, n, "%s unknown", name);
  else if (bytes % (1u << 20) == 0) snprintf(buf, n, "%s %u MiB", name, (unsigned)(bytes >> 20));
  else snprintf(buf, n, "%s %u KiB", name, (unsigned)(bytes >> 10));
}

static int cpu_describe(char *buf, int len) {
  cpu_info c;
  char id[64], l1[24], l2[24], l3[24];
  cpu_get(&c);
  if (c.arch == PFFFT_CPU_X86)
    snprintf(id, sizeof id, "x86 %s family 0x%x model 0x%x", c.vendor_text, c.id >> 8, c.id & 0xff);
  else if (c.arch == PFFFT_CPU_ARM)
    snprintf(id, sizeof id, "Arm implementer 0x%02x part 0x%03x", c.vendor, c.id);
  else
    snprintf(id, sizeof id, "unknown CPU");
  cache_text(l1, sizeof l1, "L1d", c.cache[0]);
  cache_text(l2, sizeof l2, "L2", c.cache[1]);
  cache_text(l3, sizeof l3, "L3", c.cache[2]);
  return snprintf(buf, len > 0 ? (size_t)len : 0, "%s%s, %s, %s, %s: %s", id, c.hybrid ? " (hybrid)" : "", l1, l2,
                  l3, cpu_row(&c)->name);
}

PFFFT_EXPORT int pffft_cpu_describe(char *buf, int len) { return cpu_describe(buf, len); }
PFFFT_EXPORT int pffftd_cpu_describe(char *buf, int len) { return cpu_describe(buf, len); }
