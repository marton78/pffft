/* benchmarks/platform.c -- the operating-system side of bench_fft: clocks, CPU
 * cycle counters, core selection, shared libraries, aligned memory and the
 * description of the machine. bench.h documents the interface.
 *
 *   macOS, iOS      thread_selfcounts(), QoS class instead of pinning
 *   Linux, Android  perf_event_open() on a pinned thread
 *   Windows         QueryPerformanceCounter(), no cycle counter
 *   anything else   clock_gettime() where POSIX has it, else clock(); no cycle
 *                   counter, no shared libraries
 */

/* Feature-test macros have to precede the first include. */
#if defined(__linux__) && !defined(_GNU_SOURCE)
#  define _GNU_SOURCE /* sched_setaffinity, sched_getcpu, syscall */
#endif
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#  define _DARWIN_C_SOURCE /* clock_gettime_nsec_np, sysctlbyname, RTLD_DEFAULT */
#endif
#if defined(__unix__) && !defined(__linux__) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 200112L /* clock_gettime on the BSDs and Cygwin */
#endif
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS /* getenv */
#endif

#include "bench.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  include <malloc.h>
#elif defined(__APPLE__)
#  include <TargetConditionals.h>
#  include <dlfcn.h>
#  include <pthread.h>
#  include <pthread/qos.h>
#  include <sys/sysctl.h>
#elif defined(__linux__)
#  include <dlfcn.h>
#  include <errno.h>
#  if defined(__has_include)
#    if !__has_include(<linux/perf_event.h>)
#      define BENCH_NO_PERF_EVENT /* musl systems install kernel headers separately */
#    endif
#  endif
#  if !defined(BENCH_NO_PERF_EVENT)
#    include <linux/perf_event.h>
#  endif
#  include <sched.h>
#  include <sys/ioctl.h>
#  include <sys/syscall.h>
#  include <sys/utsname.h>
#  include <unistd.h>
#  if defined(__ANDROID__)
#    include <sys/system_properties.h>
#  endif
#endif

#if defined(_M_ARM64EC) /* also defines _M_X64 */
#  define BENCH_ARCH "arm64ec"
#elif defined(__aarch64__) || defined(_M_ARM64)
#  define BENCH_ARCH "arm64"
#elif defined(__x86_64__) || defined(_M_X64)
#  define BENCH_ARCH "x86_64"
#elif defined(__i386__) || defined(_M_IX86)
#  define BENCH_ARCH "x86"
#elif defined(__arm__) || defined(_M_ARM)
#  define BENCH_ARCH "arm"
#elif defined(__riscv)
#  define BENCH_ARCH "riscv"
#elif defined(__powerpc64__)
#  define BENCH_ARCH "ppc64"
#elif defined(__wasm__)
#  define BENCH_ARCH "wasm"
#else
#  define BENCH_ARCH "unknown"
#endif

typedef struct machine {
  char os[96];
  char cpu[96];
  long cache[3]; /* bytes of L1 data, L2, L3; 0 if unknown */
} machine;

static const char *g_counter = "none";
static int g_pinned_cpu = -1;

/* Each system below provides now_ns(), read_counters(), steer(),
 * counter_open(), counter_close() and describe(). */

/* ------------------------------------------------------------------------
 * macOS, iOS
 * ---------------------------------------------------------------------- */
#if defined(__APPLE__)

/* thread_selfcounts() is exported by libsystem_kernel without a header, and
 * needs no privileges. Kind 2 fills one {instructions, cycles} pair per kind
 * of core, performance cores first. A machine with a single kind of core
 * fills only the first pair, so the second one has to start out as zero. */
typedef int (*selfcounts_fn)(int kind, void *buf, size_t nbytes);
static selfcounts_fn g_selfcounts;
static size_t g_selfcounts_bytes;

static uint64_t now_ns(void) { return clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW); }

static void read_counters(bench_snap *s) {
  uint64_t b[4] = {0, 0, 0, 0};
  if (g_selfcounts) g_selfcounts(2, b, g_selfcounts_bytes);
  s->cycles = b[1];
  s->stray_cycles = b[3];
  s->instructions = b[0] + b[2];
}

static void steer(int pin_cpu) {
  /* Threads cannot be pinned here. This class keeps the thread on the
   * performance cores most of the time; the batches that stray are told
   * apart by their efficiency-core cycles. */
  (void)pin_cpu;
  pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
}

static int counter_open(void) {
  uint64_t b[4];
  /* Looked up at run time: a system without the call has no counter, instead
   * of a program that does not link. The detour over an integer is because
   * ISO C has no conversion from an object pointer to a function pointer. */
  void *sym = dlsym(RTLD_DEFAULT, "thread_selfcounts");
  if (!sym) return 0;
  g_selfcounts = (selfcounts_fn)(uintptr_t)sym;
  for (g_selfcounts_bytes = sizeof b; g_selfcounts_bytes >= 16; g_selfcounts_bytes -= 16) {
    if (g_selfcounts(2, b, g_selfcounts_bytes) == 0) return 1;
  }
  g_selfcounts = NULL;
  return 0;
}

static void counter_close(void) { g_selfcounts = NULL; }

static const char *counter_name(void) { return "thread_selfcounts"; }

static void sysctl_text(const char *name, char *out, size_t n) {
  size_t len = n - 1;
  if (sysctlbyname(name, out, &len, NULL, 0) != 0) len = 0;
  out[len] = 0;
}

/* The cache of the performance cores if the machine has two kinds of core. */
static long sysctl_cache(const char *leaf) {
  static const char *const where[2] = {"hw.perflevel0.", "hw."};
  for (int i = 0; i < 2; i++) {
    char name[64];
    int64_t v = 0; /* the values are 4 or 8 bytes wide; all Apple CPUs are little-endian */
    size_t len = sizeof v;
    snprintf(name, sizeof name, "%s%s", where[i], leaf);
    if (sysctlbyname(name, &v, &len, NULL, 0) == 0 && v > 0) return (long)v;
  }
  return 0;
}

static void describe(machine *m) {
  char version[32];
  sysctl_text("kern.osproductversion", version, sizeof version);
  snprintf(m->os, sizeof m->os, "%s%s%s", TARGET_OS_IPHONE ? "iOS" : "macOS", version[0] ? " " : "",
           version);
  sysctl_text("machdep.cpu.brand_string", m->cpu, sizeof m->cpu);
  if (!m->cpu[0]) sysctl_text("hw.machine", m->cpu, sizeof m->cpu); /* iOS: the device model */
  m->cache[0] = sysctl_cache("l1dcachesize");
  m->cache[1] = sysctl_cache("l2cachesize");
  m->cache[2] = sysctl_cache("l3cachesize");
}

/* ------------------------------------------------------------------------
 * Linux, Android
 * ---------------------------------------------------------------------- */
#elif defined(__linux__)

static int g_perf_cycles = -1, g_perf_instructions = -1;

static uint64_t now_ns(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC_RAW, &t);
  return (uint64_t)t.tv_sec * 1000000000u + (uint64_t)t.tv_nsec;
}

static void read_counters(bench_snap *s) {
  uint64_t b[3] = {0, 0, 0}; /* number of values, cycles, instructions */
  if (g_perf_cycles >= 0 && read(g_perf_cycles, b, sizeof b) < 0) b[1] = b[2] = 0;
  s->cycles = b[1];
  s->stray_cycles = 0; /* a pinned thread never runs on another kind of core */
  s->instructions = b[2];
}

/* First line of /sys/devices/system/cpu/cpu<cpu>/<leaf>; 0 if there is none. */
static int sysfs_text(int cpu, const char *leaf, char *out, size_t n) {
  char path[128];
  FILE *f;
  int ok;
  snprintf(path, sizeof path, "/sys/devices/system/cpu/cpu%d/%s", cpu, leaf);
  f = fopen(path, "r");
  if (!f) return 0;
  ok = fgets(out, (int)n, f) != NULL;
  fclose(f);
  if (ok) out[strcspn(out, "\n")] = 0;
  return ok;
}

static long sysfs_long(int cpu, const char *leaf) {
  char text[32];
  return sysfs_text(cpu, leaf, text, sizeof text) ? strtol(text, NULL, 10) : -1;
}

/* Pinning comes before the counter is opened, because the counter is bound
 * to the pinned core: see perf_open(). */
static void steer(int pin_cpu) {
  static const char *const rank[2] = {"cpu_capacity", "cpufreq/cpuinfo_max_freq"};
  cpu_set_t allowed, one;
  int here = sched_getcpu(), cpu = pin_cpu;
  if (pin_cpu == BENCH_PIN_NONE) return;
  if (pin_cpu == BENCH_PIN_AUTO) {
    long best = -1;
    cpu = here;
    if (sched_getaffinity(0, sizeof allowed, &allowed) == 0) {
      for (int k = 0; k < 2 && best < 0; k++) {
        for (int c = 0; c < CPU_SETSIZE; c++) {
          long v = CPU_ISSET(c, &allowed) ? sysfs_long(c, rank[k]) : -1;
          /* among equals, stay where the scheduler put the thread */
          if (v > best || (v >= 0 && v == best && c == here)) {
            best = v;
            cpu = c;
          }
        }
      }
    }
  }
  if (cpu < 0 || cpu >= CPU_SETSIZE) {
    fprintf(stderr, "bench_fft: cannot pin the thread: no cpu %d\n", cpu);
    return;
  }
  CPU_ZERO(&one);
  CPU_SET(cpu, &one);
  if (sched_setaffinity(0, sizeof one, &one) == 0) g_pinned_cpu = cpu;
  else fprintf(stderr, "bench_fft: cannot pin the thread to cpu %d: %s\n", cpu, strerror(errno));
}

#  if defined(BENCH_NO_PERF_EVENT)
static int counter_open(void) { return 0; }
#  else
/* Neither glibc nor bionic wraps perf_event_open().
 *
 * The event is bound to the pinned core. On big.LITTLE every kind of core has
 * its own PMU, and an event that names no core goes to the first PMU, where
 * it stands still while the thread runs on the other kind (seen on a
 * Cortex-A78/A55 phone). */
static int perf_open(unsigned config, int group) {
  struct perf_event_attr attr;
  memset(&attr, 0, sizeof attr);
  attr.type = PERF_TYPE_HARDWARE;
  attr.size = sizeof attr;
  attr.config = config;
  attr.disabled = group < 0; /* the leader starts the group once it is complete */
  attr.exclude_kernel = 1;   /* user-space counting is what unprivileged users may do */
  attr.exclude_hv = 1;
  attr.read_format = PERF_FORMAT_GROUP;
  return (int)syscall(SYS_perf_event_open, &attr, 0 /* this thread */, g_pinned_cpu, group, 0UL);
}

static int counter_open(void) {
  /* Unpinned, a batch that moves to the other kind of core halfway is counted
   * in part, and nothing tells that it was: no counter then. */
  if (g_pinned_cpu < 0) return 0;
  g_perf_cycles = perf_open(PERF_COUNT_HW_CPU_CYCLES, -1);
  if (g_perf_cycles < 0) return 0;
  /* If this one fails the group still counts cycles, and instructions read as 0. */
  g_perf_instructions = perf_open(PERF_COUNT_HW_INSTRUCTIONS, g_perf_cycles);
  if (ioctl(g_perf_cycles, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) == 0) return 1;
  if (g_perf_instructions >= 0) close(g_perf_instructions);
  close(g_perf_cycles);
  g_perf_cycles = g_perf_instructions = -1;
  return 0;
}
#  endif

static void counter_close(void) {
  if (g_perf_instructions >= 0) close(g_perf_instructions);
  if (g_perf_cycles >= 0) close(g_perf_cycles);
  g_perf_cycles = g_perf_instructions = -1;
}

static const char *counter_name(void) { return "perf_event"; }

#  if !defined(__ANDROID__)
/* The value of the first "<key> : value" line of /proc/cpuinfo. */
static int cpuinfo_text(const char *key, char *out, size_t n) {
  char line[256];
  size_t klen = strlen(key);
  int found = 0;
  FILE *f = fopen("/proc/cpuinfo", "r");
  if (!f) return 0;
  while (!found && fgets(line, sizeof line, f)) {
    const char *colon = strchr(line, ':');
    /* "model name" must not match "model" and the other way round */
    if (!colon || strncmp(line, key, klen) != 0 || strspn(line + klen, " \t") != (size_t)(colon - line) - klen)
      continue;
    snprintf(out, n, "%s", colon + 1 + strspn(colon + 1, " \t"));
    out[strcspn(out, "\n")] = 0;
    found = 1;
  }
  fclose(f);
  return found;
}

/* arm64 /proc/cpuinfo names no model, only the implementer and part numbers
 * of the core; servers (GitHub's arm64 runners, AWS Graviton) are known by
 * them. Writes "" if there are none. */
static void arm_core_name(char *out, size_t n) {
  static const struct {
    unsigned implementer, part;
    const char *name;
  } known[] = {{0x41, 0xd0c, "Arm Neoverse N1"}, {0x41, 0xd40, "Arm Neoverse V1"}, {0x41, 0xd49, "Arm Neoverse N2"},
               {0x41, 0xd4f, "Arm Neoverse V2"}, {0x41, 0xd8e, "Arm Neoverse N3"}, {0x41, 0xd84, "Arm Neoverse V3"},
               {0x41, 0xd08, "Arm Cortex-A72"},  {0x41, 0xd0b, "Arm Cortex-A76"},  {0xc0, 0xac3, "AmpereOne"},
               {0x61, 0, "Apple"}};
  char implementer[32], part[32];
  unsigned imp, prt;
  out[0] = 0;
  if (!cpuinfo_text("CPU implementer", implementer, sizeof implementer) || !cpuinfo_text("CPU part", part, sizeof part))
    return;
  imp = (unsigned)strtoul(implementer, NULL, 0);
  prt = (unsigned)strtoul(part, NULL, 0);
  for (size_t i = 0; i < sizeof known / sizeof known[0]; i++) {
    if (known[i].implementer == imp && (known[i].part == prt || known[i].part == 0)) {
      snprintf(out, n, "%s", known[i].name);
      return;
    }
  }
  snprintf(out, n, "arm64 core, implementer 0x%x, part 0x%x", imp, prt);
}
#  endif

static void describe(machine *m) {
  int cpu = g_pinned_cpu >= 0 ? g_pinned_cpu : sched_getcpu();
#  if defined(__ANDROID__)
  char prop[PROP_VALUE_MAX];
  snprintf(m->os, sizeof m->os, "Android");
  if (__system_property_get("ro.build.version.release", prop) > 0)
    snprintf(m->os, sizeof m->os, "Android %s", prop);
  /* /proc/cpuinfo names no model on arm64; the system-on-chip is what identifies a phone */
  if (__system_property_get("ro.soc.model", prop) > 0 || __system_property_get("ro.hardware", prop) > 0)
    snprintf(m->cpu, sizeof m->cpu, "%s", prop);
#  else
  struct utsname u;
  snprintf(m->os, sizeof m->os, "Linux");
  if (uname(&u) == 0) snprintf(m->os, sizeof m->os, "%.40s %.40s", u.sysname, u.release);
  if (!cpuinfo_text("model name", m->cpu, sizeof m->cpu) && !cpuinfo_text("Hardware", m->cpu, sizeof m->cpu) &&
      !cpuinfo_text("Model", m->cpu, sizeof m->cpu))
    arm_core_name(m->cpu, sizeof m->cpu);
#  endif
  for (int i = 0; i < 8; i++) {
    char leaf[48], type[32], size[32], *unit;
    long level, bytes;
    snprintf(leaf, sizeof leaf, "cache/index%d/level", i);
    level = sysfs_long(cpu, leaf);
    snprintf(leaf, sizeof leaf, "cache/index%d/type", i);
    if (level < 1 || level > 3 || !sysfs_text(cpu, leaf, type, sizeof type) || !strcmp(type, "Instruction"))
      continue;
    snprintf(leaf, sizeof leaf, "cache/index%d/size", i);
    if (!sysfs_text(cpu, leaf, size, sizeof size)) continue;
    bytes = strtol(size, &unit, 10); /* "32K", "1024K", "8M" */
    if (*unit == 'K') bytes *= 1024;
    if (*unit == 'M') bytes *= 1024 * 1024;
    m->cache[level - 1] = bytes;
  }
}

/* ------------------------------------------------------------------------
 * Windows, and anything else: wall clock only
 * ---------------------------------------------------------------------- */
#else

#  if defined(_WIN32)

static uint64_t now_ns(void) {
  static LARGE_INTEGER freq;
  LARGE_INTEGER t;
  if (!freq.QuadPart) QueryPerformanceFrequency(&freq);
  QueryPerformanceCounter(&t);
  /* in two parts, because ticks * 1e9 overflows after half an hour of uptime */
  return (uint64_t)(t.QuadPart / freq.QuadPart) * 1000000000u +
         (uint64_t)(t.QuadPart % freq.QuadPart) * 1000000000u / (uint64_t)freq.QuadPart;
}

static void describe(machine *m) {
  HKEY key;
  DWORD len = (DWORD)sizeof m->cpu - 1; /* m is zero-filled: the text stays terminated */
  char build[16] = {0};
  DWORD blen = (DWORD)sizeof build - 1;
  snprintf(m->os, sizeof m->os, "Windows");
  /* the build number names the release (26100: Windows 11 24H2 and Server
   * 2025); ProductName says "Windows 10" also on Windows 11 */
  if (RegOpenKeyExA(HKEY_LOCAL_MACHINE, "SOFTWARE\\Microsoft\\Windows NT\\CurrentVersion", 0, KEY_READ, &key) ==
      ERROR_SUCCESS) {
    if (RegQueryValueExA(key, "CurrentBuildNumber", NULL, NULL, (LPBYTE)build, &blen) == ERROR_SUCCESS && build[0])
      snprintf(m->os, sizeof m->os, "Windows build %s", build);
    RegCloseKey(key);
  }
  if (RegOpenKeyExA(HKEY_LOCAL_MACHINE, "HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0", 0, KEY_READ,
                    &key) != ERROR_SUCCESS)
    return;
  if (RegQueryValueExA(key, "ProcessorNameString", NULL, NULL, (LPBYTE)m->cpu, &len) != ERROR_SUCCESS)
    m->cpu[0] = 0;
  RegCloseKey(key);
}

#  else

static uint64_t now_ns(void) {
#    if defined(CLOCK_MONOTONIC)
  struct timespec t;
  if (clock_gettime(CLOCK_MONOTONIC, &t) == 0) return (uint64_t)t.tv_sec * 1000000000u + (uint64_t)t.tv_nsec;
#    endif
  /* processor time, often in coarse ticks: only long batches measure well */
  return (uint64_t)((double)clock() * (1e9 / CLOCKS_PER_SEC));
}

static void describe(machine *m) { (void)m; }

#  endif

static void read_counters(bench_snap *s) { s->cycles = s->stray_cycles = s->instructions = 0; }
static void steer(int pin_cpu) { (void)pin_cpu; }
static int counter_open(void) { return 0; }
static void counter_close(void) {}
static const char *counter_name(void) { return "none"; }

#endif

/* ------------------------------------------------------------------------
 * Common to all systems
 * ---------------------------------------------------------------------- */

int bench_platform_init(int pin_cpu) {
  /* Undocumented switch for testing the wall-clock fallback on any machine. */
  const char *off = getenv("BENCH_NO_COUNTER");
  volatile unsigned spin = 0;
  bench_snap start, a, b;
  steer(pin_cpu);
  if ((off && *off && strcmp(off, "0") != 0) || !counter_open()) return 0;
  /* Work until the counter shows the thread on the preferred cores, for at
   * most 0.2 s. macOS starts a thread on an efficiency core and moves it only
   * after it has been busy for some milliseconds; calibrating before that
   * would be calibrating on the wrong core. A counter that opens but stands
   * still, as in many virtual machines, is no counter. */
  bench_snap_begin(&start);
  b = start;
  do {
    a = b;
    for (unsigned i = 0; i < 100000u; i++) spin += i;
    bench_snap_end(&b);
    if (b.cycles + b.stray_cycles == a.cycles + a.stray_cycles) {
      counter_close();
      return 0;
    }
  } while (b.stray_cycles != a.stray_cycles && b.ns - start.ns < 200000000u);
  g_counter = counter_name();
  return 1;
}

/* The clock is read on the outside and the counters on the inside, so that
 * the cycle window holds as little as possible besides the timed code. */
void bench_snap_begin(bench_snap *s) {
  s->ns = now_ns();
  read_counters(s);
}

void bench_snap_end(bench_snap *s) {
  read_counters(s);
  s->ns = now_ns();
}

static const char *trimmed(char *s) {
  size_t n = strlen(s);
  while (n > 0 && (s[n - 1] == ' ' || s[n - 1] == '\t' || s[n - 1] == '\r' || s[n - 1] == '\n')) s[--n] = 0;
  s += strspn(s, " \t");
  return *s ? s : "unknown";
}

void bench_platform_describe(FILE *f, const char *prefix) {
  static const char *const cache_key[3] = {"l1d_bytes", "l2_bytes", "l3_bytes"};
  machine m;
  memset(&m, 0, sizeof m);
  describe(&m);
  fprintf(f, "%sos: %s\n", prefix, trimmed(m.os));
  fprintf(f, "%sarch: %s\n", prefix, BENCH_ARCH);
  fprintf(f, "%scpu: %s\n", prefix, trimmed(m.cpu));
  fprintf(f, "%scounter: %s\n", prefix, g_counter);
  if (g_pinned_cpu >= 0) fprintf(f, "%spinned_cpu: %d\n", prefix, g_pinned_cpu);
  else fprintf(f, "%spinned_cpu: none\n", prefix);
  for (int i = 0; i < 3; i++) {
    if (m.cache[i] > 0) fprintf(f, "%s%s: %ld\n", prefix, cache_key[i], m.cache[i]);
  }
}

/* ------------------------------------------------------------------------
 * Aligned memory and shared libraries
 * ---------------------------------------------------------------------- */
#if defined(_WIN32)

void *bench_aligned_alloc(size_t bytes) { return _aligned_malloc(bytes, 64); }
void bench_aligned_free(void *p) { _aligned_free(p); }

void *bench_dl_open(const char *path, char *err, size_t errlen) {
  HMODULE lib = LoadLibraryA(path);
  if (!lib && err && errlen) snprintf(err, errlen, "LoadLibrary failed with error %lu", (unsigned long)GetLastError());
  return (void *)lib;
}

void *bench_dl_sym(void *lib, const char *symbol) {
  return (void *)(uintptr_t)GetProcAddress((HMODULE)lib, symbol);
}

#elif defined(__APPLE__) || defined(__linux__)

void *bench_aligned_alloc(size_t bytes) {
  void *p = NULL;
  return posix_memalign(&p, 64, bytes) == 0 ? p : NULL;
}
void bench_aligned_free(void *p) { free(p); }

void *bench_dl_open(const char *path, char *err, size_t errlen) {
  /* RTLD_LOCAL: several versions of one library are loaded side by side */
  void *lib = dlopen(path, RTLD_NOW | RTLD_LOCAL);
  if (!lib && err && errlen) {
    const char *why = dlerror();
    snprintf(err, errlen, "%s", why ? why : "dlopen failed");
  }
  return lib;
}

void *bench_dl_sym(void *lib, const char *symbol) { return dlsym(lib, symbol); }

#else

/* The original pointer is kept just below the aligned block. */
void *bench_aligned_alloc(size_t bytes) {
  void *raw = malloc(bytes + 64 + sizeof(void *));
  void **p;
  if (!raw) return NULL;
  p = (void **)(((uintptr_t)raw + sizeof(void *) + 63) & ~(uintptr_t)63);
  p[-1] = raw;
  return p;
}
void bench_aligned_free(void *p) {
  if (p) free(((void **)p)[-1]);
}

void *bench_dl_open(const char *path, char *err, size_t errlen) {
  (void)path;
  if (err && errlen) snprintf(err, errlen, "shared libraries are not supported on this system");
  return NULL;
}

void *bench_dl_sym(void *lib, const char *symbol) {
  (void)lib;
  (void)symbol;
  return NULL;
}

#endif
