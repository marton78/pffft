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
#  include <dirent.h>
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

/* --events (bench_platform_events): the names as given, and how many */
static char g_event_names[BENCH_MAX_EVENTS][96];
static int g_nevents;

/* Fills in how system-specific code opens event i, from its text; 0, or -1
 * with the reason in err. Defined per system below. */
static int event_resolve(int i, const char *text, char *err, size_t errlen);
static int event_counted(int i); /* was event i opened (bench_platform_init)? */

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

static int event_resolve(int i, const char *text, char *err, size_t errlen) {
  (void)i;
  snprintf(err, errlen, "%s: hardware events are counted only on Linux (perf_event)", text);
  return -1;
}
static int event_counted(int i) { (void)i; return 0; }

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
/* the --events: how to open each, its file and its place among the group's values */
static struct {
  uint32_t type;
  uint64_t config, config1, config2;
  int fd, slot;
} g_ev[BENCH_MAX_EVENTS];
static int g_instructions_slot = -1;

static int event_counted(int i) { return g_ev[i].fd >= 0; }

static uint64_t now_ns(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC_RAW, &t);
  return (uint64_t)t.tv_sec * 1000000000u + (uint64_t)t.tv_nsec;
}

static void read_counters(bench_snap *s) {
  /* the number of values, then cycles, and the others in the order opened */
  uint64_t b[3 + BENCH_MAX_EVENTS] = {0};
  if (g_perf_cycles >= 0 && read(g_perf_cycles, b, sizeof b) < 0) memset(b, 0, sizeof b);
  s->cycles = b[1];
  s->stray_cycles = 0; /* a pinned thread never runs on another kind of core */
  s->instructions = g_instructions_slot >= 0 ? b[1 + g_instructions_slot] : 0;
  for (int i = 0; i < g_nevents; i++) s->events[i] = g_ev[i].slot >= 0 ? b[1 + g_ev[i].slot] : 0;
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
static void counter_close(void) {}
static int event_resolve(int i, const char *text, char *err, size_t errlen) {
  (void)i;
  snprintf(err, errlen, "%s: this build has no perf_event (no linux/perf_event.h)", text);
  return -1;
}
#  else
/* Neither glibc nor bionic wraps perf_event_open().
 *
 * The event is bound to the pinned core. On big.LITTLE every kind of core has
 * its own PMU, and an event that names no core goes to the first PMU, where
 * it stands still while the thread runs on the other kind (seen on a
 * Cortex-A78/A55 phone). */
static int perf_open(uint32_t type, uint64_t config, uint64_t config1, uint64_t config2, int group) {
  struct perf_event_attr attr;
  memset(&attr, 0, sizeof attr);
  attr.type = type;
  attr.size = sizeof attr;
  attr.config = config;
  attr.config1 = config1;
  attr.config2 = config2;
  attr.disabled = group < 0; /* the leader starts the group once it is complete */
  attr.exclude_kernel = 1;   /* user-space counting is what unprivileged users may do */
  attr.exclude_hv = 1;
  attr.read_format = PERF_FORMAT_GROUP;
  return (int)syscall(SYS_perf_event_open, &attr, 0 /* this thread */, g_pinned_cpu, group, 0UL);
}

static void counter_close(void) {
  for (int i = 0; i < g_nevents; i++) {
    if (g_ev[i].fd >= 0) close(g_ev[i].fd);
    g_ev[i].fd = -1;
  }
  if (g_perf_instructions >= 0) close(g_perf_instructions);
  if (g_perf_cycles >= 0) close(g_perf_cycles);
  g_perf_cycles = g_perf_instructions = g_instructions_slot = -1;
}

static int counter_open(void) {
  int slot = 1;
  /* Unpinned, a batch that moves to the other kind of core halfway is counted
   * in part, and nothing tells that it was: no counter then. */
  if (g_pinned_cpu < 0) return 0;
  g_perf_cycles = perf_open(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES, 0, 0, -1);
  if (g_perf_cycles < 0) return 0;
  /* If this one fails the group still counts cycles, and instructions read as 0. */
  g_perf_instructions = perf_open(PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS, 0, 0, g_perf_cycles);
  if (g_perf_instructions >= 0) g_instructions_slot = slot++;
  for (int i = 0; i < g_nevents; i++) {
    g_ev[i].fd = perf_open(g_ev[i].type, g_ev[i].config, g_ev[i].config1, g_ev[i].config2, g_perf_cycles);
    if (g_ev[i].fd < 0) {
      int why = errno;
      fprintf(stderr, "bench_fft: cannot count the event %s, left out: %s%s\n", g_event_names[i], strerror(why),
              why == EINVAL ? " (unknown to this PMU, or more events than it has counters: cycles and"
                              " instructions take two)" : "");
      g_ev[i].slot = -1;
      continue;
    }
    g_ev[i].slot = slot++;
  }
  if (ioctl(g_perf_cycles, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP) == 0) return 1;
  counter_close();
  return 0;
}

/* perf's generic events, by the names perf list gives them */
#    define CACHE_EVENT(cache, result) \
      (PERF_COUNT_HW_CACHE_##cache | (PERF_COUNT_HW_CACHE_OP_READ << 8) | (PERF_COUNT_HW_CACHE_RESULT_##result << 16))
static const struct {
  const char *name;
  uint32_t type;
  uint64_t config;
} generic_events[] = {
    {"branches", PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_INSTRUCTIONS},
    {"branch-instructions", PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_INSTRUCTIONS},
    {"branch-misses", PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_MISSES},
    {"cache-references", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES},
    {"cache-misses", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES},
    {"stalled-cycles-frontend", PERF_TYPE_HARDWARE, PERF_COUNT_HW_STALLED_CYCLES_FRONTEND},
    {"stalled-cycles-backend", PERF_TYPE_HARDWARE, PERF_COUNT_HW_STALLED_CYCLES_BACKEND},
    {"ref-cycles", PERF_TYPE_HARDWARE, PERF_COUNT_HW_REF_CPU_CYCLES},
    {"L1-dcache-loads", PERF_TYPE_HW_CACHE, CACHE_EVENT(L1D, ACCESS)},
    {"L1-dcache-load-misses", PERF_TYPE_HW_CACHE, CACHE_EVENT(L1D, MISS)},
    {"L1-icache-loads", PERF_TYPE_HW_CACHE, CACHE_EVENT(L1I, ACCESS)},
    {"L1-icache-load-misses", PERF_TYPE_HW_CACHE, CACHE_EVENT(L1I, MISS)},
    {"iTLB-loads", PERF_TYPE_HW_CACHE, CACHE_EVENT(ITLB, ACCESS)},
    {"iTLB-load-misses", PERF_TYPE_HW_CACHE, CACHE_EVENT(ITLB, MISS)},
    {"branch-loads", PERF_TYPE_HW_CACHE, CACHE_EVENT(BPU, ACCESS)},
    {"branch-load-misses", PERF_TYPE_HW_CACHE, CACHE_EVENT(BPU, MISS)},
};

#    define PMU_DIR "/sys/bus/event_source/devices"

/* The first line of a file, without its line break; 0 if there is none. */
static int file_line(const char *path, char *out, size_t n) {
  FILE *f = fopen(path, "r");
  int ok;
  if (!f) return 0;
  ok = fgets(out, (int)n, f) != NULL;
  fclose(f);
  if (ok) out[strcspn(out, "\n")] = 0;
  return ok;
}

/* Puts value into the bits that PMU_DIR/<pmu>/format/<field> names, e.g.
 * "config:0-7", "config1:0-15" or "config:0-7,32-35" (the value's low bits go
 * into the first range). */
static int set_field(int i, const char *pmu, const char *field, uint64_t value, char *err, size_t errlen) {
  char path[256], fmt[128], *p;
  uint64_t *dst;
  snprintf(path, sizeof path, PMU_DIR "/%s/format/%s", pmu, field);
  if (!file_line(path, fmt, sizeof fmt)) {
    snprintf(err, errlen, "the PMU %s has no field '%s' (%s)", pmu, field, path);
    return -1;
  }
  if (!strncmp(fmt, "config2:", 8)) dst = &g_ev[i].config2, p = fmt + 8;
  else if (!strncmp(fmt, "config1:", 8)) dst = &g_ev[i].config1, p = fmt + 8;
  else if (!strncmp(fmt, "config:", 7)) dst = &g_ev[i].config, p = fmt + 7;
  else {
    snprintf(err, errlen, "%s: unknown format '%s'", path, fmt);
    return -1;
  }
  while (*p) {
    char *end;
    long lo = strtol(p, &end, 10), hi = lo;
    if (end == p) break;
    if (*end == '-') hi = strtol(end + 1, &end, 10);
    for (long b = lo; b <= hi && b < 64; b++, value >>= 1) {
      if (value & 1) *dst |= (uint64_t)1 << b;
    }
    p = end + (*end == ',');
  }
  return 0;
}

/* "field=value,field,..." of an event of pmu: a field without a value is 1, or
 * the name of one of the PMU's events, which is then read from its events/. */
static int set_terms(int i, const char *pmu, const char *terms, int depth, char *err, size_t errlen) {
  char copy[256], *t, *save = NULL;
  snprintf(copy, sizeof copy, "%s", terms);
  for (t = strtok_r(copy, ",", &save); t; t = strtok_r(NULL, ",", &save)) {
    char *eq = strchr(t, '='), path[256], text[256];
    t += strspn(t, " ");
    if (eq) {
      *eq = 0;
      if (set_field(i, pmu, t, strtoull(eq + 1, NULL, 0), err, errlen)) return -1;
      continue;
    }
    snprintf(path, sizeof path, PMU_DIR "/%s/events/%s", pmu, t);
    if (depth == 0 && file_line(path, text, sizeof text)) {
      if (set_terms(i, pmu, text, 1, err, errlen)) return -1;
    } else if (set_field(i, pmu, t, 1, err, errlen)) {
      return -1;
    }
  }
  return 0;
}

static int set_pmu(int i, const char *pmu, char *err, size_t errlen) {
  char path[256], text[32];
  snprintf(path, sizeof path, PMU_DIR "/%s/type", pmu);
  if (!file_line(path, text, sizeof text)) {
    snprintf(err, errlen, "no PMU '%s' (%s)", pmu, path);
    return -1;
  }
  g_ev[i].type = (uint32_t)strtoul(text, NULL, 10);
  return 0;
}

static int event_resolve(int i, const char *text, char *err, size_t errlen) {
  const char *slash = strchr(text, '/');
  memset(&g_ev[i], 0, sizeof g_ev[i]);
  g_ev[i].fd = -1;
  for (size_t k = 0; k < sizeof generic_events / sizeof generic_events[0]; k++) {
    if (!strcmp(text, generic_events[k].name)) {
      g_ev[i].type = generic_events[k].type;
      g_ev[i].config = generic_events[k].config;
      return 0;
    }
  }
  if (text[0] == 'r' && text[1] && strspn(text + 1, "0123456789abcdefABCDEF") == strlen(text + 1)) {
    g_ev[i].type = PERF_TYPE_RAW;
    g_ev[i].config = strtoull(text + 1, NULL, 16);
    return 0;
  }
  if (slash) { /* pmu/terms/ */
    char pmu[64], terms[256];
    size_t len = (size_t)(slash - text), tl = strlen(slash + 1);
    if (len == 0 || len >= sizeof pmu || tl == 0 || slash[tl] != '/' || tl > sizeof terms) {
      snprintf(err, errlen, "'%s': expected PMU/FIELDS/", text);
      return -1;
    }
    memcpy(pmu, text, len);
    pmu[len] = 0;
    memcpy(terms, slash + 1, tl - 1);
    terms[tl - 1] = 0;
    if (set_pmu(i, pmu, err, errlen)) return -1;
    return set_terms(i, pmu, terms, 0, err, errlen);
  }
  /* a bare name: the first PMU that has the event, the core's ("cpu" on x86) first */
  {
    static const char *const first[] = {"cpu", "cpu_core", "armv8_pmuv3_0", "armv8_pmuv3"};
    char path[512], line[256];
    for (size_t k = 0; k < sizeof first / sizeof first[0]; k++) {
      snprintf(path, sizeof path, PMU_DIR "/%s/events/%s", first[k], text);
      if (file_line(path, line, sizeof line)) {
        if (set_pmu(i, first[k], err, errlen)) return -1;
        return set_terms(i, first[k], line, 1, err, errlen);
      }
    }
    {
      DIR *d = opendir(PMU_DIR);
      struct dirent *e;
      while (d && (e = readdir(d)) != NULL) {
        if (e->d_name[0] == '.') continue;
        snprintf(path, sizeof path, PMU_DIR "/%s/events/%s", e->d_name, text);
        if (file_line(path, line, sizeof line)) {
          char pmu[256];
          int r;
          snprintf(pmu, sizeof pmu, "%s", e->d_name);
          closedir(d);
          if (set_pmu(i, pmu, err, errlen)) return -1;
          r = set_terms(i, pmu, line, 1, err, errlen);
          return r;
        }
      }
      if (d) closedir(d);
    }
  }
  snprintf(err, errlen, "unknown event '%s': not a generic name, no PMU in " PMU_DIR " has it, not rHEX", text);
  return -1;
}
#  endif

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

static int event_resolve(int i, const char *text, char *err, size_t errlen) {
  (void)i;
  snprintf(err, errlen, "%s: hardware events are counted only on Linux (perf_event)", text);
  return -1;
}
static int event_counted(int i) { (void)i; return 0; }

#endif

/* ------------------------------------------------------------------------
 * Common to all systems
 * ---------------------------------------------------------------------- */

int bench_platform_events(const char *list, char *err, size_t errlen) {
  const char *p = list;
  g_nevents = 0;
  while (*p) {
    /* the commas inside PMU/.../ belong to the event */
    size_t len = 0;
    int inside = 0;
    char text[sizeof g_event_names[0]];
    while (p[len] && (inside || p[len] != ',')) {
      if (p[len] == '/') inside = !inside && strchr(p + len + 1, '/') != NULL;
      len++;
    }
    if (len == 0 || len >= sizeof text) {
      snprintf(err, errlen, "'%.*s': expected an event of 1 to %d characters", (int)len, p, (int)sizeof text - 1);
      return -1;
    }
    if (g_nevents == BENCH_MAX_EVENTS) {
      snprintf(err, errlen, "at most %d events", BENCH_MAX_EVENTS);
      return -1;
    }
    memcpy(text, p, len);
    text[len] = 0;
    if (event_resolve(g_nevents, text, err, errlen)) return -1;
    memcpy(g_event_names[g_nevents++], text, len + 1);
    p += len + (p[len] == ',');
  }
  return 0;
}

int bench_event_count(void) { return g_nevents; }
int bench_event_counted(int i) { return i >= 0 && i < g_nevents && event_counted(i); }
const char *bench_event_name(int i) { return g_event_names[i]; }

int bench_platform_init(int pin_cpu) {
  /* Undocumented switch for testing the wall-clock fallback on any machine. */
  const char *off = getenv("BENCH_NO_COUNTER");
  volatile unsigned spin = 0;
  bench_snap start, a, b;
  steer(pin_cpu);
  if (off && *off && strcmp(off, "0") != 0) return 0;
  if (!counter_open()) {
    /* a group with the events that the PMU cannot schedule: count without them */
    if (!g_nevents) return 0;
    fprintf(stderr, "bench_fft: the counter does not open with the events: counting without them\n");
    g_nevents = 0;
    if (!counter_open()) return 0;
  }
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
      if (!g_nevents) return 0;
      /* a group that stands still may be one the PMU cannot schedule */
      fprintf(stderr, "bench_fft: the counter stands still with the events: counting without them\n");
      g_nevents = 0;
      return bench_platform_init(pin_cpu);
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
