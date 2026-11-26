/* benchmarks/bench_fft.c -- the driver of bench_fft, the pffft FFT benchmark:
 * command line, measurement loop, statistics, output.
 *
 * The cost of a transform is the fastest of many short batches, in CPU cycles
 * where the system counts them, while all engines take turns in random order
 * inside this one process. benchmarks/README.md gives the reasons. */
#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS /* fopen, gmtime */
#endif

#include "bench.h"

#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if !defined(BENCH_HAVE_F32) && !defined(BENCH_HAVE_F64)
#  error "the build must define BENCH_HAVE_F32, BENCH_HAVE_F64 or both"
#endif

/* The name of pffft's SIMD backend goes into the CSV header. Nothing else of
 * pffft is used here, and its headers do not all pass -Wpedantic. */
#ifdef BENCH_HAVE_F32
const char *pffft_simd_arch(void);
#endif
#ifdef BENCH_HAVE_F64
const char *pffftd_simd_arch(void);
#endif

#define BENCH_STR_(x) #x
#define BENCH_STR(x) BENCH_STR_(x)
#if defined(__clang__)
#  define BENCH_COMPILER "clang " __clang_version__
#elif defined(__GNUC__)
#  define BENCH_COMPILER "gcc " __VERSION__
#elif defined(_MSC_VER)
#  define BENCH_COMPILER "msvc " BENCH_STR(_MSC_VER)
#else
#  define BENCH_COMPILER "unknown"
#endif

#ifndef BENCH_TARGET_ARCH
#  define BENCH_TARGET_ARCH "" /* TARGET_C_ARCH of the build: -march or /arch */
#endif

#define MAX_SIZES 512
#define MAX_LIBS 8
#define MAX_PAIRS (1 << 24)
#define MIN_ROUNDS 10     /* --case-seconds never ends a trial before this */
/* nor before every engine has this many clean batches, unless the trial has
 * taken CLEAN_PATIENCE times --case-seconds: a batch of a large n takes tens
 * of ms, and on Apple systems, which cannot pin, all 10 of them could touch
 * an efficiency core, leaving no cycle count */
#define MIN_CLEAN 3
#define CLEAN_PATIENCE 5
#define PLACE_RANGE 16384 /* later trials shift each array by less than this */
/* An engine whose batches ran at a clock this much below that of the batches
 * around them slowed the clock itself, and its cycles are charged for it (see
 * add_result). Other engines stay within 0.5%. */
#define CLOCK_SLACK 0.02

enum { F32, F64, NPREC };
enum { SHOW_CYCLES, SHOW_NS, SHOW_MFLOPS };

static const char *const prec_name[NPREC] = {"f32", "f64"};
static const char *const xform_name[2] = {"real", "cplx"};

int bench_fftw_measure_max = 65536;

typedef struct options {
  const char *sizes;   /* NULL: the default sweep */
  const char *engines; /* NULL: all built-in engines, or all --lib engines */
  const char *csv;
  int show; /* SHOW_*, or -1: cycles if counted, else ns */
  struct {
    char label[41]; /* "pffft@" and the label fit a bench_engine name */
    const char *path;
  } lib[MAX_LIBS];
  int nlibs;
  int prec[NPREC], xform[2];
  int min_n, max_n, trials, rounds, pin, quiet, list;
  int accuracy, seeds, bits_only; /* --accuracy, --seeds, --bits-only */
  double batch_us, case_seconds;
  unsigned long seed;
  int have_seed;
} options;

/* The engines of one precision: all that exist, and those that take part. */
typedef struct registry {
  bench_engine all[BENCH_MAX_ENGINES];
  int nall, nbuiltin;
  const bench_engine *sel[BENCH_MAX_ENGINES];
  int nsel;
} registry;

/* One engine during one trial of one case. */
typedef struct slot {
  const bench_engine *e;
  bench_instance inst;
  int column; /* index among the selected engines */
  int pairs, batches, clean;
  double *ns;     /* what the clock saw of every batch, overhead included */
  double *cycles; /* the same in CPU cycles, clean batches only */
  double *ghz;    /* cycles per nanosecond of the clean batches */
  double *ratios; /* the same relative to the batches around them */
  int nratios;
  double best_cycles, best_instructions; /* of the cheapest clean batch */
  double clock;   /* the median of ratios[], 0 if unknown */
} slot;

/* One line of the CSV. A negative value was not measured. */
typedef struct result {
  const bench_engine *e;
  int prec, cplx, n, trial, pairs, batches, clean;
  double ns, ns_median, cycles, cycles_median, instructions;
  double ghz;   /* the engine's clock, cycles per nanosecond */
  double clock; /* its clock relative to the other engines' (1: the same) */
} result;

static registry g_reg[NPREC];
static result *g_results;
static int g_nresults, g_maxresults;
static int g_counted; /* is there a cycle counter? */
static int g_failed;  /* a self-check or a zero check failed */
static uint64_t g_rng;
static struct {
  uint64_t ns, cycles, instructions;
  int known;
} g_overhead; /* cost of an empty bench_snap_begin/end pair: the smallest seen */

static void die(const char *fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  fprintf(stderr, "bench_fft: ");
  vfprintf(stderr, fmt, ap);
  fprintf(stderr, "\n");
  va_end(ap);
  exit(2);
}

/* splitmix64: the C library's rand() differs between systems, and a --seed
 * should give the same order and placement everywhere. */
static uint64_t rnd(void) {
  uint64_t z = (g_rng += 0x9E3779B97F4A7C15u);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
  return z ^ (z >> 31);
}

static void shuffle(int *a, int count) {
  for (int i = count - 1; i > 0; i--) {
    int j = (int)(rnd() % (uint64_t)(i + 1)), t = a[i];
    a[i] = a[j];
    a[j] = t;
  }
}

static int compare_doubles(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

/* Sorts v, so that v[0] is the minimum afterwards. */
static double median(double *v, int n) {
  qsort(v, (size_t)n, sizeof *v, compare_doubles);
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static double scalar_at(const void *array, size_t i, int prec) {
  return prec == F32 ? (double)((const float *)array)[i] : ((const double *)array)[i];
}

/* ------------------------------------------------------------------------
 * Engines
 * ---------------------------------------------------------------------- */

static int builtin_engines(int prec, bench_engine *list, int max) {
#ifdef BENCH_HAVE_F32
  if (prec == F32) return bench_engines_f32(list, max);
#endif
#ifdef BENCH_HAVE_F64
  if (prec == F64) return bench_engines_f64(list, max);
#endif
  (void)list;
  (void)max;
  return 0;
}

static int lib_engines(int prec, const char *label, const char *path, bench_engine *list, int max, char *err,
                       size_t errlen) {
#ifdef BENCH_HAVE_F32
  if (prec == F32) return bench_lib_engines_f32(label, path, list, max, err, errlen);
#endif
#ifdef BENCH_HAVE_F64
  if (prec == F64) return bench_lib_engines_f64(label, path, list, max, err, errlen);
#endif
  (void)label, (void)path, (void)list, (void)max, (void)err, (void)errlen;
  return 0;
}

static const char *simd_name(void) {
  static char both[96];
  const char *f32 = NULL, *f64 = NULL;
#ifdef BENCH_HAVE_F32
  f32 = pffft_simd_arch();
#endif
#ifdef BENCH_HAVE_F64
  f64 = pffftd_simd_arch();
#endif
  if (!f32 || !f64) return f32 ? f32 : f64;
  if (!strcmp(f32, f64)) return f32;
  snprintf(both, sizeof both, "%s (f32), %s (f64)", f32, f64);
  return both;
}

/* Copies the next item of a comma-separated list; returns 0 at its end. */
static int next_item(const char **list, char *item, size_t size) {
  const char *s = *list, *comma;
  size_t len;
  if (!*s) return 0;
  comma = strchr(s, ',');
  len = comma ? (size_t)(comma - s) : strlen(s);
  if (len >= size) die("list item is too long: %s", s);
  memcpy(item, s, len);
  item[len] = 0;
  *list = comma ? comma + 1 : s + len;
  return 1;
}

static void select_engine(registry *r, const bench_engine *e) {
  for (int i = 0; i < r->nsel; i++) {
    if (r->sel[i] == e) return;
  }
  r->sel[r->nsel++] = e;
}

static void load_engines(const options *o) {
  char name[64], err[256];
  int loaded[MAX_LIBS] = {0}; /* engines of each --lib */
  for (int p = 0; p < NPREC; p++) {
    registry *r = &g_reg[p];
    if (!o->prec[p]) continue;
    r->nall = r->nbuiltin = builtin_engines(p, r->all, BENCH_MAX_ENGINES);
    for (int i = 0; i < o->nlibs; i++) {
      int added;
      err[0] = 0;
      added = lib_engines(p, o->lib[i].label, o->lib[i].path, r->all + r->nall, BENCH_MAX_ENGINES - r->nall, err,
                          sizeof err);
      if (added < 0) die("--lib %s=%s: %s", o->lib[i].label, o->lib[i].path, err);
      r->nall += added;
      loaded[i] += added;
    }
    if (!o->engines) {
      /* with --lib the comparison of the libraries is the point */
      for (int i = o->nlibs ? r->nbuiltin : 0; i < (o->nlibs ? r->nall : r->nbuiltin); i++)
        if (!r->all[i].on_request) select_engine(r, &r->all[i]);
    }
  }
  for (int i = 0; i < o->nlibs; i++) {
    if (!loaded[i]) die("--lib %s=%s: no pffft of the selected precisions in it", o->lib[i].label, o->lib[i].path);
  }
  if (o->engines) {
    const char *list = o->engines;
    while (next_item(&list, name, sizeof name)) {
      int found = 0;
      for (int p = 0; p < NPREC; p++) {
        for (int i = 0; i < g_reg[p].nall; i++) {
          if (strcmp(g_reg[p].all[i].name, name) != 0) continue;
          select_engine(&g_reg[p], &g_reg[p].all[i]);
          found = 1;
        }
      }
      if (!found) die("no engine '%s' (--list prints the engines)", name);
    }
  }
}

/* ------------------------------------------------------------------------
 * Measuring
 * ---------------------------------------------------------------------- */

/* A reading with the cost of taking it removed. The overhead is the smallest
 * one seen, so a real batch costs more than that; where one does not (a
 * coarse clock) the result still has to stay above zero, or the minimum over
 * the batches would be meaningless. */
static double net(double gross, uint64_t overhead) {
  if (gross > (double)overhead) return gross - (double)overhead;
  return gross > 0 ? 0.5 * gross : 0.5;
}

/* One empty begin/end pair. With a counter only pairs on the fast cores
 * count: a program starts out on a slow core, where the pair takes longer and
 * its cycles are not in the count. That is also why sampling goes on during
 * the rounds instead of ending at start-up. */
static void sample_overhead(void) {
  bench_snap a, b;
  bench_snap_begin(&a);
  bench_snap_end(&b);
  if (g_counted && (b.stray_cycles != a.stray_cycles || b.cycles == a.cycles)) return;
  if (!g_overhead.known || b.ns - a.ns < g_overhead.ns) g_overhead.ns = b.ns - a.ns;
  if (!g_overhead.known || b.cycles - a.cycles < g_overhead.cycles) g_overhead.cycles = b.cycles - a.cycles;
  if (!g_overhead.known || b.instructions - a.instructions < g_overhead.instructions)
    g_overhead.instructions = b.instructions - a.instructions;
  g_overhead.known = 1;
}

/* One pair on random input must give back the input times the engine's gain.
 * Returns the largest error relative to the largest expected value. */
static double self_check(const slot *s, int prec) {
  const size_t len = s->inst.len, bytes = len * (prec == F32 ? sizeof(float) : sizeof(double));
  double worst = 0, largest = 0;
  void *original = malloc(bytes ? bytes : 1);
  if (!original) die("out of memory");
  for (size_t i = 0; i < len; i++) {
    double v = (double)(rnd() >> 11) * (2.0 / 9007199254740992.0) - 1.0; /* [-1, 1) */
    if (prec == F32) ((float *)s->inst.in)[i] = (float)v;
    else ((double *)s->inst.in)[i] = v;
  }
  memcpy(original, s->inst.in, bytes);
  s->e->run(&s->inst, 1);
  for (size_t i = 0; i < len; i++) {
    double want = s->inst.gain * scalar_at(original, i, prec);
    double err = fabs(scalar_at(s->inst.out, i, prec) - want);
    if (fabs(want) > largest) largest = fabs(want);
    if (!(err <= worst)) worst = err;
    if (worst != worst) break; /* NaN: nothing can make it better */
  }
  free(original);
  return largest > 0 ? worst / largest : HUGE_VAL;
}

/* -0.0 counts as zero, NaN does not. */
static int all_zero(const bench_bufs *bufs, int prec) {
  const size_t count = bufs->bytes / (prec == F32 ? sizeof(float) : sizeof(double));
  for (int b = 0; b < 4; b++) {
    for (size_t i = 0; i < count; i++) {
      if (scalar_at(bufs->buf[b], i, prec) != 0) return 0;
    }
  }
  return 1;
}

/* The number of pairs that fill a batch of target_ns. Every step is judged by
 * the fastest of three attempts, so that a slow first call (cold caches, lazy
 * binding) cannot settle for a batch that is far too short. For the same
 * reason attempts that strayed to a slow core are repeated, within limits. */
static int calibrate(const slot *s, double target_ns) {
  int pairs = 1;
  double best, want;
  for (;;) {
    double fastest = 0;
    int good = 0;
    best = 0;
    for (int tries = 0; tries < 6 && good < 3; tries++) {
      bench_snap a, b;
      double ns;
      bench_snap_begin(&a);
      s->e->run(&s->inst, pairs);
      bench_snap_end(&b);
      ns = net((double)(b.ns - a.ns), g_overhead.ns);
      if (tries == 0 || ns < fastest) fastest = ns;
      if (b.stray_cycles != a.stray_cycles) continue;
      if (good++ == 0 || ns < best) best = ns;
    }
    if (!good) best = fastest;
    /* long enough to be measured well: the rest is proportion */
    if (8 * best >= target_ns || pairs >= MAX_PAIRS) break;
    pairs *= 4;
  }
  want = floor(pairs * target_ns / best + 0.5);
  return want < 1 ? 1 : want > MAX_PAIRS ? MAX_PAIRS : (int)want;
}

/* Returns the clock of the batch in cycles per nanosecond, 0 if it is not
 * clean. */
static double record(slot *s, const bench_snap *a, const bench_snap *b) {
  s->ns[s->batches++] = (double)(b->ns - a->ns);
  /* Clean: no cycles on another kind of core, which would be missing from the
   * count. A counter that stood still did not count this batch at all. */
  if (g_counted && b->stray_cycles == a->stray_cycles && b->cycles != a->cycles && b->ns > a->ns) {
    double cycles = (double)(b->cycles - a->cycles);
    if (s->clean == 0 || cycles < s->best_cycles) {
      s->best_cycles = cycles;
      s->best_instructions = (double)(b->instructions - a->instructions);
    }
    s->ghz[s->clean] = cycles / (double)(b->ns - a->ns);
    s->cycles[s->clean] = cycles;
    return s->ghz[s->clean++];
  }
  return 0;
}

/* The clock of each engine relative to the others': the median over its clean
 * batches of the batch's clock over the faster of the nearest clean batches
 * of other engines before and after it. Neighbours, so that the system
 * changing the clock during a trial does not count; the faster one, because
 * the clock recovers only some 40 us after AMX work, slowing the batch that
 * follows. who[] and ghz[] are the slot and record() of each batch in turn. */
static void clock_ratios(slot *slots, int count, const int *who, const double *ghz, int nseq) {
  for (int i = 0; i < nseq; i++) {
    double ref = 0;
    if (ghz[i] <= 0) continue;
    for (int j = i - 1; j >= 0 && ref == 0; j--)
      if (who[j] != who[i] && ghz[j] > 0) ref = ghz[j];
    for (int j = i + 1; j < nseq; j++) {
      if (who[j] == who[i] || ghz[j] <= 0) continue;
      if (ghz[j] > ref) ref = ghz[j];
      break;
    }
    if (ref > 0) {
      slot *s = &slots[who[i]];
      s->ratios[s->nratios++] = ghz[i] / ref;
    }
  }
  for (int i = 0; i < count; i++) slots[i].clock = slots[i].nratios ? median(slots[i].ratios, slots[i].nratios) : 0;
}

/* Turns the batches of a slot into the cost of one transform.
 *
 * Cycles do not count a clock that the engine itself slows: on Apple chips
 * the clock drops by 10-25% while the AMX matrix unit works (vDSP uses it for
 * f32 FFTs from n=256 or 512 on an M2), much as AVX-512 slows Intel cores.
 * The engines take turns batch by batch, so a clock that differs between
 * them is their own doing, not the system's. An engine that runs at a lower
 * clock than the others is charged its cycles at theirs, which makes it
 * compare by time. */
static void add_result(slot *s, int prec, int cplx, int n, int trial) {
  const double per_transform = 0.5 / s->pairs;
  result *r;
  if (g_nresults == g_maxresults) {
    g_maxresults = g_maxresults ? 2 * g_maxresults : 256;
    g_results = (result *)realloc(g_results, (size_t)g_maxresults * sizeof *g_results);
    if (!g_results) die("out of memory");
  }
  r = &g_results[g_nresults++];
  r->e = s->e;
  r->prec = prec;
  r->cplx = cplx;
  r->n = n;
  r->trial = trial;
  r->pairs = s->pairs;
  r->batches = s->batches;
  r->clean = s->clean;
  r->ns_median = net(median(s->ns, s->batches), g_overhead.ns) * per_transform;
  r->ns = net(s->ns[0], g_overhead.ns) * per_transform;
  r->cycles = r->cycles_median = r->instructions = r->ghz = r->clock = -1;
  if (s->clean > 0) {
    double charge = 1;
    r->ghz = median(s->ghz, s->clean);
    if (s->clock > 0) {
      r->clock = s->clock;
      if (s->clock < 1 - CLOCK_SLACK) charge = 1 / s->clock;
    }
    r->cycles_median = net(median(s->cycles, s->clean), g_overhead.cycles) * per_transform * charge;
    r->cycles = net(s->best_cycles, g_overhead.cycles) * per_transform * charge;
    if (s->best_instructions > 0) /* 0: the system counts cycles only */
      r->instructions = net(s->best_instructions, g_overhead.instructions) * per_transform;
  }
}

/* One trial of one case: every engine that can do the case is set up afresh
 * and timed against the others. dropped[] marks the selected engines that
 * failed the self-check of this case. */
static void run_trial(const options *o, const registry *reg, int prec, int cplx, int n, int trial, char *dropped) {
  const size_t scalar = prec == F32 ? sizeof(float) : sizeof(double);
  slot slots[BENCH_MAX_ENGINES];
  int order[BENCH_MAX_ENGINES];
  void *block[4], *dummy = NULL;
  bench_bufs bufs;
  bench_snap start, a, b;
  uint64_t now;
  int count = 0;

  /* Near the L1 size the cost depends by up to half a percent on where the
   * arrays and the engines' tables fall in memory. Trial 0 takes the memory
   * as it comes: arrays of exactly the size, as a program would allocate
   * them. Later trials shift the arrays and, through a block of random size
   * that is allocated first, whatever the engines allocate. */
  bufs.bytes = ((2 * (size_t)n + 64) * scalar + 63) & ~(size_t)63;
  if (trial > 0) dummy = malloc(1 + (size_t)(rnd() % (bufs.bytes + PLACE_RANGE)));
  for (int i = 0; i < 4; i++) {
    size_t shift = trial > 0 ? 64 * (size_t)(rnd() % (PLACE_RANGE / 64)) : 0;
    block[i] = bench_aligned_alloc(bufs.bytes + (trial > 0 ? PLACE_RANGE : 0));
    if (!block[i]) die("out of memory");
    bufs.buf[i] = (char *)block[i] + shift;
    memset(bufs.buf[i], 0, bufs.bytes);
  }

  /* In a random order: at small n an engine can be a few percent faster just
   * for being created after another (2.5% for f32 real n=64 on an M2, with
   * two copies of one library), presumably through where its tables fall.
   * A fixed order would hand that to the same engine in every trial. */
  for (int c = 0; c < reg->nsel; c++) order[c] = c;
  shuffle(order, reg->nsel);
  for (int k = 0; k < reg->nsel; k++) {
    int c = order[k];
    slot *s = &slots[count];
    if (dropped[c]) continue;
    memset(s, 0, sizeof *s);
    s->e = reg->sel[c];
    s->column = c;
    if (s->e->create(s->e, n, cplx, &bufs, &s->inst) == 0) count++; /* else not a case for this engine */
  }

  /* Only now, when all engines exist: memory that the check allocates and
   * frees between two engines would give the later one a different layout. */
  if (trial == 0) {
    int kept = 0;
    for (int i = 0; i < count; i++) {
      slot *s = &slots[i];
      double err = self_check(s, prec);
      if (err <= (prec == F32 ? 1e-3 : 1e-9)) {
        slots[kept++] = *s;
        continue;
      }
      fprintf(stderr, "bench_fft: %s fails the self-check for %s %s n=%d (error %.3g of the largest value): dropped\n",
              s->e->name, prec_name[prec], xform_name[cplx], n, err);
      s->e->destroy(&s->inst);
      dropped[s->column] = 1;
      g_failed = 1;
    }
    count = kept;
  }
  for (int i = 0; i < count; i++) {
    slots[i].ns = (double *)calloc((size_t)o->rounds, sizeof(double)); /* checks the product */
    slots[i].cycles = (double *)calloc((size_t)o->rounds, sizeof(double));
    slots[i].ghz = (double *)calloc((size_t)o->rounds, sizeof(double));
    slots[i].ratios = (double *)calloc((size_t)o->rounds, sizeof(double));
    if (!slots[i].ns || !slots[i].cycles || !slots[i].ghz || !slots[i].ratios) die("out of memory");
    order[i] = i;
  }

  if (count > 0) {
    /* Zeros stay zeros under any transform, so the engines can share the
     * arrays and in-place ones cannot overflow. */
    for (int i = 0; i < 4; i++) memset(bufs.buf[i], 0, bufs.bytes);
    for (int i = 0; i < count; i++) slots[i].pairs = calibrate(&slots[i], o->batch_us * 1000.0);

    int nseq = 0, *who = (int *)malloc((size_t)o->rounds * (size_t)count * sizeof(int));
    double *ghz = (double *)malloc((size_t)o->rounds * (size_t)count * sizeof(double));
    if (!who || !ghz) die("out of memory");
    bench_snap_begin(&start);
    now = start.ns;
    for (int round = 0; round < o->rounds; round++) {
      shuffle(order, count);
      for (int i = 0; i < count; i++) {
        slot *s = &slots[order[i]];
        bench_snap_begin(&a);
        s->e->run(&s->inst, s->pairs);
        bench_snap_end(&b);
        who[nseq] = order[i];
        ghz[nseq++] = record(s, &a, &b);
        now = b.ns;
      }
      sample_overhead();
      if (round + 1 >= MIN_ROUNDS && (double)(now - start.ns) >= o->case_seconds * 1e9) {
        int short_of_clean = 0;
        for (int i = 0; i < count && g_counted; i++) short_of_clean |= slots[i].clean < MIN_CLEAN;
        if (!short_of_clean || (double)(now - start.ns) >= CLEAN_PATIENCE * o->case_seconds * 1e9) break;
      }
    }

    if (!all_zero(&bufs, prec)) {
      fprintf(stderr, "bench_fft: %s %s n=%d: the arrays are not zero after timing; an engine wrote garbage\n",
              prec_name[prec], xform_name[cplx], n);
      g_failed = 1;
    }
    clock_ratios(slots, count, who, ghz, nseq);
    free(who);
    free(ghz);
    for (int i = 0; i < count; i++) {
      add_result(&slots[i], prec, cplx, n, trial);
      slots[i].e->destroy(&slots[i].inst);
      free(slots[i].ns);
      free(slots[i].cycles);
      free(slots[i].ghz);
      free(slots[i].ratios);
    }
  }
  for (int i = 0; i < 4; i++) bench_aligned_free(block[i]);
  free(dummy);
}

/* ------------------------------------------------------------------------
 * Output
 * ---------------------------------------------------------------------- */

static double mflops(const result *r) {
  return (r->cplx ? 5.0 : 2.5) * r->n * (log((double)r->n) / log(2.0)) / r->ns * 1000.0;
}

static void write_header(FILE *f, const char *prefix, const options *o, const char *date) {
  fprintf(f, "%sdate: %s\n", prefix, date);
  bench_platform_describe(f, prefix);
  fprintf(f, "%scompiler: %s\n", prefix, BENCH_COMPILER);
  fprintf(f, "%spffft_simd: %s\n", prefix, simd_name());
  fprintf(f, "%starget_arch: %s\n", prefix, BENCH_TARGET_ARCH[0] ? BENCH_TARGET_ARCH : "none");
  fprintf(f, "%strials: %d\n", prefix, o->trials);
  fprintf(f, "%srounds: %d\n", prefix, o->rounds);
  fprintf(f, "%sbatch_us: %g\n", prefix, o->batch_us);
  fprintf(f, "%sseed: %lu\n", prefix, o->seed);
}

static void write_csv(FILE *f, const options *o, const char *date) {
  fprintf(f, "# pffft-bench 1\n");
  write_header(f, "# ", o, date);
  fprintf(f, "engine,prec,xform,n,trial,pairs,batches,clean,ns,ns_median,cycles,cycles_median,instructions,mflops,ghz,clock\n");
  for (int i = 0; i < g_nresults; i++) {
    const result *r = &g_results[i];
    fprintf(f, "%s,%s,%s,%d,%d,%d,%d,%d,%.3f,%.3f,", r->e->name, prec_name[r->prec], xform_name[r->cplx], r->n,
            r->trial, r->pairs, r->batches, r->clean, r->ns, r->ns_median);
    if (r->cycles >= 0) fprintf(f, "%.3f,%.3f,", r->cycles, r->cycles_median);
    else fprintf(f, ",,");
    if (r->instructions >= 0) fprintf(f, "%.3f,", r->instructions);
    else fprintf(f, ",");
    fprintf(f, "%.1f,", mflops(r));
    if (r->ghz > 0) fprintf(f, "%.4f,", r->ghz);
    else fprintf(f, ",");
    if (r->clock > 0) fprintf(f, "%.4f\n", r->clock);
    else fprintf(f, "\n");
  }
}

static int column_width(const bench_engine *e) {
  int len = (int)strlen(e->name);
  return len < 10 ? 10 : len;
}

static void print_table_header(const registry *reg, int prec, int cplx, int show, int trials) {
  static const char *const what[3] = {"CPU cycles per transform", "nanoseconds per transform", "MFLOPS"};
  printf("\n%s %s: %s, fastest batch", prec_name[prec], xform_name[cplx], what[show]);
  if (trials > 1) printf(", median of %d trials", trials);
  printf("\n%9s ", "n");
  for (int c = 0; c < reg->nsel; c++) printf("  %*s", column_width(reg->sel[c]), reg->sel[c]->name);
  printf("\n");
}

/* The results of one case, all its trials: one table row. Returns 1 if an
 * engine was measured but has no value, because no batch was clean. */
static int print_table_row(const registry *reg, const result *r, int count, int show) {
  double *v = (double *)malloc((size_t)count * sizeof(double));
  int unclean = 0;
  if (!v) die("out of memory");
  printf("%9d%c", r->n, (r->n & (r->n - 1)) ? '*' : ' ');
  for (int c = 0; c < reg->nsel; c++) {
    int k = 0, measured = 0;
    for (int i = 0; i < count; i++) {
      double x = show == SHOW_CYCLES ? r[i].cycles : show == SHOW_NS ? r[i].ns : mflops(&r[i]);
      if (r[i].e != reg->sel[c]) continue;
      measured = 1;
      if (x >= 0) v[k++] = x;
    }
    if (k == 0) {
      printf("  %*s", column_width(reg->sel[c]), measured ? "?" : "-");
      unclean |= measured;
    } else {
      double m = median(v, k);
      printf("  %*.*f", column_width(reg->sel[c]), m < 100 ? 2 : 1, m);
    }
  }
  printf("\n");
  fflush(stdout);
  free(v);
  return unclean;
}

/* ------------------------------------------------------------------------
 * Command line
 * ---------------------------------------------------------------------- */

static void usage(FILE *f) {
  fprintf(f,
          "bench_fft [options]\n"
          "  --sizes LIST        sizes and/or keywords, comma separated: pow2; picked,\n"
          "                      28 sizes 96..1000000 with factors 3 and 5 up to 5^6;\n"
          "                      nonpow2, 3, 5, 9 and 15 times 2^k; all, pow2 and\n"
          "                      nonpow2 (default: pow2,picked)\n"
          "  --min N --max N     limits for the keyword lists (defaults 32 and 1048576)\n"
          "  --prec LIST         f32,f64 (default: every precision built in)\n"
          "  --xform LIST        real,cplx (default: both)\n"
          "  --engines LIST      engine names (default: all built-in engines but pffft-ip\n"
          "                      and pffft-un; if --lib is given, only the --lib engines)\n"
          "  --lib LABEL=PATH    add the engine pffft@LABEL from a pffft shared library\n"
          "                      (repeatable)\n"
          "  --trials T          default 1\n"
          "  --rounds R          default 100\n"
          "  --batch-us U        default 300\n"
          "  --case-seconds S    default 1.0\n"
          "  --csv FILE          write results as CSV; '-' means stdout\n"
          "  --quiet             do not print the table\n"
          "  --show WHAT         table values: cycles, ns or mflops (default: cycles if\n"
          "                      counted, else ns)\n"
          "  --seed S            default: from the clock\n"
          "  --cpu N | --no-pin  core selection where pinning exists (Linux, Android)\n"
          "  --fftw-measure-max N  default 65536\n"
          "  --list              print the engines of each precision and exit\n"
          "  --accuracy          instead of timing: the error of pffft (the --lib ones, or\n"
          "                      the one built in) against an exact FFT, and the bitwise\n"
          "                      equivalences of one build (see accuracy.c); exits 1 if\n"
          "                      a check fails\n"
          "  --seeds K           with --accuracy: inputs for at least K * 4096 points per\n"
          "                      case, at least one (default 4)\n"
          "  --bits-only         with --accuracy: only the checks and the hashes of the\n"
          "                      results, not the errors (no exact transform; fast)\n"
          "  --help\n");
}

static int int_arg(const char *option, const char *text, int lowest) {
  char *end;
  long long v = strtoll(text, &end, 10);
  if (end == text || *end || v < lowest || v > INT_MAX) die("%s: bad value '%s'", option, text);
  return (int)v;
}

static double real_arg(const char *option, const char *text) {
  char *end;
  double v = strtod(text, &end);
  if (end == text || *end || !(v >= 0 && v < HUGE_VAL)) die("%s: bad value '%s'", option, text);
  return v;
}

static void parse_flags(const char *option, const char *list, const char *const *names, int count, int *flags) {
  char item[16];
  memset(flags, 0, (size_t)count * sizeof *flags);
  while (next_item(&list, item, sizeof item)) {
    int i = 0;
    while (i < count && strcmp(item, names[i]) != 0) i++;
    if (i == count) die("%s: unknown '%s'", option, item);
    flags[i] = 1;
  }
}

static int add_series(int *sizes, int count, int odd, int lowest, int highest) {
  for (int64_t n = odd; n <= highest; n *= 2) {
    if (n < lowest) continue;
    if (count == MAX_SIZES) die("too many sizes");
    sizes[count++] = (int)n;
  }
  return count;
}

static int compare_ints(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}

/* Expands --sizes into a sorted list without duplicates. */
static int expand_sizes(const options *o, int *sizes) {
  static const int odd[5] = {1, 3, 5, 9, 15};
  /* "picked": one to six passes of radix 3 and 5 (odd parts 3 to 5^6), so that
   * the cost of each shows, the ido == 3 fast path of radix 3 (9 * 2^k: 1152,
   * 9216), and lengths people use (48000); 20480 fills the gap between 15360
   * and 32000. All are multiples of 32, so every engine can do them, real and
   * complex. */
  static const int picked[] = {96,    160,   192,   320,   480,   640,   768,    960,    1152,   1280,
                               1440,  1600,  1920,  2400,  2560,  3200,  4000,   4800,   7200,   9216,
                               15360, 20480, 32000, 48000, 72000, 100000, 300000, 1000000};
  const char *list = o->sizes ? o->sizes : "default";
  char item[32];
  int count = 0, unique = 0;
  while (next_item(&list, item, sizeof item)) {
    int first = 0, last = -1; /* range of odd[] */
    if (!strcmp(item, "pow2")) last = 0;
    else if (!strcmp(item, "nonpow2")) first = 1, last = 4;
    else if (!strcmp(item, "all")) last = 4;
    else if (!strcmp(item, "picked") || !strcmp(item, "default")) {
      if (item[0] == 'd') last = 0; /* default: pow2,picked */
      for (size_t i = 0; i < sizeof picked / sizeof picked[0]; i++) {
        if (picked[i] < o->min_n || picked[i] > o->max_n) continue;
        if (count == MAX_SIZES) die("too many sizes");
        sizes[count++] = picked[i];
      }
    } else {
      if (count == MAX_SIZES) die("too many sizes");
      sizes[count++] = int_arg("--sizes", item, 1);
    }
    for (int i = first; i <= last; i++) count = add_series(sizes, count, odd[i], o->min_n, o->max_n);
  }
  qsort(sizes, (size_t)count, sizeof *sizes, compare_ints);
  for (int i = 0; i < count; i++) {
    if (unique == 0 || sizes[i] != sizes[unique - 1]) sizes[unique++] = sizes[i];
  }
  /* The arrays and the engines' tables take up to 32 bytes per point. Where
   * size_t has 32 bits, a larger size overflows the computed sizes of the
   * allocations instead of failing to allocate. */
  if (unique > 0 && (size_t)sizes[unique - 1] > SIZE_MAX / 64)
    die("--sizes: %d is too large for this system", sizes[unique - 1]);
  return unique;
}

static void parse_options(int argc, char **argv, options *o) {
  static const char *const show_names[3] = {"cycles", "ns", "mflops"};
  memset(o, 0, sizeof *o);
  o->show = -1;
  o->min_n = 32;
  o->max_n = 1048576;
  o->trials = 1;
  o->seeds = 4;
  o->rounds = 100;
  o->batch_us = 300;
  o->case_seconds = 1.0;
  o->pin = BENCH_PIN_AUTO;
  o->xform[0] = o->xform[1] = 1;
#ifdef BENCH_HAVE_F32
  o->prec[F32] = 1;
#endif
#ifdef BENCH_HAVE_F64
  o->prec[F64] = 1;
#endif
  for (int i = 1; i < argc; i++) {
    const char *a = argv[i], *v;
    if (!strcmp(a, "--help") || !strcmp(a, "-h")) {
      usage(stdout);
      exit(0);
    }
    if (!strcmp(a, "--list")) { o->list = 1; continue; }
    if (!strcmp(a, "--accuracy")) { o->accuracy = 1; continue; }
    if (!strcmp(a, "--bits-only")) { o->bits_only = 1; continue; }
    if (!strcmp(a, "--quiet")) { o->quiet = 1; continue; }
    if (!strcmp(a, "--no-pin")) { o->pin = BENCH_PIN_NONE; continue; }
    if (i + 1 >= argc) die("%s: unknown option, or its value is missing (see --help)", a);
    v = argv[++i];
    if (!strcmp(a, "--sizes")) o->sizes = v;
    else if (!strcmp(a, "--min")) o->min_n = int_arg(a, v, 1);
    else if (!strcmp(a, "--max")) o->max_n = int_arg(a, v, 1);
    else if (!strcmp(a, "--prec")) {
      int built[NPREC];
      memcpy(built, o->prec, sizeof built);
      parse_flags(a, v, prec_name, NPREC, o->prec);
      for (int p = 0; p < NPREC; p++) {
        if (o->prec[p] && !built[p]) die("--prec: %s is not built into this program", prec_name[p]);
      }
    } else if (!strcmp(a, "--xform")) parse_flags(a, v, xform_name, 2, o->xform);
    else if (!strcmp(a, "--engines")) o->engines = v;
    else if (!strcmp(a, "--lib")) {
      const char *eq = strchr(v, '=');
      size_t len = eq ? (size_t)(eq - v) : 0;
      if (len == 0 || len >= sizeof o->lib[0].label || !eq[1] || strcspn(v, ", \t@") < len)
        die("--lib: expected LABEL=PATH, LABEL 1 to 40 characters without ',', '@' or blanks; got '%s'", v);
      if (o->nlibs == MAX_LIBS) die("--lib: at most %d libraries", MAX_LIBS);
      memcpy(o->lib[o->nlibs].label, v, len);
      for (int k = 0; k < o->nlibs; k++) {
        if (!strcmp(o->lib[k].label, o->lib[o->nlibs].label)) die("--lib: the label '%s' is given twice", o->lib[k].label);
      }
      o->lib[o->nlibs++].path = eq + 1;
    } else if (!strcmp(a, "--trials")) o->trials = int_arg(a, v, 1);
    else if (!strcmp(a, "--seeds")) o->seeds = int_arg(a, v, 1);
    else if (!strcmp(a, "--rounds")) o->rounds = int_arg(a, v, 1);
    else if (!strcmp(a, "--batch-us")) o->batch_us = real_arg(a, v);
    else if (!strcmp(a, "--case-seconds")) o->case_seconds = real_arg(a, v);
    else if (!strcmp(a, "--csv")) o->csv = v;
    else if (!strcmp(a, "--show")) {
      o->show = 0;
      while (o->show < 3 && strcmp(v, show_names[o->show]) != 0) o->show++;
      if (o->show == 3) die("%s: unknown '%s'", a, v);
    } else if (!strcmp(a, "--seed")) {
      char *end;
      o->seed = strtoul(v, &end, 10);
      if (end == v || *end) die("--seed: bad value '%s'", v);
      o->have_seed = 1;
    } else if (!strcmp(a, "--cpu")) o->pin = int_arg(a, v, 0);
    else if (!strcmp(a, "--fftw-measure-max")) bench_fftw_measure_max = int_arg(a, v, 0);
    else die("%s: unknown option (see --help)", a);
  }
}

/* bench_fft --accuracy: no timing. The libraries are loaded here, not as
 * engines, as all their functions are needed. */
static int accuracy_main(options *o, const int *sizes, int nsizes, const char *date) {
  bench_accuracy_lib libs[MAX_LIBS];
  char err[256];
  FILE *csv = NULL;
  int nlibs = o->nlibs ? o->nlibs : 1, failed;
  if (nsizes == 0) die("no sizes to run: the keyword lists are limited to --min %d .. --max %d", o->min_n, o->max_n);
  if (!o->have_seed) o->seed = 1; /* reproducible by default */
  libs[0].label = "pffft";
  libs[0].lib = NULL;
  for (int i = 0; i < o->nlibs; i++) {
    libs[i].label = o->lib[i].label;
    libs[i].lib = bench_dl_open(o->lib[i].path, err, sizeof err);
    if (!libs[i].lib) die("--lib %s=%s: %s", o->lib[i].label, o->lib[i].path, err);
  }
  if (o->csv && strcmp(o->csv, "-") != 0) {
    csv = fopen(o->csv, "w");
    if (!csv) die("cannot write %s", o->csv);
  } else if (o->csv) {
    csv = stdout;
    printf("#BEGIN-CSV\n");
  }
  if (csv) {
    fprintf(csv, "# pffft-accuracy 1\n# date: %s\n", date);
    bench_platform_describe(csv, "# ");
    fprintf(csv, "# compiler: %s\n# pffft_simd: %s\n", BENCH_COMPILER, simd_name());
    fprintf(csv, "# target_arch: %s\n", BENCH_TARGET_ARCH[0] ? BENCH_TARGET_ARCH : "none");
    fprintf(csv, "# seed: %lu\n# seeds: %d\n# bits_only: %d\n", o->seed, o->seeds, o->bits_only);
    bench_accuracy_columns(csv);
  }
  failed = bench_accuracy(csv, libs, nlibs, sizes, nsizes, o->prec, o->xform, o->seeds, o->seed, o->bits_only);
  if (csv == stdout) printf("#END-CSV\n");
  else if (csv && fclose(csv) != 0) die("cannot write %s", o->csv);
  fflush(stdout);
  fprintf(stderr, "bench_fft --accuracy: %s\n", failed ? "FAILED" : "all checks passed");
  return failed ? 1 : 0;
}

int main(int argc, char **argv) {
  static int sizes[MAX_SIZES];
  options o;
  bench_snap clock_now;
  char date[32] = "unknown";
  time_t now = time(NULL);
  struct tm *utc = gmtime(&now);
  FILE *csv = NULL;
  int nsizes, show, table, marked = 0, unclean = 0, engines = 0;
  long batches = 0, clean = 0;

  parse_options(argc, argv, &o);
  nsizes = expand_sizes(&o, sizes);
  if (utc) strftime(date, sizeof date, "%Y-%m-%dT%H:%M:%SZ", utc);
  if (o.accuracy) return accuracy_main(&o, sizes, nsizes, date);
  load_engines(&o);
  for (int p = 0; p < NPREC; p++) engines += g_reg[p].nsel;
  if (o.list) {
    for (int p = 0; p < NPREC; p++) {
      if (!o.prec[p]) continue;
      printf("%s:", prec_name[p]);
      for (int i = 0; i < g_reg[p].nall; i++) printf(" %s", g_reg[p].all[i].name);
      printf("\n");
    }
    return 0;
  }
  if (engines == 0) die("no engines to run");
  if (nsizes == 0) die("no sizes to run: the keyword lists are limited to --min %d .. --max %d", o.min_n, o.max_n);

  g_counted = bench_platform_init(o.pin);
  for (int i = 0; i < 1000; i++) sample_overhead();
  if (!o.have_seed) {
    bench_snap_begin(&clock_now);
    o.seed = (unsigned long)((clock_now.ns ^ (clock_now.ns >> 32) ^ ((uint64_t)now << 20)) & 0xFFFFFFFFu);
  }
  g_rng = o.seed;

  show = o.show >= 0 ? o.show : g_counted ? SHOW_CYCLES : SHOW_NS;
  if (show == SHOW_CYCLES && !g_counted) {
    fprintf(stderr, "bench_fft: no cycle counter on this system: showing ns\n");
    show = SHOW_NS;
  }
  table = !o.quiet && !(o.csv && !strcmp(o.csv, "-"));
  if (o.csv && strcmp(o.csv, "-") != 0) {
    csv = fopen(o.csv, "w"); /* now, not after minutes of measuring */
    if (!csv) die("cannot write %s", o.csv);
  }
  if (table) write_header(stdout, "", &o, date);

  for (int p = 0; p < NPREC; p++) {
    const registry *reg = &g_reg[p];
    for (int cplx = 0; cplx < 2 && reg->nsel > 0; cplx++) {
      int header = 0;
      if (!o.xform[cplx]) continue;
      for (int i = 0; i < nsizes; i++) {
        char dropped[BENCH_MAX_ENGINES];
        int first = g_nresults;
        memset(dropped, 0, sizeof dropped);
        for (int trial = 0; trial < o.trials; trial++) run_trial(&o, reg, p, cplx, sizes[i], trial, dropped);
        if (g_nresults == first) continue; /* no engine does this size */
        if (!table) {
          fprintf(stderr, "%s %s n=%d\n", prec_name[p], xform_name[cplx], sizes[i]);
          continue;
        }
        if (!header) print_table_header(reg, p, cplx, show, o.trials);
        header = 1;
        unclean |= print_table_row(reg, &g_results[first], g_nresults - first, show);
        marked |= (sizes[i] & (sizes[i] - 1)) != 0;
      }
    }
  }

  if (marked || unclean) printf("\n");
  if (marked) printf("* not a power of two\n");
  if (unclean) printf("? measured, but no batch was counted cleanly (all left the fast cores): see --show ns\n");
  if (o.csv && !csv) {
    printf("#BEGIN-CSV\n");
    write_csv(stdout, &o, date);
    printf("#END-CSV\n");
  } else if (csv) {
    write_csv(csv, &o, date);
    if (fclose(csv) != 0) die("cannot write %s", o.csv);
  }
  fflush(stdout);

  for (int i = 0; i < g_nresults; i++) {
    batches += g_results[i].batches;
    clean += g_results[i].clean;
  }
  if (g_counted && clean < batches)
    fprintf(stderr, "bench_fft: %ld of %ld batches left the fast cores and are not in the cycle statistics\n",
            batches - clean, batches);
  if (g_nresults == 0) {
    fprintf(stderr, "bench_fft: nothing was measured: no selected engine supports the selected sizes\n");
    return 1;
  }
  return g_failed;
}
