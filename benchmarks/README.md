# bench_fft: the FFT benchmark

`bench_fft` measures what one FFT costs, for pffft and other FFT libraries,
side by side in one process, on this computer, an Android phone attached with
adb, or an iPhone attached to a Mac. It answers three questions:

1. How fast is pffft compared with other FFT libraries?
2. How does the cost change with the length n: are powers of two much cheaper
   than other lengths, and where do the caches show?
3. Did a change make pffft faster or slower?

The other two programs here, `bench_pf_mixer_float` and `bench_pf_conv_float`,
measure the mixers (PFDSP) and the direct convolution kernels of `pf_conv`,
one build per architecture option; this file is about bench_fft only.

## Quick start

`bench_fft` is built with the benchmarks, which are on by default:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target bench_fft --config Release
build/benchmarks/bench_fft --sizes 1024,1152 --prec f32 --xform real
```
```
f32 real: CPU cycles per transform, fastest batch
        n        pffft     fftpack     kissfft   pocketfft       green        vdsp
     1024       1680.4      5384.2     12022.8      4799.9      4738.0      1219.2
     1152*      2085.2      8006.6     13214.8      6420.1           -           -
```

`*` marks lengths that are not powers of two; `-` means the engine does not
support the length. Without options bench_fft runs the default sweep: powers
of two from 32 to 2^20 and 28 picked other lengths from 96 to 1000000, both
precisions, real and complex. The picked lengths take one to six passes of
radix 3 and 5 (odd parts 3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 1125, 3125,
9375 and 15625), so that each shows its cost. On an M2, in cache, f32 real
lengths ran per n·log2(n) from 6% slower to 11% faster than the powers of two
around them with one such pass, 5-10% slower with two, 8-14% with three; out
of cache most ran faster than their neighbours (1000000: 15%). `--sizes
nonpow2` gives 3, 5, 9 and 15 times every power of two instead, `--sizes all`
adds the powers of two.
`--help` lists the options, `--list` the engines of the build. With Visual
Studio or Xcode, which build several configurations, the program is
`build/benchmarks/Release/bench_fft`.

The Python tools (Python 3.9 or later, standard library only) build bench_fft
for a target in `build_bench_cache/`, run it there and handle the results:

```sh
python3 benchmarks/run.py -o local.csv                   # this computer, default sweep
python3 benchmarks/run.py --fftw -o local.csv            # also FFTW, built by benchmarks/build_fftw.py
python3 benchmarks/run.py --target android --fftw -o android.csv
python3 benchmarks/run.py --target ios --fftw -o ios.csv
python3 benchmarks/run.py -- --sizes 1024,4096 --prec f32 --trials 3   # after --: bench_fft options
python3 benchmarks/compare.py master .                   # the working tree against master
python3 benchmarks/compare.py c14802c master --target android -o cmp.csv
python3 benchmarks/compare_many.py master a=HEAD~2 b=HEAD~1 c=HEAD --sizes 32,64,1024 -o many.csv
uv run --with matplotlib python3 benchmarks/plot.py local.csv -o local.png
uv run --with matplotlib python3 benchmarks/plot.py local.csv --engines pffft -o pffft.png
uv run --with matplotlib python3 benchmarks/plot.py local.csv android.csv -o devices.png
uv run --with matplotlib python3 benchmarks/plot.py mac-cmp.csv ios-cmp.csv --names 1.1.0,1.3.1 -o cmp.png
```

- **run.py** prints the tables and saves the CSV with `-o`. bench_fft options
  after `--` take the form `--opt value`. run.py and compare.py build for
  this computer's CPU: `TARGET_C_ARCH=native` (-march=native; with MSVC
  /arch:AVX512 or AVX2 where the CPU has them), for pffft, the other engines
  and FFTW alike; the CSV header says `target_arch`. Phones get the default
  of their toolchain.
- **compare.py** builds both revisions as shared libraries with their own
  CMake (any revision from `fde5014` on; `.` is the working tree), loads both
  into one bench_fft and prints the cost of A and B per case, and "B vs A":
  B's speed-up over A (A's cost over B's), the median over 3 launches of 2
  trials each (`--launches`, `--trials`), with its range. Positive means B is
  faster.
  `--fftw` adds FFTW as a reference. `compare.py X X` shows the noise floor.
  `-o` saves the measurements of all launches as one CSV, for plot.py.
  Each revision builds with its own defaults, so v1.0.0 and v1.1.0 get no
  NEON on arm64.
- **compare_many.py** does the same for any number of revisions, each an
  engine `pffft@LABEL` (`LABEL=REV`, or the revision itself), and prints
  each one's speed-up over the first: one commit per engine separates the
  effects of the commits of a branch in one run. A revision given twice is
  loaded twice, which shows the noise floor inside the run.
- **plot.py** draws speed against n, higher is faster: GFLOPS, counting
  5·n·log2(n) flops per complex transform and 2.5·n·log2(n) per real one
  (benchFFT's convention), so that lengths compare fairly. Where cycles were
  counted, the time per transform is its cycles at the median clock of the
  run (the subtitle says which): free of the clock's noise, and still in
  GFLOPS, which compare across devices; elsewhere it is the measured time.
  `--time` takes the measured time throughout, `--per-cycle` plots flops per
  CPU cycle. A panel per precision and transform kind; per engine a solid
  line through the powers of two and a dashed one through the other lengths,
  with a band for the range over trials. vdsp, green and ooura, which do only powers
  of two, get a dotted line there as if padded: the flops of n at the cost of
  the next power of two. Vertical lines mark where 4 arrays of n samples fill
  a cache (a rough guide), if the CSV header gives cache sizes and all files
  give the same. `--engines` picks engines. Several files are overlaid, one
  marker shape each. CSVs of `compare.py -o` give B's speed-up over A per n,
  higher is better, with the range over launches and trials, one color per
  file (device); `--names 1.1.0,1.3.1` names the revisions. It needs
  matplotlib 3.7 or later.

With everything built, the default sweep took 37 s on an Apple M2. With
`--fftw` it took 2 min 21 s there, 4 min 11 s on an Android phone (MT6878) and
2 min 39 s on an iPhone 11 Pro; compare.py took 36 s, 75 s and 49 s. Building
FFTW the first time takes about 2 minutes per target.

## How it measures, and why

- **CPU cycles, not time.** Where the system lets a program count its own
  cycles (macOS, iOS, Linux, Android), costs are in cycles per transform.
  While the data fits in the caches, cycles do not change when the clock is
  throttled; the M2's clock ranged from 0.6 to 3.4 GHz during tests, which is
  what made time-based results noisy. Time (ns) is recorded too, and is all
  there is on Windows or when the CSV header says `counter: none`.
- **The fastest of many short batches.** A batch runs forward+backward pairs
  for about 300 µs (`--batch-us`). Interference only ever adds time, so the
  fastest batch estimates the true cost. Short batches also keep a fanless
  device from heating up.
- **Everything in one process, interleaved.** The engines of a case take turns
  batch by batch in random order, for 100 rounds or about 1 s. Whatever the
  machine does hits all of them alike. Another pffft revision is one more
  engine, loaded as a shared library (`--lib LABEL=PATH`).
- **Trials.** Where the arrays and tables land in memory moves results by up
  to ±0.5% near the L1 size, and by up to ±2.5% at small n, where an engine
  can gain just from being created after another. `--trials T` repeats each
  case with fresh engines, created in a random order, and arrays shifted by
  random amounts; tables show the median. Use `--trials 3` for figures you
  publish: a single trial can spike.
- **Clean batches.** Apple systems cannot pin a thread to a core: bench_fft
  asks for the performance cores, and a batch that touched an efficiency core
  counts for time but not for cycles (bench_fft says how many did). On Linux
  and Android the thread is pinned to the fastest core (`--cpu N` picks
  another), and cycles are counted only when it is pinned. Linux allows
  counting only with `kernel.perf_event_paranoid` at 2 or lower (Debian and
  Ubuntu set more; `sudo sysctl kernel.perf_event_paranoid=2`) and, in a
  container, a seccomp profile that allows `perf_event_open`.
- **Costs that cycles do not show.** bench_fft records the clock (cycles
  per ns) of every clean batch and compares it with that of the batches of
  other engines just before and after it, so that the system changing the
  clock does not count. An engine whose median ratio is more than 2% below 1
  is charged its cycles at the others' clock, which makes it compare by time;
  the CSV columns `ghz` and `clock` show this. Other engines stay within 0.04%
  of 1 on an M2 (90th percentile), 1.5% at most. Two causes were found:
  - On Apple chips the clock drops while the AMX matrix unit works: for vDSP
    f32 FFTs from n = 256 or 512, by 7.7% on the M2 and up to 25% on the
    iPhone 11 Pro; by about 23% for a small sgemm. It stays down as long as
    AMX works (as low after 20 ms of vDSP as after 300 µs) and recovers some
    40 µs after; unlike FFTW's planning it is no one-time cost. AVX-512 slows
    Intel cores the same way.
  - Linux counts user-space cycles only. PocketFFT allocates a work array in
    every transform; at n from about 5000 the page faults that costs, in the
    kernel, made its batches look 2-41% slower in time than in cycles on an
    Android phone (MT6878).

  On the iPhone 1.2% of the cases of other engines, all at n of 262144 and
  more, were charged 2-6% for reasons not found.
- **Self-check.** Before timing, every engine must transform random data
  correctly; one that fails is dropped with a message, and bench_fft exits
  with a nonzero status.

With `compare.py master master` (3 launches of 2 trials), two copies of the
same library differed by at most 0.49% for n ≥ 128 on the M2 (0.06% on
average) and by at most 0.93% below. In earlier runs with 5 trials in one
launch it was 0.36% with 8 busy processes running, 0.62% on the iPhone and
0.75% on the Android phone (1.31% out of cache). Under the same 8 busy
processes, three runs of the old benchmark spread by 11% (pffft, f32 real,
n = 1024, MFLOPS), three runs of bench_fft by 0.001% (cycles). A real change
stands out: from c14802c to master, f32 real n = 1024 gets 6.0 to 6.3% cheaper
on the M2 and the Android phone, 9.4% on the iPhone.

What it cannot do:

- **Out of cache**, memory takes the same time at any clock, so cycles depend
  on the clock there. The CSV has both; for large n check the ns too.
- **Small n.** Up to a few hundred, single trials of identical code differ by
  up to about ±2.5%, and medians by up to 1%. Trust a difference there only
  when it is larger than the range of `compare.py X X` at that n.
- **Alone, an engine that slows the clock** has nothing to be compared with:
  run it with at least one other engine (`--engines vdsp` alone counts its
  cycles at its own clock; the `clock` column is then empty).

## Engines

| Engine | Library | Precisions | Lengths | Notes |
|---|---|---|---|---|
| `pffft` | this tree, `pffft_transform_ordered` | f32, f64 | 2^a·3^b·5^c; real: multiple of 32, complex: of 16 | |
| `pffft-ip` | this tree, `pffft_transform_ordered` in place | f32, f64 | as pffft | one array and the work array, as Ooura's FFT and FFTPACK run; only with `--engines` |
| `pffft-un` | this tree, `pffft_transform` | f32, f64 | as pffft | the spectrum in pffft's internal order: no reordering; only with `--engines` |
| `pffft@L` | a pffft shared library, `--lib L=PATH`, `pffft_transform_ordered` | those it has | as pffft | what compare.py measures |
| `fftpack` | `src/fftpack.c` | f32, f64 | any | |
| `fftw-est`, `fftw-meas` | FFTW 3.3.10 from `FFTW3_ROOT` | f32, f64 | any; `fftw-meas` up to `--fftw-measure-max` (65536) | FFTW_ESTIMATE plans from a cost model, at once, as pffft sets up; FFTW_MEASURE plans by timing candidate algorithms first (about 5 s at n = 262144 on an M2, outside the timing), an advantage no other engine takes |
| `kissfft` | submodule `kissfft` | f32 | any (real: even) | |
| `pocketfft` | submodule `pocketfft` | f32, f64 | any | |
| `green` | submodule `greenffts` | f32 | powers of two | not on Windows, where it crashes |
| `ooura` | submodule `ooura`, Takuya Ooura's split-radix `fftsg.c` | f32, f64 | powers of two | |
| `vdsp` | Apple Accelerate | f32, f64 | powers of two | Apple only; f32 from n = 256 or 512 slows the clock while it runs, and is charged for it (see above) |
| `mkl` | Intel oneMKL | f32, f64 | any | not yet run with bench_fft |

Every library is built from source, except vDSP, which is part of the
operating system, and MKL, a vendor binary. FFTW is never taken from the
system: a package manager may build it without SIMD (Homebrew's on Apple
Silicon has none), which would understate it. `benchmarks/build_fftw.py` builds it
with NEON, or SSE2, AVX and AVX2 on x86-64, plus AVX-512 where the CPU has it
(pffft has no AVX-512 code); `run.py --fftw` and `compare.py --fftw` call it
for their target. For a CMake build:

```sh
python3 benchmarks/build_fftw.py --prefix build_bench_cache/fftw-local
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPFFFT_USE_BENCH_FFTW=ON \
      -DFFTW3_ROOT=$PWD/build_bench_cache/fftw-local
```

## Devices

**Android:** `adb` on the PATH, the phone attached with USB debugging allowed,
and an Android NDK: `$ANDROID_NDK`, `$ANDROID_NDK_HOME`, or the newest in the
`ndk/` directory of `$ANDROID_HOME`, `~/Library/Android/sdk` or
`~/Android/Sdk`. `--target android:SERIAL` picks one of several phones. A run
copies its files to `/data/local/tmp/pffft-bench/` and deletes the directory
at the end, so run one at a time per phone. `--fftw` needs an arm64-v8a phone.
Cycles need perf_event access from the adb shell (the CSV header then says
`counter: perf_event`); the phone tested here (Android 16) allowed it as it
came. With `counter: none` the results are in ns.

**iOS:** a Mac with Xcode, and the iPhone attached, paired and in Developer
Mode (`xcrun devicectl list devices` shows it). `--target ios:DEVICE` takes a
UDID, identifier or name. Signing needs a development certificate in the
keychain and the wildcard profile "iOS Team Provisioning Profile: \*" that
lists the phone; Xcode makes it when it signs any app for that phone with
automatic signing. Keep the phone unlocked while it runs. The app,
`com.example.pffft.bench`, stays installed; remove it with
`xcrun devicectl device uninstall app --device UDID com.example.pffft.bench`.

## Other CPUs, on GitHub's runners

`.github/workflows/bench.yml` runs the sweep, and optionally compare.py, on
GitHub-hosted runners: by default x86_64 and arm64 Linux, Apple-silicon and
Intel macOS, x86_64 and arm64 Windows. It runs only when started by hand,
from the Actions tab or with gh, once the workflow is on the default branch:

```sh
gh workflow run bench.yml --ref my-branch -f runners=ubuntu-26.04-arm,macos-26-intel -f compare="v1.3.0 HEAD"
gh run watch
gh run download <run id> --dir bench-results
uv run --with matplotlib python3 benchmarks/plot.py bench-results/bench-ubuntu-26.04-arm/ubuntu-26.04-arm.csv -o arm.png
```

Each runner uploads an artifact `bench-<runner>` with the CSVs and what the
system says about its CPU (`<runner>-system.txt`); the workflow draws no
plots. Runners are virtual machines that others share, so expect more noise
than at home; a Linux VM may expose no cycle counter, and Windows has none
(the CSV header says `counter: none`; plots then take the measured time). On Windows
build_fftw.py uses FFTW's CMake build with MSVC, which has SSE2, AVX and
AVX2 but no AVX-512 and no NEON: FFTW is left out on Windows ARM64.

`-f compare_args=...` passes options to compare.py, `-f many="a=master
b=HEAD"` (and `many_args`) runs compare_many.py as well.

A runner label lands on one of several CPUs: ubuntu-26.04 and windows-2025
on Intel Xeon and AMD EPYC machines of several generations, some of which
come up once in 30 runs. `benchmarks/ci_dispatch.py` runs the workflow
until each CPU of a table has done its part: every run gets the table as
the input `by_cpu`, does the part of the CPU it lands on, and ends after a
few seconds on any other CPU. Give a CPU everything it should measure in
one entry (sweep, compares, compare_many); the script's docstring shows the
table:

```sh
python3 benchmarks/ci_dispatch.py table.json ci-out --ref my-branch
```

## CSV

`bench_fft --csv FILE` writes the results; with `--csv -` they go to stdout
between the lines `#BEGIN-CSV` and `#END-CSV`. `benchmarks/results.py` reads them.

```
# pffft-bench 1
# date: 2026-11-11T19:57:30Z
# os: macOS 26.5.2
...
engine,prec,xform,n,trial,pairs,batches,clean,ns,ns_median,cycles,cycles_median,instructions,mflops,ghz,clock
pffft,f32,real,1024,0,304,100,100,480.127,493.970,1681.586,1682.602,10811.535,53319.3,3.5038,0.9998
```

The `#` lines describe the machine and the run (date in UTC, os, cpu,
counter, the known cache sizes, compiler, pffft_simd, trials, seed, ...).
There is one row per engine, precision (`f32`, `f64`), transform kind
(`real`, `cplx`), `n` and trial:

| Column | Meaning |
|---|---|
| `pairs` | forward+backward pairs per batch |
| `batches`, `clean` | batches run; batches counted in cycles (all on the fast cores) |
| `ns`, `ns_median` | time per transform: fastest batch, median batch |
| `cycles`, `cycles_median` | the same in CPU cycles, over clean batches; empty without a counter or a clean batch |
| `instructions` | per transform, in the fastest clean batch |
| `mflops` | (5 complex, 2.5 real) · n · log2(n) / `ns` · 1000 |
| `ghz` | the clock of the clean batches (cycles per ns, median) |
| `clock` | that clock relative to the batches of other engines around them, median; below 0.98 the cycles were charged for it (see above); empty with one engine or without a counter |

## CMake options and tests

| Option | Default | Effect |
|---|---|---|
| `PFFFT_BUILD_BENCHMARKS` | ON if pffft is the top-level project | builds bench_fft, and the mixer and convolution benchmarks in `benchmarks/` |
| `PFFFT_USE_TYPE_FLOAT`, `PFFFT_USE_TYPE_DOUBLE` | ON | the precisions bench_fft has |
| `PFFFT_USE_FFTPACK` | ON | `fftpack` |
| `PFFFT_USE_BENCH_KISS`, `_POCKET`, `_GREEN`, `_OOURA` | ON | `kissfft`, `pocketfft`, `green`, `ooura`, if the submodule is there (`git submodule update --init`) |
| `PFFFT_USE_BENCH_VDSP` | ON (Apple only) | `vdsp` |
| `PFFFT_USE_BENCH_FFTW` | OFF | `fftw-est`, `fftw-meas`; needs `FFTW3_ROOT`, a prefix made by `build_fftw.py` |
| `PFFFT_USE_BENCH_MKL` | OFF | `mkl`; MKL's libraries and `<mkl/mkl_dfti.h>` must be on the compiler's paths |

`ctest --test-dir build -R bench_fft` (with Visual Studio or Xcode add
`-C Release`) runs the two smoke tests, `bench_fft_smoke` (table, n = 64, 96,
128) and `bench_fft_csv` (CSV, n = 64, 96): every engine's self-check and a
short timing, under a second each.
