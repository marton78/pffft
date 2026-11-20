# pffft benchmark results

Measured with [bench_fft](../benchmarks/README.md) in November 2026: pffft
1.3.1 against the FFT libraries it is usually compared with, on eight
machines, and against pffft 1.1.0, which is as fast as Julien Pommier's
original, on three. Every chart has a panel per precision (f32, f64) and
transform kind (real, complex) and covers 44 lengths from 32 to 2^20: the
powers of two and 28 other lengths with factors 3 and 5.

## Summary

pffft's speed relative to each library: the geometric mean over the 176
cases of a machine (2 precisions × 2 kinds × 44 lengths) of the other
library's time over pffft's. Above 1, pffft is faster; "-" means the library
is not available there.

| Machine | FFTW (measure) | FFTW (estimate) | vDSP | PocketFFT | KISS FFT | FFTPACK | Green FFT |
|---|---|---|---|---|---|---|---|
| Apple M2, macOS | 0.86 | 1.17 | 0.80 | 1.80 | 4.51 | 2.26 | 2.31 |
| iPhone 11 Pro (A13), iOS | 0.93 | 1.26 | 0.93 | 1.96 | 4.65 | 2.18 | 2.57 |
| MediaTek MT6878, Android | 1.09 | 1.37 | - | 1.96 | 4.26 | 2.17 | 2.02 |
| Arm Neoverse N2, Linux | 1.17 | 1.60 | - | 1.81 | 4.55 | 2.63 | 4.23 |
| Microsoft Cobalt 100 (Neoverse N2), Windows | - | - | - | 2.27 | 3.97 | 2.42 | - |
| Intel Xeon 6975P-C, Linux | 0.55 | 0.84 | - | 1.51 | 5.59 | 2.12 | 3.34 |
| AMD EPYC 7763, Windows | 0.70 | 1.60 | - | 2.93 | 5.09 | 2.81 | - |
| Intel Core i7-8700B, macOS (VM) | 0.76 | 1.28 | 0.63 | 1.97 | 4.26 | 2.16 | 2.43 |

KISS FFT and Green FFT are measured in single precision only.

- **Portable libraries.** pffft was faster than PocketFFT, KISS FFT, FFTPACK
  and Green FFT in 1397 of the 1408 cases. The 11 exceptions are complex
  transforms on the two x86 machines with macOS and Linux, 10 of them at n
  of 72000 and more, out of cache; there pffft was at most 8.7% slower.
- **FFTW** plans a transform in one of several modes; two are measured here.
  `FFTW_ESTIMATE` plans at once, as pffft and the other libraries do.
  `FFTW_MEASURE` times candidate algorithms on the machine and keeps the
  fastest, which takes up to a few seconds per length. The benchmark leaves
  that planning time out, so the comparison favours `FFTW_MEASURE`: in an
  application its cost pays off only over many transforms of one length.
  - *Real transforms.* pffft was faster than `FFTW_ESTIMATE` on every
    machine (geometric means from 1.02 on the Xeon 6975P-C to 1.75 on the
    Neoverse N2). On arm64 it was also faster than `FFTW_MEASURE` on the
    iPhone, the MediaTek phone and the Neoverse N2 (1.10, 1.27 and 1.39)
    and level with it on the M2 (1.00).
  - *Complex transforms.* pffft was faster than `FFTW_ESTIMATE` everywhere
    but on the Xeon 6975P-C (1.08 to 1.53; 0.69 there). `FFTW_MEASURE` was
    faster on Apple's chips (pffft 0.73 and 0.78) and level on the MediaTek
    phone and the Neoverse N2 (0.93 and 0.99).
  - *x86.* `FFTW_MEASURE` was faster on the x86 machines (pffft 0.55 to
    0.76 over all cases): FFTW has AVX2 and AVX-512 code, while pffft's
    single precision uses 4-wide SSE there (its double precision uses AVX).
- **vDSP**, Apple's FFT, was faster in single precision. On Apple's chips it
  hands single-precision FFTs from n = 256 or 512 to the AMX matrix unit, a
  coprocessor that the cores of a cluster share: one thread at a time
  benefits, while pffft's NEON code runs on every core. While AMX works, the
  chip lowers its clock, by 7.7% on the M2 and up to 25% on the iPhone 11
  Pro, for as long as AMX is busy; bench_fft charges vDSP for that (see
  [How it was measured](#how-it-was-measured)). In double precision, where
  vDSP did not lower the clock, pffft was level or faster on Apple's chips
  (geometric means 0.99 to 1.32 per kind) and slower on the Intel Mac
  (0.74 real, 0.64 complex).

## Performance improvements since Pommier's pffft

![pffft 1.3.1 vs 1.1.0 on Apple M2, iPhone 11 Pro and MediaTek MT6878](benchmarks/pffft-1.3.1-vs-1.1.0.png)

pffft 1.1.0 (here commit c14802c) predates the performance work of this
fork. Its single-precision FFT is Julien Pommier's original code, so its
speed is the speed of Pommier's pffft (which has no double precision). Since
then pffft has gained, among others, fused multiply-adds (NEON and x86 FMA),
NEON double precision with native intrinsics instead of an emulation, faster
NEON transposes, fast paths in the radix-3 and radix-4 real butterflies, a
radix-8 pass for complex n = 32, and code layout that keeps the radix-3 and
radix-5 real kernels out of line.

Both versions run in one process and take turns, batch by batch
(benchmarks/compare.py). Averaged over the lengths, 1.3.1 was 11 to 13%
faster below n = 4096 (11.1% on the M2, 12.6% on the iPhone, 12.0% on the
MediaTek phone) and 5 to 8% faster from 4096 up (5.5%, 8.1% and 5.3%), up
to 24% at single lengths. In cache (n up to 16384) it was faster at every
length, by 2% at least; out of cache the gain shrinks, and on the MediaTek
phone it is within the noise at many lengths from 2^16 up.

## Machines

How to read the charts: each shows one machine, with the speed in GFLOPS on
a log scale (higher is faster) over the transform length n.

- **Solid lines** are the powers of two.
- **Dashed lines** are the other lengths, which the library transforms as
  they are.
- **Dotted lines** are the other lengths for a library that only does powers
  of two (vDSP, Green FFT): the flops of n in the time of the next power of
  two, as if the input were padded with zeros to it.
- **Bands** around a line show the range over the 3 trials.
- **Grey vertical lines** mark the lengths at which four arrays of n
  samples, about the working set of a transform, fill the L1 and L2 caches
  (and L3, where the system reports one).

Where cycles were counted, the subtitle gives the clock at which they were
converted to time.

### Apple M2, macOS
![FFT speed on Apple M2](benchmarks/apple-m2-macos.png)

### iPhone 11 Pro (A13), iOS
![FFT speed on iPhone 11 Pro](benchmarks/iphone-11-pro-ios.png)

### MediaTek MT6878, Android
![FFT speed on MediaTek MT6878](benchmarks/mediatek-mt6878-android.png)

### Arm Neoverse N2, Linux (GitHub runner)
![FFT speed on Arm Neoverse N2](benchmarks/neoverse-n2-linux.png)

### Microsoft Cobalt 100, Windows (GitHub runner)
FFTW is left out: its CMake build, the one Windows can use, has no NEON.
Green FFT crashes on Windows.

![FFT speed on Microsoft Cobalt 100](benchmarks/cobalt-100-windows.png)

### Intel Xeon 6975P-C, Linux (GitHub runner)
![FFT speed on Intel Xeon 6975P-C](benchmarks/xeon-6975p-c-linux.png)

### AMD EPYC 7763, Windows (GitHub runner)
FFTW is built with its CMake build here, with AVX2 (the CPU has no AVX-512).

![FFT speed on AMD EPYC 7763](benchmarks/epyc-7763-windows.png)

### Intel Core i7-8700B, macOS (GitHub runner, VM)
![FFT speed on Intel Core i7-8700B](benchmarks/core-i7-8700b-macos.png)

## How it was measured

- **Speed** is in GFLOPS, counting 5·n·log2(n) flops per complex transform
  and 2.5·n·log2(n) per real one (benchFFT's convention), for one transform:
  half a forward and backward pair, out of place, ordered output
  (pffft_transform_ordered). It is the fastest of many short batches,
  median of 3 trials.
- **Time.** Where the system lets a program count its CPU cycles (the M2,
  the iPhone, the MediaTek phone, the Neoverse N2 and the Xeon 6975P-C), a
  transform's time is its cycles at the median clock of the run, which
  keeps the clock's noise out. Elsewhere (Windows, macOS virtual machines)
  it is the measured time.
- **Clock.** Cycles are fair only if every library runs at the same clock.
  bench_fft compares the clock of each batch with that of the neighbouring
  batches of the other libraries; a library that lowers the clock by more
  than 2%, such as vDSP through AMX, is charged for the cycles it costs the
  others: its count is divided by the ratio of the clocks.
  [benchmarks/README.md](../benchmarks/README.md) explains the method and
  its limits.
- **Builds.** Every library was compiled from source with the same flags,
  the best for the CPU: -march=native with Clang and GCC (the M2, the Linux
  and macOS runners), /arch:AVX2 with MSVC on the EPYC, MSVC's default on
  the Cobalt 100. The phones were built with their toolchain's defaults
  (NEON). FFTW 3.3.10 was built with its SIMD code (NEON, AVX, AVX2,
  AVX-512 where the CPU has it). vDSP is the system's Accelerate.
- **Versions.** 1.3.1 is pffft as of commit e392851; 1.1.0 is commit
  c14802c. The phones and the M2 ran the benchmark locally; the other
  machines are GitHub-hosted runners (shared virtual machines), measured by
  .github/workflows/bench.yml.

To repeat a measurement:

```sh
python3 benchmarks/run.py --fftw -o local.csv            # this computer
python3 benchmarks/compare.py c14802c master -o cmp.csv  # two revisions
uv run --with matplotlib python3 benchmarks/plot.py local.csv -o local.png
```
