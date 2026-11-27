# Code layout

How the speed of pffft depends on where its code lies in memory, how to tell
that from the effect of a change, and what is known about it per CPU. Which
check a loss needs is in AGENTS.md (Thresholds that depend on the CPU, "Code
layout").

## Why

The same instructions can run at different speeds at different addresses.
The decoders and op caches of a CPU work on 32- and 64-byte blocks of code,
so where a loop starts in such a block decides how many cycles each of its
iterations takes to fetch. The branch predictors index their tables by
address bits. A change to pffft moves all the code after it, so a comparison
of two revisions mixes two things: what the change does, and where it put the
code.

For short transforms the second part can dominate. pffft's master, built
eight times with its code moved by a pad (below), and nothing else changed:

| machine | compiler | band: max over the cases (case) | median over the cases |
|---|---|---|---|
| AMD EPYC 7763 (Zen 3) | MSVC 19.52 | 16% (f64 real 64) | 1.4% |
| AMD EPYC 7763 (Zen 3) | GCC 15.2 | 12% (f32 real 128) | 0.5% |
| Arm Neoverse N2 | GCC 15.2 | 6% (f32 real 64) | 0.2% |
| Intel Core i7-8700B | Clang 21 | 5% (f32 complex 32) | 0.5% |
| AMD EPYC 9V45 (Zen 4) | GCC 15.2 | 4% (f64 real 256) | 0.4% |
| Microsoft Cobalt 100 | MSVC 19.52 | 4% (f32 real 64) | 0.2% |
| Raspberry Pi 2 (Cortex-A7) | GCC 8 | 2% (f64 complex 2400) | 0.3% |
| AMD EPYC 9V74 | GCC 15.2 | 2% (f64 complex 96) | 0.2% |
| Apple M2 | Clang 21 | 2% (f64 real 64) | 0.0% |
| Intel Xeon 6975P-C | GCC 15.2 | 1% (f64 real 512) | 0.2% |
| Intel Xeon Platinum 8573C | MSVC 19.52 | 1% (f64 real 128) | 0.2% |
| iPhone 11 Pro (A13) | Clang 21 | 0.2% (f32 complex 512) | 0.0% |

The band of a case is how far the speed of the same code moves over eight
layouts (max - min, counting at each layout only what two copies of the
library agree on), for n = 32 to 2400, both precisions, real and complex
(November 2026, `build_bench_cache/layout/table.md`). Longer transforms,
which spend their time on memory, move far less.

So on the EPYC 7763, the N2 and the i7, a comparison of two builds says little
about a short transform below 5 to 15%: a loss of #140 on the 7763 of 9.4% at
f64 real 128 lies inside master's own band of 7%, and the "-19.6% at f32 real
64" of #135 inside 8%. On the recent Intel and Apple cores a single comparison
is fine. And the other way round: a loss that one calls layout without
evidence may be real. Three of them were (below).

## The tool

```sh
python3 benchmarks/compare_many.py master . --pads spread --sizes 32,64,96,128,160,256,480,512,1024 -o cmp.csv
python3 benchmarks/compare_many.py master m2=master --pads spread ...    # the noise at matched layouts
python3 benchmarks/compare_many.py --report cmp.csv                      # the report of a saved run
python3 benchmarks/layout.py diff cmp.csv.libs/master+0.dylib cmp.csv.libs/.+0.dylib   # what the change changed
```

`compare_many.py --pads LIST` builds every revision once per pad and loads
them all into one bench_fft, which measures them batch by batch in random
order (as compare.py does for two). `--pads spread` is 0, 176, 400, 720, 928,
1264, 1600 and 1888 bytes: every multiple of 16 modulo 128 once, and the
bits above (128 to 1024) varied too. Pads of 4096 or more would add nothing,
as the loader puts a library at an address that is a multiple of the page.
Two revisions times eight pads are the 16 libraries that bench_fft loads at
most.

**The pad** comes from `benchmarks/layout_pad.cmake`, which compare_many.py
gives to the revision's own CMake as its project include, so that it works on
any revision without touching its sources:

- GCC and Clang: a header forced into the first C file of the library, whose
  top-level asm puts the pad's bytes at the start of its code section.
  Everything after it, in this object and those linked after it, moves by the
  pad.
- MSVC: one more object, first in the library, with a function of that many
  bytes of no-ops in section `.text$a`, which the linker puts in front of all
  objects' code (`.text$mn`). x64 functions start with a 2-byte instruction
  (66 90), so that the function is pad - 3 one-byte no-ops and a return.

Both had to be code inside the code section. A pad of data went to a section
of its own with MSVC, and nothing moved; a pad as an object of its own was
rounded by GCC to the alignment of the next object's code section (32 or 64
bytes), so that the kernels never moved within the 64-byte blocks.

With a pad of 0 the library is the plain build. The objects of pffft are the
same at every pad: `layout.py diff` finds the functions equal, but for the
addresses in the jump tables that MSVC keeps in the code section after some
functions, which it cannot tell from code.
The run reads where the functions went (nm, or the map that MSVC's linker
writes) and prints how far each pad moved them, warns if a pad moved nothing
or repeated another pad's layout, and says per revision how many offsets
modulo 64 the kernels of the short transforms took over the pads, and
whether the two revisions' kernels took the same ones.

**What the pads sample** depends on the compiler:

- MSVC and Clang move every function by the pad: each kernel takes four
  offsets modulo 64.
- GCC aligns functions and loops to 32 or 64 bytes for many CPUs where that
  skips few bytes (`.p2align 5,,15`), and so takes up part of each pad. On the
  N2 most kernels took two offsets, with GCC 15.2 on x86 three quarters of
  them only one: there the pads hardly sample the layouts.
  `--cmake '*=-DPFFFT_LAYOUT_ALIGN=16'` builds every label with functions and
  loops aligned to 16 bytes exactly, where the pads move all kernels on the
  N2; on x86 GCC 15.2 still aligns some loops to 64 bytes that way. It is a
  build of its own, so give it to every label.
- Two different revisions cannot be put at the same layout: the code that
  changed moves the kernels after it by other amounts (on Arm by multiples of
  4 bytes, which pads of 16 do not undo). The pads sample each revision's
  layouts, and the band says how much layout moves each one.

**The report** gives per case each revision's speed-up over the first one at
the same pad (median over the trials), the range over the pads, the verdict,
and each revision's band:

- `slower` / `faster`: at every pad;
- `slower*` / `faster*`: at every pad, and by more than both bands;
- `layout`: the sign depends on the pad.

Call a loss **real** if it is `slower*`, larger than the noise of the A/A run
at that case (the largest difference between two copies of master at the same
pad), and the kernels were sampled: on the machine's compiler they took two or
more offsets, the same ones in both revisions. Call it **layout** if the
verdict says so. Anything else is undetermined: more pads, more launches, or
another compiler (MSVC or Clang on x86).

**Noise beyond the address.** On an EPYC 9V45, a virtual machine, the
library file that the build wrote ran 14% slower than its byte-identical copy
at the same pad in every trial of three launches, and in another run another
one 40%, also with the load order shuffled. With every library loaded from a
new copy of its file in each launch, which `--pads` now does (`--shuffle`),
no library was slow throughout and the A/A noise there fell to 3%. Why is
open; on such a machine, look at the A/A copies before trusting a
difference.

**Counters.** On Linux, `--events LIST` (`bench_fft --events`) counts up to
four more hardware events per transform with the cycles, and the report
prints them per pad next to the instructions: perf's generic names
(`branch-misses`, `L1-icache-load-misses`, ...), a PMU's own events
(`stall_frontend` on Arm), the fields of an event from the CPU vendor's list
(`cpu/event=0x79,umask=0x08/`: IDQ.DSB_UOPS on Intel), or raw codes. A layout
case keeps its instruction count.

`-o FILE` also saves the functions' addresses (`FILE.layout.json`) and the
libraries (`FILE.libs/`, with MSVC's map), so that `layout.py loops LIB` (the
loop heads modulo 64 and 32) and `layout.py diff A B` (which functions two
builds differ in, with the addresses taken out) also work on another
computer, e.g. on the runners' artifacts.

## What it found

The first round (November 2026) checked the losses of #133 to #140 that had
been put down to layout:

| case | single comparison | matched layouts | verdict |
|---|---|---|---|
| #140, N2, f32 real 64 | -4.7%, same instructions | -4.3 to +4.4% (band 5.6%) | layout |
| #140, EPYC 7763 MSVC, f64 real 128 | -9.4% | -13.1 to +4.8% | layout |
| #135, EPYC 7763 MSVC, f32 real 64 | -19.6% | master's band there: 8% | layout |
| #139, EPYC 7763 MSVC, f32 real 96 | -8.9% | -13.7 to +5.0% | layout |
| #140, i7-8700B, f32 real 96 / 128 / 160 | -7.5 / -8.7 / -5.7% | -7.4..-9.0 / -5.5..-9.1 / -4.4..-6.0% | **real** |
| #140, i7-8700B, f64 real 96 / 128 / 160 | -7.4 / -6.4 / -5.5% | -2.7 to -7.9% (64 to 160) | **real** |
| #140, Raspberry Pi 2, f64 real 32 / 64 / 128 | | -1.7..-2.5 / -2.4..-4.4 / -1.8..-3.3% | **real** |
| #133, Apple M2, f32 real 15360 | -2.9% | -2.9 to -3.8% at 16 layouts, bands 0.3% | **real** |

The real ones have causes in the code: #140 changed the real drivers
(`rfftb1_ps` and `rfftf1_ps`, `pffft_transform_internal`), #133 the radix-3
and radix-5 passes that 15360 runs (`layout.py diff`).

**What decides the band.** With MSVC and Clang, the pads 0/1600, 176/1264,
400/720 and 928/1888 put every function at the same offset modulo 64. On the
EPYC 7763 the speed agrees within those pairs (f32 real 64: band 9%, within
the pairs 0.6%; f32 real 96: 16% and 3%), on the i7 too: there the offset
modulo 64 decides, which is the op cache's line. On the N2 and the Cobalt it
does not (N2 f32 real 64: band 6%, within the pairs 4%): bits above 64 bytes
matter, the branch predictors perhaps. The N2's counters show the same
instructions at every pad, no branch mispredictions, and no wasted
speculation; the slow layouts issue fewer instructions per cycle without the
front-end stalls that the PMU counts. On the Xeon 6975P-C, Intel's events
(IDQ.DSB_UOPS against IDQ.MITE_UOPS) show the short transforms delivered from
the decoded-uop cache at every pad, which is why its band is 1%.

## What to do with it

- **Losses lists** (pr-policy): on the machines with a wide band (EPYC 7763,
  N2, i7, 9V45), give each loss of a short transform the matched-layout
  verdict, measured with MSVC or Clang where possible; on the others a single
  comparison does.
- **Pinning the alignment** (ticket 320): on the EPYC 7763 with MSVC, one of
  the four offsets modulo 64 that the pads tried is the fastest or within 2%
  of it at every length, and master as built lies up to 10% from it (f32
  real 96; f64 real 96 6%, f64 real 128 5%):

  | code moved by (modulo 64) | mean loss against each case's best | worst case |
  |---|---|---|
  | 0 (master as built) | 1.0% | 10.5% (f32 real 96) |
  | 16 | 0.4% | 2.3% |
  | 32 | 1.2% | 9.7% |
  | 48 | 2.1% | 19.2% (f64 real 64) |

  Aligning the hot kernels to 64 bytes would make that layout stay where it
  is when unrelated code changes (GCC and Clang: an attribute or
  `-falign-functions=64`; MSVC has no option for it, but clang-cl does).
  Which kernel carries the effect is the next thing to find. On the N2 and the
  Cobalt pinning to 64 bytes would not help.
