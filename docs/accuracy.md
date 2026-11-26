# Accuracy checks

How pffft's accuracy is checked during development: the tool, what it
measures, and how its tolerances were set. Which check a change needs, and on
which platform, is in AGENTS.md (Accuracy checks).

A change that alters the arithmetic (FMA placement, another radix order, other
twiddles) need not give master's bits, but it must not be less accurate; a
change that only moves data must give master's bits. `accuracy.py` checks
either:

```sh
python3 benchmarks/accuracy.py master .                  # this computer
python3 benchmarks/accuracy.py master . --bits           # the same bits as master (fast)
python3 benchmarks/accuracy.py master . --target rosetta # x86 SSE/AVX/FMA under Rosetta (Apple arm64)
python3 benchmarks/accuracy.py master . -v               # every length, not only the flagged ones
```

It builds each revision twice as a shared library, as configured by default
and with `PFFFT_SPLIT_MAX=1024` (every array of more than 1 KiB split into
blocks), and runs `bench_fft --accuracy` on all four with the same random
inputs (uniform in [-1, 1), rounded to the precision): by default the powers
of two and 3, 5, 9 and 15 times them from 16 to 2^22, and the picked lengths,
with as many inputs per length as make 16384 points (one from n = 16384 up).
Per library and case it measures, against an exact transform,

- `fwd`: the relative RMS error of `pffft_transform_ordered` forward on a
  random signal, pooled over the inputs: sqrt(Σ|y − ref|² / Σ|ref|²);
- `bwd`: the same backward, on a random spectrum (any values in the ordered
  layout are a valid one), against the unscaled inverse transform;
- `rt`: the round trip, backward(forward(x)) / n against x, shown but not
  judged: the backward pass uses the forward pass's twiddles conjugated, so
  their errors largely cancel in it;

and it checks bit for bit, within each build, the ordered transform in place
against out of place (both directions; also the unordered forward), and the
ordered transform against `pffft_transform` + `pffft_zreorder` (both
directions). It fails if

- one of those checks fails, or an error of one input exceeds 2e-6 (f32) or
  5e-15 (f64), about ten times what pffft gives;
- the split build of a revision differs from its default build in one bit;
- B's error over A's exceeds 1.05 at one case (`--case-tol`), or its
  geometric mean over the lengths of a precision and kind exceeds 1.01
  (`--class-tol`), for `fwd` or `bwd`;
- with `--bits`, B differs from A in one bit (then no errors are measured).

## The reference

The exact transform is a mixed-radix FFT in double-double arithmetic (about
106 bits; accuracy.c), with twiddles from an exact integer reduction of the
angle and sine and cosine from their series. Against a DFT computed with
mpmath at 60 digits its error is about 1e-32, far below pffft's in both
precisions. Another FFT library in double precision would not do as the
reference for f64: its error is the size of pffft's and correlated with it
(PocketFFT and pffft both descend from FFTPACK), so the distance between the
two rewards resembling it rather than being accurate.

## Calibration

Pooled over this many points the ratios are stable. On the M2, 3c480ca
(explicit FMA in the radix-3/5 codelets) against its parent changed the bits
of half the cases, but no ratio by more than 0.7% (and all class means were
slightly better). #128 (complex radix-8 last pass) left f64 complex 2.9% more
accurate on average, and f32 complex 2400 5.4% less accurate. #133 (real
radix-8 and complex radix-6 last passes) cost no accuracy in real transforms
but made complex ones with a radix-6 pass (3 · 2^k and 9 · 2^k) up to 34%
less accurate (f64) and 21% (f32), through the error of the twiddles, which
the new factor order exposed. The gate flags both. `accuracy.py X X` must pass
with every case bit for bit equal. A full run took 1 minute 50 on the M2 and
3 minutes 30 under Rosetta, builds included, a `--bits` run 30 seconds.

## bench_fft --accuracy

`bench_fft --accuracy` alone checks the pffft built into it (or `--lib`
libraries): `--seeds K` (at least K · 4096 points per case), `--bits-only`
(no errors, only the checks and the hashes), `--sizes`, `--prec`, `--xform`
as for timing, and `--csv` for one row per library and case:

```
engine,prec,xform,n,inputs,fwd,bwd,rt,fwd_max,bwd_max,rt_max,checks,hash_fwd,hash_bwd,hash_un
pffft,f32,real,128,128,9.0162e-08,8.9543e-08,1.1926e-07,1.0764e-07,1.0468e-07,1.4266e-07,,04c6bc42f58dfdfe,c7a6d58f8640c60f,93ba9f37cb01e4e4
```

`_max` is the largest error of one input; `checks` lists the failed checks;
the hashes cover the ordered forward, the ordered backward and the unordered
forward results of all inputs.


`ctest` runs `test_accuracy` and `test_accuracy_split` (tests/test_accuracy.c:
the same measures with accuracy.c, built with the tests rather than the
benchmarks, the second against the library that splits every array of more
than 1 KiB) on lengths of every kind of pass up to 4096, in under a second.
CI runs them on every backend it builds, emulators included.
