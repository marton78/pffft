# Performance Testing How-To

Practical guide to measuring and comparing pffft's performance, covering both
the automated `perf.py ab` workflow and the direct/manual technique used to
validate it. Read the caveats section before trusting a "neutral" verdict.

## Quick start: the `ab` workflow

```bash
# Compare HEAD against the accepted chain head (or --base <rev> explicitly)
./bench/perf.py ab --label my-change --sizes short --max-len 2048

# Fold a genuinely-faster result into bench_chain.json (human-gated; refuses
# on a "slower" verdict unless --force)
./bench/perf.py accept

# Show the accepted chain plus the last ab run's verdicts
./bench/perf.py status

# Render evolution charts (curves + waterfall) over the accepted chain
./bench/perf.py evolution
```

`ab` builds two worktrees (`.perf/wt-base`, `.perf/wt-var`), interleaves short
timed invocations of both binaries per size, and reports a Mann-Whitney U /
Hodges-Lehmann verdict per `(algo, xform)` group. Key flags: `--variant`/
`--base` (git revisions, default `HEAD` / the chain head), `--prec flt|dbl|
both`, `--sizes short|pow2|nonpow2|all`, `--target local` (repeatable; also
`ssh://host`, `adb[:serial]`, `ios[:udid]`), `--invocations`/`--runs`.

Requirements: Python 3 with `numpy`, `scipy` (`pip3 install numpy scipy`);
`matplotlib` too for chart rendering.

## Caveats — read before trusting a verdict

**1. Interleaved short runs are noisier than they look.** `ab` alternates
tiny (`--runs` reps, a few seconds each) invocations of the base and variant
binaries. This is convenient (no need to rebuild/rerun by hand) but far more
susceptible to scheduler jitter, thermal state transitions, and turbo-boost
ramp-up/down between binaries than one sustained run of each. A genuine
2-4% throughput difference can register as "neutral" under `--sizes short
--invocations 3 --runs 3` and show up clearly under a direct, sustained
comparison (see below). If a change is architecturally significant (a full
rewrite, a fused-multiply-add path, anything touching the FFT inner loop),
**cross-check with the direct method** before accepting a "neutral"
verdict as final.

**2. `vdsp`/`fftpack`/other competitor libraries are not your code.**
`bench_pffft` benchmarks every compiled-in library every run. A "slower"
flag on `vdsp` (Apple's Accelerate, never touched by a pffft-only commit) is
noise — sub-percent shifts are statistically significant but practically
meaningless once timings are this reproducible. Only trust verdicts on the
algo you actually changed (`pffft`/`pffftu`).

**3. On battery, Mac performance can look worse than it should — the
iPhone won't.** If the host machine throttles under battery/thermal
pressure mid-session, only *local* numbers degrade; a connected iOS/Android
device's numbers are unaffected. A sudden Mac-only regression mid-series is
a power/thermal artifact, not a code regression.

**4. Sequential order confounds unless you check for it.** If you benchmark
several checkpoints back-to-back, a monotonic drift (things keep getting
faster, or slower, checkpoint after checkpoint regardless of what changed)
usually means the machine is warming up/cooling down, not that the code is
improving/regressing. Re-run an early checkpoint last; if it comes back to
its original number, the drift wasn't thermal.

**5. `git worktree` "HEAD" means the worktree's OWN head, not yours.**
If you ever touch `bench/targets.py::sync_worktree`, remember: a linked
worktree has its own detached-HEAD state. `git -C <worktree> reset --hard
HEAD` is a no-op relative to the worktree, not "catch up with the main
checkout." `sync_worktree` resolves refs against `REPO_ROOT` first for
exactly this reason — don't reintroduce the bug.

## Direct measurement (bypasses `ab`'s interleaving entirely)

For a trustworthy, low-noise comparison across several versions — or to
sanity-check an `ab` verdict — build each version once and run it in one
sustained block instead of interleaving:

```bash
# One worktree per version you want to compare
git worktree add --detach /tmp/ckpt-master   <sha-or-rev>
git worktree add --detach /tmp/ckpt-mychange HEAD

# Build each (same flags perf.py's own build() uses — disable competitor
# libs you don't need, keeps builds fast)
for wt in /tmp/ckpt-master /tmp/ckpt-mychange; do
  cmake -S "$wt" -B "$wt/build" -DCMAKE_BUILD_TYPE=Release \
    -DPFFFT_USE_TYPE_FLOAT=ON -DPFFFT_USE_TYPE_DOUBLE=ON -DPFFFT_USE_SIMD=ON \
    -DPFFFT_USE_BENCH_FFTW=OFF -DPFFFT_USE_BENCH_GREEN=OFF -DPFFFT_USE_BENCH_KISS=OFF \
    -DPFFFT_USE_BENCH_POCKET=OFF -DPFFFT_USE_BENCH_MKL=OFF -DPFFFT_USE_FFTPACK=OFF \
    -DPFFFT_USE_BENCH_FFTS=OFF -DPFFFT_USE_BENCH_AVFFT=OFF \
    -DPFFFT_BUILD_TESTS=OFF -DPFFFT_BUILD_EXAMPLES=OFF -DPFFFT_BUILD_BENCHMARKS=ON -Wno-dev
  cmake --build "$wt/build" --target bench_pffft_float bench_pffft_double --parallel
done

# One sustained run per version per precision, 20 reps for a real distribution
for label in master mychange; do
  wt=/tmp/ckpt-$label
  tree=$(git -C "$wt" rev-parse HEAD^{tree})
  bin=$wt/build/benchmarks
  "$bin/bench_pffft_float"  --size 256,1024,4096 --runs 20 --samples /tmp/$label.csv \
    --meta label=$label --meta git_tree=$tree --meta target=local --meta host=$(hostname)
  "$bin/bench_pffft_double" --size 256,1024,4096 --runs 20 --samples /tmp/$label.csv \
    --meta label=$label --meta git_tree=$tree --meta target=local --meta host=$(hostname)
done

# Clean up worktrees when done
git worktree remove /tmp/ckpt-master --force
git worktree remove /tmp/ckpt-mychange --force
```

Each `bench_pffft_{float,double}` invocation appends to the same samples
file (append-guarded: mismatched provenance on the same file refuses and
tells you exactly which key differs). `--size` restricts to the sizes you
name (comma-separated; non-power-of-two sizes additionally need
`--non-pow2`). `--runs N` repeats the whole size/algo sweep N times, giving
N MFLOPS data points per `(algo, prec, xform, size)` group — the raw
material for a distribution plot.

Read the resulting CSVs directly if you just want numbers:

```python
import sys; sys.path.insert(0, "bench")
from samples import read_samples
prov, rows = read_samples("/tmp/mychange.csv")
vals = sorted(s.mflops for s in rows
             if s.algo == "pffft" and s.prec == "flt" and s.xform == "real" and s.size == 1024)
print("median:", vals[len(vals)//2])
```

## Plotting the evolution: `bench/plot_evolution.py`

Renders one panel per requested FFT size, each showing a box or violin plot
of the MFLOPS distribution per version (one input file = one version,
identified by its `label` provenance key):

```bash
./bench/plot_evolution.py /tmp/master.csv /tmp/mychange.csv \
    --order master,mychange --sizes 256,1024,4096 \
    --algo pffft --prec flt --xform real --kind box \
    --out bench_results/evolution_pffft_flt_real.png
```

- `files`: one or more samples CSVs or globs (e.g. `.perf/samples/*.csv`).
  Multiple files sharing a `label` are concatenated into one version.
- `--order`: comma-separated version order for the x-axis; defaults to the
  order versions are first seen among the input files.
- `--sizes`: comma-separated FFT sizes, one panel each.
- `--algo`/`--prec`/`--xform`: which series to plot (default `pffft`/`flt`/
  `real`) — rerun with different values to see other combinations; the
  script intentionally does one `(algo, prec, xform)` selection per call
  rather than faceting everything into one busy figure.
- `--kind box|violin`: box plots show quartiles/outliers cleanly; violins
  show the full density shape (useful when a distribution is bimodal, e.g.
  a benchmark that sometimes gets scheduled off a performance core).

A box/violin plot is also a good bug detector: if every version's box looks
suspiciously identical regardless of what changed, don't assume "the
change is neutral" — check whether the underlying metric is actually
varying at all (see caveat #1's root cause: this exact tool caught a bug
where the comparison metric was accidentally constant by design).

## Device targets (SSH / Android / iOS)

`--target ssh://host`, `--target adb[:serial]`, `--target ios[:udid]` (all
repeatable, combine with `--target local`) drive remote benchmarking. Gotchas
specific to each, learned the hard way:

- **iOS (`ios-deploy`)**: needs `--debug` to launch reliably (plain
  `--noninteractive` install-and-run doesn't always start the process) and
  every forwarded `--meta` value must avoid spaces and shell metacharacters
  — `ios-deploy`'s own argument re-splitting isn't shell-quote-safe. Both are
  handled automatically by `IosTarget`; if you see "installs but never
  launches" or an indefinite hang after install, check for a **stale
  `ios-deploy` process from a previous run still holding the device**
  (`ps aux | grep ios-deploy`, kill it) before assuming it's a new bug.
  Keep the phone screen unlocked during a run.
- **Android (`adb`)**: needs the NDK (`~/Library/Android/sdk/ndk` or
  `$ANDROID_NDK`) and a connected/authorized device (`adb devices`).
- **SSH (Raspberry Pi, etc.)**: never launch an unbounded build (e.g.
  `make -jN` compiling something heavy like cmake from source) on a
  resource-constrained remote device without a way to kill it remotely —
  a saturated Pi can become unreachable over SSH (port 22 times out during
  banner exchange even though it still answers ping) and need a physical
  power cycle to recover.

## Recommended process for evaluating a real optimization

1. Fast iteration: `./bench/perf.py ab --label X --sizes short --max-len 2048`
   on `--target local`. Good for "did I break anything / rough direction."
2. Before believing a "neutral" verdict on something that should matter
   (a rewrite, an FMA fusion, anything in the FFT inner loop): run the
   **direct measurement** method above for the specific sizes/precision you
   care about, 15-20 reps, and eyeball it with `plot_evolution.py`.
3. Only `accept` into `bench_chain.json` once you trust the number — accept
   is a human decision, not a formality; it refuses `slower` verdicts for a
   reason.
4. For a device you don't personally have running 24/7 (Pi, phone), budget
   real wall-clock time — a device sweep with several sizes and invocations
   easily takes minutes, and a `--prec both` sweep roughly doubles that.
