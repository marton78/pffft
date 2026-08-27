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

**3. Thermal state matters most on fanless hardware, and it's not just a
"battery" thing.** A MacBook Air has no fan — sustained load saturates its
passive heatsink within a few minutes regardless of AC vs battery, producing
a real cold-start-then-plateau pattern (confirmed: identical behavior with
the charger connected). Re-running the very first checkpoint late in a
session reproduced the same ~35% "improvement" with zero code change,
proving it was thermal, not any commit. Phones show a different character —
not a clean monotonic drift, but real session-to-session bounce (one
iPhone 11 Pro session moved by up to 6% between two measurements of
*byte-identical compiled code*) — so don't assume "it's a phone, it won't
throttle." Mitigate with `bench/collect_series.py --warmup-runs N`
(discards N sustained reps before recording, on every requested target) and
physically elevate the device for airflow; see "Reducing thermal noise"
below.

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

**6. A significant p-value is not the same as a real effect.** With 20
sustained reps, `bench/stats.py::compare`'s Mann-Whitney U test flags
`p < 0.01` almost every time — even between two runs of *the exact same
binary*. Measured proof: comparing two checkpoints on Android where zero
intervening commit touched the tested sizes, 11 of 12 `(prec, xform, size)`
groups came back "significant" at shifts as small as 0.08%. The variance in
these benchmarks is small enough that MWU's null hypothesis ("these two
samples are from literally the same distribution") is nearly always false
in a trivial sense, regardless of whether the code changed. Don't accept a
verdict on `p < alpha` alone — see "Objective significance" below for the
noise-floor criterion this project actually uses.

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

## Objective significance: a noise floor, not just a p-value

Caveat #6 above shows why `p < 0.01` alone can't tell you a change is real.
What actually distinguishes signal from noise in this project:

1. **Measure an empirical noise floor per platform first.** Take two
   samples files where you *know* nothing relevant changed between them
   (e.g. two checkpoints where the diff between them doesn't touch the
   `(prec, xform, size)` group you're inspecting) and run `compare()` on
   them. The resulting `|shift_pct|` values are pure noise — their max
   (not median; be conservative) is your floor. Measured so far: Android
   ~0.5%, Mac (once past the cold-start plateau) ~0.4%, iPhone ~6.5% —
   these differ by more than 10x, so a single blanket threshold across
   platforms is wrong.
2. **Require the effect to clear the floor, not just beat alpha.** A
   commit's median Hodges-Lehmann shift has to exceed the platform's own
   noise floor by a real margin (this project uses the raw floor as a strict
   pass/fail line; a stricter version would require 3-5x it) before it
   counts as a real effect — this is minimum-effect-size / equivalence
   testing in spirit (TOST), just without the extra machinery: instead of
   asking "is the shift different from exactly zero" (always yes, per
   caveat #6), ask "is the shift bigger than what this platform's own
   same-code noise looks like."
3. **Negative-control your own commit's scope.** If a commit touches only
   `pf_neon_double.h`, its `flt` groups should sit *inside* the noise floor.
   If they don't, that's the tell the "effect" you're seeing on the touched
   groups is session noise too, not the commit's doing — this caught every
   false-positive iOS swing during this project's own investigation (`flt`
   moving in lockstep with `dbl` on a double-only commit).
4. **Require cross-platform sign agreement** (same direction on 2 of 3
   platforms, or an explicit code-scope reason one platform is exempt, e.g.
   an ARMv7-only guard on arm64-only test hardware) before trusting a
   single-platform result.

## `bench/collect_series.py`: sustained multi-platform direct measurement

A driver around `bench/targets.py`'s cross-build machinery for exactly the
"direct measurement" method above, across `local`/`adb`/`ios` in one call,
one size per invocation (`ios-deploy`'s `--args` re-splitting silently
strips commas from a joined `--size` list — this bit us once):

```bash
./bench/collect_series.py --ref <sha-or-rev> --label my-checkpoint \
    --target local --target adb --target ios \
    --sizes 256,1024,4096 --runs 20 --warmup-runs 3 \
    --out .perf/series
```

- `--warmup-runs N`: runs N discarded sustained reps per target/precision
  *before* the recorded ones, to reach thermal steady state first (see
  caveat #3). Applies to every `--target` given, not just `local`.
- `--wt`: worktree to build in (default `.perf/wt-seq`), reused
  incrementally across calls with different `--ref` — unlike `sync_worktree`
  used by `ab`, this does a plain detached checkout, not `git clean -ffd`,
  so `build/`, `build-android/`, `build-ios/` stay incremental across
  checkpoints (a full Xcode rebuild per checkpoint is otherwise the
  dominant cost).
- Output: one CSV per `(label, target, prec)` under `--out`, directly
  consumable by `plot_evolution.py` or `plot_summary.py`.

## `bench/plot_summary.py`: reader-facing before/after plots

`plot_evolution.py`'s many-column box plot is the right tool for *deriving*
a verdict (the engineering audit trail — did this specific commit move the
numbers?) but the wrong one for *communicating* it: a reader has to know
what a box plot and 10 commit labels mean before the picture says anything.
`plot_summary.py` renders exactly two bars — baseline vs final — with a
shaded noise-floor band (from the previous section) so "the change is
bigger than measurement noise" is visible without reading any of this file:

```bash
./bench/plot_summary.py baseline.csv final.csv \
    --baseline-label master --final-label final \
    --platform "Android (Motorola Edge 50 Neo)" \
    --sizes 256,1024,4096 --prec flt --xform real \
    --noise-floor-pct 0.6 \
    --out bench_results/summary-adb-flt-real.png
```

Each panel's title states the percent change and a plain-language verdict
(`faster` / `slower` / `within noise`) computed directly from the
noise-floor comparison — a bar that doesn't clear the shaded band is
labeled `within noise` even if the raw percentage looks nonzero. Use
`--caption` for a one-line footnote when the baseline needed a caveat (e.g.
substituting a warm checkpoint for a cold-start-confounded `master` on
fanless hardware — see caveat #3). Recommended structure for a public
README: lead with 2-3 of these summary plots, one line explaining the gray
band, and link the full per-commit evolution plots as the "how we know"
appendix rather than including them inline.

## Reducing thermal/session noise before measuring

- **Software (most reliable, works regardless of physical setup):**
  `--warmup-runs` on `collect_series.py` (above) — start every *recorded*
  measurement from the same already-warmed-up state instead of fighting the
  warm-up transient.
- **Elevate the device for airflow.** A phone or laptop resting flat on a
  desk traps heat against its underside; propping it on something a few cm
  tall lets air circulate underneath. Free, zero risk.
- **Avoid direct cold-pack contact for condensation reasons.** A chilled
  gel pad (the reusable kind used for bruises) cools faster than air but
  risks condensation if it's colder than the room's dew point — moisture
  can wick into a phone's mic/speaker/USB-C openings or, worse, a laptop's
  fan-less intake vents near the logic board. A passive aluminum stand or
  cooling plate (no chilling) is the safer middle ground; a clip-on phone
  cooling fan (a real, inexpensive product for mobile gaming) is the
  purpose-built solution if you want active cooling without condensation
  risk.
- **Fanless hardware (e.g. MacBook Air) throttles fastest and hardest** —
  it has no active cooling at all, so `--warmup-runs` matters most there.

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
   (a rewrite, an FMA fusion, anything in the FFT inner loop): run
   `bench/collect_series.py` (direct, sustained, multi-platform in one
   call, `--warmup-runs` to control for thermal state) for the specific
   sizes/precision you care about, 15-20 reps, and check it against your
   platforms' noise floors (see "Objective significance") rather than
   eyeballing `p < alpha` alone.
3. Only `accept` into `bench_chain.json` once the effect clears the noise
   floor and agrees in direction across platforms — accept is a human
   decision, not a formality; it refuses `slower` verdicts for a reason.
4. Render `plot_summary.py` before/after plots for anything going in a
   README or PR description; keep `plot_evolution.py`'s full per-commit
   view as the audit trail, not the headline artifact.
5. For a device you don't personally have running 24/7 (Pi, phone), budget
   real wall-clock time — a device sweep with several sizes and invocations
   easily takes minutes, and a `--prec both` sweep roughly doubles that.
