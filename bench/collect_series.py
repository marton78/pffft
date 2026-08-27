#!/usr/bin/env python3
"""Sustained direct MFLOPS measurement of ONE checkpoint on several targets.

Unlike `perf.py ab` (which interleaves short A/B invocations), this builds a
checkpoint once per target and runs one long sustained sweep per precision,
appending to a per-(label, target, prec) samples file. Those files are the
input for `bench/plot_evolution.py` (one figure per platform, one panel per
size, one box per variant).

Deliberately does NOT use targets.sync_worktree(): that runs `git clean -ffd`,
which nukes the worktree's build/ build-android/ build-ios/ directories and
forces a full (Xcode!) rebuild at every checkpoint. A plain detached checkout
keeps every build tree incremental.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bench"))

from samples import read_samples                                  # noqa: E402
from targets import build, make_target, sanitize_meta_value       # noqa: E402


def sh(*cmd: str) -> str:
    return subprocess.run(cmd, capture_output=True, check=True, text=True).stdout


def checkout(wt: Path, ref: str) -> str:
    """Detached-checkout `ref` (resolved in REPO) in `wt`; return tree hash."""
    sha = sh("git", "-C", str(REPO), "rev-parse", f"{ref}^{{commit}}").strip()
    sh("git", "-C", str(wt), "checkout", "--detach", "--force", sha)
    return sh("git", "-C", str(wt), "rev-parse", "HEAD^{tree}").strip()


def safe_name(s: str) -> str:
    return s.replace("/", "_").replace(":", "_")


def run_one(tgt, bindir: Path, label: str, tree: str, sizes: list[int],
            runs: int, out_dir: Path, precs: list[str],
            max_len: int | None, warmup_steady_sec: float = 0.0) -> list[str]:
    """One target, both precisions, sustained. Returns log lines.

    `warmup_steady_sec`, when > 0, is passed straight to the binary as
    --warmup-steady: it runs an in-process adaptive warmup (burst the
    largest requested size until throughput stabilizes, capped at this
    many seconds) immediately before the recorded reps, in the SAME
    invocation/app-launch -- avoiding both a second iOS app relaunch and
    guesswork over how many reps a fixed discard count should be.

    One size per invocation (matches perf.py's cmd_ab): ios-deploy's
    app-arg safety filter strips commas from --size lists, silently
    collapsing a joined list to a single garbage size on iOS.
    """
    sizes_meta = sanitize_meta_value(",".join(str(n) for n in sizes))
    logs: list[str] = []
    for prec in precs:
        exe = tgt.binary_path(bindir, prec)
        path = out_dir / f"{label}-{safe_name(tgt.name)}-{prec}.csv"
        meta = [f"label={label}", f"git_tree={tree}",
                f"sizes={sizes_meta}", f"runs={runs}",
                *tgt.meta_list(),
                f"datetime={time.strftime('%Y-%m-%dT%H:%MZ', time.gmtime())}"]
        for size in sizes:
            extra = [] if size & (size - 1) == 0 else ["--non-pow2"]
            cmd = [str(exe),
                   "--size", str(size),
                   "--runs", str(runs),
                   *(["--max-len", str(max_len)] if max_len else []),
                   *(["--warmup-steady", str(warmup_steady_sec)] if warmup_steady_sec > 0 else []),
                   "--samples", str(path), *extra,
                   *[x for m in meta for x in ("--meta", m)]]
            t0 = time.time()
            for attempt in range(3):
                n_lines = sum(1 for _ in tgt.run(cmd))
                if n_lines > 0:
                    break
                # Empty capture (e.g. iOS console-scrape flake): if this was
                # the invocation that just created the local file, it wrote
                # a 0-byte file that would satisfy exists() on the NEXT
                # invocation and silently switch it to header-less append.
                # Never leave that behind -- retry, and only keep the file
                # if this attempt (or a later one) actually produced rows.
                if path.exists() and path.stat().st_size == 0:
                    path.unlink()
                if attempt == 2:
                    raise RuntimeError(
                        f"{tgt.name} {prec} n={size}: empty capture "
                        f"after 3 attempts (device/console flake)")
                time.sleep(3)
            logs.append(f"  {tgt.name:<24} {prec} n={size:<6} "
                        f"{time.time() - t0:6.1f}s -> {path.name}")
    return logs


def summarize(path: Path, algo: str = "pffft") -> str:
    import statistics
    prov, rows = read_samples(path)
    out = []
    keys = sorted({(r.xform, r.size) for r in rows if r.algo == algo})
    for xform, size in keys:
        v = [r.mflops for r in rows
             if r.algo == algo and r.xform == xform and r.size == size]
        if v:
            out.append(f"{xform}/{size}={statistics.median(v):.0f}")
    return " ".join(out)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--target", action="append", required=True)
    ap.add_argument("--sizes", default="256,1024,4096")
    ap.add_argument("--runs", type=int, default=20)
    ap.add_argument("--prec", default="both", choices=["flt", "dbl", "both"])
    ap.add_argument("--wt", default=".perf/wt-seq")
    ap.add_argument("--out", default=".perf/series")
    ap.add_argument("--max-len", type=int, default=None)
    ap.add_argument("--warmup-runs", type=int, default=0,
                    help="extra sustained reps run BEFORE the recorded ones, "
                    "discarded, to reach thermal steady state first on every "
                    "requested target (mitigates the cold-start plateau seen "
                    "on a fanless Mac, and reduces session-to-session drift "
                    "on phones too); ignored if --warmup-steady is given")
    ap.add_argument("--warmup-steady", type=float, default=0.0,
                    help="in-process adaptive warmup cap in seconds, passed "
                    "to the binary as --warmup-steady (see bench_pffft.c): "
                    "bursts the largest requested size until MFLOPS "
                    "stabilizes, in the SAME invocation as the recorded "
                    "reps. Preferred over --warmup-runs: detects actual "
                    "steady state instead of guessing a rep count, and "
                    "halves iOS app-relaunch overhead")
    ap.add_argument("--serial-measure", action="store_true",
                    help="measure targets one after another (no host contention)")
    a = ap.parse_args()

    wt = (REPO / a.wt) if not Path(a.wt).is_absolute() else Path(a.wt)
    out_dir = (REPO / a.out) if not Path(a.out).is_absolute() else Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    sizes = [int(x) for x in a.sizes.split(",")]
    precs = ["flt", "dbl"] if a.prec == "both" else [a.prec]

    tree = checkout(wt, a.ref)
    print(f"== {a.label}: {a.ref} tree={tree[:8]}", flush=True)

    bindir = build(wt)                       # host build (also iOS/adb anchor)
    tgts = [make_target(t) for t in a.target]

    # Cross-builds first, serialized: they are host-CPU heavy and must not
    # overlap a measurement.
    for tgt in tgts:
        for prec in precs:
            tgt.binary_path(bindir, prec)
    if a.warmup_steady:
        print(f"-- warmup: adaptive, <= {a.warmup_steady:.0f}s per target/prec/size "
              f"(in-process, see --warmup-steady)", flush=True)
    elif a.warmup_runs:
        print(f"-- warmup: {a.warmup_runs} discarded reps per target/prec "
              f"(thermal/power steady-state, not recorded)", flush=True)
        warm_dir = out_dir / ".warmup-scratch"
        warm_dir.mkdir(parents=True, exist_ok=True)
        for tgt in tgts:
            for prec in precs:
                exe = tgt.binary_path(bindir, prec)
                for size in sizes:
                    extra = [] if size & (size - 1) == 0 else ["--non-pow2"]
                    cmd = [str(exe), "--size", str(size),
                          "--runs", str(a.warmup_runs),
                          "--samples", str(warm_dir / f"warmup-{tgt.name}-{prec}.csv"),
                          *extra]
                    for _ in tgt.run(cmd):
                        pass
        import shutil
        shutil.rmtree(warm_dir, ignore_errors=True)
    print(f"-- built for: {', '.join(t.name for t in tgts)}", flush=True)

    fails: list[str] = []
    if a.serial_measure or len(tgts) == 1:
        for tgt in tgts:
            try:
                print("\n".join(run_one(tgt, bindir, a.label, tree, sizes,
                                        a.runs, out_dir, precs, a.max_len,
                                        a.warmup_steady)),
                      flush=True)
            except Exception as e:                       # noqa: BLE001
                fails.append(f"{tgt.name}: {e}")
    else:
        with cf.ThreadPoolExecutor(max_workers=len(tgts)) as ex:
            futs = {ex.submit(run_one, tgt, bindir, a.label, tree, sizes,
                              a.runs, out_dir, precs, a.max_len,
                              a.warmup_steady): tgt
                    for tgt in tgts}
            for fut in cf.as_completed(futs):
                tgt = futs[fut]
                try:
                    print("\n".join(fut.result()), flush=True)
                except Exception as e:                   # noqa: BLE001
                    fails.append(f"{tgt.name}: {e}")

    print(f"\n== medians ({a.label}) ==")
    for tgt in tgts:
        for prec in precs:
            p = out_dir / f"{a.label}-{safe_name(tgt.name)}-{prec}.csv"
            if p.exists():
                print(f"  {tgt.name:<24} {prec} pffft {summarize(p)}")
    for f in fails:
        print(f"FAILED {f}", file=sys.stderr)
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
