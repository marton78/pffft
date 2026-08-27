#!/usr/bin/env python3
"""Like collect_series.py, but for benchmarking commits from a branch that
predates the samples/CLI benchmark harness (e.g. `perf`, `master`).

The harness itself (benchmarks/bench_pffft.c's --samples/--runs/--meta CLI,
bench/collect_series.py, bench/targets.py, bench/plot_evolution.py) lives
only on `bench-update`, which forked off the *same* base commit as `perf`
(fd0d9c5) but never received perf's optimization commits, and vice versa.
A plain `git checkout <perf-commit>` therefore lands on a tree whose
benchmarks/bench_pffft.c predates --samples entirely.

Fix: keep the worktree pinned to bench-update's harness (CMakeLists.txt,
cmake/, benchmarks/, bench/, include/) and overlay ONLY src/ (the actual
FFT implementation) from the requested ref. Verified safe for this sweep:
bench-update makes zero changes under src/ relative to the fork point, and
every perf-branch commit in this sweep touches only files under src/ (the
one exception, a docs-only commit, touches nothing under src/ either).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bench"))

from collect_series import run_one, summarize          # noqa: E402
from targets import build, make_target, sanitize_meta_value  # noqa: E402

HARNESS_REF = "bench-update"


def sh(*cmd: str, cwd: Path | None = None) -> str:
    return subprocess.run(cmd, capture_output=True, check=True, text=True,
                          cwd=cwd).stdout


def overlay_checkout(wt: Path, harness_ref: str, src_ref: str) -> str:
    """Pin `wt` to harness_ref, then overlay src/ from src_ref.

    Returns src_ref's own tree hash (resolved in REPO) as the provenance
    `git_tree` value -- it is what actually varies between labels here.
    """
    harness_sha = sh("git", "-C", str(REPO), "rev-parse",
                     f"{harness_ref}^{{commit}}").strip()
    src_sha = sh("git", "-C", str(REPO), "rev-parse",
                f"{src_ref}^{{commit}}").strip()
    sh("git", "-C", str(wt), "checkout", "--detach", "--force", harness_sha)
    sh("git", "-C", str(wt), "checkout", src_sha, "--", "src")
    dirty = sh("git", "-C", str(wt), "status", "--porcelain", "--", "src")
    if not dirty and src_sha != harness_sha:
        # Identical src/ content is plausible for a docs-only commit; not
        # an error, just means this label is byte-identical to its parent.
        pass
    return sh("git", "-C", str(REPO), "rev-parse", f"{src_ref}^{{tree}}").strip()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True,
                    help="perf/master commit-ish supplying src/")
    ap.add_argument("--harness-ref", default=HARNESS_REF,
                    help="branch supplying the benchmark harness (default: %(default)s)")
    ap.add_argument("--label", required=True)
    ap.add_argument("--target", action="append", required=True)
    ap.add_argument("--sizes", default="256,1024,4096")
    ap.add_argument("--runs", type=int, default=20)
    ap.add_argument("--prec", default="both", choices=["flt", "dbl", "both"])
    ap.add_argument("--wt", default=".perf/wt-seq")
    ap.add_argument("--out", default=".perf/series")
    ap.add_argument("--max-len", type=int, default=None)
    ap.add_argument("--warmup-runs", type=int, default=0)
    ap.add_argument("--warmup-steady", type=float, default=0.0,
                    help="in-process adaptive warmup cap in seconds (see "
                    "bench_pffft.c --warmup-steady); preferred over "
                    "--warmup-runs, ignored if both given")
    ap.add_argument("--serial-measure", action="store_true")
    a = ap.parse_args()

    wt = (REPO / a.wt) if not Path(a.wt).is_absolute() else Path(a.wt)
    out_dir = (REPO / a.out) if not Path(a.out).is_absolute() else Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    sizes = [int(x) for x in a.sizes.split(",")]
    precs = ["flt", "dbl"] if a.prec == "both" else [a.prec]

    tree = overlay_checkout(wt, a.harness_ref, a.ref)
    print(f"== {a.label}: src={a.ref} harness={a.harness_ref} tree={tree[:8]}",
         flush=True)

    bindir = build(wt)
    tgts = [make_target(t) for t in a.target]

    for tgt in tgts:
        for prec in precs:
            tgt.binary_path(bindir, prec)

    if a.warmup_steady:
        print(f"-- warmup: adaptive, <= {a.warmup_steady:.0f}s per target/prec/size "
              f"(in-process, see --warmup-steady)", flush=True)
    elif a.warmup_runs:
        print(f"-- warmup: {a.warmup_runs} discarded reps per target/prec",
             flush=True)
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
    for tgt in tgts:
        try:
            print("\n".join(run_one(tgt, bindir, a.label, tree, sizes,
                                    a.runs, out_dir, precs, a.max_len,
                                    a.warmup_steady)),
                 flush=True)
        except Exception as e:                          # noqa: BLE001
            fails.append(f"{tgt.name}: {e}")

    print(f"\n== medians ({a.label}) ==")
    for tgt in tgts:
        for prec in precs:
            p = out_dir / f"{a.label}-{tgt.name.replace('/', '_').replace(':', '_')}-{prec}.csv"
            if p.exists():
                print(f"  {tgt.name:<24} {prec} pffft {summarize(p)}")
    for f in fails:
        print(f"FAILED {f}", file=sys.stderr)
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
