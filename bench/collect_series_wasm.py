#!/usr/bin/env python3
"""Direct sustained MFLOPS measurement of a WASM build, run under Node.

Same harness-overlay trick as collect_series_overlay.py (the samples/CLI
benchmark harness lives only on bench-update; perf-branch commits predate
it), but targets Emscripten instead of a native Target. The WASM commit
(e836207) additionally touches CMakeLists.txt/cmake/target_optimizations.cmake
to add its dedicated WASM SIMD backend, so --overlay-extra lets the caller
pull those in too (the pre-wasm checkpoint doesn't need them: bench-update's
existing NEON-emulation Emscripten path already matches what that commit's
own tree would produce, since it makes zero cmake changes itself).

No Target class here: emcc's default Node output only needs stdout capture
(``--samples -``), mirrored into a local CSV exactly like SshTarget/
AdbTarget/IosTarget in bench/targets.py -- no NODERAWFS/real-FS wiring
needed.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "bench"))

from targets import BASE_CMAKE_FLAGS, sanitize_meta_value  # noqa: E402

HARNESS_REF = "bench-update"


def sh(*cmd: str, cwd: Path | None = None) -> str:
    return subprocess.run(cmd, capture_output=True, check=True, text=True,
                          cwd=cwd).stdout


def overlay_checkout(wt: Path, harness_ref: str, src_ref: str,
                     extra_paths: list[str]) -> str:
    harness_sha = sh("git", "-C", str(REPO), "rev-parse",
                     f"{harness_ref}^{{commit}}").strip()
    src_sha = sh("git", "-C", str(REPO), "rev-parse",
                f"{src_ref}^{{commit}}").strip()
    sh("git", "-C", str(wt), "checkout", "--detach", "--force", harness_sha)
    sh("git", "-C", str(wt), "checkout", src_sha, "--", "src", *extra_paths)
    return sh("git", "-C", str(REPO), "rev-parse", f"{src_ref}^{{tree}}").strip()


def build(wt: Path) -> Path:
    bdir = wt / "build-wasm"
    if not (bdir / "CMakeCache.txt").exists():
        # Emscripten's musl-based libc headers gate strdup/strndup behind
        # POSIX feature-test macros that -std=c99 alone doesn't imply
        # (unlike Apple's libc, which the native build relies on).
        subprocess.run(["emcmake", "cmake", "-S", str(wt), "-B", str(bdir),
                        *BASE_CMAKE_FLAGS,
                        "-DCMAKE_C_FLAGS=-D_DEFAULT_SOURCE"], check=True)
    subprocess.run(["emmake", "cmake", "--build", str(bdir),
                    "--target", "bench_pffft_float", "bench_pffft_double",
                    "-j", "8"], check=True)
    return bdir / "benchmarks"


def binary_path(bindir: Path, prec: str) -> Path:
    exe = {"flt": "bench_pffft_float", "dbl": "bench_pffft_double"}[prec]
    return bindir / f"{exe}.js"


def run_node(js_path: Path, args: list[str]) -> list[str]:
    proc = subprocess.run(["node", str(js_path), *args],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"node {js_path} failed ({proc.returncode}): "
                           f"{proc.stderr[-4000:]}")
    return [ln for ln in proc.stdout.splitlines() if ln]


def write_samples(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        rows = [ln for ln in lines if not ln.startswith(("#", "algo,"))]
        with open(path, "a") as fh:
            fh.write("\n".join(rows) + ("\n" if rows else ""))
    else:
        path.write_text("\n".join(lines) + ("\n" if lines else ""))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--harness-ref", default=HARNESS_REF)
    ap.add_argument("--overlay-extra", default="",
                    help="comma-separated extra paths to overlay from --ref "
                         "besides src/ (e.g. CMakeLists.txt,cmake)")
    ap.add_argument("--label", required=True)
    ap.add_argument("--sizes", default="256,1024,4096")
    ap.add_argument("--runs", type=int, default=20)
    ap.add_argument("--warmup-runs", type=int, default=0)
    ap.add_argument("--prec", default="both", choices=["flt", "dbl", "both"])
    ap.add_argument("--wt", default=".perf/wt-wasm")
    ap.add_argument("--out", default=".perf/series")
    a = ap.parse_args()

    wt = (REPO / a.wt) if not Path(a.wt).is_absolute() else Path(a.wt)
    out_dir = (REPO / a.out) if not Path(a.out).is_absolute() else Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    sizes = [int(x) for x in a.sizes.split(",")]
    precs = ["flt", "dbl"] if a.prec == "both" else [a.prec]
    extra = [p for p in a.overlay_extra.split(",") if p]

    tree = overlay_checkout(wt, a.harness_ref, a.ref, extra)
    print(f"== {a.label}: src={a.ref} harness={a.harness_ref} "
         f"extra={extra} tree={tree[:8]}", flush=True)

    bindir = build(wt)
    print("-- built for: wasm-node", flush=True)

    if a.warmup_runs:
        print(f"-- warmup: {a.warmup_runs} discarded reps per prec/size",
             flush=True)
        for prec in precs:
            js = binary_path(bindir, prec)
            for size in sizes:
                extra_flag = [] if size & (size - 1) == 0 else ["--non-pow2"]
                run_node(js, ["--size", str(size), "--runs", str(a.warmup_runs),
                             "--samples", "-", *extra_flag])

    sizes_meta = sanitize_meta_value(",".join(str(n) for n in sizes))
    for prec in precs:
        js = binary_path(bindir, prec)
        path = out_dir / f"{a.label}-wasm-node-{prec}.csv"
        meta = [f"label={a.label}", f"git_tree={tree}",
               f"sizes={sizes_meta}", f"runs={a.runs}",
               "target=wasm-node", "host=node+emcc",
               f"datetime={time.strftime('%Y-%m-%dT%H:%MZ', time.gmtime())}"]
        for size in sizes:
            extra_flag = [] if size & (size - 1) == 0 else ["--non-pow2"]
            t0 = time.time()
            lines = run_node(js, ["--size", str(size), "--runs", str(a.runs),
                                  "--samples", "-", *extra_flag,
                                  *[x for m in meta for x in ("--meta", m)]])
            write_samples(path, lines)
            print(f"  wasm-node               {prec} n={size:<6} "
                 f"{time.time() - t0:6.1f}s -> {path.name}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
