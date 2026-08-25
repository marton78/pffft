"""Worktree management and Target abstraction for the benchmark harness."""

from __future__ import annotations

import os
import platform
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator

REPO_ROOT = Path(__file__).resolve().parent.parent

BASE_CMAKE_FLAGS = [
    "-DCMAKE_BUILD_TYPE=Release",
    "-DPFFFT_USE_TYPE_FLOAT=ON", "-DPFFFT_USE_TYPE_DOUBLE=ON",
    "-DPFFFT_USE_SIMD=ON",
    "-DPFFFT_USE_BENCH_FFTW=OFF", "-DPFFFT_USE_BENCH_GREEN=OFF",
    "-DPFFFT_USE_BENCH_KISS=OFF", "-DPFFFT_USE_BENCH_POCKET=OFF",
    "-DPFFFT_USE_BENCH_MKL=OFF", "-DPFFFT_USE_FFTPACK=OFF",
    "-DPFFFT_USE_BENCH_FFTS=OFF", "-DPFFFT_USE_BENCH_AVFFT=OFF",
    "-DPFFFT_BUILD_TESTS=OFF", "-DPFFFT_BUILD_EXAMPLES=OFF",
    "-DPFFFT_BUILD_BENCHMARKS=ON", "-Wno-dev",
]


def sh(*cmd: str) -> str:
    """Run a command, return its captured stdout."""
    return subprocess.run(cmd, capture_output=True, check=True,
                          text=True).stdout


def ensure_worktrees() -> tuple[Path, Path]:
    """Return paths to .perf/wt-base and .perf/wt-var, creating them once."""
    base = REPO_ROOT / ".perf" / "wt-base"
    var = REPO_ROOT / ".perf" / "wt-var"
    if not (base / ".git").exists():
        sh("git", "-C", str(REPO_ROOT), "worktree", "add",
           str(base), "HEAD")
    if not (var / ".git").exists():
        sh("git", "-C", str(REPO_ROOT), "worktree", "add", str(var), "HEAD")
    return base, var


def sync_worktree(wt: Path, ref: str) -> str:
    """Point a worktree at ref (commit-ish or raw tree) and return its
    content-identity tree hash."""
    sh("git", "-C", str(REPO_ROOT), "fetch", "--all", "--quiet")
    kind = subprocess.run(["git", "-C", str(REPO_ROOT), "cat-file", "-t", ref],
                          capture_output=True, text=True)
    if kind.returncode == 0 and kind.stdout.strip() == "tree":
        # `git reset --hard` needs a commit; check out a bare tree via
        # read-tree instead (chain heads are stored as tree hashes).
        sh("git", "-C", str(wt), "reset", "--hard", "HEAD")
        sh("git", "-C", str(wt), "clean", "-ffd", "-e", ".perf")
        sh("git", "-C", str(wt), "read-tree", "-u", "--reset", ref)
        return ref
    sh("git", "-C", str(wt), "reset", "--hard", ref)
    # -ffd: remove nested git dirs too; -e .perf: keep local benchmark state
    sh("git", "-C", str(wt), "clean", "-ffd", "-e", ".perf")
    return sh("git", "-C", str(wt), "rev-parse", "HEAD^{tree}").strip()


def build(wt: Path,
          targets: tuple[str, ...] = ("bench_pffft_float",
                                      "bench_pffft_double")) -> Path:
    """Configure (once) and incrementally build targets; return bin dir."""
    bdir = wt / "build"
    if not (bdir / "CMakeCache.txt").exists():
        subprocess.run(["cmake", "-S", str(wt), "-B", str(bdir),
                        *BASE_CMAKE_FLAGS], check=True)
    subprocess.run(["cmake", "--build", str(bdir),
                    "--target", *targets,
                    "-j", str(os.cpu_count())], check=True)
    return bdir / "benchmarks"


class Target(ABC):
    """A machine that can run benchmark binaries."""

    name: str

    @abstractmethod
    def meta(self) -> dict[str, str]:
        """Environment metadata (must include target= and host=)."""

    @abstractmethod
    def binary_path(self, wt_build: Path, prec: str) -> Path:
        """Path to the benchmark binary for prec in ("flt", "dbl")."""

    def run(self, cmd: list[str],
            cwd: Path | None = None) -> Iterator[str]:
        """Run cmd on this target, yielding stdout lines."""
        proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                                stderr=subprocess.DEVNULL, text=True)
        try:
            for line in proc.stdout:
                yield line.rstrip("\n")
            ret = proc.wait()
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()
        if ret != 0:
            raise RuntimeError(
                f"command failed ({ret}): {' '.join(cmd)}")

    def meta_list(self) -> list[str]:
        """Metadata as k=v strings."""
        return [f"{k}={v}" for k, v in self.meta().items()]


class LocalTarget(Target):
    name = "local"

    def meta(self) -> dict[str, str]:
        return {"target": self.name, "host": platform.node()}

    def binary_path(self, wt_build: Path, prec: str) -> Path:
        exe = {"flt": "bench_pffft_float",
               "dbl": "bench_pffft_double"}[prec]
        return wt_build / exe


def make_target(spec: str) -> Target:
    """Dispatch a target spec to a Target instance."""
    if spec == "local":
        return LocalTarget()
    if spec.startswith("ssh://"):
        raise NotImplementedError("ssh targets arrive in Task 8")
    if spec.startswith("adb"):
        raise NotImplementedError("adb targets arrive in Task 9")
    if spec.startswith("ios"):
        raise NotImplementedError("ios targets arrive in Task 10")
    raise ValueError(f"unknown target spec: {spec!r}")
