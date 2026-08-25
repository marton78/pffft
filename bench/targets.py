"""Worktree management and Target abstraction for the benchmark harness."""

from __future__ import annotations

import os
import platform
import re
import shlex
import subprocess
import sys
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator

REPO_ROOT = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(REPO_ROOT))
from cross_build_android import (detect_device_march, find_ndk,          # noqa: E402
                                 pffft_android_cmake_argv)

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

# Sources rsynced to a remote target: everything cmake needs to build the
# benchmark binary (competitor libs are OFF, so no submodules are required).
REMOTE_SOURCES = ("src", "include", "CMakeLists.txt", "cmake",
                  "benchmarks/bench_pffft.c")
REMOTE_EXCLUDES = (".git", ".perf", "build")
REMOTE_ROOT = "fftbench"           # under $HOME on the remote machine

def sh(*cmd: str) -> str:
    """Run a command, return its captured stdout."""
    return subprocess.run(cmd, capture_output=True, check=True,
                          text=True).stdout


def clean_value(value: str) -> str:
    """Sanitize a provenance value: no newlines or commas."""
    return re.sub(r"[,\r\n]+", ";", str(value))


def stream_stdout(argv: list[str],
                  cwd: Path | None = None) -> Iterator[str]:
    """Run argv, yield rstripped stdout lines; raise on nonzero exit."""
    proc = subprocess.Popen(argv, cwd=cwd, stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)
    ret: int | None
    try:
        for line in proc.stdout:
            yield line.rstrip("\n")
        ret = proc.wait()
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait()
    if ret != 0:
        raise RuntimeError(f"command failed ({ret}): {' '.join(argv)}")


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
    content-identity tree hash.

    `ref` is resolved against REPO_ROOT before touching `wt`. `wt` is a
    linked worktree with its OWN HEAD file: `git -C wt reset --hard HEAD`
    resets it to wt's own frozen detached HEAD, not REPO_ROOT's current
    branch tip -- passing symbolic refs straight through silently makes
    the default `--variant HEAD` never advance past the commit the
    worktree was first created at. Resolve first, always via REPO_ROOT.
    """
    sh("git", "-C", str(REPO_ROOT), "fetch", "--all", "--quiet")
    kind = subprocess.run(["git", "-C", str(REPO_ROOT), "cat-file", "-t", ref],
                          capture_output=True, text=True)
    if kind.returncode == 0 and kind.stdout.strip() == "tree":
        # `git reset --hard` needs a commit; check out a bare tree via
        # read-tree instead (chain heads are stored as tree hashes).
        # A full tree hash needs no REPO_ROOT-relative resolution.
        sh("git", "-C", str(wt), "reset", "--hard", "HEAD")
        sh("git", "-C", str(wt), "clean", "-ffd", "-e", ".perf")
        sh("git", "-C", str(wt), "read-tree", "-u", "--reset", ref)
        return ref
    resolved = sh("git", "-C", str(REPO_ROOT), "rev-parse", ref).strip()
    sh("git", "-C", str(wt), "reset", "--hard", resolved)
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
        yield from stream_stdout(cmd, cwd=cwd)

    def meta_list(self) -> list[str]:
        """Metadata as k=v strings (values sanitized)."""
        return [f"{k}={clean_value(v)}" for k, v in self.meta().items()]


class LocalTarget(Target):
    name = "local"

    def meta(self) -> dict[str, str]:
        return {"target": self.name, "host": platform.node()}

    def binary_path(self, wt_build: Path, prec: str) -> Path:
        exe = {"flt": "bench_pffft_float",
               "dbl": "bench_pffft_double"}[prec]
        return wt_build / exe


class SshTarget(Target):
    """A remote machine reached over ssh/rsync (e.g. a Raspberry Pi).

    Sources are rsynced to ``~/fftbench/<worktree-name>/`` and built
    natively there with cmake (incremental across runs); benchmark
    invocations execute remotely with ``--samples -`` so the CSV streams
    back over stdout and is written locally into SAMPLES_DIR.
    """

    def __init__(self, host: str, user: str | None = None):
        self.host = host
        self.user = user or None
        self._remote_name: str | None = None
        self._prepared: set[str] = set()

    @property
    def name(self) -> str:
        return self.host

    @property
    def dest(self) -> str:
        """ssh/rsync destination: [user@]host."""
        return f"{self.user}@{self.host}" if self.user else self.host

    def ssh_argv(self, *remote_cmd: str) -> list[str]:
        return ["ssh", self.dest, *remote_cmd]

    def sync_argv(self, wt_name: str) -> list[str]:
        """rsync command pushing one worktree's sources to the device."""
        return ["rsync", "-az", "--delete",
                *(x for e in REMOTE_EXCLUDES for x in ("--exclude", e)),
                "-R", *REMOTE_SOURCES,
                f"{self.dest}:{REMOTE_ROOT}/{wt_name}/"]

    def build_script(self, wt_name: str) -> str:
        """Remote shell snippet configuring + building one arm's tree."""
        rdir = f"$HOME/{REMOTE_ROOT}/{wt_name}"
        return (f"cmake -S {rdir} -B {rdir}/build "
                f"{' '.join(BASE_CMAKE_FLAGS)} && "
                f"cmake --build {rdir}/build --target bench_pffft_float "
                f"bench_pffft_double -j$(nproc 2>/dev/null || echo 4)")

    def remote_host(self) -> str:
        """Actual hostname of the device (cached); falls back to spec."""
        if self._remote_name is None:
            try:
                self._remote_name = clean_value(
                    sh(*self.ssh_argv("uname", "-n")).strip()) or self.host
            except Exception:
                self._remote_name = self.host
        return self._remote_name

    def _prepare(self, wt_build: Path) -> None:
        """Sync sources and build remotely (once per worktree per process)."""
        wt_name = wt_build.parent.parent.name      # wt-base / wt-var
        if wt_name in self._prepared:
            return
        subprocess.run(self.ssh_argv(f"mkdir -p {REMOTE_ROOT}/{wt_name}"),
                       check=True)
        subprocess.run(self.sync_argv(wt_name),
                       cwd=wt_build.parent.parent, check=True)
        subprocess.run(self.ssh_argv(self.build_script(wt_name)), check=True)
        self._prepared.add(wt_name)

    def meta(self) -> dict[str, str]:
        return {"target": f"ssh://{self.dest}", "host": self.remote_host()}

    def binary_path(self, wt_build: Path, prec: str) -> Path:
        self._prepare(wt_build)
        exe = {"flt": "bench_pffft_float",
               "dbl": "bench_pffft_double"}[prec]
        # Path RELATIVE to the remote home dir: ssh lands the shell there,
        # matching how the rsync destination is resolved. Kept free of
        # $-expansions so every token can be shell-quoted in run().
        return (Path(REMOTE_ROOT) / wt_build.parent.parent.name
                / "build" / "benchmarks" / exe)

    def run(self, cmd: list[str],
            cwd: Path | None = None) -> Iterator[str]:
        """Run cmd remotely; capture --samples CSV from stdout to disk.

        Every token is shell-quoted before ssh hands the joined string to
        the remote shell, so meta values containing spaces survive intact.
        """
        cmd = list(cmd)
        samples_path: Path | None = None
        if "--samples" in cmd:
            i = cmd.index("--samples")
            if i + 1 >= len(cmd):
                raise ValueError("dangling --samples argument")
            samples_path = Path(cmd[i + 1])
            cmd[i + 1] = "-"       # binary streams pure CSV to stdout
        lines: list[str] = []
        for line in stream_stdout(self.ssh_argv(*(shlex.quote(t)
                                                  for t in cmd))):
            lines.append(line)
            yield line
        if samples_path is not None:
            samples_path.parent.mkdir(parents=True, exist_ok=True)
            if samples_path.exists():
                # perf.py reuses one file per (label, target, tree) across
                # invocations: append data rows only, never a second header.
                rows = [ln for ln in lines
                        if ln and not ln.startswith(("#", "algo,"))]
                with open(samples_path, "a") as fh:
                    fh.write("\n".join(rows) + ("\n" if rows else ""))
            else:
                samples_path.write_text(
                    "\n".join(lines) + ("\n" if lines else ""))


class AdbTarget(Target):
    """An Android device driven over adb (optionally a specific serial).

    Binaries are cross-compiled on this host with the Android NDK toolchain,
    reusing pffft_android_cmake_argv() from cross_build_android.py (same
    flags as that script's main()); the results are pushed to
    /data/local/tmp and run there with --samples - so stdout streams back
    and is captured into a local sample file, mirroring SshTarget.
    """

    DEVICE_DIR = "/data/local/tmp"

    def __init__(self, serial: str | None = None, abi: str = "arm64-v8a",
                 api: int = 24):
        self.serial = serial
        self.abi = abi
        self.api = api
        self._prepared: set[str] = set()

    @property
    def name(self) -> str:
        return f"adb:{self.serial}" if self.serial else "adb"

    def adb(self, *args: str) -> list[str]:
        argv = ["adb"]
        if self.serial:
            argv += ["-s", self.serial]
        return argv + list(args)

    def _getprop(self, prop: str) -> str:
        return sh(*self.adb("shell", "getprop", prop)).strip()

    def meta(self) -> dict[str, str]:
        host = " ".join(part for part in (
            self._getprop("ro.product.model"),
            self._getprop("ro.board.platform")) if part)
        return {"target": self.name, "host": clean_value(host)}

    def binary_path(self, wt_build: Path, prec: str) -> Path:
        self._prepare(wt_build)
        exe = {"flt": "bench_pffft_float",
               "dbl": "bench_pffft_double"}[prec]
        # On-device path; resolved by the device shell in run().
        return Path(self.DEVICE_DIR) / wt_build.parent.parent.name / exe

    def _prepare(self, wt_build: Path) -> None:
        """Cross-compile one worktree with the NDK and push the binaries."""
        wt_name = wt_build.parent.parent.name      # wt-base / wt-var
        if wt_name in self._prepared:
            return
        ndk = find_ndk()
        if ndk is None or not (ndk / "build" / "cmake"
                               / "android.toolchain.cmake").is_file():
            raise RuntimeError(
                "Android NDK not found; install it (e.g. under "
                "~/Library/Android/sdk/ndk) or export ANDROID_NDK")
        wt = wt_build.parent.parent
        bdir = wt / "build-android"
        march = detect_device_march(self.adb(), self.abi)
        if march is None:
            march = {"arm64-v8a": "armv8-a",
                     "armeabi-v7a": "armv7-a"}.get(self.abi, "none")
        if not (bdir / "CMakeCache.txt").exists():
            subprocess.run(pffft_android_cmake_argv(
                wt, bdir, ndk / "build" / "cmake" / "android.toolchain.cmake",
                self.abi, self.api, march), check=True)
        subprocess.run(["cmake", "--build", str(bdir), "--config", "Release",
                        "--", f"-j{os.cpu_count() or 4}"], check=True)
        devdir = f"{self.DEVICE_DIR}/{wt_name}"
        subprocess.run(self.adb("shell", f"mkdir -p {devdir}"), check=True)
        for exe in ("bench_pffft_float", "bench_pffft_double"):
            exe_path = bdir / "benchmarks" / exe
            if exe_path.is_file():
                subprocess.run(self.adb("push", str(exe_path), devdir),
                               check=True, capture_output=True)
                subprocess.run(self.adb("shell", f"chmod +x {devdir}/{exe}"),
                               check=True)
        self._prepared.add(wt_name)

    def device_shell_cmd(self, cmd: list[str]) -> str:
        """Quote a command list into one string for the device shell."""
        return " ".join(shlex.quote(str(a)) for a in cmd)

    def run(self, cmd: list[str],
            cwd: Path | None = None) -> Iterator[str]:
        """Run cmd on the device; capture --samples CSV from stdout."""
        cmd = list(cmd)
        samples_path: Path | None = None
        if "--samples" in cmd:
            i = cmd.index("--samples")
            if i + 1 >= len(cmd):
                raise ValueError("dangling --samples argument")
            samples_path = Path(cmd[i + 1])
            cmd[i + 1] = "-"       # binary streams pure CSV to stdout
        argv = self.adb("shell", self.device_shell_cmd(cmd))
        lines: list[str] = []
        for line in stream_stdout(argv):
            lines.append(line)
            yield line
        if samples_path is not None:
            samples_path.parent.mkdir(parents=True, exist_ok=True)
            if samples_path.exists():
                # Accumulate: the C writer appends per invocation, so an
                # existing file only receives fresh DATA rows — never a
                # second header, which would truncate or corrupt it.
                rows = [ln for ln in lines
                        if not ln.startswith("#") and not ln.startswith("algo,")]
                with open(samples_path, "a") as f:
                    f.write("\n".join(rows) + ("\n" if rows else ""))
            else:
                samples_path.write_text(
                    "\n".join(lines) + ("\n" if lines else ""))


def make_target(spec: str) -> Target:
    """Dispatch a target spec to a Target instance."""
    if spec == "local":
        return LocalTarget()
    if spec.startswith("ssh://"):
        rest = spec[len("ssh://"):].rstrip("/")
        user, sep, host = rest.rpartition("@")
        if not sep:
            user, host = None, rest
        if not host:
            raise ValueError(f"missing host in target spec: {spec!r}")
        return SshTarget(host, user)
    if spec == "adb":
        return AdbTarget()
    if spec.startswith("adb:"):
        serial = spec[len("adb:"):]
        if not serial:
            raise ValueError(f"missing serial in target spec: {spec!r}")
        return AdbTarget(serial)
    if spec.startswith("ios"):
        raise NotImplementedError("ios targets arrive in Task 10")
    raise ValueError(f"unknown target spec: {spec!r}")
