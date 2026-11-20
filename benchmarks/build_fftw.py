#!/usr/bin/env python3
"""Build FFTW 3.3.10 as static libraries, float and double, with SIMD.

    python3 benchmarks/build_fftw.py --prefix DIR
    cmake ... -DPFFFT_USE_BENCH_FFTW=ON -DFFTW3_ROOT=DIR

Package-manager builds of FFTW may lack SIMD (Homebrew's on Apple Silicon
plans with scalar codelets), which would understate FFTW in the benchmark.
Targets call build() with their own compiler, --host triplet and flags.
Needs a POSIX shell and make; on Windows, CMake and Visual Studio instead,
for x64 only: FFTW's CMake build, which Windows needs, has SSE2, AVX and
AVX2 but neither AVX-512 nor NEON.
"""
import argparse
import hashlib
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

VERSION = "3.3.10"
URL = f"https://www.fftw.org/fftw-{VERSION}.tar.gz"
# fftw.org publishes only an MD5 (8ccbf6a5ea78a16dbc3e1306e234cc5c, matched);
# this SHA-256 equals the one in conda-forge's recipe for 3.3.10.
SHA256 = "56c932549852cddcfafdab3820b0200c7742675be92179e59e6215b340e26467"
CACHE = Path(__file__).resolve().parent.parent / "build_bench_cache"
FLOAT_ONLY = {"--enable-sse", "--enable-altivec"}  # configure rejects them in double
# For every aarch64 target. CNTVCT_EL0 is FFTW's only timer on aarch64, Apple
# included (its mach_absolute_time branch is dead code); without a timer
# FFTW_MEASURE silently plans from estimates, and plans worse than
# FFTW_ESTIMATE. User space may read CNTVCT_EL0 on Linux, Android and Apple.
ARM64_ARGS = ["--enable-neon", "--enable-armv8-cntvct-el0"]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch():
    path = CACHE / f"fftw-{VERSION}.tar.gz"
    if path.is_file() and sha256(path) == SHA256:
        return path
    CACHE.mkdir(exist_ok=True)
    part = path.with_name(f"{path.name}.{os.getpid()}")
    print(f"build_fftw: downloading {URL}", file=sys.stderr)
    try:
        urllib.request.urlretrieve(URL, part)
        if sha256(part) != SHA256:
            sys.exit(f"build_fftw: {URL} does not have SHA-256 {SHA256}")
        part.replace(path)
    finally:
        part.unlink(missing_ok=True)  # after a failed or interrupted download
    return path


def run(cmd, cwd, env, log):
    with open(log, "a") as f:
        if subprocess.run([str(c) for c in cmd], cwd=cwd, env=env, stdout=f,
                          stderr=subprocess.STDOUT).returncode:
            tail = Path(log).read_text(errors="replace").splitlines()[-40:]
            sys.exit("\n".join(tail) + f"\nbuild_fftw: '{' '.join(map(str, cmd))}' failed, see {log}")


def unpack(tarball, into):
    with tarfile.open(tarball) as t:
        t.extractall(into, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))


# configure's SIMD flags and the options of FFTW's CMake build that do the same
CMAKE_SIMD = {"--enable-sse2": "ENABLE_SSE2", "--enable-avx": "ENABLE_AVX", "--enable-avx2": "ENABLE_AVX2"}


def build_cmake(prefix, configure_args):
    """build() on Windows, where configure cannot run: FFTW's own CMake build
    with the default generator (Visual Studio, MSVC). It knows SSE2, AVX and
    AVX2 only; other SIMD flags are dropped, and on ARM64, where it has no
    NEON, it refuses, as FFTW would be scalar there."""
    if platform.machine().lower() in ("arm64", "aarch64"):
        sys.exit("build_fftw: FFTW's CMake build, the one Windows can use, has no NEON: "
                 "FFTW would run scalar on Windows ARM64, so it is not built there")
    simd = [f"-D{CMAKE_SIMD[a]}=ON" for a in configure_args if a in CMAKE_SIMD]
    dropped = [a for a in configure_args if a not in CMAKE_SIMD]
    common = [f"-DCMAKE_INSTALL_PREFIX={prefix}", "-DBUILD_SHARED_LIBS=OFF", "-DBUILD_TESTS=OFF",
              "-DDISABLE_FORTRAN=ON", "-DCMAKE_BUILD_TYPE=Release",
              "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"] + simd  # FFTW asks for CMake 3.0
    args = {"float": common + ["-DENABLE_FLOAT=ON"], "double": common}
    stamp = prefix / "build_fftw.txt"
    recipe = "".join(f"cmake {' '.join(a)}\n" for a in args.values())
    if stamp.is_file() and stamp.read_text() == recipe:
        return prefix
    if dropped:
        print(f"build_fftw: FFTW's CMake build has no option for {' '.join(dropped)}: left out",
              file=sys.stderr)
    work = prefix / "build"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    log = work / "log.txt"
    unpack(fetch(), work)
    for prec, a in args.items():
        print(f"build_fftw: {prec}: cmake {' '.join(a)}", file=sys.stderr)
        run(["cmake", "-S", work / f"fftw-{VERSION}", "-B", work / prec, *a], work, None, log)
        run(["cmake", "--build", work / prec, "--config", "Release", "--target", "install",
             "--parallel", str(os.cpu_count() or 4)], work, None, log)
    shutil.rmtree(work)
    stamp.write_text(recipe)
    return prefix


def build(prefix, *, cc="cc", cflags="-O3", host=None, configure_args=()):
    """Build libfftw3f.a and libfftw3.a and fftw3.h into prefix (lib/, include/)
    unless they are already there; return prefix. host is the autoconf --host
    triplet for cross builds, None for native. configure_args carries the SIMD
    flags, chosen by the caller (the target). On Windows (native only) FFTW's
    CMake build makes fftw3f.lib and fftw3.lib instead; cc and cflags do not
    apply there."""
    prefix = Path(prefix).resolve()
    if os.name == "nt" and not host:
        return build_cmake(prefix, configure_args)
    common = [f"--prefix={prefix}", "--disable-shared", "--disable-fortran",
              "--disable-doc"] + ([f"--host={host}"] if host else [])
    args = {"float": common + ["--enable-float", *configure_args],
            "double": common + [a for a in configure_args if a not in FLOAT_ONLY]}
    # A prefix built with other settings is rebuilt, never measured unnoticed.
    stamp = prefix / "build_fftw.txt"
    recipe = "".join(f"CC={cc} CFLAGS={cflags} configure {' '.join(a)}\n"
                     for a in args.values())
    if stamp.is_file() and stamp.read_text() == recipe:
        return prefix
    env = dict(os.environ, CC=cc, CFLAGS=cflags)
    # The NDK has no <host>-ar, and the build machine's ar cannot index ELF
    # objects on macOS: use the archiver that comes with the compiler.
    tools = Path(shutil.which(cc) or cc).parent
    for var, tool in (("AR", "llvm-ar"), ("RANLIB", "llvm-ranlib")):
        if host and (tools / tool).is_file():
            env[var] = str(tools / tool)
    work = prefix / "build"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    log = work / "log.txt"
    unpack(fetch(), work)
    for prec, a in args.items():
        print(f"build_fftw: {prec}: configure {' '.join(a)}", file=sys.stderr)
        (work / prec).mkdir()
        run([work / f"fftw-{VERSION}" / "configure", *a], work / prec, env, log)
        run(["make", f"-j{os.cpu_count() or 4}", "install"], work / prec, env, log)
    shutil.rmtree(work)
    stamp.write_text(recipe)
    return prefix


def native_simd_args():
    """SIMD configure flags for the machine this runs on."""
    machine = platform.machine().lower()
    if machine in ("arm64", "aarch64"):
        # FFTW's config.guess of 2018 calls Apple Silicon "arm", and configure
        # allows NEON in double precision only on "aarch64".
        return ARM64_ARGS + ["--build=aarch64-apple-darwin"] * (sys.platform == "darwin")
    if machine in ("x86_64", "amd64"):
        # FFTW checks the CPU at run time before it uses AVX or AVX2 codelets.
        args = ["--enable-sse2", "--enable-avx", "--enable-avx2"]
        if sys.platform == "darwin":
            avx512 = subprocess.run(["sysctl", "-n", "hw.optional.avx512f"],
                                    capture_output=True, text=True).stdout.strip() == "1"
        elif os.name == "nt":
            avx512 = False  # FFTW's CMake build, the one Windows can use, has no AVX-512
        else:
            cpuinfo = Path("/proc/cpuinfo")
            avx512 = cpuinfo.is_file() and "avx512f" in cpuinfo.read_text().split()
        return args + ["--enable-avx512"] * avx512
    return []


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--prefix", required=True, type=Path)
    print(build(parser.parse_args().prefix, configure_args=native_simd_args()))
