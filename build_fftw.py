#!/usr/bin/env python3
"""Build FFTW 3 (float + double) from source with optimal flags for the host CPU.

Produces a local install prefix containing include/fftw3.h and
lib/libfftw3{f,}.a, which can be passed to the pffft benchmark build as
  cmake -DPFFFT_USE_BENCH_FFTW=ON -DFFTW3_ROOT=<prefix> ...

Usage:
    python3 build_fftw.py [--prefix <dir>] [--version <x.y.z>]

Options:
    --prefix <dir>   Install destination (default: ./fftw_build/prefix)
    --version <ver>  FFTW version to download (default: 3.3.10)

Key configure flags chosen automatically per platform:
  ARM64 (Apple Silicon, Linux arm64):
    --enable-neon
    --enable-armv8-cntvct-el0  (direct hardware cycle counter; discovered by
                                 andrej5elin — see
                                 github.com/andrej5elin/howto_fftw_apple_silicon)
  x86-64:
    --enable-sse2 / --enable-avx / --enable-avx2 / --enable-avx512
    based on what the CPU actually supports

CFLAGS uses -march=native so the compiler can also exploit any other
micro-architectural details of the host machine.

Not needed on Android — use cross_build_android.py --fftw instead.
"""

import argparse
import os
import platform
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path


FFTW_DEFAULT_VERSION = "3.3.10"
FFTW_URL_TEMPLATE = "https://www.fftw.org/fftw-{version}.tar.gz"


# ── helpers ───────────────────────────────────────────────────────────────────

def run(cmd, **kwargs):
    print("$", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kwargs)


def banner(text):
    line = "─" * 56
    print(f"\n{line}")
    for t in text:
        print(f"  {t}")
    print(line)


# ── CPU feature detection ─────────────────────────────────────────────────────

def detect_fftw_flags():
    """Return (float_flags, double_flags) for optimal FFTW ./configure.

    NEON codelets in FFTW are single-precision only, so --enable-neon must
    not be passed to the double-precision build.
    """
    machine = platform.machine().lower()

    if machine in ("arm64", "aarch64"):
        # --enable-armv8-cntvct-el0: use the ARM virtual counter (CNTVCT_EL0)
        # directly from userspace for high-resolution plan timing. Without it
        # FFTW cannot find a cycle counter and falls back to ESTIMATE mode
        # (no measurement), producing poor plans.
        # Discovered by andrej5elin:
        #   https://github.com/andrej5elin/howto_fftw_apple_silicon
        shared = ["--enable-armv8-cntvct-el0"]
        return shared + ["--enable-neon"], shared   # neon: float only

    if machine in ("x86_64", "amd64"):
        flags = _x86_fftw_flags()
        return flags, flags   # all x86 SIMD flags apply to both precisions

    print(f"WARNING: unrecognised machine '{machine}', using no SIMD flags",
          file=sys.stderr)
    return [], []


def _x86_fftw_flags():
    """Detect x86 SIMD capabilities and return corresponding FFTW flags."""
    flags = ["--enable-sse2"]   # baseline for x86-64

    # Try to read /proc/cpuinfo (Linux) or use sysctl (macOS)
    cpu_features = _read_x86_features()

    ordered = [
        ("avx",     "--enable-avx"),
        ("avx2",    "--enable-avx2"),
        ("avx512f", "--enable-avx512"),
    ]
    for feat, flag in ordered:
        if feat in cpu_features:
            flags.append(flag)

    return flags


def _read_x86_features():
    """Return a set of lower-case CPU feature strings for the host x86 CPU."""
    # Linux
    proc = Path("/proc/cpuinfo")
    if proc.exists():
        text = proc.read_text()
        for line in text.splitlines():
            if line.startswith("flags"):
                return set(line.split(":")[1].lower().split())
        return set()

    # macOS
    if platform.system() == "Darwin":
        feats = set()
        checks = {
            "hw.optional.avx1_0": "avx",
            "hw.optional.avx2_0": "avx2",
            "hw.optional.avx512f": "avx512f",
        }
        for key, feat in checks.items():
            r = subprocess.run(["sysctl", "-n", key],
                               capture_output=True, text=True)
            if r.stdout.strip() == "1":
                feats.add(feat)
        return feats

    return set()


# ── download + extract ────────────────────────────────────────────────────────

def fetch_fftw(version, src_dir):
    """Download and extract fftw-<version>.tar.gz into src_dir if needed."""
    url = FFTW_URL_TEMPLATE.format(version=version)
    tarball = src_dir / f"fftw-{version}.tar.gz"
    tree    = src_dir / f"fftw-{version}"

    src_dir.mkdir(parents=True, exist_ok=True)

    if not tarball.is_file():
        banner([f"Downloading FFTW {version}", url])
        urllib.request.urlretrieve(url, tarball)
        print(f"  Saved: {tarball}")

    if not tree.is_dir():
        print(f"Extracting {tarball.name} ...")
        with tarfile.open(tarball, "r:gz") as tf:
            tf.extractall(src_dir)

    return tree


# ── build ─────────────────────────────────────────────────────────────────────

def build_fftw(version, prefix, work_dir):
    src_dir  = work_dir / "src"
    fftw_src = fetch_fftw(version, src_dir)
    prefix.mkdir(parents=True, exist_ok=True)

    float_flags, double_flags = detect_fftw_flags()
    cpu_count = os.cpu_count() or 4

    # -march=native lets the compiler exploit every host micro-arch detail
    env = os.environ.copy()
    env["CFLAGS"] = "-O3 -fomit-frame-pointer -fstrict-aliasing -march=native"

    configure = str(fftw_src / "configure")
    base_args = [
        configure,
        "--disable-shared",
        "--enable-static",
        f"--prefix={prefix}",
    ]

    for label, prec_flags, extra in [
        ("float",  float_flags,  ["--enable-float"]),
        ("double", double_flags, []),
    ]:
        banner([f"Building FFTW {version} — {label} precision",
                f"flags: {' '.join(prec_flags)}"])
        bdir = work_dir / f"build_{label}"
        bdir.mkdir(parents=True, exist_ok=True)
        run(base_args + prec_flags + extra, cwd=bdir, env=env)
        run(["make", "-j", str(cpu_count)], cwd=bdir)
        run(["make", "install"], cwd=bdir)

    print(f"\nFFTW {version} installed to: {prefix}")
    return prefix


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--prefix", default=None,
                        help="Install prefix (default: ./fftw_build/prefix)")
    parser.add_argument("--version", default=FFTW_DEFAULT_VERSION,
                        help=f"FFTW version (default: {FFTW_DEFAULT_VERSION})")
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    work_dir   = script_dir / "fftw_build"
    prefix     = Path(args.prefix) if args.prefix else work_dir / "prefix"

    float_flags, double_flags = detect_fftw_flags()
    print(f"Host:    {platform.system()} {platform.machine()}")
    print(f"Prefix:  {prefix}")
    print(f"Float:   {' '.join(float_flags)}")
    print(f"Double:  {' '.join(double_flags)}")

    if platform.system() == "Windows":
        print("ERROR: build_fftw.py requires make/autoconf and does not support "
              "Windows natively.\n  Use WSL or a Linux/macOS host.", file=sys.stderr)
        sys.exit(1)

    build_fftw(args.version, prefix, work_dir)

    print("\nTo use with the pffft benchmark:")
    print(f"  cmake ... -DPFFFT_USE_BENCH_FFTW=ON -DFFTW3_ROOT={prefix}")


if __name__ == "__main__":
    main()
