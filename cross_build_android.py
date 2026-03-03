#!/usr/bin/env python3
"""Build pffft benchmarks for Android and optionally run them on a connected device.

Usage:
    python cross_build_android.py <id> [options] [extra cmake options]

Arguments:
    id              Short identifier used in directory names (e.g. "arm64-clang")

Options:
    --abi <ABI>     Android ABI to build for (default: arm64-v8a)
                    Other options: armeabi-v7a, x86_64, x86
    --api <level>   Android API level (default: 21)
    --arch <march>  -march value passed to TARGET_C_ARCH / TARGET_CXX_ARCH.
                    If omitted and a device is connected, auto-detected from
                    /proc/cpuinfo (e.g. armv8.2-a+dotprod+fp16).
                    Falls back to armv8-a / armv7-a for arm ABIs when no
                    device is available.
    --fftw          Cross-compile FFTW 3.3.10 and include it in the benchmark
                    (requires autoconf/make; not supported on Windows)
    --no-run        Build only; do not push or run on device
    --output-dir    Local directory to store pulled benchmark CSVs
                    (default: bench_results_android_<id>)
    --ndk <path>    Override NDK root directory
    --serial <s>    ADB device serial (passed as -s <s> to adb)

Environment variables:
    ANDROID_NDK     NDK root (overrides auto-detection)
    ANDROID_SDK     SDK root for auto-detecting NDK

Examples:
    python cross_build_android.py arm64
    python cross_build_android.py arm64 --fftw
    python cross_build_android.py arm64 --no-run
    python cross_build_android.py arm64 --output-dir my_results -DPFFFT_USE_BENCH_GREEN=OFF
"""

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from datetime import date
from pathlib import Path


FFTW_VERSION = "3.3.10"
FFTW_URL = f"https://www.fftw.org/fftw-{FFTW_VERSION}.tar.gz"

# autoconf --host triple and extra configure flags per ABI.
# --enable-armv8-cntvct-el0: read the ARM virtual counter directly from
#   userspace (CNTVCT_EL0 register) for high-resolution plan timing.
#   Without this FFTW falls back to ESTIMATE mode (no measurement at all)
#   or a slow gettimeofday timer. Accessible on Android since Linux 4.x.
#   Discovered by andrej5elin: https://github.com/andrej5elin/howto_fftw_apple_silicon
# --enable-neon: enable NEON SIMD codelets (float and double).
# Tuple: (host_triple, float_extra_flags, double_extra_flags)
# FFTW NEON codelets are single-precision only — --enable-neon must not be
# passed to the double-precision build or configure will error out.
ABI_CONFIGURE = {
    "arm64-v8a":   ("aarch64-linux-android",
                    ["--enable-armv8-cntvct-el0", "--enable-neon"],  # float
                    ["--enable-armv8-cntvct-el0"]),                   # double
    "armeabi-v7a": ("armv7a-linux-androideabi",
                    ["--with-slow-timer", "--enable-neon"],
                    ["--with-slow-timer"]),
    "x86_64":      ("x86_64-linux-android",
                    ["--with-slow-timer"], ["--with-slow-timer"]),
    "x86":         ("i686-linux-android",
                    ["--with-slow-timer"], ["--with-slow-timer"]),
}


# ── helpers ───────────────────────────────────────────────────────────────────

def banner(text):
    line = "═" * 56
    print(f"\n{line}")
    for t in text:
        print(f"  {t}")
    print(line)


def run(cmd, check=True, **kwargs):
    """Run a command, printing it first, streaming output live."""
    print("$", " ".join(str(c) for c in cmd))
    return subprocess.run(cmd, check=check, **kwargs)


def adb_output(adb_cmd, shell_cmd):
    """Run an adb shell command and return stdout stripped of CR/LF."""
    result = subprocess.run(
        adb_cmd + ["shell", shell_cmd],
        capture_output=True, text=True,
    )
    return result.stdout.strip().replace("\r", "")


# ── device capability detection ───────────────────────────────────────────────

def detect_device_march(adb_cmd, abi):
    """Inspect /proc/cpuinfo on the connected device and return a -march string
    that uses every ISA extension the device actually supports.

    Only meaningful for ARM ABIs; returns None for x86*.
    """
    if not abi.startswith("arm"):
        return None

    result = subprocess.run(
        adb_cmd + ["shell", "grep -m1 Features /proc/cpuinfo"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        return None

    feats = result.stdout.lower()

    # Feature flag → clang -march extension name
    # Listed roughly in ARMv8.x generation order.
    ext_map = [
        ("asimdrdm", "rdm"),        # ARMv8.1 NEON rounding-doubles
        ("asimddp",  "dotprod"),    # ARMv8.2 NEON 8-bit dot product
        ("asimdhp",  "fp16"),       # ARMv8.2 NEON half-precision arithmetic
        ("lrcpc",    "rcpc"),       # ARMv8.3 RCpc load-acquire
    ]

    # Determine base march from highest confirmed generation
    if "asimddp" in feats or "asimdhp" in feats:
        base = "armv8.2-a"
    elif "asimdrdm" in feats:
        base = "armv8.1-a"
    else:
        base = "armv8-a"

    extras = [ext for flag, ext in ext_map if flag in feats]
    march = base + "".join(f"+{e}" for e in extras)
    print(f"  Detected device march: {march}  (features: {feats.split('features')[1].strip().lstrip(': ')})")
    return march


# ── NDK detection ─────────────────────────────────────────────────────────────

def find_ndk():
    """Return the path to the newest NDK installation, or None."""
    if os.environ.get("ANDROID_NDK"):
        return Path(os.environ["ANDROID_NDK"])

    home = Path.home()
    sdk_candidates = [
        Path(os.environ["ANDROID_SDK"]) if os.environ.get("ANDROID_SDK") else None,
        home / "Library" / "Android" / "sdk",          # macOS
        home / "Android" / "Sdk",                       # Linux
        Path("C:/Users") / os.environ.get("USERNAME", "") / "AppData" / "Local" / "Android" / "Sdk",  # Windows
        Path("/opt/android-sdk"),
    ]

    for sdk in sdk_candidates:
        if sdk is None or not sdk.is_dir():
            continue
        ndk_dir = sdk / "ndk"
        if not ndk_dir.is_dir():
            continue
        versions = sorted(ndk_dir.iterdir(), key=lambda p: [
            int(x) if x.isdigit() else x for x in p.name.split(".")
        ])
        if versions:
            return versions[-1]

    return None


def ndk_host_tag():
    """Return the NDK prebuilt host directory name for the current OS."""
    system = platform.system()
    if system == "Darwin":
        return "darwin-x86_64"   # NDK ships x86_64 binaries even on Apple Silicon
    if system == "Linux":
        return "linux-x86_64"
    if system == "Windows":
        return "windows-x86_64"
    raise RuntimeError(f"Unsupported host OS: {system}")


# ── FFTW cross-compilation ────────────────────────────────────────────────────

def build_fftw(ndk_root, abi, api, march, build_dir, cpu_count):
    """Download and cross-compile FFTW (float + double) for the given ABI.

    Returns the path to the install prefix containing include/ and lib/.
    Not supported on Windows (requires autoconf/make).
    """
    if platform.system() == "Windows":
        print("ERROR: --fftw cross-compilation requires autoconf/make and is not supported on Windows.", file=sys.stderr)
        print("  Use WSL or a Linux/macOS host instead.", file=sys.stderr)
        sys.exit(1)

    abi_cfg = ABI_CONFIGURE.get(abi)
    if abi_cfg is None:
        print(f"ERROR: --fftw: no autoconf host triple known for ABI '{abi}'", file=sys.stderr)
        sys.exit(1)
    host_triple, float_flags, double_flags = abi_cfg

    toolbin = ndk_root / "toolchains" / "llvm" / "prebuilt" / ndk_host_tag() / "bin"
    cc      = str(toolbin / f"{host_triple}{api}-clang")
    ar      = str(toolbin / "llvm-ar")
    ranlib  = str(toolbin / "llvm-ranlib")

    if not Path(cc).is_file():
        print(f"ERROR: NDK clang not found: {cc}", file=sys.stderr)
        print(f"  Try a higher --api value or a newer NDK.", file=sys.stderr)
        sys.exit(1)

    fftw_src_dir  = build_dir / "fftw_src"
    fftw_prefix   = build_dir / "fftw_prefix"
    tarball       = fftw_src_dir / f"fftw-{FFTW_VERSION}.tar.gz"
    fftw_src_dir.mkdir(parents=True, exist_ok=True)
    fftw_prefix.mkdir(parents=True, exist_ok=True)

    # Download
    if not tarball.is_file():
        banner([f"Downloading FFTW {FFTW_VERSION}"])
        print(f"  {FFTW_URL}")
        urllib.request.urlretrieve(FFTW_URL, tarball)
        print(f"  Saved to {tarball}")

    # Extract (idempotent)
    fftw_tree = fftw_src_dir / f"fftw-{FFTW_VERSION}"
    if not fftw_tree.is_dir():
        print(f"Extracting {tarball.name} ...")
        with tarfile.open(tarball, "r:gz") as tf:
            tf.extractall(fftw_src_dir)

    base_configure_args = [
        str(fftw_tree / "configure"),
        f"--host={host_triple}",
        "--disable-shared",
        "--enable-static",
        f"--prefix={fftw_prefix}",
    ]

    env = os.environ.copy()
    # Override FFTW's default CFLAGS which include -mtune=native; that flag
    # targets the *host* machine when cross-compiling, not the Android device.
    # We use -march instead, which controls both code generation and scheduling.
    march_flag = f"-march={march}" if march and march != "none" else ""
    env.update({
        "CC":     cc,
        "AR":     ar,
        "RANLIB": ranlib,
        "CFLAGS": f"-O3 -fomit-frame-pointer -fstrict-aliasing {march_flag}".strip(),
    })

    # Build float (fftw3f) and double (fftw3) from separate build dirs.
    # NEON codelets are float-only — float_flags and double_flags differ.
    for label, prec_flags, extra in [
        ("float",  float_flags,  ["--enable-float"]),
        ("double", double_flags, []),
    ]:
        banner([f"Building FFTW {FFTW_VERSION} — {label} precision"])
        bdir = build_dir / f"fftw_build_{label}"
        bdir.mkdir(parents=True, exist_ok=True)
        run(base_configure_args + prec_flags + extra, cwd=bdir, env=env)
        run(["make", "-j", str(cpu_count)], cwd=bdir)
        run(["make", "install"], cwd=bdir)

    return fftw_prefix


# ── argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("id", help="Short build identifier (used in directory names)")
    parser.add_argument("--abi", default="arm64-v8a",
                        help="Android ABI (default: arm64-v8a)")
    parser.add_argument("--api", default="21",
                        help="Android API level (default: 21)")
    parser.add_argument("--arch", default=None,
                        help="-march value for TARGET_C/CXX_ARCH (auto-detected from ABI if omitted)")
    parser.add_argument("--fftw", action="store_true",
                        help="Cross-compile FFTW and include it in the benchmark")
    parser.add_argument("--no-run", action="store_true",
                        help="Build only; skip push and run on device")
    parser.add_argument("--output-dir", default=None,
                        help="Local directory for pulled CSV results")
    parser.add_argument("--ndk", default=None,
                        help="Override NDK root directory")
    parser.add_argument("--serial", default=None,
                        help="ADB device serial number")

    args, extra_cmake = parser.parse_known_args()
    return args, extra_cmake


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args, extra_cmake = parse_args()

    script_dir = Path(__file__).resolve().parent

    # ── resolve NDK ───────────────────────────────────────────────────────────

    ndk_root = Path(args.ndk) if args.ndk else find_ndk()

    if ndk_root is None or not (ndk_root / "build" / "cmake" / "android.toolchain.cmake").is_file():
        print("ERROR: Android NDK not found. Set ANDROID_NDK or pass --ndk <path>.", file=sys.stderr)
        sys.exit(1)

    toolchain = ndk_root / "build" / "cmake" / "android.toolchain.cmake"
    print(f"Using NDK: {ndk_root}")

    # ── resolve adb command ───────────────────────────────────────────────────
    # (needed early for march auto-detection)

    adb_base = ["adb"]
    if args.serial:
        adb_base = ["adb", "-s", args.serial]

    # ── resolve march ─────────────────────────────────────────────────────────

    march = args.arch
    if march is None and not args.no_run:
        # Try to auto-detect from the connected device before building
        result = subprocess.run(adb_base + ["devices"], capture_output=True, text=True)
        device_available = any(
            line.endswith("device")
            for line in result.stdout.splitlines()[1:]
            if line.strip()
        )
        if device_available:
            print("\nAuto-detecting device CPU capabilities...")
            march = detect_device_march(adb_base, args.abi)
    if march is None:
        march = {"arm64-v8a": "armv8-a", "armeabi-v7a": "armv7-a"}.get(args.abi, "none")
        print(f"  Using fallback march: {march}")

    # ── directories ───────────────────────────────────────────────────────────

    build_dir  = script_dir / f"build_android_{args.id}"
    output_dir = Path(args.output_dir) if args.output_dir else script_dir / f"bench_results_android_{args.id}"
    device_dir = f"/data/local/tmp/pffft_bench_{args.id}"
    device_out = f"{device_dir}/results"

    cpu_count = os.cpu_count() or 4

    # ── optionally cross-compile FFTW ─────────────────────────────────────────

    fftw_cmake_args = []
    if args.fftw:
        fftw_prefix = build_fftw(ndk_root, args.abi, args.api, march, build_dir, cpu_count)
        fftw_cmake_args = [
            "-DPFFFT_USE_BENCH_FFTW=ON",
            f"-DFFTW3_ROOT={fftw_prefix}",
        ]
    else:
        fftw_cmake_args = ["-DPFFFT_USE_BENCH_FFTW=OFF"]

    # ── build pffft + benchmarks ──────────────────────────────────────────────

    banner([
        "Building pffft for Android",
        f"ABI={args.abi}  API={args.api}  march={march}  id={args.id}",
        f"FFTW={'yes' if args.fftw else 'no'}",
    ])

    # Preserve the FFTW build dirs; wipe only the pffft cmake build dir
    pffft_build_dir = build_dir / "pffft"
    if pffft_build_dir.exists():
        shutil.rmtree(pffft_build_dir)
    pffft_build_dir.mkdir(parents=True)

    run([
        "cmake", "-S", str(script_dir), "-B", str(pffft_build_dir),
        f"-DCMAKE_TOOLCHAIN_FILE={toolchain}",
        f"-DANDROID_ABI={args.abi}",
        f"-DANDROID_PLATFORM=android-{args.api}",
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DTARGET_C_ARCH={march}",
        f"-DTARGET_CXX_ARCH={march}",
        "-DPFFFT_BUILD_TESTS=OFF",
        "-DPFFFT_BUILD_BENCHMARKS=ON",
        "-DPFFFT_BUILD_EXAMPLES=OFF",
        "-DPFFFT_USE_BENCH_MKL=OFF",
        "-DPFFFT_USE_BENCH_FFTS=OFF",
    ] + fftw_cmake_args + extra_cmake)

    run([
        "cmake", "--build", str(pffft_build_dir),
        "--config", "Release",
        "--", f"-j{cpu_count}",
    ])

    print(f"\nBuild complete: {pffft_build_dir}")

    if args.no_run:
        print("Skipping device run (--no-run).")
        return

    # ── check device ──────────────────────────────────────────────────────────

    result = subprocess.run(adb_base + ["devices"], capture_output=True, text=True)
    connected = any(
        line.endswith("device")
        for line in result.stdout.splitlines()[1:]
        if line.strip()
    )
    if not connected:
        print("ERROR: No Android device connected (or adb not authorized).", file=sys.stderr)
        print("  Connect a device and run: adb devices", file=sys.stderr)
        sys.exit(1)

    # ── push binaries ─────────────────────────────────────────────────────────

    banner(["Pushing binaries to device"])

    run(adb_base + ["shell", f"mkdir -p {device_dir}"])

    pushed = []
    for exe in ("bench_pffft_float", "bench_pffft_double"):
        exe_path = pffft_build_dir / "benchmarks" / exe
        if not exe_path.is_file():
            print(f"WARNING: {exe} not found (was it built?)", file=sys.stderr)
            continue
        run(adb_base + ["push", str(exe_path), f"{device_dir}/{exe}"])
        run(adb_base + ["shell", f"chmod +x {device_dir}/{exe}"])
        print(f"Pushed: {exe}")
        pushed.append(exe)

    # ── collect device info ───────────────────────────────────────────────────

    output_dir.mkdir(parents=True, exist_ok=True)

    print("\nCollecting device info...")
    model       = adb_output(adb_base, "getprop ro.product.model")
    android_ver = adb_output(adb_base, "getprop ro.build.version.release")
    api_level   = adb_output(adb_base, "getprop ro.build.version.sdk")
    abi_list    = adb_output(adb_base, "getprop ro.product.cpu.abilist")
    soc         = adb_output(adb_base, "getprop ro.hardware")
    ndk_ver     = ndk_root.name

    device_info_lines = [
        f"Device: {model}",
        f"Android: {android_ver}  (API {api_level})",
        f"ABI: {abi_list}",
        f"SOC: {soc}",
        f"Arch: {args.abi}  march={march}",
        f"NDK: {ndk_ver}",
        f"FFTW: {FFTW_VERSION if args.fftw else 'no'}",
        f"Date: {date.today().isoformat()}",
    ]
    for line in device_info_lines:
        print(line)
    (output_dir / "device_info.txt").write_text("\n".join(device_info_lines) + "\n")

    info = {
        "cpu":      f"{model} ({soc})",
        "arch":     args.abi,
        "os":       f"Android {android_ver}",
        "compiler": f"NDK {ndk_ver} clang",
        "date":     date.today().isoformat(),
    }
    info_path = output_dir / "info.json"
    info_path.write_text(json.dumps(info, indent=2) + "\n")
    print(f"Wrote {info_path}")

    # ── run benchmarks ────────────────────────────────────────────────────────

    banner(["Running benchmarks on device"])

    run(adb_base + ["shell", f"mkdir -p {device_out}"])

    for exe in pushed:
        print(f"\n── Running {exe} ──")
        result = run(
            adb_base + ["shell", f"cd {device_dir} && ./{exe} --output-dir {device_out}"],
            check=False,
        )
        if result.returncode != 0:
            print(f"WARNING: {exe} exited with status {result.returncode}", file=sys.stderr)

    # ── pull results ──────────────────────────────────────────────────────────

    banner([f"Pulling results to {output_dir}"])

    ls_result = subprocess.run(
        adb_base + ["shell", f"ls {device_out}/*.csv 2>/dev/null"],
        capture_output=True, text=True,
    )
    csv_files = [p.strip().replace("\r", "") for p in ls_result.stdout.splitlines() if p.strip()]

    for csv in csv_files:
        run(adb_base + ["pull", csv, str(output_dir) + "/"])
        print(f"Pulled: {Path(csv).name}")

    # ── clean up device ───────────────────────────────────────────────────────

    subprocess.run(adb_base + ["shell", f"rm -rf {device_dir}"], check=False)

    print(f"\nDone. Results in: {output_dir}")
    print(f"\nTo generate charts:")
    print(f'  python3 bench/make_charts.py "{output_dir}"')


if __name__ == "__main__":
    main()
