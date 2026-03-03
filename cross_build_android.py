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
    --arch <march>  -march value passed to TARGET_C_ARCH / TARGET_CXX_ARCH
                    (default: armv8-a for arm64-v8a, armv7-a for armeabi-v7a, none for x86*)
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
    python cross_build_android.py arm64 --no-run
    python cross_build_android.py arm64 --output-dir my_results -DPFFFT_USE_BENCH_GREEN=OFF
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path


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


# ── NDK detection ─────────────────────────────────────────────────────────────

def find_ndk():
    """Return the path to the newest NDK installation, or None."""
    # 1. Explicit env var
    if os.environ.get("ANDROID_NDK"):
        return Path(os.environ["ANDROID_NDK"])

    # 2. Well-known SDK locations
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


# ── argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        # Let unknown args pass through as extra cmake flags
        add_help=True,
    )
    parser.add_argument("id", help="Short build identifier (used in directory names)")
    parser.add_argument("--abi", default="arm64-v8a",
                        help="Android ABI (default: arm64-v8a)")
    parser.add_argument("--api", default="21",
                        help="Android API level (default: 21)")
    parser.add_argument("--arch", default=None,
                        help="-march value for TARGET_C/CXX_ARCH (auto-detected from ABI if omitted)")
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

    # ── resolve march default ─────────────────────────────────────────────────

    march = args.arch
    if march is None:
        march = {"arm64-v8a": "armv8-a", "armeabi-v7a": "armv7-a"}.get(args.abi, "none")

    # ── directories ───────────────────────────────────────────────────────────

    build_dir  = script_dir / f"build_android_{args.id}"
    output_dir = Path(args.output_dir) if args.output_dir else script_dir / f"bench_results_android_{args.id}"
    device_dir = f"/data/local/tmp/pffft_bench_{args.id}"
    device_out = f"{device_dir}/results"

    # ── build ─────────────────────────────────────────────────────────────────

    banner([
        "Building pffft for Android",
        f"ABI={args.abi}  API={args.api}  march={march}  id={args.id}",
    ])

    if build_dir.exists():
        shutil.rmtree(build_dir)
    build_dir.mkdir(parents=True)

    cpu_count = os.cpu_count() or 4

    run([
        "cmake", "-S", str(script_dir), "-B", str(build_dir),
        f"-DCMAKE_TOOLCHAIN_FILE={toolchain}",
        f"-DANDROID_ABI={args.abi}",
        f"-DANDROID_PLATFORM=android-{args.api}",
        "-DCMAKE_BUILD_TYPE=Release",
        f"-DTARGET_C_ARCH={march}",
        f"-DTARGET_CXX_ARCH={march}",
        "-DPFFFT_BUILD_TESTS=OFF",
        "-DPFFFT_BUILD_BENCHMARKS=ON",
        "-DPFFFT_BUILD_EXAMPLES=OFF",
        "-DPFFFT_USE_BENCH_FFTW=OFF",
        "-DPFFFT_USE_BENCH_MKL=OFF",
        "-DPFFFT_USE_BENCH_FFTS=OFF",
    ] + extra_cmake)

    run([
        "cmake", "--build", str(build_dir),
        "--config", "Release",
        "--", f"-j{cpu_count}",
    ])

    print(f"\nBuild complete: {build_dir}")

    if args.no_run:
        print("Skipping device run (--no-run).")
        return

    # ── check device ──────────────────────────────────────────────────────────

    adb_base = ["adb"]
    if args.serial:
        adb_base = ["adb", "-s", args.serial]

    result = subprocess.run(adb_base + ["devices"], capture_output=True, text=True)
    # "device" at end of a line (not "offline", not just the header)
    connected = any(
        line.endswith("device")
        for line in result.stdout.splitlines()[1:]  # skip "List of devices attached"
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
        exe_path = build_dir / "benchmarks" / exe
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
