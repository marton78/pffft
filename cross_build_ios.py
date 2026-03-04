#!/usr/bin/env python3
"""Build pffft benchmarks for iOS and optionally run them on a connected device.

Usage:
    python cross_build_ios.py <id> [options] [extra cmake options]

Arguments:
    id              Short identifier used in directory names (e.g. "arm64")

Options:
    --fftw          Cross-compile FFTW 3.3.10 and include it in the benchmark
                    (requires autoconf/make)
    --no-run        Build only; do not deploy or run on device
    --output-dir    Local directory to store pulled benchmark CSVs
                    (default: bench_results_ios_<id>)
    --sdk <path>    Override iOS SDK root directory
    --deployment-target <ver>
                    iOS deployment target (default: 15)
    --sign <identity>
                    Code signing identity for device deployment
                    (auto-detected from keychain if omitted)
    --team-id <id>  Apple development team ID (auto-detected from cert)
    --bundle-id <id>
                    App bundle ID prefix (default: com.example.pffft)

Environment variables:
    IOS_SDK         iOS SDK root (overrides auto-detection)

Device deployment prerequisites:
    1. An Apple Developer account (free or paid) signed into Xcode
    2. A code signing identity (auto-detected, or pass --sign)
    3. The target device MUST be registered in the Apple Developer portal
       under the same team as the signing identity. If ios-deploy fails
       with 0xe8000067, the device UDID is not in the provisioning profile.
       Register it at https://developer.apple.com/account/resources/devices
    4. ios-deploy installed: brew install ios-deploy

Examples:
    python cross_build_ios.py arm64
    python cross_build_ios.py arm64 --fftw
    python cross_build_ios.py arm64 --no-run
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


# ── iOS SDK detection ─────────────────────────────────────────────────────────

def find_ios_sdk():
    """Return path to the iOS SDK (arm64 variant), or None if not found."""
    if os.environ.get("IOS_SDK"):
        return Path(os.environ["IOS_SDK"])

    try:
        result = subprocess.run(
            ["xcrun", "--sdk", "iphoneos", "--show-sdk-path"],
            capture_output=True, text=True, check=True,
        )
        sdk_path = result.stdout.strip()
        if sdk_path and Path(sdk_path).is_dir():
            return Path(sdk_path)
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

    # Try common Xcode locations
    xcode_sdks = Path("/Applications/Xcode.app/Contents/Developer"
                      "/Platforms/iPhoneOS.platform/Developer/SDKs")
    if xcode_sdks.is_dir():
        sdk_dirs = sorted(
            [d for d in xcode_sdks.iterdir() if d.is_dir()],
            key=lambda p: p.name, reverse=True,
        )
        if sdk_dirs:
            return sdk_dirs[0]

    return None


# ── device detection ──────────────────────────────────────────────────────────

def list_connected_devices():
    """List connected iOS devices using ios-deploy.

    Returns list of dicts with 'id' and 'name' keys, or empty list.
    """
    try:
        result = subprocess.run(
            ["ios-deploy", "-c", "--timeout", "3"],
            capture_output=True, text=True, check=False,
        )
        if result.returncode != 0:
            return []

        devices = []
        for line in result.stdout.splitlines():
            # Format varies; look for lines with device identifiers
            if "Found" in line and "connected" in line:
                # e.g. "Found 00008110-... (N104AP, iPhone 15 Pro, ...)"
                parts = line.split("Found ", 1)
                if len(parts) == 2:
                    rest = parts[1]
                    dev_id = rest.split()[0] if rest.split() else "unknown"
                    # Extract name from parenthetical
                    name = "iOS device"
                    if "(" in rest and ")" in rest:
                        paren = rest[rest.index("(") + 1:rest.rindex(")")]
                        name = paren
                    devices.append({"id": dev_id, "name": name})
        return devices
    except FileNotFoundError:
        return []


def collect_device_info(device_id):
    """Collect device info via ideviceinfo if available.

    Returns dict with device properties.
    """
    info = {
        "device_id": device_id,
        "timestamp": date.today().isoformat(),
    }

    try:
        result = subprocess.run(
            ["ideviceinfo", "-u", device_id],
            capture_output=True, text=True, check=False,
        )
        if result.returncode == 0:
            for line in result.stdout.splitlines():
                if ": " in line:
                    key, val = line.split(": ", 1)
                    if key == "ProductName":
                        info["product"] = val
                    elif key == "ProductVersion":
                        info["ios_version"] = val
                    elif key == "DeviceName":
                        info["device_name"] = val
                    elif key == "HardwareModel":
                        info["hardware"] = val
                    elif key == "ProductType":
                        info["model"] = val
    except FileNotFoundError:
        print("  Note: Install libimobiledevice for device info: brew install libimobiledevice",
              file=sys.stderr)

    return info


# ── FFTW cross-compilation ────────────────────────────────────────────────────

# iOS arm64 FFTW configure flags, mirroring Android's ABI_CONFIGURE.
# --enable-neon: NEON SIMD codelets (float only — do NOT pass for double).
# iOS uses mach_absolute_time for timing, which FFTW auto-detects on Darwin.
# Unlike Android, we do NOT need --enable-armv8-cntvct-el0 (Linux-only).
IOS_FFTW_FLOAT_FLAGS  = ["--enable-neon"]
IOS_FFTW_DOUBLE_FLAGS = []


def build_fftw_ios(sdk_root, deployment_target, build_dir, cpu_count):
    """Download and cross-compile FFTW (float + double) for iOS arm64.

    Returns the path to the install prefix containing include/ and lib/.
    """
    try:
        result = subprocess.run(
            ["xcrun", "-f", "clang"],
            capture_output=True, text=True, check=True,
        )
        clang = result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: xcrun clang not found. Ensure Xcode is installed.", file=sys.stderr)
        sys.exit(1)

    fftw_src_dir = build_dir / "fftw_src"
    fftw_prefix  = build_dir / "fftw_prefix"
    tarball      = fftw_src_dir / f"fftw-{FFTW_VERSION}.tar.gz"
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

    sysroot = str(sdk_root)
    deployment_flag = f"-miphoneos-version-min={deployment_target}"

    base_configure_args = [
        str(fftw_tree / "configure"),
        "--host=aarch64-apple-darwin",
        "--disable-shared",
        "--enable-static",
        f"--prefix={fftw_prefix}",
    ]

    env = os.environ.copy()
    cflags = f"-O3 -fomit-frame-pointer -fstrict-aliasing {deployment_flag} -isysroot {sysroot} -arch arm64"
    env.update({
        "CC":      clang,
        "AR":      "libtool",
        "ARFLAGS": "-static -o",
        "RANLIB":  "ranlib",
        "CFLAGS":  cflags,
        "LDFLAGS": f"{deployment_flag} -isysroot {sysroot} -arch arm64",
    })

    for label, prec_flags, extra in [
        ("float",  IOS_FFTW_FLOAT_FLAGS,  ["--enable-float"]),
        ("double", IOS_FFTW_DOUBLE_FLAGS, []),
    ]:
        banner([f"Building FFTW {FFTW_VERSION} — {label} precision (iOS arm64)"])
        bdir = build_dir / f"fftw_build_{label}"
        bdir.mkdir(parents=True, exist_ok=True)
        run(base_configure_args + prec_flags + extra, cwd=bdir, env=env)
        run(["make", "-j", str(cpu_count)], cwd=bdir)
        run(["make", "install"], cwd=bdir)

    return fftw_prefix


# ── pffft CMake build ─────────────────────────────────────────────────────────

def build_pffft_ios(script_dir, build_dir, deployment_target,
                    fftw_cmake_args, extra_cmake, cpu_count,
                    signing_identity=None, team_id=None, bundle_id=None):
    """Build pffft benchmarks for iOS using CMake + Xcode.

    When signing_identity and team_id are provided, builds signed .app bundles
    via the Xcode generator and xcodebuild (which manages provisioning profiles
    automatically).  Otherwise falls back to Unix Makefiles (unsigned, --no-run).

    Returns path to the pffft build directory.
    """
    pffft_build_dir = build_dir / "pffft"
    if pffft_build_dir.exists():
        shutil.rmtree(pffft_build_dir)
    pffft_build_dir.mkdir(parents=True)

    use_xcode = signing_identity is not None and team_id is not None

    banner([
        "Building pffft benchmarks for iOS arm64",
        f"Deployment target: {deployment_target}",
        f"Generator: {'Xcode (signed)' if use_xcode else 'Unix Makefiles (unsigned)'}",
    ])

    cmake_args = [
        "cmake", "-S", str(script_dir), "-B", str(pffft_build_dir),
        "-DCMAKE_SYSTEM_NAME=iOS",
        "-DCMAKE_SYSTEM_PROCESSOR=arm64",
        f"-DCMAKE_OSX_DEPLOYMENT_TARGET={deployment_target}",
        "-DCMAKE_OSX_ARCHITECTURES=arm64",
        # Enable NEON — all iOS arm64 devices support it
        "-DTARGET_C_ARCH=armv8-a",
        "-DTARGET_CXX_ARCH=armv8-a",
        "-DPFFFT_BUILD_TESTS=OFF",
        "-DPFFFT_BUILD_BENCHMARKS=ON",
        "-DPFFFT_BUILD_EXAMPLES=OFF",
        "-DPFFFT_USE_BENCH_MKL=OFF",
        "-DPFFFT_USE_BENCH_FFTS=OFF",
    ] + fftw_cmake_args + extra_cmake

    if use_xcode:
        cmake_args += [
            "-G", "Xcode",
            f"-DCMAKE_XCODE_ATTRIBUTE_DEVELOPMENT_TEAM={team_id}",
            "-DCMAKE_XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY=Apple Development",
            "-DCMAKE_XCODE_ATTRIBUTE_CODE_SIGN_STYLE=Automatic",
        ]
        if bundle_id:
            cmake_args += [
                f"-DCMAKE_XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER={bundle_id}",
            ]
    else:
        cmake_args += ["-DCMAKE_BUILD_TYPE=Release"]

    run(cmake_args)

    build_cmd = [
        "cmake", "--build", str(pffft_build_dir),
        "--config", "Release",
    ]
    if use_xcode:
        # xcodebuild handles parallelism and signing
        build_cmd += [
            "--",
            "-jobs", str(cpu_count),
            "-allowProvisioningUpdates",
            f"DEVELOPMENT_TEAM={team_id}",
            "CODE_SIGN_IDENTITY=Apple Development",
            "CODE_SIGN_STYLE=Automatic",
        ]
    else:
        build_cmd += ["--", f"-j{cpu_count}"]

    run(build_cmd)

    print(f"\nBuild complete: {pffft_build_dir}")
    return pffft_build_dir


# ── code signing ──────────────────────────────────────────────────────────────

def find_signing_identity():
    """Find a valid iOS development signing identity.

    Returns the identity string (e.g. "Apple Development: ...") or None.
    """
    try:
        result = subprocess.run(
            ["security", "find-identity", "-v", "-p", "codesigning"],
            capture_output=True, text=True, check=True,
        )
        for line in result.stdout.splitlines():
            for prefix in ("Apple Development", "iPhone Developer"):
                if prefix in line:
                    start = line.index('"') + 1
                    end = line.index('"', start)
                    return line[start:end]
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError):
        pass
    return None


def find_team_id(signing_identity):
    """Extract the development team ID (OU field) from a signing certificate.

    Returns the team ID string or None.
    """
    try:
        # Find the certificate matching the identity and extract OU
        result = subprocess.run(
            ["security", "find-certificate", "-c", signing_identity, "-p"],
            capture_output=True, text=True, check=True,
        )
        result2 = subprocess.run(
            ["openssl", "x509", "-noout", "-subject"],
            input=result.stdout, capture_output=True, text=True, check=True,
        )
        # Parse "subject=...OU=XXXXXXXXXX,..."
        for part in result2.stdout.replace(" ", "").split(","):
            if part.startswith("OU="):
                return part[3:]
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass
    return None


# ── deployment via ios-deploy ─────────────────────────────────────────────────

def deploy_and_run(app_bundle, device_id):
    """Deploy app bundle to device via ios-deploy and run it.

    Returns True on success.
    """
    banner([f"Deploying {app_bundle.name} to device"])

    result = run(
        ["ios-deploy", "--bundle", str(app_bundle), "--id", device_id,
         "--noninteractive", "--justlaunch"],
        check=False,
    )
    if result.returncode != 0:
        print(f"WARNING: ios-deploy failed for {app_bundle.name}", file=sys.stderr)
        print("  If error 0xe8000067: device UDID not in provisioning profile.", file=sys.stderr)
        print("  Register at https://developer.apple.com/account/resources/devices", file=sys.stderr)
        return False

    return True


# ── argument parsing ──────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("id", help="Short build identifier (used in directory names)")
    parser.add_argument("--fftw", action="store_true",
                        help="Cross-compile FFTW and include it in the benchmark")
    parser.add_argument("--no-run", action="store_true",
                        help="Build only; skip deployment and run on device")
    parser.add_argument("--output-dir", default=None,
                        help="Local directory for pulled CSV results")
    parser.add_argument("--sdk", default=None,
                        help="Override iOS SDK root directory")
    parser.add_argument("--deployment-target", default="15",
                        help="iOS deployment target (default: 15)")
    parser.add_argument("--sign", default=None,
                        help="Code signing identity (auto-detected if omitted)")
    parser.add_argument("--team-id", default=None,
                        help="Apple development team ID (auto-detected from cert if omitted)")
    parser.add_argument("--bundle-id", default="com.example.pffft",
                        help="App bundle ID prefix (default: com.example.pffft)")

    args, extra_cmake = parser.parse_known_args()
    return args, extra_cmake


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args, extra_cmake = parse_args()
    script_dir = Path(__file__).resolve().parent

    if platform.system() != "Darwin":
        print("ERROR: iOS cross-compilation is only supported on macOS.", file=sys.stderr)
        sys.exit(1)

    # ── resolve iOS SDK ───────────────────────────────────────────────────────

    sdk_root = Path(args.sdk) if args.sdk else find_ios_sdk()

    if sdk_root is None or not sdk_root.is_dir():
        print("ERROR: iOS SDK not found. Install Xcode or set IOS_SDK.", file=sys.stderr)
        sys.exit(1)

    print(f"Using iOS SDK: {sdk_root}")

    # ── resolve device ────────────────────────────────────────────────────────

    device = None
    if not args.no_run:
        devices = list_connected_devices()

        if not devices:
            print("ERROR: No iOS device connected.", file=sys.stderr)
            print("  Connect an iOS 15+ device and ensure ios-deploy can access it.", file=sys.stderr)
            print("  Install ios-deploy: brew install ios-deploy", file=sys.stderr)
            sys.exit(1)

        if len(devices) > 1:
            print(f"Multiple devices found: {[d['name'] for d in devices]}", file=sys.stderr)
            print("  Using first device.", file=sys.stderr)

        device = devices[0]
        print(f"Target device: {device['name']} ({device['id']})")

    # ── resolve signing identity and team ────────────────────────────────────

    signing_identity = None
    team_id = None
    if not args.no_run:
        signing_identity = args.sign or find_signing_identity()
        if signing_identity is None:
            print("ERROR: No code signing identity found.", file=sys.stderr)
            print("  iOS requires signed binaries for device deployment.", file=sys.stderr)
            print("  Options:", file=sys.stderr)
            print("    1. Install Xcode and sign in with an Apple ID", file=sys.stderr)
            print("    2. Pass --sign <identity> explicitly", file=sys.stderr)
            print("  List identities: security find-identity -v -p codesigning", file=sys.stderr)
            sys.exit(1)
        print(f"Signing identity: {signing_identity}")

        team_id = args.team_id or find_team_id(signing_identity)
        if team_id is None:
            print("ERROR: Could not determine development team ID.", file=sys.stderr)
            print("  Pass --team-id <TEAM_ID> explicitly.", file=sys.stderr)
            sys.exit(1)
        print(f"Development team: {team_id}")

    banner([
        "iOS pffft benchmark cross-compiler",
        f"ID={args.id}  Target={args.deployment_target}  FFTW={'yes' if args.fftw else 'no'}",
    ])

    # ── directories ───────────────────────────────────────────────────────────

    build_dir  = script_dir / f"build_ios_{args.id}"
    output_dir = Path(args.output_dir) if args.output_dir else script_dir / f"bench_results_ios_{args.id}"

    cpu_count = os.cpu_count() or 4

    # ── optionally cross-compile FFTW ─────────────────────────────────────────

    fftw_cmake_args = []
    if args.fftw:
        fftw_prefix = build_fftw_ios(sdk_root, args.deployment_target, build_dir, cpu_count)
        fftw_cmake_args = [
            "-DPFFFT_USE_BENCH_FFTW=ON",
            f"-DFFTW3_ROOT={fftw_prefix}",
        ]
    else:
        fftw_cmake_args = ["-DPFFFT_USE_BENCH_FFTW=OFF"]

    # ── build pffft + benchmarks ──────────────────────────────────────────────

    pffft_build_dir = build_pffft_ios(
        script_dir, build_dir, args.deployment_target,
        fftw_cmake_args, extra_cmake, cpu_count,
        signing_identity=signing_identity, team_id=team_id,
        bundle_id=args.bundle_id,
    )

    if args.no_run:
        print("Skipping device deployment (--no-run).")
        return

    # ── collect device info ───────────────────────────────────────────────────

    output_dir.mkdir(parents=True, exist_ok=True)

    print("\nCollecting device info...")
    device_info = collect_device_info(device["id"])

    device_info_lines = [
        f"Device: {device_info.get('device_name', device['name'])}",
        f"Model: {device_info.get('model', 'unknown')}",
        f"iOS: {device_info.get('ios_version', 'unknown')}",
        f"Hardware: {device_info.get('hardware', 'unknown')}",
        f"Arch: arm64",
        f"FFTW: {FFTW_VERSION if args.fftw else 'no'}",
        f"Date: {device_info.get('timestamp')}",
    ]
    for line in device_info_lines:
        print(line)
    (output_dir / "device_info.txt").write_text("\n".join(device_info_lines) + "\n")

    info_json = {
        "cpu":      f"{device_info.get('device_name', device['name'])} ({device_info.get('hardware', 'unknown')})",
        "arch":     "arm64",
        "os":       f"iOS {device_info.get('ios_version', 'unknown')}",
        "compiler": "Apple Clang (Xcode)",
        "date":     device_info.get("timestamp", date.today().isoformat()),
    }
    info_path = output_dir / "info.json"
    info_path.write_text(json.dumps(info_json, indent=2) + "\n")
    print(f"Wrote {info_path}")

    # ── deploy and run benchmarks ─────────────────────────────────────────────

    banner(["Running benchmarks on device"])

    for exe_name in ("bench_pffft_float", "bench_pffft_double"):
        # Xcode puts signed .app bundles under Release-iphoneos/
        app_bundle = pffft_build_dir / "benchmarks" / "Release-iphoneos" / f"{exe_name}.app"
        if not app_bundle.is_dir():
            # Fallback: CMake default (Unix Makefiles) path
            app_bundle = pffft_build_dir / "benchmarks" / f"{exe_name}.app"
        if not app_bundle.is_dir():
            print(f"WARNING: {exe_name}.app not found (was it built?)", file=sys.stderr)
            continue

        deploy_and_run(app_bundle, device["id"])

    # ── result collection notes ───────────────────────────────────────────────

    banner(["Result Collection"])

    print("iOS app sandbox prevents direct file pull (unlike Android's adb pull).")
    print("\nOptions to retrieve benchmark CSV results:")
    print("  1. Xcode > Window > Devices and Simulators > download app container")
    print("  2. ideviceinstaller: ideviceinstaller -o copy -u <device-id> ...")
    print("  3. SSH on jailbroken device: scp root@<ip>:/.../Documents/*.csv .")
    print(f"\nOnce CSVs are in {output_dir}/, generate charts:")
    print(f'  python3 bench/make_charts.py "{output_dir}"')


if __name__ == "__main__":
    main()
