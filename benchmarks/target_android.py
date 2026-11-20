"""The Android target: bench_fft built with the NDK, run on a device over adb.

    python3 benchmarks/run.py --target android[:SERIAL] -- --sizes 1024
    python3 benchmarks/compare.py master . --target android

Needs adb on the PATH and an Android NDK: $ANDROID_NDK or $ANDROID_NDK_HOME,
else the newest one in the ndk/ directory of $ANDROID_HOME or of the SDK's
default location. The ABI is the device's (ro.product.cpu.abi), the API level
24. A run copies its binaries into /data/local/tmp/pffft-bench/ and deletes
that directory when it ends, so one device serves one run at a time.
"""
from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path

import targets
from targets import Error

API = 24
DEVICE_DIR = "/data/local/tmp/pffft-bench"
TOOLCHAIN = Path("build", "cmake", "android.toolchain.cmake")


def find_ndk() -> Path:
    for var in ("ANDROID_NDK", "ANDROID_NDK_HOME"):
        if os.environ.get(var):
            ndk = Path(os.environ[var])
            if not (ndk / TOOLCHAIN).is_file():
                raise Error(f"${var}={ndk} is not an Android NDK: it has no {TOOLCHAIN}")
            return ndk
    sdks = [os.environ.get("ANDROID_HOME"), Path.home() / "Library/Android/sdk", Path.home() / "Android/Sdk"]
    ndks = [p for sdk in sdks if sdk for p in Path(sdk, "ndk").glob("*") if (p / TOOLCHAIN).is_file()]
    if not ndks:
        raise Error("no Android NDK found: set ANDROID_NDK, or install one with the SDK manager")
    return max(ndks, key=lambda p: [int(x) for x in p.name.split(".") if x.isdigit()])


class AndroidTarget(targets.Target):
    """An Android device, the only one attached or the one with this serial."""

    def __init__(self, serial: str = None):
        self.adb = ["adb"] + (["-s", serial] if serial else [])
        self.abi = self.call("shell", "getprop ro.product.cpu.abi")
        self.name = f"android-{self.abi}"
        self.ndk = find_ndk()
        self.cmake_args = [f"-DCMAKE_TOOLCHAIN_FILE={self.ndk / TOOLCHAIN}",
                           f"-DANDROID_ABI={self.abi}", f"-DANDROID_PLATFORM={API}"]

    def call(self, *args: str, check: bool = True) -> str:
        return targets.call(*self.adb, *args, check=check)

    def fftw_build_args(self):
        if self.abi != "arm64-v8a":
            raise Error(f"--fftw on Android is set up for arm64-v8a only, not {self.abi}")
        import build_fftw
        host = "aarch64-linux-android"
        cc = next(self.ndk.glob(f"toolchains/llvm/prebuilt/*/bin/{host}{API}-clang"), None)
        if cc is None:
            raise Error(f"{self.ndk} has no {host}{API}-clang")
        return {"cc": str(cc), "host": host, "configure_args": build_fftw.ARM64_ARGS}

    def run(self, harness, bench_args, libs=None):
        files = {f"{DEVICE_DIR}/bench_fft": harness}
        device_libs = {}
        for label, path in (libs or {}).items():
            device_libs[label] = f"{DEVICE_DIR}/libpffft-{label}.so"
            files[device_libs[label]] = path
        cmd = [f"{DEVICE_DIR}/bench_fft", *bench_args, *targets.lib_options(device_libs), "--csv", "-"]
        try:
            self.call("shell", f"mkdir -p {DEVICE_DIR}")
            for remote, local in files.items():
                self.call("push", str(local), remote)
            # bench_fft's stderr (progress, warnings) reaches the user directly
            p = subprocess.run(self.adb + ["shell", shlex.join(cmd)], stdin=subprocess.DEVNULL,
                               stdout=subprocess.PIPE, text=True)
        finally:
            self.call("shell", f"rm -rf {DEVICE_DIR}", check=False)
        return targets.console_csv(p.stdout, p.returncode, cmd)
