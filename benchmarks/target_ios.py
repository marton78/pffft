"""The iOS target: an iPhone attached to this Mac, which needs Xcode.

    python3 benchmarks/run.py --target ios[:DEVICE] -- --sizes 1024
    python3 benchmarks/compare.py master . --target ios

DEVICE is the device's UDID, CoreDevice identifier or name (xcrun devicectl
list devices); without it the only iOS device attached is used. bench_fft and
the pffft libraries are built here with the iOS SDK. Each run packs them into
one app, signs it, installs it with devicectl and launches it with --console,
which connects the app's standard streams to ours. The libraries go into the
app's Frameworks/ and bench_fft gets them as
@executable_path/Frameworks/<label>.dylib, a path dlopen resolves itself. Keep
the device unlocked while it runs.

Signing needs a development identity in the keychain and a wildcard
provisioning profile (application identifier TEAM.*) for its certificate that
lists the device: "iOS Team Provisioning Profile: *", which Xcode makes when it
signs an app automatically for a device registered with the team.
"""
from __future__ import annotations

import hashlib
import json
import plistlib
import shutil
import subprocess
import tempfile
from pathlib import Path

import targets
from targets import Error, call

BUNDLE_ID = "com.example.pffft.bench"
MIN_IOS = "15.0"
PROFILE_DIRS = (Path.home() / "Library/Developer/Xcode/UserData/Provisioning Profiles",  # Xcode 16 on
                Path.home() / "Library/MobileDevice/Provisioning Profiles")  # older Xcode


def devicectl(*args: str) -> dict:
    """Run xcrun devicectl; the result part of its JSON output (the only
    output of it meant for programs)."""
    with tempfile.TemporaryDirectory() as tmp:
        call("xcrun", "devicectl", *args, "--quiet", "--json-output", f"{tmp}/out.json")
        return json.loads(Path(tmp, "out.json").read_text())["result"]


def signing(udid: str):
    """(identity, profile file, application identifier) for the app on device
    udid, from the wildcard profile that lists the device and expires last;
    the identity is the SHA-1 of one of its certificates in the keychain."""
    found = call("security", "find-identity", "-v", "-p", "codesigning")
    identities = {line.split()[1] for line in found.splitlines() if line.strip()[:1].isdigit()}
    best = None
    for path in (f for d in PROFILE_DIRS if d.is_dir() for f in d.glob("*.mobileprovision")):
        p = subprocess.run(["security", "cms", "-D", "-i", str(path)], stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL)
        profile = plistlib.loads(p.stdout) if p.returncode == 0 else {}
        certs = {hashlib.sha1(c).hexdigest().upper() for c in profile.get("DeveloperCertificates", ())}
        if (str(profile.get("Entitlements", {}).get("application-identifier")).endswith(".*")
                and udid in profile.get("ProvisionedDevices", ()) and certs & identities
                and (best is None or profile["ExpirationDate"] > best[0]["ExpirationDate"])):
            best = profile, path, min(certs & identities)
    if best is None:
        raise Error(f"no wildcard provisioning profile that lists device {udid} and has a certificate "
                    f"of a code-signing identity in the keychain, in {' or '.join(map(str, PROFILE_DIRS))}: "
                    "let Xcode sign an app for the device automatically, which makes one")
    profile, path, identity = best
    return identity, path, profile["Entitlements"]["application-identifier"][:-1] + BUNDLE_ID


class IosTarget(targets.Target):
    """An iOS device attached to this Mac, the only one or the one named."""

    name = "ios"
    cmake_args = ["-DCMAKE_SYSTEM_NAME=iOS", "-DCMAKE_OSX_ARCHITECTURES=arm64",
                  "-DCMAKE_SYSTEM_PROCESSOR=arm64",  # older revisions enable NEON only if it is set
                  f"-DCMAKE_OSX_DEPLOYMENT_TARGET={MIN_IOS}",
                  "-DCMAKE_MACOSX_BUNDLE=OFF"]  # a plain executable: run() makes the app

    def __init__(self, device: str = None):
        devices = [d for d in devicectl("list", "devices")["devices"]
                   if d["hardwareProperties"].get("platform") == "iOS"
                   and d["connectionProperties"].get("transportType")  # attached by cable or network
                   and (device is None or device in (d["identifier"], d["hardwareProperties"].get("udid"),
                                                      d["deviceProperties"].get("name")))]
        if not devices:
            raise Error(f"no iOS device {repr(device) + ' ' if device else ''}attached (xcrun devicectl list devices)")
        if len(devices) > 1:
            raise Error("several iOS devices attached, choose one with --target ios:<UDID>: "
                        + ", ".join(d["hardwareProperties"]["udid"] for d in devices))
        self.udid = devices[0]["hardwareProperties"]["udid"]

    def fftw_build_args(self):
        import build_fftw
        sdk = call("xcrun", "--sdk", "iphoneos", "--show-sdk-path")
        return {"cc": call("xcrun", "--sdk", "iphoneos", "-f", "clang"), "host": "aarch64-apple-darwin",
                "cflags": f"-O3 -arch arm64 -isysroot {sdk} -miphoneos-version-min={MIN_IOS}",
                "configure_args": build_fftw.ARM64_ARGS}

    def make_app(self, harness: Path, libs: dict) -> Path:
        """bench_fft.app with harness and libs (label -> library), signed for
        the device."""
        identity, profile, app_id = signing(self.udid)
        app = targets.CACHE / self.name / "bench_fft.app"
        shutil.rmtree(app, ignore_errors=True)
        (app / "Frameworks").mkdir(parents=True)
        shutil.copy2(harness, app / "bench_fft")
        shutil.copy2(profile, app / "embedded.mobileprovision")
        info = {"CFBundleIdentifier": BUNDLE_ID, "CFBundleExecutable": "bench_fft", "CFBundleName": "bench_fft",
                "CFBundleVersion": "1"}
        (app / "Info.plist").write_bytes(plistlib.dumps(info))
        sign = ["codesign", "--force", "--timestamp=none", "--sign", identity]
        for label, lib in libs.items():
            shutil.copy2(lib, app / "Frameworks" / f"{label}.dylib")
            call(*sign, str(app / "Frameworks" / f"{label}.dylib"))
        ent = app.with_name("entitlements.plist")
        ent.write_bytes(plistlib.dumps({"application-identifier": app_id}))
        call(*sign, "--entitlements", str(ent), str(app))
        return app

    def bench_pids(self) -> set:
        found = devicectl("device", "info", "processes", "--device", self.udid)["runningProcesses"]
        return {p["processIdentifier"] for p in found
                if p.get("executable", "").endswith("/bench_fft.app/bench_fft")}

    def run(self, harness, bench_args, libs=None):
        libs = libs or {}
        app = self.make_app(harness, libs)
        targets.log(f"installing {app} on the device")
        devicectl("device", "install", "app", "--device", self.udid, str(app))
        paths = {label: f"@executable_path/Frameworks/{label}.dylib" for label in libs}
        # Launched normally (in the foreground or not), bench_fft is killed
        # after about 20 s, as an app without a user interface; started
        # stopped, as if for a debugger, and resumed here, it may run on.
        cmd = ["xcrun", "devicectl", "device", "process", "launch", "--device", self.udid, "--quiet",
               "--console", "--start-stopped", "--terminate-existing", BUNDLE_ID,
               *bench_args, *targets.lib_options(paths), "--csv", "-"]
        old = self.bench_pids()
        # bench_fft's stderr (progress, warnings) reaches the user directly
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.DEVNULL, text=True)
        try:
            for _ in range(30):  # a query takes about a second
                new = self.bench_pids() - old
                if new or p.poll() is not None:
                    break
            else:
                raise Error("bench_fft did not start on the device")
            if new:
                devicectl("device", "process", "resume", "--device", self.udid, "--pid", str(max(new)))
            out = p.communicate()[0]
        finally:
            if p.poll() is None:
                p.kill()
        return targets.console_csv(out, p.returncode, cmd)
