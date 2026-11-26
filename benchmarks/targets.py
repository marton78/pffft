"""Build bench_fft and pffft shared libraries for a machine, and run them there.

A target is a machine that runs bench_fft: this computer (LocalTarget), this
Mac running x86-64 code under Rosetta (RosettaTarget, for results only), an
Android phone (target_android.py) or an iPhone (target_ios.py). Everything is
built on this computer, in <repo>/build_bench_cache/ (ignored by git):

    src/<commit>/                    a revision, exported with git archive
    <target>/harness-worktree/       bench_fft from the working tree
    <target>/harness-fftw-worktree/  the same, with FFTW (fftw=True)
    <target>/lib-<commit>/           pffft of a revision as a shared library
    <target>/lib-worktree/           the same, from the working tree
    <target>/lib-split-<commit>/     the same, splitting every array of more
                                     than 1 KiB (accuracy.py)
    fftw-<target>/                   FFTW built by build_fftw.py for the target

A revision's export and build are kept by commit hash; the working tree ('.')
is built in place and incrementally, never cached by hash.

A remote target is a Target subclass in its own module. It sets cmake_args
(its toolchain) and fftw_build_args(), implements run(), and overrides more
only where it must, reusing the helpers here: call(), export_rev(), build_dir(),
cmake_build(), find_file(), lib_options(), console_csv() (which finds the CSV
in console output that mixes in other lines). make_target() imports it only
when asked for it:

    target_android.AndroidTarget(serial)   for "android" / "android:<serial>"
    target_ios.IosTarget(udid)             for "ios" / "ios:<udid>"

where serial and udid are None when the spec names none.
"""
from __future__ import annotations

import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import results

REPO = Path(__file__).resolve().parent.parent
CACHE = REPO / "build_bench_cache"

# Both revisions of a comparison, and the harness, get the same flags. The
# architecture is pinned to pffft's default (the toolchain's own target) so
# that a revision with another default cannot differ from the others.
COMMON_ARGS = ["-DCMAKE_BUILD_TYPE=Release", "-DTARGET_C_ARCH=none"]
HARNESS_ARGS = COMMON_ARGS + ["-DPFFFT_BUILD_BENCHMARKS=ON", "-DPFFFT_BUILD_TESTS=OFF",
                              "-DPFFFT_BUILD_EXAMPLES=OFF", "-DPFFFT_USE_BENCH_FFTW=OFF"]
LIB_ARGS = COMMON_ARGS + ["-DBUILD_SHARED_LIBS=ON", "-DPFFFT_BUILD_BENCHMARKS=OFF",
                          "-DPFFFT_BUILD_TESTS=OFF", "-DPFFFT_BUILD_EXAMPLES=OFF"]
LIB_NAMES = ("libpffft.dylib", "libpffft.so", "pffft.dll", "libpffft.dll")


class Error(Exception):
    """A failure the tools report in one message, without a traceback."""


def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


class Target:
    """A machine that runs bench_fft. build_harness() and build_lib() build on
    this computer with CMake, adding cmake_args; a target whose binaries need
    more (packaging, signing, copying to the device) overrides them."""

    name: str  # names its directory in the build cache, e.g. "local"
    cmake_args: list = []  # for every build: toolchain file, ABI, ...

    def fftw_build_args(self) -> dict:
        """Keyword arguments of build_fftw.build() for this target: cc, cflags,
        host, configure_args (SIMD flags)."""
        raise NotImplementedError

    def build_harness(self, src: Path, cmake_args=(), fftw: bool = False) -> object:
        """Build bench_fft from source tree src for this target, incrementally.
        With fftw, build FFTW for the target first and link it. Returns what
        run() needs to start it: here the file on this computer."""
        args = HARNESS_ARGS + self.cmake_args + list(cmake_args)
        if fftw:
            import build_fftw  # only --fftw needs it
            prefix = build_fftw.build(CACHE / f"fftw-{self.name}", **self.fftw_build_args())
            args += [f"-DFFTW3_ROOT={prefix}", "-DPFFFT_USE_BENCH_FFTW=ON"]
        build = build_dir(self.name, "harness-fftw" if fftw else "harness", src)
        cmake_build(src, build, args, "bench_fft")
        return find_file(build, ("bench_fft", "bench_fft.exe"))

    def build_lib(self, src: Path, label: str, cmake_args=(), variant: str = "") -> Path:
        """Build pffft from source tree src as a shared library for this target,
        with that tree's own CMake; return the library file on this computer.
        label is a short name, unique within one run and safe in file names: a
        remote target can use it to name the copy on the device. A variant
        (e.g. "split") is built with cmake_args added, in a directory of its
        own."""
        build = build_dir(self.name, f"lib{'-' + variant if variant else ''}", src)
        try:
            cmake_build(src, build, LIB_ARGS + self.cmake_args + list(cmake_args), "PFFFT")
            return find_file(build, LIB_NAMES)
        except Error as e:
            if src == REPO:
                raise
            raise Error(f"{e}\nif {src.name} is older than commit fde5014 (February 2026): pffft builds "
                        "as a shared library with its own CMake only from that commit on") from None

    def run(self, harness: object, bench_args: list, libs: dict = None) -> str:
        """Run bench_fft with bench_args and --csv -, plus --lib label=<path
        there> for each label -> library file in libs. Its stderr (progress,
        warnings) goes to the user as it comes. Returns the CSV text; raises
        Error if bench_fft produced none."""
        raise NotImplementedError


def make_target(spec: str) -> Target:
    """'local', 'rosetta[:arch]', 'android[:serial]' or 'ios[:udid]'."""
    kind, _, device = spec.partition(":")
    if kind == "local" and not device:
        return LocalTarget()
    if kind == "rosetta":
        return RosettaTarget(device or "haswell")
    if kind == "android":
        import target_android
        return target_android.AndroidTarget(device or None)
    if kind == "ios":
        import target_ios
        return target_ios.IosTarget(device or None)
    raise Error(f"unknown target '{spec}': expected local, rosetta[:arch], android[:serial] or ios[:udid]")


def call(*cmd, check: bool = True) -> str:
    """Run cmd without input; its stdout, stripped, or Error with its message."""
    try:
        p = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True)
    except FileNotFoundError:
        raise Error(f"{cmd[0]} not found: is it on the PATH?") from None
    if p.returncode and check:
        raise Error(f"{' '.join(map(str, cmd))}: {(p.stderr or p.stdout).strip()}")
    return p.stdout.strip()


def export_rev(rev: str):
    """Return (commit hash, source tree) of a revision; ('.', REPO) for '.'.
    A revision is exported once into CACHE/src/<commit>. Submodules are not
    exported: building the library does not need them."""
    if rev == ".":
        return ".", REPO
    commit = call("git", "-C", str(REPO), "rev-parse", "--verify", "--end-of-options", rev + "^{commit}")
    dst = CACHE / "src" / commit
    if not dst.is_dir():
        dst.parent.mkdir(parents=True, exist_ok=True)
        tmp = Path(tempfile.mkdtemp(prefix=commit + ".", dir=dst.parent))
        try:
            tar = subprocess.run(["git", "-C", str(REPO), "archive", "--format=tar", commit],
                                 stdout=subprocess.PIPE, check=True).stdout
            with tarfile.open(fileobj=io.BytesIO(tar)) as t:
                # the "data" filter exists from Python 3.12 and in later 3.9 to 3.11 releases
                t.extractall(tmp, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))
            os.rename(tmp, dst)
        except OSError:
            if not dst.is_dir():  # else another run exported it meanwhile
                raise
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    return commit, dst


def build_dir(target: str, kind: str, src: Path) -> Path:
    """CACHE/<target>/<kind>-<commit>, or <kind>-worktree for the working tree."""
    return CACHE / target / f"{kind}-{'worktree' if src == REPO else src.name}"


def cmake_build(src: Path, build: Path, args, cmake_target: str) -> None:
    """Configure and build one CMake target, incrementally. Configuring again
    each time is cheap and picks up changed arguments. The tools' output is
    shown only if they fail."""
    first = not (build / "CMakeCache.txt").exists()
    # Ninja where there is one, but not on Windows: there it would take whatever
    # compiler the PATH has; the default generator is Visual Studio, with MSVC
    gen = ["-G", "Ninja"] if first and shutil.which("ninja") and os.name != "nt" else []
    log(f"building {cmake_target} in {build}")
    for cmd in (["cmake", "-S", str(src), "-B", str(build), *gen, *args],
                ["cmake", "--build", str(build), "--target", cmake_target, "--config", "Release",
                 "--parallel", str(os.cpu_count() or 4)]):
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if p.returncode != 0:
            raise Error(f"{p.stdout}\nfailed ({p.returncode}): {' '.join(cmd)}")


def find_file(build: Path, names) -> Path:
    """The first of names in build or the subdirectories where generators put
    their output."""
    for sub in ("", "benchmarks", "Release", "benchmarks/Release"):
        for name in names:
            if (build / sub / name).is_file():
                return build / sub / name
    raise Error(f"none of {', '.join(names)} found in {build}")


def lib_options(libs: dict) -> list:
    """--lib options for bench_fft; libs maps label -> path as bench_fft sees it."""
    return [a for label, path in (libs or {}).items() for a in ("--lib", f"{label}={path}")]


def local_arch() -> str:
    """TARGET_C_ARCH for the best code for this computer's CPU: native for
    GCC and Clang. MSVC has no such flag: AVX512 or AVX2 where the CPU has
    them (x64), else its default."""
    if os.name != "nt":
        return "native"
    import ctypes
    import platform
    if platform.machine().lower() in ("amd64", "x86_64"):
        present = ctypes.windll.kernel32.IsProcessorFeaturePresent
        if present(41):  # PF_AVX512F_INSTRUCTIONS_AVAILABLE
            return "AVX512"
        if present(40):  # PF_AVX2_INSTRUCTIONS_AVAILABLE
            return "AVX2"
    return "none"


class LocalTarget(Target):
    """This computer. Everything is built for its CPU (-march=native or the
    best /arch), pffft, the other engines and FFTW alike: a benchmark that
    only runs here may use all of it, and a comparison should give every
    library the same chance."""

    name = "local"

    def __init__(self):
        self.arch = local_arch()
        self.cmake_args = [f"-DTARGET_C_ARCH={self.arch}"]  # after COMMON_ARGS: wins

    def fftw_build_args(self):
        import build_fftw
        cflags = "-O3 -march=native" if self.arch == "native" else "-O3"
        return {"cflags": cflags, "configure_args": build_fftw.native_simd_args()}

    def run(self, harness, bench_args, libs=None):
        cmd = [str(harness), *bench_args, *lib_options(libs), "--csv", "-"]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
        return console_csv(p.stdout, p.returncode, cmd)


class RosettaTarget(Target):
    """x86-64 code on an Apple arm64 Mac, run by Rosetta 2 with AVX, AVX2 and
    FMA advertised (ROSETTA_ADVERTISE_AVX=1, macOS 15 on), built for -march=arch
    (default haswell: AVX2 and FMA). For checking results, as accuracy.py does
    for the x86 SIMD code: timings under translation say nothing about an x86
    CPU."""

    def __init__(self, arch: str):
        import platform
        if sys.platform != "darwin" or platform.machine() != "arm64":
            raise Error("the rosetta target needs an Apple arm64 Mac")
        self.name = f"rosetta-{arch}"
        self.cmake_args = ["-DCMAKE_OSX_ARCHITECTURES=x86_64", "-DCMAKE_SYSTEM_NAME=Darwin",
                           "-DCMAKE_SYSTEM_PROCESSOR=x86_64",  # older revisions enable SSE only if it is set
                           f"-DTARGET_C_ARCH={arch}"]

    def fftw_build_args(self):
        raise Error("no FFTW for the rosetta target")

    def run(self, harness, bench_args, libs=None):
        cmd = [str(harness), *bench_args, *lib_options(libs), "--csv", "-"]
        env = dict(os.environ, ROSETTA_ADVERTISE_AVX="1")
        p = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, env=env)
        return console_csv(p.stdout, p.returncode, cmd)


def console_csv(output: str, status: int, cmd) -> str:
    """The CSV in bench_fft's console output. bench_fft exits nonzero also
    after measuring, when an engine failed its self-check: then the CSV is
    still returned, as its stderr has told the user."""
    try:
        csv = results.extract_csv(output)
    except ValueError:
        sys.stderr.write(output)  # --help, --list, or why it failed
        raise Error(f"bench_fft {f'failed ({status}) and ' if status else ''}wrote no CSV: "
                    f"{' '.join(map(str, cmd))}") from None
    if status != 0:
        log(f"warning: bench_fft exited with status {status}; see its messages above")
    return csv
