import argparse
import json
import re
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import targets
from perf import cmd_ab
from samples import read_samples
from targets import LocalTarget, SshTarget, make_target, sanitize_meta_value
from targets import AdbTarget, IosTarget


# ---- provenance value sanitizer ------------------------------------------------

def test_sanitize():
    assert sanitize_meta_value("iPhone 11 Pro,\nv13") == "iPhone_11_Pro__v13"


def test_sanitize_passes_plain_values_through():
    assert sanitize_meta_value("raspberrypi") == "raspberrypi"
    assert sanitize_meta_value("iOS 17.5") == "iOS_17.5"


# ---- sync_worktree: real git, no mocks (this is the class of bug that slips
#      through when every other test mocks sync_worktree away) -----------------

def _git(cwd, *args):
    subprocess.run(["git", "-C", str(cwd), *args], check=True,
                   capture_output=True, text=True)


def test_sync_worktree_tracks_repo_root_head_advancing(tmp_path, monkeypatch):
    """Regression: sync_worktree(wt, "HEAD") once did `git -C wt reset --hard
    HEAD`, which resolves HEAD inside the LINKED WORKTREE's own detached
    HEAD -- a no-op once the worktree exists. New commits on REPO_ROOT's
    branch after worktree creation were silently never picked up by the
    default `--variant HEAD`."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "config", "user.email", "t@example.com")
    _git(repo, "config", "user.name", "t")
    (repo / "a.txt").write_text("v1\n")
    _git(repo, "add", "a.txt")
    _git(repo, "commit", "-q", "-m", "v1")

    monkeypatch.setattr(targets, "REPO_ROOT", repo)
    wt = tmp_path / "wt"
    _git(repo, "worktree", "add", str(wt), "HEAD")

    # Advance REPO_ROOT's branch with a NEW commit, created AFTER the
    # worktree already exists (exactly the cherry-pick-then-rerun-ab flow).
    (repo / "a.txt").write_text("v2\n")
    _git(repo, "commit", "-aq", "-m", "v2")

    targets.sync_worktree(wt, "HEAD")
    assert (wt / "a.txt").read_text() == "v2\n", \
        "sync_worktree('HEAD') must follow REPO_ROOT's current branch tip"

# ---- spec parsing / dispatch -------------------------------------------------

def test_make_target_local():
    t = make_target("local")
    assert isinstance(t, LocalTarget)
    assert t.name == "local"
    assert t.meta()["target"] == "local"


def test_ssh_spec_parses_user_and_host():
    t = make_target("ssh://pi@raspberrypi")
    assert isinstance(t, SshTarget)
    assert t.host == "raspberrypi" and t.user == "pi"
    assert t.dest == "pi@raspberrypi"


def test_ssh_spec_without_user_and_trailing_slash():
    t = make_target("ssh://raspi/")
    assert t.host == "raspi" and t.user is None
    assert t.name == "raspi"


@pytest.mark.parametrize("spec", ["ftp://x", "ssh://", "ssh://user@"])
def test_bad_specs_rejected(spec):
    with pytest.raises(ValueError):
        make_target(spec)


# ---- provenance meta ----------------------------------------------------------

def test_remote_host_from_uname(monkeypatch):
    calls = []

    def fake_sh(*cmd):
        calls.append(cmd)
        return "octopi\n"

    monkeypatch.setattr(targets, "sh", fake_sh)
    t = SshTarget("raspi", "pi")
    assert t.meta() == {"target": "ssh://pi@raspi", "host": "octopi"}
    # cached: sh must not be consulted again
    monkeypatch.setattr(
        targets, "sh",
        lambda *cmd: (_ for _ in ()).throw(AssertionError("not cached")))
    assert t.remote_host() == "octopi"


def test_meta_list_sanitizes_spaces_newlines_and_commas(monkeypatch):
    monkeypatch.setattr(targets, "sh", lambda *cmd: "my,host\nv7\n")
    assert SshTarget("h").meta_list() == \
        ["target=ssh://h", "host=my_host_v7"]


def test_remote_host_falls_back_to_spec_on_failure(monkeypatch):
    def boom(*cmd):
        raise RuntimeError("ssh down")

    monkeypatch.setattr(targets, "sh", boom)
    assert SshTarget("raspi").remote_host() == "raspi"


# ---- command construction (no execution) --------------------------------------

def test_sync_argv_shape():
    argv = make_target("ssh://pi@raspi").sync_argv("wt-base")
    assert argv[0] == "rsync" and "-az" in argv and "--delete" in argv
    assert argv[-1] == "pi@raspi:fftbench/wt-base/"
    i = argv.index("-R")
    assert argv[i + 1:i + 6] == ["src", "include", "CMakeLists.txt", "cmake",
                                 "benchmarks/bench_pffft.c"]
    pairs = set(zip(argv, argv[1:]))
    for e in (".git", ".perf", "build"):
        assert ("--exclude", e) in pairs


def test_build_script_shape():
    script = make_target("ssh://raspi").build_script("wt-var")
    assert "$HOME/fftbench/wt-var" in script
    assert "cmake -S $HOME/fftbench/wt-var -B $HOME/fftbench/wt-var/build" \
        in script
    assert "--target bench_pffft_float bench_pffft_double" in script


def test_binary_path_is_remote_path_under_fftbench(monkeypatch):
    monkeypatch.setattr(subprocess, "run",
                        lambda cmd, **kw:
                        subprocess.CompletedProcess(cmd, 0))
    wt_build = Path("/repo/.perf/wt-base/build/benchmarks")
    p = SshTarget("raspi").binary_path(wt_build, "dbl")
    assert str(p) == "fftbench/wt-base/build/benchmarks/" \
                     "bench_pffft_double"


# ---- run(): stdout capture -> local CSV ----------------------------------------

class FakePopen:
    """Minimal Popen stand-in emitting canned stdout lines."""

    ret = 0
    last_argv = None

    def __init__(self, argv, **kw):
        self.argv = list(argv)
        FakePopen.last_argv = self.argv
        self.stdout = iter(["# pffft-bench-samples v2\n",
                            "# host=fake-pi\n",
                            "algo,prec,xform,size,sample_ms,n_iter\n",
                            "pffft,flt,real,64,1.5,100000\n"])

    def wait(self):
        return self.ret

    def poll(self):
        return self.ret

    def kill(self):
        pass


@pytest.fixture()
def fake_ssh(monkeypatch):
    monkeypatch.setattr(subprocess, "Popen", FakePopen)


@pytest.fixture()
def fake_popen(monkeypatch):
    monkeypatch.setattr(subprocess, "Popen", FakePopen)


def test_run_captures_samples_stdout_to_local_file(tmp_path, fake_ssh):
    out = tmp_path / "samples" / "x.csv"
    cmd = ["fftbench/wt-var/build/benchmarks/bench_pffft_float",
           "--size", "64", "--runs", "2", "--samples", str(out),
           "--meta", "host=octopi"]
    lines = list(make_target("ssh://pi@raspi").run(cmd))
    # ssh invocation received --samples - (pure CSV over stdout)
    a = FakePopen.last_argv
    assert a[:2] == ["ssh", "pi@raspi"]
    assert a[a.index("--samples") + 1] == "-"
    # all lines yielded AND written verbatim to the local file
    assert lines[0] == "# pffft-bench-samples v2"
    prov, rows = read_samples(out)
    assert prov["host"] == "fake-pi"
    assert rows[0].sample_ms == 1.5


def test_run_without_samples_arg_yields_only(fake_ssh):
    lines = list(SshTarget("raspi").run(["uname", "-a"]))
    assert lines[0] == "# pffft-bench-samples v2"
    assert lines[-1] == "pffft,flt,real,64,1.5,100000"


def test_run_raises_on_failure(tmp_path, fake_ssh, monkeypatch):
    monkeypatch.setattr(FakePopen, "ret", 3)
    with pytest.raises(RuntimeError, match=r"failed \(3\)"):
        list(SshTarget("raspi").run(["false"]))
    assert not (tmp_path / "never.csv").exists()


def test_run_shell_quotes_every_remote_token(fake_ssh):
    """Meta values with spaces must survive the remote shell intact."""
    cmd = ["fftbench/wt-var/build/benchmarks/bench_pffft_float",
           "--size", "64", "--runs", "2",
           "--samples", "/tmp/x.csv",
           "--meta", "host=motorola edge 50 neo",
           "--meta", "compiler=arm-linux-gnueabihf-gcc (9.0)"]
    list(SshTarget("raspi").run(cmd))
    a = FakePopen.last_argv
    assert a[:2] == ["ssh", "raspi"]
    # unquoting the joined remote command recovers argv exactly
    # (--samples path is swapped for "-" on the wire)
    expect = list(cmd)
    expect[expect.index("--samples") + 1] = "-"
    assert shlex.split(" ".join(a[2:])) == expect


def test_run_appends_rows_to_existing_samples_file(tmp_path, fake_ssh):
    """perf.py reuses one file across invocations; headers never duplicate."""
    out = tmp_path / "samples" / "acc.csv"
    tgt = SshTarget("raspi")
    cmd = ["fftbench/wt-var/build/benchmarks/bench_pffft_float",
           "--size", "64", "--runs", "2", "--samples", str(out)]
    list(tgt.run(cmd))                      # first run creates the file
    list(tgt.run(cmd))                      # second run appends rows only
    prov, rows = read_samples(out)
    assert prov["host"] == "fake-pi"
    assert [r.rep for r in rows] == [0, 1]  # accumulated, not overwritten
    text = out.read_text()
    assert text.count("# pffft-bench-samples v2") == 1
    assert text.count("algo,prec,xform,size,sample_ms,n_iter") == 1


# ---- cmd_ab drives SshTarget unchanged (fully faked, no live run) ---------------

class AbFakePopen(FakePopen):
    """Emits plausible CSV for remote bench runs; stubs every other process."""

    @staticmethod
    def _lines(argv):
        joined = " ".join(argv)
        if argv[0] == "ssh":
            if "uname" in joined:
                return ["octopi\n"]
            if "cmake" in joined or "mkdir" in joined:
                return []
            return AbFakePopen._bench_csv(argv[2:])
        if argv[:1] == ["rsync"]:
            return []
        return None                       # cc --version etc. -> communicate()

    @staticmethod
    def _bench_csv(cmd):
        size = int(cmd[cmd.index("--size") + 1])
        prec = "dbl" if "double" in cmd[0] else "flt"
        arm = Path(cmd[0]).parts[1]       # fftbench/<arm>/build/...
        ms = {"wt-base": 1.50, "wt-var": 1.45}[arm]
        metas, it = [], iter(range(len(cmd)))
        for i in it:
            if cmd[i] == "--meta":
                metas.append(f"# {cmd[i + 1]}")
                next(it, None)
        rows = [f"pffft,{prec},{xf},{size},{ms + j * 0.01},100000"
                for xf in ("real", "cplx") for j in range(3)]
        return ["# pffft-bench-samples v2", *metas,
                "algo,prec,xform,size,sample_ms,n_iter", *rows]


    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False
    def __init__(self, argv, stdout=None, stderr=None, stdin=None,
                 cwd=None, text=None, **kw):
        self.argv = self.args = list(argv)
        self._out = self._lines(list(argv))
        self.stdout = iter(self._out or []) if stdout is not None else None

    def communicate(self, input=None, timeout=None):
        if self.argv[0] == "cc":
            return "arm-linux-gnueabihf-gcc (fake) 9.0\n", ""
        if self._out is None:
            return "", ""
        return "\n".join(l.rstrip("\n") for l in self._out) + "\n", ""


def test_cmd_ab_with_ssh_target(tmp_path, monkeypatch):
    recorded = []
    real_run = subprocess.run

    def fake_run(cmd, **kw):
        if cmd[0] in ("rsync", "ssh"):
            recorded.append(list(cmd))
            return subprocess.CompletedProcess(cmd, 0, stdout="octopi\n")
        return real_run(cmd, **kw)

    monkeypatch.setattr(subprocess, "Popen", AbFakePopen)
    monkeypatch.setattr(subprocess, "run", fake_run)

    wts = {}
    for name, tree in (("wt-base", "aaaa1111"), ("wt-var", "bbbb2222")):
        bindir = tmp_path / name / "build" / "benchmarks"
        bindir.mkdir(parents=True)
        wts[name] = (bindir, tree)

    monkeypatch.setattr("perf.ensure_worktrees",
                        lambda: (tmp_path / "wt-base", tmp_path / "wt-var"))
    (tmp_path / "perf").mkdir()
    monkeypatch.setattr("perf.sync_worktree",
                        lambda wt, ref: wts[wt.name][1])
    monkeypatch.setattr("perf.build", lambda wt: wts[wt.name][0])
    monkeypatch.setattr("perf._dirty", lambda ref: False)
    monkeypatch.setattr("perf.SAMPLES_DIR", tmp_path / "samples")
    monkeypatch.setattr("perf.PERF_DIR", tmp_path / "perf")

    args = argparse.Namespace(
        label="sshcheck", base_label="base", base="HEAD", variant="HEAD",
        runs=2, invocations=1, prec="flt", sizes="short", max_len=64,
        target=["ssh://pi@raspi"])
    cmd_ab(args)

    # one rsync per arm pushed to ~/fftbench/<worktree>/ on the Pi
    rsyncs = sorted(c[-1] for c in recorded if c[0] == "rsync")
    assert [r.split(":")[1] for r in rsyncs] == \
        ["fftbench/wt-base/", "fftbench/wt-var/"]
    # one remote cmake build per arm
    builds = [c for c in recorded
              if c[0] == "ssh" and "cmake -S" in c[2]]
    arms = {re.search(r"fftbench/([\w-]+)", c[2]).group(1) for c in builds}
    assert arms == {"wt-base", "wt-var"}
    # sample files land locally under SAMPLES_DIR, named by REMOTE host
    files = sorted((tmp_path / "samples").glob("*.csv"))
    assert {f.name for f in files} == \
        {"base-raspi-aaaa1111.csv", "sshcheck-raspi-bbbb2222.csv"}
    prov, rows = read_samples(files[-1])
    assert prov["host"] == "octopi"            # REMOTE provenance, not Mac's
    assert prov["target"] == "ssh://pi@raspi"
    # analysis ran end-to-end over the captured files
    ab = json.loads((tmp_path / "perf" / "last_ab.json").read_text())
    assert ab["targets"] == ["raspi"]
    assert len(ab["comparisons"]) > 0



# ---- adb: spec parsing / dispatch ----------------------------------------------

def test_adb_spec_without_serial():
    t = make_target("adb")
    assert isinstance(t, AdbTarget)
    assert t.name == "adb"
    assert t.adb("devices") == ["adb", "devices"]


def test_adb_spec_with_serial():
    t = make_target("adb:ZY22KQZVPM")
    assert t.serial == "ZY22KQZVPM"
    assert t.name == "adb:ZY22KQZVPM"
    assert t.adb("shell", "x") == ["adb", "-s", "ZY22KQZVPM", "shell", "x"]


def test_adb_empty_serial_rejected():
    with pytest.raises(ValueError):
        make_target("adb:")


# ---- adb: provenance meta -------------------------------------------------------

def test_adb_meta_from_device_props(monkeypatch):
    props = {}
    def fake_sh(*cmd):
        val = {"ro.product.model": "motorola edge 50 neo",
               "ro.board.platform": "mt6878"}[cmd[-1]]
        props[cmd[-1]] = val
        return val + "\n"
    monkeypatch.setattr(targets, "sh", fake_sh)
    t = make_target("adb")
    assert t.meta() == {"target": "adb",
                        "host": "motorola_edge_50_neo_mt6878"}
    assert set(props) == {"ro.product.model", "ro.board.platform"}


# ---- adb: build reuses cross_build_android's shared cmake factory ---------------

def test_adb_build_uses_shared_cmake_argv(tmp_path, monkeypatch):
    ndk = tmp_path / "ndk"
    toolchain = ndk / "build" / "cmake" / "android.toolchain.cmake"
    toolchain.parent.mkdir(parents=True)
    toolchain.write_text("")
    monkeypatch.setattr(targets, "find_ndk", lambda: ndk)
    monkeypatch.setattr(targets, "detect_device_march",
                        lambda adb, abi: "armv8.2-a+dotprod")

    ran = []
    def fake_run(argv, check=False, **kw):
        ran.append([str(a) for a in argv])
        class R: returncode = 0
        return R()
    monkeypatch.setattr(targets.subprocess, "run", fake_run)

    wt = tmp_path / "wt-var"
    wt_build = wt / "build" / "benchmarks"     # shape perf.py hands us
    wt_build.mkdir(parents=True)
    bdir = wt / "build-android" / "benchmarks"
    bdir.mkdir(parents=True)
    for f in ("bench_pffft_float", "bench_pffft_double"):
        (bdir / f).write_text("")              # pretend the build produced them
    t = make_target("adb:SER")
    exe = t.binary_path(wt_build, "flt")

    # import reuse: the factory is cross_build_android's, not a copy
    import cross_build_android
    assert targets.pffft_android_cmake_argv \
        is cross_build_android.pffft_android_cmake_argv

    configure = next(c for c in ran if c[:2] == ["cmake", "-S"])
    assert f"-DCMAKE_TOOLCHAIN_FILE={toolchain}" in configure
    assert "-DANDROID_ABI=arm64-v8a" in configure
    assert "-DTARGET_C_ARCH=armv8.2-a+dotprod" in configure
    build = next(c for c in ran if c[:2] == ["cmake", "--build"])
    assert str(wt / "build-android") in build
    pushes = [c for c in ran if "push" in c]
    assert any(str(wt / "build-android" / "benchmarks" /
                   "bench_pffft_float") in " ".join(c) for c in pushes)
    # binary_path resolves to the on-device location under DEVICE_DIR
    assert str(exe) == "/data/local/tmp/wt-var/bench_pffft_float"


def test_adb_prepare_raises_without_ndk(monkeypatch, tmp_path):
    monkeypatch.setattr(targets, "find_ndk", lambda: None)
    wt_build = tmp_path / "wt-base" / "build" / "benchmarks"
    with pytest.raises(RuntimeError, match="NDK"):
        AdbTarget()._prepare(wt_build)


# ---- adb: run(): stdout capture -> local CSV -------------------------------------

def test_adb_run_quotes_and_captures_samples(tmp_path, fake_popen):
    out = tmp_path / "samples" / "y.csv"
    cmd = ["/data/local/tmp/wt-base/bench_pffft_float",
           "--size", "64", "--runs", "2",
           "--samples", str(out),
           "--meta", "host=moto g(9)"]
    lines = list(make_target("adb").run(cmd))
    assert lines[0] == "# pffft-bench-samples v2"
    assert lines[-1] == "pffft,flt,real,64,1.5,100000"
    # --samples was rewritten to "-" and values quoted for the device shell
    joined = " ".join(FakePopen.last_argv)
    assert "--samples -" in joined
    assert "'host=moto g(9)'" in joined
    rows = read_samples(out)[1]
    assert rows[0].sample_ms == 1.5


def test_adb_run_appends_data_rows_to_existing_file(tmp_path, monkeypatch):
    out = tmp_path / "samples" / "acc.csv"
    cmd = ["/data/local/tmp/wt-base/bench_pffft_float",
           "--size", "64", "--runs", "2",
           "--samples", str(out)]
    t = make_target("adb")

    def popen_emitting(n_data_rows):
        class P(FakePopen):
            def __init__(self, argv, **kw):
                super().__init__(argv, **kw)
                self.stdout = iter(
                    ["# pffft-bench-samples v2\n",
                     "# host=fake-pi\n",
                     "algo,prec,xform,size,sample_ms,n_iter\n"]
                    + [f"pffft,flt,real,64,{1.5 + i},100000\n"
                       for i in range(n_data_rows)])
        return P

    monkeypatch.setattr(subprocess, "Popen", popen_emitting(1))
    list(t.run(cmd))                       # creates file: header + 1 row
    monkeypatch.setattr(subprocess, "Popen", popen_emitting(2))
    list(t.run(cmd))                       # appends data rows only

    assert open(out).read().count("# pffft-bench-samples v2") == 1
    rows = read_samples(out)[1]
    assert len(rows) == 3                  # 1 + 2 accumulated
    groups = {}
    for r in rows:
        k = (r.algo, r.prec, r.xform, r.size)
        groups.setdefault(k, []).append(r)
    assert [len(v) for v in groups.values()] == [3]
    # rep indices for the single accumulated group are 0, 1, 2 in file order
    (rows_,) = groups.values()
    assert [(r.algo, r.prec, r.xform, r.size) for r in rows_] \
        == [(r.algo, r.prec, r.xform, r.size) for r in rows_[:1]] * 3
    assert [r.sample_ms for r in rows_] == [1.5, 1.5, 2.5]


# ---- ios: spec parsing / dispatch ----------------------------------------------

def test_ios_spec_without_udid():
    t = make_target("ios")
    assert isinstance(t, IosTarget)
    assert t.name == "ios"


def test_ios_spec_with_udid():
    t = make_target("ios:00008030-00116DEC0CF0802E")
    assert isinstance(t, IosTarget)
    assert t.udid == "00008030-00116DEC0CF0802E"
    assert t.name == "ios:00008030-00116DEC0CF0802E"


def test_ios_empty_udid_rejected():
    with pytest.raises(ValueError):
        make_target("ios:")


# ---- ios: provenance meta -------------------------------------------------------
def test_ios_bare_spec_name_stable_across_udid_resolution(monkeypatch):
    """Auto-detected UDID must never leak into name (perf.py derives
    sample-file paths from name before meta()/run() resolve the device)."""
    monkeypatch.setattr(targets, "list_connected_devices",
                        lambda: [{"id": "AUTO-UDID", "name": "iPhone"}])
    t = make_target("ios")
    assert t.name == "ios"
    assert t._resolve_udid() == "AUTO-UDID"
    assert t.name == "ios"          # unchanged after resolution


def test_ios_meta_condenses_device_info_through_sanitizer(monkeypatch):
    monkeypatch.setattr(targets, "collect_device_info", lambda udid: {
        "device_name": "iPhone 11 Pro",
        "hardware": "D421AP",
        "ios_version": "17.5\n"})
    assert make_target("ios:F00D").meta() == \
        {"target": "ios:F00D", "host": "iPhone_11_Pro_D421AP_iOS_17.5_"}


# ---- ios: build reuses cross_build_ios's shared functions ------------------------

def test_ios_prepare_uses_shared_signing_and_build(tmp_path, monkeypatch):
    calls = {}

    def fake_build(script_dir, build_dir, deployment_target, fftw_cmake,
                   extra_cmake, cpus, **kw):
        calls.update(script_dir=script_dir, build_dir=build_dir,
                     identity=kw.get("signing_identity"),
                     team_id=kw.get("team_id"),
                     profile=kw.get("provisioning_profile"))

    monkeypatch.setattr(targets, "find_ios_sdk", lambda: Path("/sdk"))
    monkeypatch.setattr(targets, "find_signing_identity",
                        lambda: "Apple Development: x")
    monkeypatch.setattr(targets, "find_team_id", lambda ident: "TEAM1234")
    monkeypatch.setattr(targets, "build_pffft_ios", fake_build)
    wt = tmp_path / "wt-var"
    (wt / "build" / "benchmarks").mkdir(parents=True)
    IosTarget("UDID")._prepare(wt / "build" / "benchmarks")
    assert calls["script_dir"] == wt   # build from THIS arm's worktree, not REPO_ROOT
    assert calls["build_dir"] == wt / "build-ios"
    assert calls["identity"] == "Apple Development: x"
    assert calls["team_id"] == "TEAM1234"
    assert calls["profile"] is None      # automatic provisioning


def test_ios_prepare_builds_each_arm_from_its_own_worktree(tmp_path, monkeypatch):
    """Regression: _prepare once built every arm from REPO_ROOT, so an A/B
    comparison silently compiled and ran the IDENTICAL binary for both arms.
    Each worktree (wt-base/wt-var) must be built as its own cmake source."""
    seen_script_dirs = []

    def fake_build(script_dir, build_dir, deployment_target, fftw_cmake,
                   extra_cmake, cpus, **kw):
        seen_script_dirs.append(script_dir)

    monkeypatch.setattr(targets, "find_ios_sdk", lambda: Path("/sdk"))
    monkeypatch.setattr(targets, "find_signing_identity",
                        lambda: "Apple Development: x")
    monkeypatch.setattr(targets, "find_team_id", lambda ident: "TEAM1234")
    monkeypatch.setattr(targets, "build_pffft_ios", fake_build)
    base = tmp_path / "wt-base"
    var = tmp_path / "wt-var"
    (base / "build" / "benchmarks").mkdir(parents=True)
    (var / "build" / "benchmarks").mkdir(parents=True)
    tgt = IosTarget("UDID")
    tgt._prepare(base / "build" / "benchmarks")
    tgt._prepare(var / "build" / "benchmarks")
    assert seen_script_dirs == [base, var]
    assert seen_script_dirs[0] != seen_script_dirs[1]


def test_ios_prepare_raises_without_sdk(monkeypatch, tmp_path):
    monkeypatch.setattr(targets, "find_ios_sdk", lambda: None)
    wt = tmp_path / "wt-base" / "build" / "benchmarks"
    wt.mkdir(parents=True)
    with pytest.raises(RuntimeError):
        IosTarget()._prepare(wt)


# ---- ios: binary_path resolves the signed .app bundle ----------------------------

def test_ios_binary_path_is_app_bundle(tmp_path):
    wt = tmp_path / "wt-var"
    bundle = (wt / "build-ios" / "pffft" / "benchmarks"
              / "Release-iphoneos" / "bench_pffft_float.app")
    bundle.mkdir(parents=True)

    class T(IosTarget):
        def _prepare(self, wt_build):
            pass

    assert T("UDID").binary_path(wt / "build" / "benchmarks", "flt") == bundle


# ---- ios: run() deploys via deploy_and_run and captures stdout locally -------------

IOS_STDOUT = (
    "\x1b[1G ios-deploy install chatter\n"
    "[100%] Installed package\n"
    "# pffft-bench-samples v2\n"
    "# host=fake-phone\n"
    "algo,prec,xform,size,sample_ms,n_iter\n"
    "pffft,flt,real,64,1.5,100000\n"
    "(lldb) quit\n")


def test_ios_run_captures_samples_stdout_to_local_file(tmp_path, monkeypatch):
    out = tmp_path / "samples" / "z.csv"
    cmd = ["wt-base/build-ios/pffft/benchmarks/Release-iphoneos/"
           "bench_pffft_float.app",
           "--size", "64", "--runs", "2",
           "--samples", str(out),
           "--meta", "host=iPhone"]
    recorded = {}

    def fake_deploy(bundle, udid, app_args=None):
        recorded.update(bundle=bundle, udid=udid, args=list(app_args))
        return IOS_STDOUT

    monkeypatch.setattr(targets, "deploy_and_run", fake_deploy)
    lines = list(make_target("ios:F00D").run(cmd))
    # ios-deploy chatter stripped; only the benchmark's own CSV survives
    assert lines[0] == "# pffft-bench-samples v2"
    assert lines[-1] == "pffft,flt,real,64,1.5,100000"
    # app launched through --args with --samples rewritten to "-"
    i = recorded["args"].index("--samples")
    assert recorded["args"][i + 1] == "-"
    assert recorded["udid"] == "F00D"
    rows = read_samples(out)[1]
    assert rows[0].sample_ms == 1.5


def test_ios_run_appends_data_rows_to_existing_file(tmp_path, monkeypatch):
    """Accumulation semantics match Ssh/Adb: data rows only, never a header."""
    out = tmp_path / "samples" / "acc.csv"
    (tmp_path / "samples").mkdir()
    out.write_text("# pffft-bench-samples v2\n"
                   "# host=x\n"
                   "algo,prec,xform,size,sample_ms,n_iter\n"
                   "pffft,flt,real,64,1.5,100000\n")
    cmd = ["b.app", "--size", "64", "--runs", "2", "--samples", str(out)]
    monkeypatch.setattr(
        targets, "deploy_and_run",
        lambda bundle, udid, app_args=None:
            IOS_STDOUT.replace("fake-phone", "x"))
    list(IosTarget("U").run(cmd))
    assert open(out).read().count("# pffft-bench-samples v2") == 1
    assert [r.sample_ms for r in read_samples(out)[1]] == [1.5, 1.5]


def test_ios_run_raises_on_deploy_failure(tmp_path, monkeypatch):
    out = tmp_path / "never.csv"
    cmd = ["b.app", "--size", "64", "--samples", str(out)]
    monkeypatch.setattr(targets, "deploy_and_run",
                        lambda bundle, udid, app_args=None: None)
    with pytest.raises(RuntimeError):
        list(IosTarget("U").run(cmd))
    assert not out.exists()


# ---- ios: shared argv builder lives in cross_build_ios -----------------------------
def test_ios_deploy_argv_sanitizes_unsafe_chars_in_app_args():
    """ios-deploy's --args re-splitter breaks on spaces/shell metacharacters
    even when shell-quoted (verified against a real device); every token's
    unsafe characters must become '_' so it always survives as one word."""
    from cross_build_ios import ios_deploy_argv
    argv = ios_deploy_argv(Path("/a/B.app"), "UDID",
                           ["--runs", "2",
                            "--meta", "compiler=Apple clang (x)"])
    assert argv[:7] == ["ios-deploy", "--bundle", "/a/B.app",
                        "--id", "UDID", "--noninteractive", "--debug"]
    assert argv[7] == "--args"
    assert argv[8] == "--runs 2 --meta compiler=Apple_clang__x_"
    assert " " not in "compiler=Apple_clang__x_"