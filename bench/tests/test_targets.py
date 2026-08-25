import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import targets
from samples import read_samples
from targets import AdbTarget, make_target


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
def fake_popen(monkeypatch):
    monkeypatch.setattr(subprocess, "Popen", FakePopen)


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
                        "host": "motorola edge 50 neo mt6878"}
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
