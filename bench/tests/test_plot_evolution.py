import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_evolution import load_versions, select_mflops


def _write(path, label, size, sample_ms_list):
    rows = "\n".join(
        f"pffft,flt,real,{size},{ms},100000" for ms in sample_ms_list)
    path.write_text(
        "# pffft-bench-samples v2\n"
        f"# label={label}\n"
        "algo,prec,xform,size,sample_ms,n_iter\n"
        f"{rows}\n")


def test_load_versions_groups_by_label(tmp_path):
    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    _write(a, "base", 256, [150.0, 151.0])
    _write(b, "opt", 256, [140.0])

    versions = load_versions([str(a), str(b)])
    assert set(versions) == {"base", "opt"}
    assert len(versions["base"]) == 2
    assert len(versions["opt"]) == 1


def test_load_versions_uses_filename_stem_without_label(tmp_path):
    p = tmp_path / "unlabeled.csv"
    p.write_text(
        "# pffft-bench-samples v2\n"
        "algo,prec,xform,size,sample_ms,n_iter\n"
        "pffft,flt,real,64,150.0,100000\n")
    versions = load_versions([str(p)])
    assert list(versions) == ["unlabeled"]


def test_load_versions_extends_same_label_across_files(tmp_path):
    a = tmp_path / "a1.csv"
    b = tmp_path / "a2.csv"
    _write(a, "same", 256, [150.0])
    _write(b, "same", 256, [151.0])
    versions = load_versions([str(a), str(b)])
    assert len(versions["same"]) == 2


def test_select_mflops_filters_exactly(tmp_path):
    p = tmp_path / "s.csv"
    _write(p, "x", 256, [150.0, 150.0])
    versions = load_versions([str(p)])
    vals = select_mflops(versions["x"], "pffft", "flt", "real", 256)
    assert len(vals) == 2
    assert select_mflops(versions["x"], "pffft", "flt", "real", 999) == []
    assert select_mflops(versions["x"], "vdsp", "flt", "real", 256) == []
