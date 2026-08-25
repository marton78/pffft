import sys, textwrap
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from samples import SCHEMA_MAGIC, read_samples, provenance_diff, group_key

HEADER = textwrap.dedent("""\
    # pffft-bench-samples v2
    # label=base
    # git_tree=aaa
    # datetime=2026-08-25T09:00Z
    # simd_arch=4xNEON
    algo,prec,xform,size,sample_ms,n_iter
    pffft,flt,real,64,152.0,412300
    pffft,flt,real,64,150.0,412300
    pffftu,flt,real,64,180.0,412300
    pffft,flt,cplx,64,290.0,206150
""")

def write(tmp_path):
    p = tmp_path / "s.csv"
    p.write_text(HEADER)
    return p

def test_reads_rows_and_rep(tmp_path):
    prov, rows = read_samples(write(tmp_path))
    assert prov["label"] == "base" and prov["simd_arch"] == "4xNEON"
    assert len(rows) == 4
    assert rows[0].rep == 0 and rows[1].rep == 1 and rows[2].rep == 0
    assert group_key(rows[3]) == ("pffft", "flt", "cplx", 64)

def test_mflops_matches_reference_formula(tmp_path):
    import math
    prov, rows = read_samples(write(tmp_path))
    r = rows[0]
    expected = 2 * r.n_iter * 2.5 * 64 * (math.log(64) / math.log(2)) / 1e6 / (152.0 / 1000)
    assert abs(r.mflops - expected) < 1e-9
    rc = rows[3]
    exp_c = 2 * rc.n_iter * 5 * 64 * (math.log(64) / math.log(2)) / 1e6 / (290.0 / 1000)
    assert abs(rc.mflops - exp_c) < 1e-9

def test_provenance_diff_ignores_datetime(tmp_path):
    a = {"label": "x", "datetime": "1", "host": "h"}
    b = {"label": "x", "datetime": "2", "host": "H"}
    d = provenance_diff(a, b)
    assert d == {"host": ("h", "H")}
