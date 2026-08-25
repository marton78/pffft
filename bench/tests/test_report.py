import sys, statistics
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from samples import read_samples
import report

def make(tmp_path, name, text):
    p = tmp_path / name
    p.write_text("\n".join(line[4:] for line in text.strip("\n").splitlines()) + "\n")
    return p

# Two synthetic sessions, same (pffft,flt,real,host-a) group, overlapping sizes.
A = """\
    # pffft-bench-samples v2
    # label=base
    # git_tree=aaa111
    # host=host-a
    # datetime=2026-08-25T09:00Z
    algo,prec,xform,size,sample_ms,n_iter
    pffft,flt,real,64,152.0,412300
    pffft,flt,real,64,150.0,412300
    pffft,flt,real,256,700.0,412300
    pffft,flt,cplx,64,290.0,206150
"""
B = """\
    # pffft-bench-samples v2
    # label=opt
    # git_tree=bbb222
    # host=host-a
    # datetime=2026-08-26T10:00Z
    algo,prec,xform,size,sample_ms,n_iter
    pffft,flt,real,64,140.0,412300
    pffft,flt,real,128,300.0,412300
"""

def median_mflops(rows):
    return statistics.median(r.mflops for r in rows)

def test_round_trip_cells_match_medians(tmp_path):
    fa = make(tmp_path, "a.csv", A)
    fb = make(tmp_path, "b.csv", B)
    written = report.render(tmp_path / "out", [fa, fb])
    dest = tmp_path / "out" / "pffft-flt-real-host-a.csv"
    assert dest in written
    lines = dest.read_text().splitlines()
    assert lines[0] == "label,host,tree,datetime,64,128,256"
    base, opt = lines[1], lines[2]
    assert base == "base,host-a,aaa111,2026-08-25T09:00Z,%0.1f,,%0.1f" % (
        median_mflops(read_samples(fa)[1][:2]), median_mflops(read_samples(fa)[1][2:3]))
    assert opt == "opt,host-a,bbb222,2026-08-26T10:00Z,%0.1f,%0.1f," % (
        median_mflops(read_samples(fb)[1][:1]), median_mflops(read_samples(fb)[1][1:2]))

def test_group_split_by_host_and_missing_provenance(tmp_path):
    fa = make(tmp_path, "a.csv", A)
    # no provenance at all -> empty strings, separate group from host-a file
    fb = make(tmp_path, "b.csv", B.replace("    # host=host-a\n", ""))
    written = report.render(tmp_path / "out", [fa, fb])
    names = {w.name for w in written}
    assert names == {"pffft-flt-real-host-a.csv",
                     "pffft-flt-cplx-host-a.csv",
                     "pffft-flt-real-unknown.csv"}
    row = (tmp_path / "out" / "pffft-flt-real-unknown.csv").read_text().splitlines()[1]
    assert row.startswith("opt,,bbb222,")


C1 = """\
    # pffft-bench-samples v2
    # label=space-host
    # git_tree=ccc333
    # host=host a
    # datetime=2026-08-25T09:00Z
    algo,prec,xform,size,sample_ms,n_iter
    pffft,flt,real,64,100.0,412300
"""
C2 = """\
    # pffft-bench-samples v2
    # label=base
    # git_tree=aaa111
    # host=host-a
    # datetime=2026-08-25T09:00Z
    algo,prec,xform,size,sample_ms,n_iter
    pffft,flt,real,64,152.0,412300
"""

def test_sanitized_name_collision_keeps_both_tables(tmp_path):
    f1 = make(tmp_path, "c1.csv", C1)
    f2 = make(tmp_path, "c2.csv", C2)
    written = report.render(tmp_path / "out", [f1, f2])
    out = tmp_path / "out"
    names = {p.name for p in written}
    assert names == {"pffft-flt-real-host-a.csv", "pffft-flt-real-host-a-2.csv"}
    assert len(list(out.iterdir())) == 2
    r1 = (out / "pffft-flt-real-host-a.csv").read_text().splitlines()[1]
    r2 = (out / "pffft-flt-real-host-a-2.csv").read_text().splitlines()[1]
    assert r1.startswith("space-host,host a,ccc333,")   # sorts first, keeps plain name
    assert r2.startswith("base,host-a,aaa111,")         # collides -> "-2" suffix

def test_cli_dirs_and_out_dir(tmp_path):
    d = tmp_path / "samples"
    d.mkdir()
    make(d, "a.csv", A)
    make(d, "b.csv", B)
    out = tmp_path / "tables"
    rc = report.main(["--out-dir", str(out), str(d)])
    assert rc == 0
    assert (out / "pffft-flt-real-host-a.csv").exists()
