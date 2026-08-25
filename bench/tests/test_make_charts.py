import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import make_charts as mc


def samples_file(tmp_path, name='s.csv', label='base', tree='aaa',
                 target='local', rows=(
                     'pffft,flt,real,64,152.0,412300',
                     'pffft,flt,real,64,150.0,412300',
                     'pffft,flt,real,64,154.0,412300',
                     'pffftu,flt,real,64,180.0,412300')):
    header = (
        '# pffft-bench-samples v2\n'
        f'# label={label}\n'
        f'# git_tree={tree}\n'
        f'# target={target}\n'
        'algo,prec,xform,size,sample_ms,n_iter\n'
    )
    p = tmp_path / name
    p.write_text(header + '\n'.join(rows) + '\n')
    return p


def test_read_samples_file_medians_over_reps(tmp_path):
    p = samples_file(tmp_path)
    data = mc.read_samples_file(p)
    assert set(data) == {'pffft', 'pffftu'}
    sizes, medians = data['pffft'][('flt', 'real')]
    assert sizes == [64]
    assert medians == [5208.0]  # derived MFLOPS at median sample_ms=152


def test_read_samples_file_rejects_non_magic(tmp_path):
    p = tmp_path / 'junk.csv'
    p.write_text('garbage,not,samples\n')
    assert mc.read_samples_file(p) is None
    assert mc.read_provenance(p) is None


def test_scan_directory_skips_non_magic_and_merges(tmp_path):
    samples_file(tmp_path, name='a.csv', label='l1')
    samples_file(tmp_path, name='b.csv', label='l2', tree='bbb')
    (tmp_path / 'junk.csv').write_text('nope\n')
    prov, panels = mc.scan_directory(tmp_path)
    assert prov['label'] == 'l1'  # first readable file's header
    # same algo from both files: one raw series per contributing file
    series = panels[('flt', 'real')]
    pffft = [s for s in series if s[0] == 'pffft']
    assert len(pffft) == 2
    assert all(s[2] == [64] and s[3] == [5208.0] for s in pffft)
    # merge_dirs folds them into a single per-size-median series
    _, mpanels = mc.merge_dirs([tmp_path])
    merged = [s for s in mpanels[('flt', 'real')] if s[0] == 'pffft']
    assert merged == [('pffft', 'default', [64], [5208.0])]


def test_draw_combined_single_panel():
    panels = [(('flt', 'real'), [('pffft', 'default', 'PFFFT',
                                  [32, 64], [100.0, 200.0], 0)])]
    fig = mc.draw_combined(panels, {'label': 'x'}, 'title')
    assert len(fig.axes) == 1
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_draw_combined_two_panels():
    panels = [(k, [('pffft', 'default', 'PFFFT', [32, 64], [100.0, 200.0], 0)])
              for k in (('flt', 'real'), ('flt', 'cplx'))]
    fig = mc.draw_combined(panels, None, 'title')
    assert len(fig.axes) == 2
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_filter_head_files_by_tree(tmp_path):
    head = samples_file(tmp_path, name='head.csv', label='opt', tree='aaa')
    base = samples_file(tmp_path, name='base.csv', label='base', tree='bbb')
    kept = mc.filter_head_files([base, head], 'aaa', 'opt')
    assert kept == [head]


def test_filter_head_files_label_fallback(tmp_path):
    # no git_tree provenance on either side -> match by label
    p = tmp_path / 's.csv'
    p.write_text(textwrap.dedent("""\
        # pffft-bench-samples v2
        # label=opt
        algo,prec,xform,size,sample_ms,n_iter
        pffft,flt,real,64,152.0,412300
        """))
    other = tmp_path / 'o.csv'
    other.write_text(textwrap.dedent("""\
        # pffft-bench-samples v2
        # label=base
        algo,prec,xform,size,sample_ms,n_iter
        pffft,flt,real,64,152.0,412300
        """))
    kept = mc.filter_head_files([other, p], '', 'opt')
    assert kept == [p]
