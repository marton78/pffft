import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from perf import balance_order, pick_sizes, collect_group, _series
from stats import compare


def test_balance_order_is_balanced():
    for seed in range(20):
        seq = balance_order(8, random.Random(seed))
        assert seq.count("A") == 4 and seq.count("B") == 4

def test_balance_order_even_counts_only():
    # Contract: sequence length is 2*(n//2) — callers always pass an even n.
    seq = balance_order(10, random.Random(1))
    assert len(seq) == 10
    assert set(seq) <= {"A", "B"}

def test_balance_order_shuffles_deterministically():
    assert balance_order(8, random.Random(42)) == balance_order(8, random.Random(42))


def test_pick_sizes_short():
    s = pick_sizes("short")
    assert 32 in s and 2048 in s and 96 in s and max(s) <= 2048
    assert {32, 64, 128, 256, 512, 1024, 2048} <= set(s)

def test_pick_sizes_short_is_sorted_union_of_pow2_and_nonpow2():
    s = pick_sizes("short")
    assert s == sorted(s)
    assert set(SHORT := [32, 64, 128, 256, 512, 1024, 2048]) <= set(s)
    assert {96, 192, 480} <= set(s)


def test_pick_sizes_pow2_spans_32_to_2_21():
    p = pick_sizes("pow2")
    assert p[0] == 32 and p[-1] == 32 << 16
    assert all(n & (n - 1) == 0 for n in p)
    assert all(32 <= n <= (1 << 21) for n in p)


def test_pick_sizes_nonpow2_matches_c_binary_active_list():
    npow2 = pick_sizes("nonpow2")
    assert 96 in npow2 and 15360 in npow2
    assert all(n & (n - 1) for n in npow2)


def test_pick_sizes_all_is_sorted_union():
    assert pick_sizes("all") == sorted(
        set(pick_sizes("pow2")) | set(pick_sizes("nonpow2")))


def test_pick_sizes_unknown_mode_raises():
    with pytest.raises(ValueError):
        pick_sizes("bogus")


# ---- collect_group/_series feed compare() the right metric --------------------
# Regression: collect_group once collected sample_ms, the C binary's fixed
# ~150ms calibrated timing window (see bench_pffft.c) that stays roughly
# constant regardless of code speed -- the actual throughput lives in n_iter
# (derived MFLOPS). Comparing sample_ms directly meant every ab verdict came
# back "neutral" no matter how much faster or slower the variant really was.

def _write_samples(path, label, n_iter_values, sample_ms=150.0, size=1024):
    rows = "\n".join(
        f"pffft,flt,real,{size},{sample_ms},{n}" for n in n_iter_values)
    content = ("# pffft-bench-samples v2\n"
              f"# label={label}\n"
               "algo,prec,xform,size,sample_ms,n_iter\n"
              f"{rows}\n")
    path.write_text(content)


def test_collect_group_compares_throughput_not_calibration_window(tmp_path):
    base = tmp_path / "base.csv"
    var = tmp_path / "var.csv"
    # IDENTICAL sample_ms (both hit the same ~150ms calibration window) but
    # the variant does clearly MORE iterations per window => genuinely faster.
    _write_samples(base, "base", [100_000 + i * 10 for i in range(20)])
    _write_samples(var, "var", [120_000 + i * 10 for i in range(20)])

    sb = _series(collect_group(base)[0], "flt")
    sv = _series(collect_group(var)[0], "flt")
    b, v = sb[("pffft", "real")][1024], sv[("pffft", "real")][1024]
    c = compare(b, v)
    assert c.faster is True, (
        "collect_group must feed compare() derived throughput (MFLOPS), "
        "not the near-constant calibrated sample_ms -- a real speedup must "
        "be detectable even when every sample_ms is ~identical")
    assert c.shift_pct < -5
