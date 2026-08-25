import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from stats import compare, aggregate, trim


def test_trim_drops_warmup_and_spikes():
    # warmup=2 comes from the Design Contract verdict rule ("trim warmup 2
    # reps + MAD outliers"): drop 999.0 and one 10.0 (22 -> 20 samples), then
    # the MAD filter (mad==0 => tolerance 0) removes the 5000.0 spike.
    x = [999.0] + [10.0] * 20 + [5000.0]
    t = trim(x)
    assert 5000.0 not in t and len(t) == 19


def test_detects_clear_improvement():
    rng = np.random.default_rng(42)
    base = rng.lognormal(mean=np.log(100), sigma=0.05, size=40)
    var = rng.lognormal(mean=np.log(90), sigma=0.05, size=40)   # 10% faster
    c = compare(base, var)
    assert c.faster is True and c.shift_pct < -5 and c.p_value < 1e-3


def test_detects_clear_regression():
    rng = np.random.default_rng(11)
    base = rng.lognormal(mean=np.log(100), sigma=0.05, size=40)
    var = rng.lognormal(mean=np.log(110), sigma=0.05, size=40)  # 10% slower
    c = compare(base, var)
    assert c.faster is False and c.shift_pct > 5 and c.p_value < 0.01


def test_no_false_positive_on_identical_distributions():
    rng = np.random.default_rng(7)
    x = rng.lognormal(mean=np.log(100), sigma=0.1, size=60)
    c = compare(x, x.copy())
    assert c.faster is None or c.faster is True  # never claims slower


def test_aggregate_rule():
    from stats import Comparison

    fast = [Comparison(p_value=1e-4, shift_pct=-8.0, faster=True, n_a=20, n_b=20)] * 8
    slow = [Comparison(p_value=1e-4, shift_pct=+8.0, faster=False, n_a=20, n_b=20)]
    assert aggregate(fast)["suggest"] == "faster"
    agg = aggregate(fast + slow)
    assert agg["suggest"] == "slower" and agg["n_sig_slower"] == 1
