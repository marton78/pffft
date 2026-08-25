"""A/B runtime statistics: MAD trim + Mann-Whitney U + Hodges-Lehmann shift."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
from scipy.stats import mannwhitneyu

ALPHA = 0.01
WARMUP = 2              # Design Contract verdict rule: trim warmup 2 reps
TRIM_MAD = 5.0
SIG_FRAC_FASTER = 0.7   # share of groups that must improve for "faster"


def trim(x: Sequence[float], warmup: int = WARMUP,
         mad_mult: float = TRIM_MAD) -> np.ndarray:
    """Drop warmup reps, then MAD outliers (median +/- mad_mult scaled-MAD).

    When the scaled MAD is 0 (degenerate, e.g. constant runtimes plus one
    spike) the threshold is 0, so only samples equal to the median survive.
    """
    a = np.asarray(x, dtype=float)[warmup:]
    med = np.median(a)
    mad = np.median(np.abs(a - med)) * 1.4826
    return a[np.abs(a - med) <= mad_mult * mad]


def _hl_location(x: np.ndarray) -> float:
    """Hodges-Lehmann location estimate: median of pairwise averages."""
    return float(np.median((x[:, None] + x[None, :]).ravel() / 2.0))


def hodges_lehmann_shift(a: np.ndarray, b: np.ndarray) -> float:
    """Shift of b relative to a as fraction of |a|: median of (ai+bj)/2 pairs
    vs a's HL location. Negative => b (variant) has smaller runtimes."""
    pair_med = float(np.median(((a[:, None] + b[None, :]) / 2.0).ravel()))
    loc_a = _hl_location(a)
    return (pair_med - loc_a) / abs(loc_a) if loc_a != 0 else 0.0


@dataclass
class Comparison:
    p_value: float
    shift_pct: float
    faster: bool | None     # None = not significant
    n_a: int
    n_b: int


def compare(base: Sequence[float], variant: Sequence[float],
            alpha: float = ALPHA) -> Comparison:
    """Compare trimmed base vs variant runtimes (lower is better).
    Variant is significantly faster iff p < alpha AND Hodges-Lehmann
    shift < 0; significantly slower iff p < alpha AND shift > 0;
    `faster` is None when p >= alpha. The two-sided test flags both
    improvements and pessimizations.
    """
    a, b = trim(base), trim(variant)
    if len(a) < 3 or len(b) < 3:
        return Comparison(1.0, 0.0, None, len(a), len(b))
    # H0: distributions equal; two-sided so significant slowdowns are flagged too
    res = mannwhitneyu(b, a, alternative="two-sided")
    shift = hodges_lehmann_shift(a, b)
    if res.pvalue >= alpha:
        faster = None
    else:
        faster = shift < 0     # negative shift = variant runtimes smaller = faster
    return Comparison(float(res.pvalue), shift * 100.0, faster,
                      len(a), len(b))


def aggregate(comparisons: Iterable[Comparison]) -> dict:
    comps = list(comparisons)
    nf = sum(c.faster is True for c in comps)
    ns = sum(c.faster is False for c in comps)
    if ns > 0:
        suggest = "slower"
    elif nf >= SIG_FRAC_FASTER * max(len(comps), 1):
        suggest = "faster"
    else:
        suggest = "neutral"
    return {"n_groups": len(comps), "n_sig_faster": nf, "n_sig_slower": ns,
            "suggest": suggest}
