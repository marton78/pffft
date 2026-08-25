"""Reader for pffft-bench-samples v2 files (see docs plan 2026-08-25).

MFLOPS is DERIVED here, never trusted from the file, so formula bugs are
fixable post hoc:  mflops = 2*n_iter*(cplx?5:2.5)*N*log2(N)/1e6/(sample_ms/1000)
This reproduces bench_record() in benchmarks/bench_pffft.c.
"""
from __future__ import annotations
import math
import os
import sys
from dataclasses import dataclass
from typing import BinaryIO, Iterable, TextIO

SCHEMA_MAGIC = "# pffft-bench-samples v2"
VOLATILE_KEYS = frozenset({"datetime"})
COLUMNS = ("algo", "prec", "xform", "size", "sample_ms", "n_iter")


@dataclass
class Sample:
    algo: str
    prec: str
    xform: str
    size: int
    sample_ms: float
    n_iter: int
    rep: int = -1

    @property
    def mflops(self) -> float:
        flops_per_fft = 5.0 if self.xform == "cplx" else 2.5
        return (2 * self.n_iter * flops_per_fft * self.size
                * (math.log(self.size) / math.log(2)) / 1e6
                / (self.sample_ms / 1000))


def group_key(s: Sample):
    return (s.algo, s.prec, s.xform, s.size)


def read_samples(path_or_file) -> tuple[dict[str, str], list[Sample]]:
    """Returns (provenance, rows). Raises ValueError on bad magic/columns."""
    fh = open(os.fspath(path_or_file)) if isinstance(path_or_file, (str, bytes, os.PathLike)) else path_or_file
    with fh:
        first = fh.readline().rstrip("\n")
        if first != SCHEMA_MAGIC:
            raise ValueError(f"bad schema line: {first!r}")
        prov: dict[str, str] = {}
        rows: list[Sample] = []
        counts: dict[tuple, int] = {}
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if "=" not in line:
                    continue
                k, _, v = line[2:].partition("=")
                prov[k] = v
                continue
            parts = line.split(",")
            if parts == list(COLUMNS):
                continue
            if len(parts) != len(COLUMNS):
                raise ValueError(f"bad row: {line!r}")
            s = Sample(parts[0], parts[1], parts[2], int(parts[3]),
                       float(parts[4]), int(parts[5]))
            key = group_key(s)
            s.rep = counts.get(key, 0)
            counts[key] = s.rep + 1
            rows.append(s)
    return prov, rows


def provenance_diff(a: dict[str, str], b: dict[str, str]) -> dict[str, tuple[str, str]]:
    keys = (set(a) | set(b)) - VOLATILE_KEYS
    return {k: (a.get(k), b.get(k)) for k in sorted(keys) if a.get(k) != b.get(k)}
