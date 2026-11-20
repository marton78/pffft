"""Read the CSV that bench_fft writes, aggregate over trials, format tables."""
from __future__ import annotations

import statistics
from pathlib import Path

INTS = ("n", "trial", "pairs", "batches", "clean")
FLOATS = ("ns", "ns_median", "cycles", "cycles_median", "instructions", "mflops", "ghz", "clock")
UNITS = {"cycles": "CPU cycles per transform", "ns": "nanoseconds per transform", "mflops": "MFLOPS"}


def extract_csv(text: str) -> str:
    """The lines between #BEGIN-CSV and #END-CSV of console output, which may
    mix in other lines (and carriage returns, from device consoles). Raises
    ValueError if the block is missing or unfinished."""
    lines = [line.rstrip("\r") for line in text.splitlines()]
    try:
        begin = lines.index("#BEGIN-CSV")
        end = lines.index("#END-CSV", begin)
    except ValueError:
        raise ValueError("no complete #BEGIN-CSV ... #END-CSV block in the output") from None
    return "\n".join(lines[begin + 1:end]) + "\n"


def read_csv(source):
    """Read bench_fft CSV from a file (a Path, or a str without line breaks)
    or from text; text may be console output with the CSV between #BEGIN-CSV
    and #END-CSV. Returns (meta, rows): meta maps the header's keys to their
    values (strings); each row is a dict with numbers converted and None for
    empty fields."""
    text = Path(source).read_text() if isinstance(source, Path) or "\n" not in source else source
    if "#BEGIN-CSV" in text:
        text = extract_csv(text)
    meta, rows, columns = {}, [], None
    for number, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if line.startswith("# pffft-bench "):
            if line != "# pffft-bench 1":
                raise ValueError(f"line {number}: unknown format '{line}'")
            meta["format"] = "1"
        elif line.startswith("#"):
            key, sep, value = line[1:].partition(":")
            if sep:
                meta[key.strip()] = value.strip()
        elif line and columns is None:
            columns = line.split(",")
            if "format" not in meta or columns[0] != "engine":
                raise ValueError(f"line {number}: not a pffft-bench CSV")
        elif line:
            fields = line.split(",")
            if len(fields) != len(columns):
                raise ValueError(f"line {number}: {len(fields)} fields, expected {len(columns)}")
            row = dict(zip(columns, fields))
            for key in INTS:
                row[key] = int(row[key])
            for key in FLOATS:  # ghz and clock are missing from older files
                row[key] = float(row[key]) if row.get(key) else None
            rows.append(row)
    if columns is None:
        raise ValueError("no pffft-bench CSV found")
    return meta, rows


def write_csv(path: Path, meta, rows) -> None:
    """Write meta and rows as read_csv returns them, as a pffft-bench CSV."""
    def cell(v):
        return "" if v is None else str(v)
    lines = ["# pffft-bench 1"] + [f"# {k}: {v}" for k, v in meta.items() if k != "format"]
    columns = list(rows[0]) if rows else ["engine"]
    lines.append(",".join(columns))
    lines += [",".join(cell(r[c]) for c in columns) for r in rows]
    Path(path).write_text("\n".join(lines) + "\n")


def by_case(rows, field):
    """{(engine, prec, xform, n): [value of each trial]}, in the order of the
    rows; trials without the value are left out."""
    cases = {}
    for r in rows:
        if r[field] is not None:
            cases.setdefault((r["engine"], r["prec"], r["xform"], r["n"]), []).append(r[field])
    return cases


def summary(values):
    """(median, min, max)"""
    return statistics.median(values), min(values), max(values)


def fmt(x) -> str:
    return "-" if x is None else f"{x:.2f}" if abs(x) < 100 else f"{x:.1f}"


def is_pow2(n: int) -> bool:
    return n & (n - 1) == 0


def table(title: str, header, rows) -> str:
    """Right-aligned columns under a title line."""
    widths = [max(len(str(r[i])) for r in [header, *rows]) for i in range(len(header))]
    lines = [title] + ["  ".join(str(c).rjust(w) for c, w in zip(r, widths)) for r in [header, *rows]]
    return "\n".join(lines) + "\n"


def engine_tables(rows, field: str) -> str:
    """One table per (precision, transform kind): a row per size, a column per
    engine, the median over trials; non-powers of two marked with '*'."""
    cases = by_case(rows, field)
    engines = list(dict.fromkeys(r["engine"] for r in rows))
    blocks, starred = [], False
    for prec, xform in dict.fromkeys((r["prec"], r["xform"]) for r in rows):
        sizes = sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)})
        cols = [e for e in engines if any((e, prec, xform, n) in cases for n in sizes)]
        if not cols:
            continue
        trials = max(len(v) for k, v in cases.items() if k[1:3] == (prec, xform))
        title = f"{prec} {xform}: {UNITS[field]}, fastest batch"
        title += f", median of {trials} trials" if trials > 1 else ""
        body = []
        for n in sizes:
            vals = [cases.get((e, prec, xform, n)) for e in cols]
            body.append([f"{n}{' ' if is_pow2(n) else '*'}"] + [fmt(v and statistics.median(v)) for v in vals])
            starred |= not is_pow2(n)
        blocks.append(table(title, ["n "] + cols, body))
    return "\n".join(blocks) + ("\n* not a power of two\n" if starred else "")
