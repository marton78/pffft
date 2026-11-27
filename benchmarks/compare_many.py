#!/usr/bin/env python3
"""Time many pffft revisions side by side in one bench_fft process.

    python3 benchmarks/compare_many.py --report FILE
    python3 benchmarks/compare_many.py REV [REV ...] [--tune LABEL=SPEC ...] [--sizes LIST]
                                  [--prec LIST] [--xform LIST] [--launches 3] [--trials 2]
                                  [--pads LIST] [--events LIST]
                                  [--target T] [-o FILE]

A REV is a git revision, optionally LABEL=REV. Each distinct commit is built once
as a shared library; a revision given twice is loaded twice (as a copy), which
measures the noise floor inside the run. The engines are pffft@LABEL (LABEL
defaults to the revision as given). Prints, per case, each revision's speed-up
over the first one (median and range over launches x trials).

--tune LABEL=SPEC makes the setups of pffft@LABEL with the library's default
tuning changed by SPEC (bench_fft --help: e.g. split=8m,prefetch=never). A
label with a --tune loads the library of its commit itself, not a copy, so
that tunings of one commit compare the same code: e.g.
"base=HEAD s8m=HEAD --tune s8m=split=8m".

--cmake LABEL="ARGS" builds the library of LABEL with these CMake arguments
added, e.g. --cmake out="-DPFFFT_DISABLE_PREFETCH=ON"; --cmake "*=ARGS" adds
them to every label (e.g. -DTARGET_C_ARCH=x86-64-v3 for a generic build).
Labels with the same commit and arguments share one build.

--pads LIST compares the revisions at matched code layouts. Every revision is
built once per pad: that many bytes come before the code of the library
(layout_pad.cmake, which works on any revision), so that the same objects lie
at other addresses. The engines are pffft@LABEL+PAD. The report gives, per
case, each revision's speed-up over the first one at the same pad, the range
over the pads, and whether it is slower at every pad, faster at every pad or
depends on the pad (a layout case); and the band: how much the first
revision's own speed moves with the pad. Use multiples of 16 bytes (the
alignment of functions). --pads spread is 0,176,400,720,928,1264,1600,1888:
each pad modulo 128 once (16 apart), and the bits above it (128 to 1024)
varied too, for the predictors and caches that index by them; a pad of 4096
or more would add nothing, as the loader puts a library at an address that
is a multiple of the page. How far the code really moves depends on the
compiler (layout_pad.cmake): the run prints it, and per revision how many
offsets modulo 64 the kernels of the short transforms took over the pads. The libraries' function addresses, which show that each
pad moved the code, are printed and, with -o FILE, saved to FILE.layout.json,
and the libraries to FILE.libs/ (layout.py loops and diff read them).
bench_fft loads at most 16 libraries: revisions x pads <= 16.

--shuffle (always on with --pads) loads the libraries in another order in
each launch, each from a new copy of its file. On an EPYC 9V45 (a VM) the
library file that the build wrote ran 14%, and in another run 40%, slower
than its byte-identical copy at the same pad, in every trial of three
launches, also with the load order shuffled; with new copies per launch no
library was slow throughout. The cause is open.
"""
from __future__ import annotations

import argparse
import hashlib
import shlex
import shutil
import statistics
import sys
from pathlib import Path

import json
import random

import layout
import results
import targets

DEFAULT_SIZES = "1024,2048,4000,4096,8192,15360,16384,32768,48000,65536,100000,131072,262144"
MAX_LIBS = 16  # bench_fft's MAX_LIBS
PAD_CMAKE = targets.REPO / "benchmarks" / "layout_pad.cmake"
SPREAD = "0,176,400,720,928,1264,1600,1888"  # --pads spread: every 16 bytes modulo 128, bits 7-10 varied


def pad_args(pad: int) -> list:
    """CMake arguments that build any revision with pad bytes before its code."""
    return [f"-DCMAKE_PROJECT_PRETTY_FAST_FFT_INCLUDE={PAD_CMAKE.as_posix()}", f"-DPFFFT_LAYOUT_PAD={pad}"]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("revs", nargs="*", metavar="REV")
    ap.add_argument("--report", metavar="FILE", type=Path,
                    help="print the report of a CSV that compare_many.py -o saved, instead of measuring")
    ap.add_argument("--tune", action="append", default=[], metavar="LABEL=SPEC",
                    help="tuning of the setups of pffft@LABEL (repeatable)")
    ap.add_argument("--cmake", action="append", default=[], metavar="LABEL=ARGS",
                    help="CMake arguments for the library of LABEL, or of every label with * (repeatable)")
    ap.add_argument("--pads", metavar="LIST", help="build every revision at each of these pads, "
                    "in bytes, and compare the revisions pad by pad (e.g. 0,16,32,48, or spread)")
    ap.add_argument("--events", metavar="LIST", help="bench_fft --events: hardware events to count too (Linux)")
    ap.add_argument("--shuffle", action="store_true",
                    help="load the libraries in another order, each from a new copy, in each launch (always with --pads)")
    ap.add_argument("--sizes", default=DEFAULT_SIZES)
    ap.add_argument("--prec")
    ap.add_argument("--xform", default="real")
    ap.add_argument("--launches", type=int, default=3)
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("-o", metavar="FILE", type=Path)
    ap.add_argument("--target", default="local", help="local, android[:serial] or ios[:udid] (default: local)")
    args = ap.parse_args()
    if args.report:
        print(saved_report(args.report))
        return
    if not args.revs:
        ap.error("give at least one REV, or --report FILE")

    tuned = {t.partition("=")[0] for t in args.tune}
    cmake = {}
    for c in args.cmake:
        label, sep, extra = c.partition("=")
        if not sep:
            ap.error(f"--cmake: expected LABEL=ARGS, got '{c}'")
        cmake.setdefault(label, []).extend(shlex.split(extra))
    pads = [None]
    if args.pads:
        try:
            pads = sorted({int(p) for p in (SPREAD if args.pads == "spread" else args.pads).split(",")})
        except ValueError:
            ap.error(f"--pads: expected numbers of bytes, got '{args.pads}'")
        if min(pads) < 0:
            ap.error("--pads: a pad cannot be negative")
    if len(args.revs) * len(pads) > MAX_LIBS:
        ap.error(f"{len(args.revs)} revisions x {len(pads)} pads: bench_fft loads at most {MAX_LIBS} libraries")
    target = targets.make_target(args.target)
    harness = target.build_harness(targets.REPO)
    libs, commits, seen, built, bases = {}, {}, {}, {}, []
    for spec in args.revs:
        base, _, rev = spec.rpartition("=")
        base = base or rev
        if base in bases:
            base += "-copy"
        bases.append(base)
        commit, src = targets.export_rev(rev)
        for pad in pads:
            label = base if pad is None else f"{base}+{pad}"
            extra = cmake.get("*", []) + cmake.get(base, []) + ([] if pad is None else pad_args(pad))
            # a build directory per set of arguments; the variant names it
            variant = hashlib.sha1(" ".join(extra).encode()).hexdigest()[:10] if extra else ""
            path = target.build_lib(src, label.replace("/", "_"), extra, variant)
            built[label] = path
            if path in seen.values() and base not in tuned:  # same library: load a copy, as for another revision
                copy = path.with_name(f"copy-{label.replace('/', '_')}-{path.name}")
                shutil.copy2(path, copy)
                path = copy
            seen[label] = path
            libs[label] = path
            commits[label] = commit
    tunes = [f"{b}+{p}={spec}" if p is not None else f"{b}={spec}"
             for b, _, spec in (t.partition("=") for t in args.tune) for p in pads]
    places = layout_report(built, bases, pads, commits) if args.pads else None
    engines = [f"pffft@{label}" for label in libs]
    bench_args = ["--sizes", args.sizes, "--trials", str(args.trials), "--engines", ",".join(engines)]
    for t in tunes:
        bench_args += ["--tune", t]
    for option in ("prec", "xform", "events"):
        if getattr(args, option):
            bench_args += [f"--{option}", getattr(args, option)]
    rows, meta = [], {}
    for launch in range(args.launches):
        order = list(libs.items())
        if args.pads or args.shuffle:  # see --shuffle
            random.Random(launch).shuffle(order)
            order = fresh_copies(order, launch)
        meta, launch_rows = results.read_csv(target.run(harness, bench_args, dict(order)))
        for r in launch_rows:
            r["trial"] += launch * args.trials
        rows += launch_rows
    meta.update(trials=args.launches * args.trials, launches=args.launches,
                revisions=" ".join(f"{l}={c[:10]}" for l, c in commits.items()))
    if args.tune:
        meta.update(tunes=" ".join(args.tune))
    if args.cmake:
        meta.update(cmake=" ; ".join(args.cmake))
    if args.pads:
        meta.update(pads=args.pads)
    if args.pads or args.shuffle:
        meta.update(load_order="shuffled per launch, each library from a new copy")
    if args.o:
        results.write_csv(args.o, meta, rows)
        if places:
            Path(str(args.o) + ".layout.json").write_text(json.dumps(places, indent=1) + "\n")
            save_libs(built, Path(str(args.o) + ".libs"))
    print(f"{meta.get('cpu')}, {meta.get('os')}, counter: {meta.get('counter')}")
    print(report_pads(rows, bases, pads) if args.pads else report(rows, engines))
    if args.pads and rows:
        counted = [c for c in rows[0] if c.startswith("ev:")]
        if any(r["instructions"] is not None for r in rows):
            counted.insert(0, "instructions")
        for field in counted:
            print(report_counts(rows, bases, pads, field))


def saved_report(path: Path) -> str:
    """The report of a saved run: the revisions and pads from its header,
    the layout lines from FILE.layout.json if it is there."""
    meta, rows = results.read_csv(path)
    labels = [x.partition("=")[0] for x in meta.get("revisions", "").split()]
    out = [f"{meta.get('cpu')}, {meta.get('os')}, counter: {meta.get('counter')}"]
    if not meta.get("pads"):
        engines = list(dict.fromkeys(r["engine"] for r in rows))
        return "\n".join(out) + "\n" + report(rows, engines)
    bases = list(dict.fromkeys(l.rpartition("+")[0] for l in labels))
    pads = sorted({int(l.rpartition("+")[2]) for l in labels})
    lay = Path(str(path) + ".layout.json")
    if lay.is_file():
        places = json.loads(lay.read_text())
        for b in bases:
            out.append(f"layout: {b}: {coverage(places, b, pads)}")
            if b != bases[0]:
                out.append(f"layout: {b}: {overlap(places, bases[0], b, pads)}")
    out.append(report_pads(rows, bases, pads))
    counted = [c for c in rows[0] if c.startswith("ev:")] if rows else []
    if any(r["instructions"] is not None for r in rows):
        counted.insert(0, "instructions")
    out += [report_counts(rows, bases, pads, f) for f in counted]
    return "\n".join(out)


def fresh_copies(order, launch):
    """[(label, path)] with each library copied to a new file for this launch
    (see --shuffle)."""
    out = []
    for label, path in order:
        d = path.parent / f"launch-{launch}"
        d.mkdir(exist_ok=True)
        copy = d / path.name
        copy.unlink(missing_ok=True)  # a new file, with new pages in memory
        shutil.copy2(path, copy)
        out.append((label, copy))
    return out


def save_libs(built, where: Path) -> None:
    """Copies of the libraries as LABEL.<suffix>, and with MSVC their maps as
    LABEL.map, for layout.py (loops, diff) on another computer."""
    where.mkdir(parents=True, exist_ok=True)
    for label, path in built.items():
        name = label.replace("/", "_")
        shutil.copy2(path, where / f"{name}{path.suffix}")
        if path.suffix.lower() == ".dll":
            try:
                shutil.copy2(layout.find_map(path), where / f"{name}.map")
            except layout.Error:
                pass


def layout_report(built, bases, pads, commits) -> dict:
    """Read the function addresses of every library, say how far each pad
    moved the code against the smallest pad of the same revision (the groups
    of functions that moved by one amount, the largest first), warn if a pad
    moved nothing or gave the layout of another pad, and return {label:
    {commit, pad, shift, symbols}} for the layout file."""
    places = {}
    for base in bases:
        ref, seen = None, {}
        for pad in pads:
            label = f"{base}+{pad}"
            try:
                syms = layout.text_symbols(built[label])
            except layout.Error as e:
                targets.log(f"layout: {label}: no addresses ({e})")
                places[label] = {"commit": commits[label], "pad": pad, "shift": None, "symbols": {}}
                continue
            if ref is None:
                ref = (pad, syms)
            moved = layout.shift(ref[1], syms)
            groups = ", ".join(f"{d:+} ({k})" for d, k in moved["groups"].items())
            note = ""
            if pad != ref[0] and not any(moved["groups"]):
                note = "; WARNING: the pad moved nothing"
            key = tuple(sorted(moved["groups"].items()))
            if pad != ref[0] and key in seen:
                note += f"; WARNING: the same layout as {base}+{seen[key]}"
            seen.setdefault(key, pad)
            targets.log(f"layout: {label}: functions moved against {base}+{ref[0]} by {groups}{note}")
            places[label] = {"commit": commits[label], "pad": pad, "shift": moved, "symbols": syms}
        targets.log(f"layout: {base}: {coverage(places, base, pads)}")
        if base != bases[0]:
            targets.log(f"layout: {base}: {overlap(places, bases[0], base, pads)}")
    return places


def overlap(places, a, b, pads) -> str:
    """For how many kernels revision b took the same offsets modulo 64 over
    the pads as revision a. Only then do the pads compare the two at the same
    layouts; with other offsets, "slower at every pad" can still be layout
    (code before the kernel that changed size moves it by other amounts)."""
    sa = [places.get(f"{a}+{p}", {}).get("symbols") for p in pads]
    sb = [places.get(f"{b}+{p}", {}).get("symbols") for p in pads]
    if not all(sa) or not all(sb):
        return "no kernel offsets"
    names = [n for n in sa[0] if n.split("#")[0] in KERNELS and all(n in s for s in sa + sb)]
    same = [n for n in names if sorted(s[n] % 64 for s in sa) == sorted(s[n] % 64 for s in sb)]
    differ = [n for n in names if n not in same]
    return (f"kernels at the same offsets modulo 64 over the pads as {a}: {len(same)} of {len(names)}"
            + (f"; others: {' '.join(differ[:12])}{' ...' if len(differ) > 12 else ''}" if differ else ""))


# the kernels of the short transforms: drivers, passes, reorderings (f32; f64 as NAME#2 or pffftd_)
KERNELS = ("pffft_transform_internal", "pffftd_transform_internal", "rfftf1_ps", "rfftb1_ps", "cfftf1_ps",
           "radf2_ps", "radf3_ps", "radf4_ps", "radf5_ps", "radf6_ps", "radf8_ps",
           "radb2_ps", "radb3_ps", "radb4_ps", "radb5_ps", "radb6_ps", "radb8_ps",
           "passf2_ps", "passf3_ps", "passf4_ps", "passf5_ps", "passf6_ps", "passf8_ps",
           "real_finalize_reorder", "real_reorder_preprocess", "cplx_finalize_reorder", "cplx_reorder_preprocess")


def coverage(places, base, pads) -> str:
    """How many offsets modulo 64 the kernels took over the pads: a pad that
    moves a kernel only by multiples of 64 does not vary where its loops lie
    in the 64-byte blocks."""
    syms = [places.get(f"{base}+{p}", {}).get("symbols") for p in pads]
    if not all(syms):
        return "no kernel offsets"
    count = {}
    for name in syms[0]:
        if name.split("#")[0] in KERNELS and all(name in s for s in syms):
            count[name] = len({s[name] % 64 for s in syms})
    if not count:
        return "no kernels found"
    hist = {}
    for k in count.values():
        hist[k] = hist.get(k, 0) + 1
    fixed = [n for n, k in count.items() if k == 1]
    return ("kernels by number of offsets modulo 64 over the pads: "
            + ", ".join(f"{k}: {hist[k]}" for k in sorted(hist, reverse=True))
            + (f"; the same at every pad: {' '.join(fixed)}" if fixed else ""))


def report_pads(rows, bases, pads) -> str:
    """Per (prec, xform, n): the speed-up of each revision over the first one
    at the same pad, in % (median over trials), the range over the pads and
    the verdict; and the band, each revision's speed at each pad over its
    speed at the first pad."""
    out = []
    base = bases[0]
    for prec, xform in dict.fromkeys((r["prec"], r["xform"]) for r in rows):
        out.append(f"\n{prec} {xform}: speed-up over {base} at the same pad, % (median over trials); >0: faster")
        out.append(f"  band: each revision at each pad over itself at +{pads[0]}, max - min;"
                   " slower*/faster*: at every pad, and by more than both bands")
        head = "n".rjust(7) + "  " + "revision".ljust(14) + "".join(f"+{p}".rjust(7) for p in pads)
        out.append(head + "     min     max  verdict" + "       band")
        for n in sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)}):
            same = [r for r in rows if (r["prec"], r["xform"], r["n"]) == (prec, xform, n)]
            field = "cycles" if all(r["cycles"] is not None for r in same) else "ns"
            by = {}
            for r in same:
                by.setdefault(r["engine"], {})[r["trial"]] = r[field]

            def gain(a, b):  # speed-up of engine b over engine a, median over the trials of both
                ta, tb = by.get(f"pffft@{a}", {}), by.get(f"pffft@{b}", {})
                d = [(ta[t] / tb[t] - 1) * 100 for t in ta if t in tb]
                return statistics.median(d) if d else None

            base_band = [gain(f"{base}+{pads[0]}", f"{base}+{p}") for p in pads]
            for i, rev in enumerate(bases):
                band = [gain(f"{rev}+{pads[0]}", f"{rev}+{p}") for p in pads]
                band_txt = rng(band)
                if i == 0:
                    vals, verdict, tail = band, "(band)", ""
                else:
                    vals = [gain(f"{base}+{p}", f"{rev}+{p}") for p in pads]
                    verdict = judge(vals, max(spread(base_band), spread(band)))
                    tail = band_txt
                cells = "".join(("-" if v is None else f"{v:+.1f}").rjust(7) for v in vals)
                known = [v for v in vals if v is not None]
                span = (f"{min(known):+.2f}".rjust(8) + f"{max(known):+.2f}".rjust(8)) if known else " " * 16
                out.append(str(n).rjust(7) + "  " + rev[:14].ljust(14) + cells + span + "  " + verdict.ljust(10)
                           + tail.rjust(11) + ("" if field == "cycles" else "  (ns)"))
    return "\n".join(out)


def report_counts(rows, bases, pads, field) -> str:
    """Per (prec, xform, n) and revision: field (instructions or an event)
    per transform at each pad, the median over trials."""
    out = []
    for prec, xform in dict.fromkeys((r["prec"], r["xform"]) for r in rows):
        out.append(f"\n{prec} {xform}: {field.removeprefix('ev:')} per transform (median over trials)")
        out.append("n".rjust(7) + "  " + "revision".ljust(14) + "".join(f"+{p}".rjust(10) for p in pads))
        for n in sorted({r["n"] for r in rows if (r["prec"], r["xform"]) == (prec, xform)}):
            for rev in bases:
                cells = ""
                for p in pads:
                    v = [r[field] for r in rows if (r["prec"], r["xform"], r["n"], r["engine"]) ==
                         (prec, xform, n, f"pffft@{rev}+{p}") and r[field] is not None]
                    cells += (f"{statistics.median(v):.1f}" if v else "-").rjust(10)
                out.append(str(n).rjust(7) + "  " + rev[:14].ljust(14) + cells)
    return "\n".join(out)


def spread(vals) -> float:
    known = [v for v in vals if v is not None]
    return max(known) - min(known) if known else 0.0


def rng(vals) -> str:
    known = [v for v in vals if v is not None]
    return f"{spread(known):.2f}" if known else "-"


def judge(vals, band=None) -> str:
    """slower: at every pad; faster: at every pad; layout: the sign depends on
    the pad. A * when the effect at every pad is larger than band, the wider
    of the two revisions' bands: the pads sample the layouts of each revision
    (a changed revision moves its kernels by other amounts than the pads do,
    so that no pad puts both at the same offsets), and the band says how much
    layout moves each one's speed; an effect beyond it is not layout."""
    known = [v for v in vals if v is not None]
    if not known:
        return "-"
    if max(known) < 0:
        return "slower" + ("*" if band is not None and max(known) < -band else "")
    if min(known) > 0:
        return "faster" + ("*" if band is not None and min(known) > band else "")
    return "layout"


if __name__ == "__main__":
    try:
        main()
    except targets.Error as e:
        sys.exit(f"compare_many.py: {e}")
