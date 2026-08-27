#!/usr/bin/env python3
"""Profile-guided optimization harness for PFFFT. See docs/superpowers/plans/."""
from __future__ import annotations

import argparse
import json
import random
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from samples import read_samples, group_key
from stats import SIG_FRAC_FASTER, aggregate, compare
from targets import (BASE_CMAKE_FLAGS, REPO_ROOT, build, ensure_worktrees,
                     make_target, sync_worktree)

PERF_DIR = REPO_ROOT / ".perf"
SAMPLES_DIR = PERF_DIR / "samples"
CHAIN_FILE = REPO_ROOT / "bench_chain.json"

SHORT_POW2 = [32, 64, 128, 256, 512, 1024, 2048]
SHORT_NPOW2 = [96, 192, 480]


def balance_order(n: int, rng: random.Random) -> list[str]:
    """Balanced shuffled A/B sequence: half A, half B, Fisher-Yates shuffle."""
    seq = ["A"] * (n // 2) + ["B"] * (n // 2)
    rng.shuffle(seq)
    return seq


def pick_sizes(mode: str) -> list[int]:
    pow2 = [32 << e for e in range(17)]           # 32 .. 2^21 (C binary's max)
    npow2 = [96, 160, 192, 320, 480, 640, 768, 960, 1152, 1280, 1440, 1600,
             1920, 2400, 2560, 3200, 4000, 4800, 7200, 9216, 15360]
    if mode == "short":  return sorted(set(SHORT_POW2 + SHORT_NPOW2))
    if mode == "pow2":   return pow2
    if mode == "nonpow2": return npow2
    if mode == "all":    return sorted(set(pow2) | set(npow2))
    raise ValueError(f"--sizes: unknown mode {mode!r}")


def collect_group(path: Path) -> tuple[dict[tuple, list[float]], dict[str, str]]:
    """-> ({(algo,prec,xform,size): [cost,...]}, provenance)

    `cost` = -mflops, NOT sample_ms. sample_ms is a fixed calibrated
    timing window (~150ms by design; the binary picks n_iter to fill it,
    see bench_pffft.c) and is therefore roughly CONSTANT regardless of
    code speed -- the entire throughput signal lives in n_iter, i.e. the
    derived MFLOPS. Comparing sample_ms directly (the original bug here)
    compares two near-identical calibration windows and can never detect
    a real performance difference; every verdict comes back "neutral"
    regardless of the actual code. Negating mflops keeps the existing
    "lower cost = faster" convention (Comparison/aggregate/HL shift) working
    unmodified: higher throughput => lower (more negative) cost.
    """
    prov, rows = read_samples(path)
    out: dict[tuple, list[float]] = {}
    for r in rows:
        out.setdefault(group_key(r), []).append(-r.mflops)
    return out, prov


def _dirty(ref: str) -> bool:
    """True iff the working tree differs from ref (diff-index --quiet, inverted)."""
    return subprocess.run(["git", "diff-index", "--quiet", ref, "--"],
                          cwd=REPO_ROOT).returncode != 0


def sample_file(label: str, target_name: str, tree: str) -> Path:
    return SAMPLES_DIR / f"{label}-{target_name.replace('/', '_')}-{tree[:8]}.csv"


def load_chain() -> dict | None:
    """Load bench_chain.json, or None when absent/corrupt."""
    try:
        return json.loads(CHAIN_FILE.read_text())
    except (OSError, ValueError):
        return None


def save_chain(chain: dict) -> None:
    CHAIN_FILE.write_text(json.dumps(chain, indent=2) + "\n")


def _chain_head(chain: dict | None) -> dict | None:
    if not chain:
        return None
    steps = chain.get("steps") or []
    if steps:
        return steps[-1]
    return chain.get("base") or None


def chain_head_tree() -> str | None:
    """Tree hash of the current chain head (last step, else base), if any."""
    head = _chain_head(load_chain())
    return head.get("tree") if head else None


def chain_head_label() -> str | None:
    """Label of the current chain head (last step, else base), if any."""
    head = _chain_head(load_chain())
    return head.get("label") if head else None


def counts_at_size(path: Path, prec: str, size: int) -> dict[tuple[str, str], int]:
    """Rows per (algo, xform) already recorded for one (prec, size) in an arm."""
    counts: dict[tuple[str, str], int] = {}
    if not path.exists():
        return counts
    _, rows = read_samples(path)
    for r in rows:
        if r.prec == prec and r.size == size:
            counts[(r.algo, r.xform)] = counts.get((r.algo, r.xform), 0) + 1
    return counts


def _missing_invocations(counts: dict[tuple[str, str], int],
                         need_rows: int, runs: int) -> int:
    """Invocations still required: worst deficit across (algo, xform) groups."""
    if not counts:
        return -(-need_rows // runs)
    worst_deficit = max(need_rows - c for c in counts.values())
    return -(-worst_deficit // runs)


def _series(groups: dict[tuple, list[float]],
            prec: str) -> dict[tuple[str, str], dict[int, list[float]]]:
    """{(algo, xform): {size: [cost]}} for one precision (see collect_group)."""
    out: dict[tuple[str, str], dict[int, list[float]]] = {}
    for (algo, p, xform, size), cost in groups.items():
        if p == prec:
            out.setdefault((algo, xform), {})[size] = cost
    return out


def cmd_ab(args):
    rng = random.Random(time.time_ns())
    targets = [make_target(t) for t in (args.target or ["local"])]
    base_wt, var_wt = ensure_worktrees()
    base_ref = args.base
    if base_ref is None:
        head_tree = chain_head_tree()
        if head_tree:
            # No explicit --base: continue the accepted-optimization chain.
            base_ref = head_tree
            own = subprocess.run(["git", "rev-parse", "HEAD^{tree}"],
                                 cwd=REPO_ROOT, capture_output=True, text=True)
            own_tree = own.stdout.strip() if own.returncode == 0 else None
            if own_tree != head_tree:
                print(f"base: {chain_head_label()} @ {head_tree[:8]} "
                      f"(chain head)", flush=True)
    base_ref = base_ref or "HEAD"
    var_ref = args.variant or "HEAD"
    label = args.label
    base_label = args.base_label or chain_head_label() or f"{label}-base"
    base_tree = sync_worktree(base_wt, base_ref)
    var_tree = sync_worktree(var_wt, var_ref)
    var_dirty = _dirty(var_ref)
    base_bin = build(base_wt)
    var_bin = build(var_wt)
    sizes = [n for n in pick_sizes(args.sizes) if n <= args.max_len]

    compiler = subprocess.run(["cc", "--version"], capture_output=True,
                              text=True).stdout.splitlines()[0].strip()
    common_meta = [
        f"compiler={compiler}",
        f"flags={' '.join(BASE_CMAKE_FLAGS)}",
        f"sizes={args.sizes}",
        f"datetime={time.strftime('%Y-%m-%dT%H:%MZ', time.gmtime())}",
    ]

    arms = {
        "A": (base_wt, base_bin, base_tree, base_label),
        "B": (var_wt, var_bin, var_tree, label),
    }
    precs = ["flt", "dbl"] if args.prec == "both" else [args.prec]

    for tgt in targets:
        paths = {w: sample_file(lab, tgt.name, tree)
                 for w, (_, _, tree, lab) in arms.items()}
        for prec in precs:
            for size in sizes:
                # Resume support: only issue invocations whose rows are missing.
                missing = {w: _missing_invocations(
                               counts_at_size(paths[w], prec, size),
                               args.invocations * args.runs, args.runs)
                           for w in ("A", "B")}
                order = ["A"] * missing["A"] + ["B"] * missing["B"]
                rng.shuffle(order)
                for who in order:
                    _, bindir, tree, lab = arms[who]
                    exe = tgt.binary_path(bindir, prec)
                    meta = [
                        f"label={lab}", f"git_tree={tree}",
                        f"dirty={'true' if who == 'B' and var_dirty else 'false'}",
                        *tgt.meta_list(), *common_meta,
                    ]
                    # NB: bench_pffft matches flags with exact strcmp, so
                    # values must be separate argv tokens, never --opt=value.
                    cmd = [
                        str(exe), "--size", str(size), "--runs", str(args.runs),
                        "--samples", str(paths[who]),
                        *([] if size & (size - 1) == 0 else ["--non-pow2"]),
                        *[x for m in meta for x in ("--meta", m)],
                    ]
                    print(f"[{who}] {tgt.name} {prec} n={size} "
                          f"-> {paths[who].name} ({missing[who]} left)",
                          flush=True)
                    for _ in tgt.run(cmd):
                        pass

    # ---- analyze ----
    files = {"base": [], "var": []}
    comparisons: list[dict] = []
    group_aggs: dict[str, dict] = {}

    print("\n== A/B verdict (negative shift = variant faster) ==")
    print(f"{'algo':<10} {'xform':<5} {'shift%':>9}  {'verdict':<8} p_min")
    for tgt in targets:
        fb = sample_file(base_label, tgt.name, base_tree)
        fv = sample_file(label, tgt.name, var_tree)
        files["base"].append(str(fb))
        files["var"].append(str(fv))
        for p in precs:
            if not fb.exists() or not fv.exists():
                print(f"({tgt.name} {p}: missing sample file(s), skipped)")
                continue
            sb = _series(collect_group(fb)[0], p)
            sv = _series(collect_group(fv)[0], p)
            for gx in sorted(set(sb) & set(sv)):
                comps, done_sizes = [], []
                for size in sorted(set(sb[gx]) & set(sv[gx])):
                    b, v = sb[gx][size], sv[gx][size]
                    n = min(len(b), len(v))     # align reps by index
                    if n:
                        comps.append(compare(b[:n], v[:n]))
                        done_sizes.append(size)
                if not comps:
                    continue
                shifts = [c.shift_pct for c in comps]
                agg = aggregate(comps)
                p_min = min(c.p_value for c in comps)
                print(f"{gx[0]:<10} {gx[1]:<5} {float(np.median(shifts)):>+8.2f}%  "
                      f"{agg['suggest']:<8} {p_min:.4f}")
                for size, c in zip(done_sizes, comps):
                    comparisons.append({"target": tgt.name, "prec": p,
                                        "algo": gx[0], "xform": gx[1],
                                        "size": size, **asdict(c)})
                group_aggs[f"{tgt.name}/{p}/{gx[0]}/{gx[1]}"] = {
                    **agg, "median_shift_pct": float(np.median(shifts))}

    # Top-level verdict at GROUP grain (Design Contract): per-(algo, xform)
    # group suggests vote; any slower => slower, else >=70% faster => faster.
    group_suggests = [a["suggest"] for a in group_aggs.values()]
    if any(s == "slower" for s in group_suggests):
        suggest = "slower"
    elif (group_suggests
          and sum(s == "faster" for s in group_suggests)
          >= SIG_FRAC_FASTER * len(group_suggests)):
        suggest = "faster"
    else:
        suggest = "neutral"
    final = {"n_groups": len(group_suggests),
             "n_sig_faster": sum(s == "faster" for s in group_suggests),
             "n_sig_slower": sum(s == "slower" for s in group_suggests),
             "suggest": suggest}
    result = {
        "label": label, "base_label": base_label,
        "base_tree": base_tree, "var_tree": var_tree,
        "files": files, "targets": [t.name for t in targets],
        "suggest": suggest,
        "comparisons": comparisons, "groups": group_aggs,
        "aggregate": final,
    }
    out = PERF_DIR / "last_ab.json"
    out.write_text(json.dumps(result, indent=2))
    print(f"\nsuggest: {suggest} "
          f"({final['n_sig_faster']} groups faster, "
          f"{final['n_sig_slower']} groups slower of {final['n_groups']})")
    print(f"wrote {out}")


def _step_verdicts(last: dict) -> list[dict]:
    """Per-group verdict entries for a chain step, from last_ab.json."""
    p_min: dict[str, float] = {}
    for c in last.get("comparisons", []):
        k = f"{c['target']}/{c['prec']}/{c['algo']}/{c['xform']}"
        p_min[k] = min(p_min.get(k, 1.0), c["p_value"])
    out = []
    for key, agg in sorted(last.get("groups", {}).items()):
        target, prec, algo, xform = key.split("/")
        out.append({"target": target, "prec": prec, "algo": algo,
                    "xform": xform,
                    "sizes_total": agg["n_groups"],
                    "sizes_sig_faster": agg["n_sig_faster"],
                    "sizes_sig_slower": agg["n_sig_slower"],
                    "shift_pct": agg.get("median_shift_pct"),
                    "p_min": p_min.get(key)})
    return out


def cmd_accept(args):
    """Fold the last `ab` result into bench_chain.json (human-gated)."""
    last_file = PERF_DIR / "last_ab.json"
    try:
        last = json.loads(last_file.read_text())
    except (OSError, ValueError) as e:
        raise SystemExit(f"accept: cannot read {last_file} ({e}); "
                         "run `perf.py ab` first")
    label, tree = last["label"], last["var_tree"]
    suggest = last.get("suggest")
    groups = last.get("groups", {})
    print(f"label: {label} @ {tree[:8]}   suggest: {suggest}")
    print(f"{'group':<40} {'shift%':>9} {'faster':>7} {'slower':>7}")
    for key, agg in sorted(groups.items()):
        print(f"{key:<40} {agg.get('median_shift_pct', 0.0):>+8.2f}% "
              f"{agg['n_sig_faster']:>7} {agg['n_sig_slower']:>7}")
    if suggest == "slower" and not args.force:
        raise SystemExit("accept: refusing: A/B verdict is 'slower' "
                         "(pass --force to override)")
    if not args.yes:
        reply = input(f"accept '{label}' @ {tree[:8]} into "
                      f"{CHAIN_FILE.name}? [y/N] ")
        if reply.strip().lower() not in ("y", "yes"):
            raise SystemExit("aborted")
    chain = load_chain()
    if not chain or not isinstance(chain.get("base"), dict):
        chain = {"base": {"label": last.get("base_label"),
                          "tree": last.get("base_tree")},
                 "steps": []}
    def _chain_path(p):
        """Store PERF_DIR-relative when possible (make_charts resolves both)."""
        try:
            return str(Path(p).relative_to(PERF_DIR))
        except ValueError:
            return str(p)
    sample_files = sorted({_chain_path(p)
                           for p in last.get("files", {}).get("base", [])
                           + last.get("files", {}).get("var", [])})
    chain.setdefault("steps", []).append({
        "label": label, "tree": tree, "verdicts": _step_verdicts(last),
        "accepted": True, "sample_files": sample_files})
    save_chain(chain)
    print(f"accepted '{label}' @ {tree[:8]} -> {CHAIN_FILE}")
    print(f"chain head is now {label} @ {tree[:8]}")


def cmd_status(_args=None):
    """Pretty-print the chain base, accepted steps, and current head."""
    chain = load_chain()
    if not chain:
        print(f"no chain manifest at {CHAIN_FILE}")
        return
    base = chain.get("base") or {}
    print(f"base: {base.get('label')} @ {(base.get('tree') or '')[:8]}")
    steps = chain.get("steps") or []
    if not steps:
        print("steps: (none)")
    for s in steps:
        vs = s.get("verdicts") or []
        shifts = [v["shift_pct"] for v in vs if v.get("shift_pct") is not None]
        med = f"{np.median(shifts):+.2f}%" if shifts else "n/a"
        nf = sum(v.get("sizes_sig_faster", 0) for v in vs)
        ns = sum(v.get("sizes_sig_slower", 0) for v in vs)
        print(f"  step {s.get('label')} @ {(s.get('tree') or '')[:8]}  "
              f"median shift {med}, {nf} sizes sig-faster / {ns} sig-slower "
              f"across {len(vs)} verdict(s), "
              f"{len(s.get('sample_files') or [])} sample file(s)")
    head_t, head_l = chain_head_tree(), chain_head_label()
    print(f"head: {head_l} @ {head_t[:8] if head_t else '?'}")


def cmd_evolution(_args=None):
    """Render evolution charts via make_charts.py --evolution.

    On failure (usually an empty chain or no sample data yet), fall back
    to printing a chain summary.
    """
    script = Path(__file__).resolve().parent / "make_charts.py"
    r = subprocess.run([sys.executable, str(script), "--evolution",
                        CHAIN_FILE.name], cwd=REPO_ROOT)
    if r.returncode != 0:
        print("\nevolution: chart generation failed (usually an empty chain "
              "or no sample data yet); printing chain summary instead\n")
        cmd_status()


def main():
    ap = argparse.ArgumentParser(prog="perf.py")
    sub = ap.add_subparsers(dest="cmd", required=True)
    ab = sub.add_parser("ab")
    ab.add_argument("--label", required=True)
    ab.add_argument("--base-label", default=None,
                    help="default: chain head label, else '<label>-base'")
    ab.add_argument("--variant", default="HEAD")
    ab.add_argument("--base", default=None)
    ab.add_argument("--runs", type=int, default=4)
    ab.add_argument("--invocations", type=int, default=5)
    ab.add_argument("--prec", choices=["flt", "dbl", "both"], default="flt")
    ab.add_argument("--sizes", choices=["short", "pow2", "nonpow2", "all"],
                    default="short")
    ab.add_argument("--max-len", type=int, default=1 << 30)
    ab.add_argument("--target", action="append", default=None)
    ab.set_defaults(fn=cmd_ab)
    acc = sub.add_parser(
        "accept", help="fold the last `ab` result into bench_chain.json")
    acc.add_argument("--yes", action="store_true",
                     help="skip interactive confirmation")
    acc.add_argument("--force", action="store_true",
                     help="accept even when the verdict suggests 'slower'")
    acc.set_defaults(fn=cmd_accept)
    st = sub.add_parser("status",
                        help="print chain base, accepted steps, and head")
    st.set_defaults(fn=cmd_status)
    evo = sub.add_parser(
        "evolution",
        help="render evolution charts via make_charts.py --evolution; falls "
             "back to a chain summary when chart generation fails")
    evo.set_defaults(fn=cmd_evolution)
    args = ap.parse_args()
    SAMPLES_DIR.mkdir(parents=True, exist_ok=True)
    args.fn(args)


if __name__ == "__main__":
    main()
