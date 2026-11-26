#!/usr/bin/env python3
"""Run .github/workflows/bench.yml until every entry of a table has run on its CPU.

    python3 benchmarks/ci_dispatch.py TABLE.json OUT --ref BRANCH [--repo OWNER/NAME] [-k 3] [--waves 30]

GitHub's runner labels land on several CPUs: ubuntu-26.04 on Intel Xeon and AMD
EPYC machines of several generations, windows-2025 likewise. A run gets, as
the workflow's input by_cpu, what every CPU still missing on its label should
do, and does the part of the CPU it lands on: a run on a wanted CPU is never
wasted, one on another CPU ends after a few seconds. Several runs per label go
out at once (one per missing CPU, at most -k). Give each CPU everything it
should measure in one entry: some CPUs come up once in 30 runs.

TABLE.json:

    {"entries": {
      "xeon6975": {"runner": "ubuntu-26.04", "cpu": "6975P", "cycles": true,
                   "job": {"sweep": true, "fftw": true,
                           "compares": [{"revs": "master HEAD", "name": "master"}],
                           "compare_args": "--sizes pow2", "many": "a=master b=HEAD", "many_args": ""}},
      "epyc7763": {"runner": "windows-2025", "cpu": "7763",
                   "job": {"compares": [{"revs": "master HEAD", "name": "master"}]}}}}

"cpu" is a regex for the model name, "cycles" asks for a cycle counter. Runs
are repeated, in waves, until every entry has its files: <runner>.csv for the
sweep, <runner>-cmp-<name>.csv per compare, <runner>-many.csv. OUT/status.json
records the run and files of each entry, OUT/<entry>/run-<id>/ holds the
artifacts. Started again with the same OUT, it goes on where it stopped.
Needs gh, logged in, and bench.yml on the repository's default branch.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import time


def sh(*args: str) -> str:
    return subprocess.run(args, capture_output=True, text=True, check=True).stdout


def run_ids(repo: str) -> list[int]:
    return [r["databaseId"] for r in json.loads(sh("gh", "run", "list", "-R", repo, "--workflow", "bench.yml",
                                                     "-L", "40", "--json", "databaseId"))]


def expected(runner: str, job: dict) -> list[str]:
    """The files that an entry's job produces."""
    files = [f"{runner}.csv"] if job.get("sweep") else []
    files += [f"{runner}-cmp-{c['name']}.csv" for c in job.get("compares", [])]
    return files + ([f"{runner}-many.csv"] if job.get("many") else [])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("table", help="the table of entries (JSON, see above)")
    ap.add_argument("out", help="directory for status.json and the artifacts")
    ap.add_argument("--ref", required=True, help="the branch whose bench.yml and code run")
    ap.add_argument("--repo", help="OWNER/NAME (default: gh's repository of this checkout)")
    ap.add_argument("-k", type=int, default=3, help="runs per label at once, at most (default: 3)")
    ap.add_argument("--waves", type=int, default=30, help="give up after this many waves (default: 30)")
    args = ap.parse_args()
    repo = args.repo or sh("gh", "repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner").strip()
    entries = json.load(open(args.table))["entries"]
    os.makedirs(args.out, exist_ok=True)
    status_path = os.path.join(args.out, "status.json")
    status = json.load(open(status_path)) if os.path.exists(status_path) else {}

    for wave in range(1, args.waves + 1):
        pending = {n: e for n, e in entries.items() if not status.get(n, {}).get("ok")}
        if not pending:
            break
        print(f"wave {wave}: {sorted(pending)}", flush=True)
        runs = []
        for runner in sorted({e["runner"] for e in pending.values()}):
            mine = [e for e in pending.values() if e["runner"] == runner]
            by_cpu = json.dumps({e["cpu"]: e["job"] for e in mine})
            for _ in range(min(len(mine), args.k)):
                before = set(run_ids(repo))
                sh("gh", "workflow", "run", "bench.yml", "-R", repo, "--ref", args.ref, "-f", f"runners={runner}",
                   "-f", "sweep=false", "-f", "fftw=false", "-f", f"by_cpu={by_cpu}")
                rid = None
                for _ in range(40):  # the run appears a few seconds after the dispatch
                    time.sleep(3)
                    new = [i for i in run_ids(repo) if i not in before]
                    if new:
                        rid = max(new)
                        break
                runs.append((runner, rid))
                print(f"  {runner}: run {rid}", flush=True)
        for _, rid in runs:
            if rid:
                subprocess.run(["gh", "run", "watch", str(rid), "-R", repo, "--interval", "60"], capture_output=True)
        for runner, rid in runs:
            if not rid:
                continue
            tmp = os.path.join(args.out, "_runs", f"run-{rid}")
            shutil.rmtree(tmp, ignore_errors=True)
            os.makedirs(tmp)
            subprocess.run(["gh", "run", "download", str(rid), "-R", repo, "--dir", tmp], capture_output=True)
            system = glob.glob(f"{tmp}/*/*-system.txt")
            m = re.search(r"(model name|brand_string|^Name)\s*:?\s*(.*)", open(system[0]).read(), re.I | re.M) \
                if system else None
            cpu = m.group(2).strip() if m else ""
            entry_file = glob.glob(f"{tmp}/*/*-entry.txt")
            key = open(entry_file[0]).read().strip() if entry_file else None
            name = next((n for n, e in entries.items() if e["runner"] == runner and e["cpu"] == key), None)
            if name is None or status.get(name, {}).get("ok"):
                print(f"  {runner}: run {rid} on {cpu!r}: " + ("nothing to do" if name is None else f"{name} again"),
                      flush=True)
                shutil.rmtree(tmp, ignore_errors=True)
                continue
            entry = entries[name]
            files = [glob.glob(f"{tmp}/*/{f}") for f in expected(runner, entry["job"])]
            ok = all(files) and (not entry.get("cycles") or
                                 all("perf_event" in open(f[0]).read(3000) for f in files))
            dst = os.path.join(args.out, name, f"run-{rid}")
            shutil.rmtree(dst, ignore_errors=True)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(tmp, dst)
            status[name] = {"run": rid, "cpu": cpu, "ok": bool(ok), "files": [f[0].replace(tmp, dst) for f in files if f]}
            print(f"  {name}: run {rid} on {cpu!r}: ok={bool(ok)}", flush=True)
            with open(status_path, "w") as f:
                json.dump(status, f, indent=1)
    print("done:", json.dumps({n: status.get(n, {}).get("ok", False) for n in entries}), flush=True)


if __name__ == "__main__":
    main()
