#!/usr/bin/env python3
"""fmcw_radar's own budget report of its default config, kept as data (row A of the table).

The table's row A is not computed here: its numbers are what fmcw_radar's
``tools/budget_report.py configs/radar_60g_sensing.yaml --json`` writes at the
pinned commit.  That report is not a file of fmcw_radar's repository, so it
cannot be vendored like one; this tool regenerates it from the pinned commit
(exported with ``git archive`` into a temporary directory, so the sibling's
working tree is not touched) and either records it or compares it with the
recorded copy:

    python tools/upstream_report.py record      # write src/jcas60/data/upstream/fmcw_radar_budget.json
    python tools/upstream_report.py check       # regenerate and compare; exit 1 on a difference
    python tools/upstream_report.py check --strict   # a missing sibling checkout is a failure, not a skip

Numbers are compared to a relative 1e-9 (the report is floating point from
numpy / scipy; another library version may move the last bits), everything
else exactly.  The sibling is found like ``tools/vendor_check.py`` finds it.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tarfile
import tempfile
from io import BytesIO
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RECORD = REPO / "src" / "jcas60" / "data" / "upstream" / "fmcw_radar_budget.json"
UPSTREAM = {"repo": "fmcw_radar", "sha": "c303566180c2ef22bd2bf6ffaa1f08e59174e110",
            "config": "configs/radar_60g_sensing.yaml",
            "command": "python tools/budget_report.py configs/radar_60g_sensing.yaml --json <out>"}
ROOTS = ("/home/user", "/home/user/alexyudragonsword-collab")
RTOL = 1e-9


def sibling(explicit: str | None = None) -> Path | None:
    cands = [explicit or os.environ.get("JCAS60_SIBLING_fmcw_radar")] if (explicit or os.environ.get("JCAS60_SIBLING_fmcw_radar")) \
        else [os.path.join(r, UPSTREAM["repo"]) for r in ROOTS]
    return next((Path(c) for c in cands if c and (Path(c) / ".git").exists()), None)


def regenerate(root: Path) -> dict:
    sha = UPSTREAM["sha"]
    if subprocess.run(["git", "-C", str(root), "cat-file", "-e", f"{sha}^{{commit}}"], capture_output=True).returncode:
        subprocess.run(["git", "-C", str(root), "fetch", "--quiet", "--depth=1", "origin", sha], capture_output=True)
    tar = subprocess.run(["git", "-C", str(root), "archive", sha], capture_output=True, check=True).stdout
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(fileobj=BytesIO(tar)) as t:
            t.extractall(tmp, filter="data")
        out = Path(tmp) / "report.json"
        subprocess.run([sys.executable, "tools/budget_report.py", UPSTREAM["config"], "-o", str(Path(tmp) / "report.md"),
                        "--json", str(out)], cwd=tmp, check=True, capture_output=True)
        return json.loads(out.read_text(encoding="utf-8"))


def differences(a, b, path: str = "") -> list[str]:
    if isinstance(a, dict) and isinstance(b, dict):
        out = [f"{path}/{k}: only in one" for k in sorted(set(a) ^ set(b))]
        for k in sorted(set(a) & set(b)):
            out += differences(a[k], b[k], f"{path}/{k}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: lengths {len(a)} and {len(b)}"]
        return [d for i, (x, y) in enumerate(zip(a, b)) for d in differences(x, y, f"{path}[{i}]")]
    if isinstance(a, float) or isinstance(b, float):
        if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
            return [] if math.isclose(a, b, rel_tol=RTOL, abs_tol=1e-300) or (math.isnan(a) and math.isnan(b)) \
                else [f"{path}: {a!r} != {b!r}"]
    return [] if a == b else [f"{path}: {a!r} != {b!r}"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("record", "check"))
    ap.add_argument("--sibling", default=None, help="path of the fmcw_radar checkout")
    ap.add_argument("--strict", action="store_true", help="fail when the sibling checkout is absent")
    args = ap.parse_args(argv)
    root = sibling(args.sibling)
    if root is None:
        print("fmcw_radar checkout not found: the recorded report cannot be verified here")
        return 1 if (args.strict or args.mode == "record") else 0
    report = regenerate(root)
    if args.mode == "record":
        RECORD.parent.mkdir(parents=True, exist_ok=True)
        RECORD.write_text(json.dumps({"_provenance": UPSTREAM, "report": report}, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {RECORD.relative_to(REPO)}")
        return 0
    have = json.loads(RECORD.read_text(encoding="utf-8"))
    diffs = ([] if have["_provenance"] == UPSTREAM else ["_provenance differs from this tool's pin"]) \
        + differences(have["report"], json.loads(json.dumps(report)))
    for d in diffs:
        print(d)
    print(f"fmcw_radar@{UPSTREAM['sha'][:7]} budget report: " + ("identical to the recorded copy" if not diffs else f"{len(diffs)} difference(s)"))
    return 1 if diffs else 0


if __name__ == "__main__":
    raise SystemExit(main())
