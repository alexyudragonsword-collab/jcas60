#!/usr/bin/env python3
"""Do the examples still print what ``examples/expected`` holds?

The docs and the commit messages quote what the examples print; a change that
moves a number without failing a test would leave them stale.  Each example's
full output is kept in ``examples/expected/<name>.txt``:

    python tools/example_drift.py compare     # run every example, list each line that differs; exit 1 if any does
    python tools/example_drift.py accept      # after the docs are fixed: overwrite the expected outputs

(The siblings' tool also looks each moved number up in the docs; with one
example this library compares the lines and leaves the lookup to ``grep``.)
"""
from __future__ import annotations

import argparse
import difflib
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
EXAMPLES = REPO / "examples"
EXPECTED = EXAMPLES / "expected"


def run(example: Path) -> str:
    env = {**os.environ, "PYTHONPATH": str(REPO / "src")}
    return subprocess.run([sys.executable, str(example)], capture_output=True, text=True, check=True, env=env).stdout


def examples() -> list[Path]:
    return sorted(p for p in EXAMPLES.glob("[0-9][0-9]_*.py"))


def compare(out=sys.stdout) -> int:
    bad = 0
    for ex in examples():
        want = EXPECTED / f"{ex.stem}.txt"
        got = run(ex)
        old = want.read_text(encoding="utf-8") if want.exists() else ""
        if got != old:
            bad += 1
            diff = list(difflib.unified_diff(old.splitlines(), got.splitlines(), f"expected/{want.name}", "now", lineterm="", n=0))
            print(f"{ex.name}: {sum(1 for d in diff if d[:1] in '+-' and d[:3] not in ('+++', '---'))} line(s) differ", file=out)
            for d in diff[:80]:
                print("  " + d[:240], file=out)
        else:
            print(f"{ex.name}: as expected", file=out)
    return 1 if bad else 0


def accept() -> int:
    for ex in examples():
        (EXPECTED / f"{ex.stem}.txt").write_text(run(ex), encoding="utf-8")
        print(f"wrote examples/expected/{ex.stem}.txt")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=("compare", "accept"))
    return compare() if ap.parse_args(argv).mode == "compare" else accept()


if __name__ == "__main__":
    raise SystemExit(main())
