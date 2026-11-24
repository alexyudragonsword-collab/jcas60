"""Milestone 0: FMCW against sensing with the 802.11ay Golay training field, on one set of front-end numbers.

Prints the comparison table ``jcas60.table`` builds: row A is fmcw_radar's
default radar as it reports itself; rows B and C are the ay single-carrier
L-CEF used for sensing, monostatic and bistatic, each under both readings of
47 CFR 15.255 (communication limit, radar limit).  Every number carries its
tag; the closed form and the waveform simulation are both run and their
difference is a column (section 3).

    python examples/00_comparison_table.py

The output is kept in ``examples/expected/00_comparison_table.txt``;
``tools/example_drift.py`` compares the two.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jcas60 import ENGINE_BANNER, table  # noqa: E402


def main() -> int:
    text = table.render()
    assert text.startswith(ENGINE_BANNER)
    sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
