#!/usr/bin/env python3
"""Example 03: a receive array (milestone 1, second step).

Example 02 put N power amplifiers on the transmit side.  Here N receive
chains sit behind N elements, each with its own ADC, all on one LO, and the
beams are formed afterwards:

1. Three ways to spend N elements -- receive only, transmit swept, both --
   under both regulatory readings and at the same air time.
2. The simulated chains combined: the SNR, and how well one frame tells the
   target's direction (against the Cramer-Rao bound).
3. The transmit leakage in every chain.

    python examples/03_receive_array.py

About 55 s.  ``tools/example_drift.py`` compares the output with
``examples/expected/03_receive_array.txt``.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jcas60 import rxarray  # noqa: E402

if __name__ == "__main__":
    sys.stdout.write(rxarray.render())
