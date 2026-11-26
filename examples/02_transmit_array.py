#!/usr/bin/env python3
"""Example 02: the comparison table's question asked again with a transmit array (milestone 1, first step).

One antenna radiates 22 dBm of EIRP against a communication limit of 40 dBm.
N elements, each behind its own PA and steered in phase:

1. EIRP, back-off and the rows' SNR under both regulatory readings.
2. The beam's width, how many beams cover the field of view, and what a
   sweep of them costs in air time.
3. Phase steering over a waveform 2.2 GHz wide: the squint.
4. The link on both roads at the array that reaches the communication limit.
5. The monostatic row: N power amplifiers beside the receiver.

    python examples/02_transmit_array.py

About 45 s.  ``tools/example_drift.py`` compares the output with
``examples/expected/02_transmit_array.txt``.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jcas60 import arrays  # noqa: E402

if __name__ == "__main__":
    sys.stdout.write(arrays.render())
