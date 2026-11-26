#!/usr/bin/env python3
"""Example 01: eight questions the comparison table raises (the ROADMAP's list, in its order).

1. The receiver's compression is computed at the ADC's two samples a chip:
   what is it in continuous time, and how much rides on the limiter's knee?
2. The receiver's gain setting as a variable: which setting leaves the most
   SNR with the leakage there?
3. A target beyond the clean window: where does it land and how strong?
4. What does a fast target lose inside one packet?
5. The ADC that nothing dithers, as a statistic over the target's delay, and
   with a dither.
6. The leakage's delay as a variable, and a near reflector as a second path.
7. Every row on one LO.
8. The bistatic rows' leftover carrier offset.

    python examples/01_receiver_and_window.py

About 130 s.  ``tools/example_drift.py`` compares the output with
``examples/expected/01_receiver_and_window.txt``.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jcas60 import studies  # noqa: E402

if __name__ == "__main__":
    sys.stdout.write(studies.render())
