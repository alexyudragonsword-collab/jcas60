#!/usr/bin/env python3
"""Example 01: three questions the comparison table raises (ROADMAP items 1-3).

1. The receiver's compression is computed at the ADC's two samples a chip:
   what is it in continuous time, and how much rides on the limiter's knee?
2. The receiver's gain setting as a variable: which setting leaves the most
   SNR with the leakage there?
3. A target beyond the clean window: where does it land and how strong?

    python examples/01_receiver_and_window.py

About a minute.  ``tools/example_drift.py`` compares the output with
``examples/expected/01_receiver_and_window.txt``.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from jcas60 import studies  # noqa: E402

if __name__ == "__main__":
    sys.stdout.write(studies.render())
