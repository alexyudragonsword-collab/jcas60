"""Six files of fmcw_radar (see ``jcas60.vendor``): detection, CFAR, windows,
range-Doppler reading, the regulatory check, and its example 00 as printed.

``C0`` is here for ``processing/rd.py``: upstream it comes from
``config/schema.py`` (``C0 = 299_792_458.0``), which is not copied.
"""
from pathlib import Path

C0 = 299_792_458.0

HERE = Path(__file__).resolve().parent
EXAMPLE_00 = HERE / "examples" / "expected" / "00_single_target_rd_map.txt"
