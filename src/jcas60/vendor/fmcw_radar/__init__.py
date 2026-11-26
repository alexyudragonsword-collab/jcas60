"""Eight files of fmcw_radar (see ``jcas60.vendor``): detection, CFAR and its
statistics, windows, range-Doppler reading, the regulatory check, the limiter,
and its example 00 as printed.

``C0`` is the value fmcw_radar's ``config/schema.py`` holds
(``C0 = 299_792_458.0``; that file is not copied).  It was here for
``processing/rd.py``, which since fmcw_radar#44 no longer needs it at load;
``jcas60.physconst`` still re-exports it from here.
"""
from pathlib import Path

C0 = 299_792_458.0

HERE = Path(__file__).resolve().parent
EXAMPLE_00 = HERE / "examples" / "expected" / "00_single_target_rd_map.txt"
