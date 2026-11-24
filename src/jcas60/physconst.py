"""Physical constants and the dB helpers.  The only module allowed to hold
multi-digit literals (``tests/test_laws.py::test_law1_no_stray_numbers``)."""
from __future__ import annotations

import math

from .vendor.fmcw_radar import C0  # noqa: F401  (re-exported: the value fmcw_radar's processing uses)

K_BOLTZMANN = 1.380649e-23      # J/K, exact (SI 2019); fmcw_radar's config schema has the same
T0_K = 290.0                    # K, the noise-figure reference temperature


def db(x: float) -> float:
    return 10.0 * math.log10(x)


def undb(x: float) -> float:
    return 10.0 ** (x / 10.0)


#: k T0 in dBm/Hz (-173.98)
KT0_DBM_HZ = db(K_BOLTZMANN * T0_K) + 30.0
