"""The single-carrier Golay training field used for sensing: the L-CEF.

Built by mmw60_trx's ``ay/preamble.py`` (vendored; mmw60_trx#89) from two
spec_sources entries, which this module reads and hands over as arguments:

* ``dmg.golay_sequences`` -- Ga128 / Gb128 as tabulated (and the recursion's
  D_k / W_k, which give the same sequences: a test holds the two equal);
* ``dmg.cef_golay_layout`` -- Gu512 and Gv512 as four signed blocks each and
  the SC order Gu512, Gv512, Gv128; the pi/2 rotation is the builder's.

Layout (chips; n = 0 is the first chip of the CEF)::

    n = -128 .. -1     the cyclic prefix: Gu512's last 128 chips = -Ga128 (the STF's last block)
    n =    0 .. 511    Gu512
    n =  512 .. 1023   Gv512        (Gu512's last 128 = Gv512's last 128: Gv512's cyclic prefix)
    n = 1024 .. 1151   Gv128 = Gv512's first 128: the cyclic postfix

every chip times exp(j pi n / 2).  For a path delayed by 0 .. 128 chips each
512-chip window therefore holds a cyclic shift of its sequence, and the two
cyclic correlations sum to 1024 at the path's delay and to zero elsewhere.
Inside that zone each sequence is already clean on its own (it is made of
Ga128 / Gb128); as a pair the two are periodically complementary at every
lag.  That is exact for a path that does not move.

What stays here is the reading of the entries and the correlation windows
the waveform simulation uses (it correlates at two samples per chip, which
the vendored chip-rate ``correlate`` does not).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from . import sources as S
from .vendor.mmw60.ay import preamble as P


def tabulated(name: str) -> np.ndarray:
    """``ga128`` ... ``gb32`` as the entry tabulates them, in transmit order."""
    return np.asarray(S.raw(f"dmg.golay_sequences#{name}"), dtype=float)


def from_recursion(n: int) -> tuple[np.ndarray, np.ndarray]:
    """(Ga_n, Gb_n) from the entry's D_k / W_k: the recursion read backwards, as the entry says."""
    return P.pair_from_recursion(S.raw(f"dmg.golay_sequences#d_{n}"), S.raw(f"dmg.golay_sequences#w_{n}"))


@dataclass(frozen=True)
class Field:
    chips: np.ndarray          # the transmitted chips, rotated: prefix, Gu512, Gv512, Gv128
    n0: int                    # index in ``chips`` of the CEF's first chip (the prefix length)
    windows: tuple             # ((start, reference), ...): each correlation window's start in ``chips`` and its rotated sequence
    zero_sidelobe_chips: int   # delays 0 .. this many chips keep the cyclic property

    @property
    def n_correlated(self) -> int:
        return sum(len(r) for _, r in self.windows)


@lru_cache(maxsize=1)
def l_cef() -> Field:
    ga, gb = tabulated("ga128"), tabulated("gb128")
    gu, gv = P.gu_gv(ga, gb, S.raw("dmg.cef_golay_layout#gu512"), S.raw("dmg.cef_golay_layout#gv512"))
    order = list(S.raw("dmg.cef_golay_layout#sc_order"))
    if order != ["Gu512", "Gv512", "Gv128"]:
        raise ValueError("dmg.cef_golay_layout: the SC order is not Gu512, Gv512, Gv128")
    seg = len(ga)
    if not np.array_equal(P.cyclic_prefix(gu, gv, order, rotated=False), -ga):
        raise ValueError("the CEF's cyclic prefix is not -Ga128: the layout entry and the sequences disagree")
    field = P.cef(gu, gv, order)
    if field.size != int(S.num("dmg.cef_golay_layout#cef_chips").value):
        raise ValueError("the CEF is not the length dmg.cef_golay_layout gives")
    chips = np.concatenate([P.cyclic_prefix(gu, gv, order), field])
    chips.setflags(write=False)
    wins = tuple((seg + k * len(gu), chips[seg + k * len(gu): seg + (k + 1) * len(gu)]) for k in range(2))
    return Field(chips=chips, n0=seg, windows=wins, zero_sidelobe_chips=seg)
