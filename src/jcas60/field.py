"""The single-carrier Golay training field used for sensing: the L-CEF.

Built from two spec_sources entries and nothing else:

* ``dmg.golay_sequences`` -- Ga128 / Gb128 as tabulated (and the recursion's
  D_k / W_k, which the vendored ``golay_pair`` turns into the same sequences:
  a test holds the two equal);
* ``dmg.cef_golay_layout`` -- Gu512 and Gv512 as four signed segments each,
  the SC order Gu512, Gv512, Gv128, and the pi/2 rotation.

Layout (chips; n = 0 is the first chip of the CEF)::

    n = -128 .. -1     the STF's last 128 chips = -Ga128 = Gu512's last 128: its cyclic prefix
    n =    0 .. 511    Gu512
    n =  512 .. 1023   Gv512        (Gu512's last 128 = Gv512's last 128: Gv512's cyclic prefix)
    n = 1024 .. 1151   Gv128 = Gv512's first 128: the cyclic postfix

every chip times exp(j pi n / 2).  For a path delayed by 0 .. 128 chips each
512-chip window therefore holds a cyclic shift of its sequence, and the two
cyclic correlations sum to 1024 at the path's delay and to zero elsewhere:
the pair is *periodically* complementary (spec_sources checks this; so does
``tests/test_field.py``).  That is exact for a path that does not move.

This is a communication waveform's construction and belongs in mmw60_trx
(mmw60_trx#89); it lives here until that library builds it.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from . import sources as S
from .vendor.mmw60.ay.golay import golay_pair


def tabulated(name: str) -> np.ndarray:
    """``ga128`` ... ``gb32`` as the entry tabulates them, in transmit order."""
    return np.asarray(S.raw(f"dmg.golay_sequences#{name}"), dtype=float)


def from_recursion(n: int) -> tuple[np.ndarray, np.ndarray]:
    """(Ga_n, Gb_n) from the entry's D_k / W_k through the vendored recursion,
    reversed as the entry says (``Ga(n) = A_K(N - 1 - n)``)."""
    a, b = golay_pair(S.raw(f"dmg.golay_sequences#d_{n}"), S.raw(f"dmg.golay_sequences#w_{n}"))
    return a[::-1].copy(), b[::-1].copy()


def _concat(segments, ga, gb) -> np.ndarray:
    parts = {"Ga128": ga, "Gb128": gb}
    return np.concatenate([(-1.0 if s[0] == "-" else 1.0) * parts[s.lstrip("+-")] for s in segments])


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
    gu = _concat(S.raw("dmg.cef_golay_layout#gu512"), ga, gb)
    gv = _concat(S.raw("dmg.cef_golay_layout#gv512"), ga, gb)
    if list(S.raw("dmg.cef_golay_layout#sc_order")) != ["Gu512", "Gv512", "Gv128"]:
        raise ValueError("dmg.cef_golay_layout: the SC order is not Gu512, Gv512, Gv128")
    seg = len(ga)
    prefix, postfix = gu[-seg:], gv[:seg]
    if not (np.array_equal(prefix, -ga) and np.array_equal(gv[-seg:], prefix)):
        raise ValueError("the CEF's cyclic prefixes are not -Ga128: the layout entry and the sequences disagree")
    flat = np.concatenate([prefix, gu, gv, postfix])
    if flat.size - seg != int(S.num("dmg.cef_golay_layout#cef_chips").value):
        raise ValueError("the CEF is not the length dmg.cef_golay_layout gives")
    n = np.arange(-seg, flat.size - seg)
    chips = flat * np.exp(0.5j * np.pi * n)                       # the entry's rotation exp(j pi n / 2)
    chips.setflags(write=False)
    wins = tuple((seg + k * len(gu), chips[seg + k * len(gu): seg + (k + 1) * len(gu)]) for k in range(2))
    return Field(chips=chips, n0=seg, windows=wins, zero_sidelobe_chips=seg)
