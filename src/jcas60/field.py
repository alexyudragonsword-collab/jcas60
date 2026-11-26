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
    """``ga128``, ``gb128``, ``ga64``, ``gb64`` as the entry tabulates them, in transmit order."""
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


def trn_r_subfields_per_unit() -> int:
    """How many subfields one unit of a TRN-R packet's field has, from the entry's text (``L_T = 0, L_P = 0, L_M = <n>``)."""
    import re
    text = str(S.raw("edmg.trn_field#trn_r_packet"))
    got = dict((k, int(v)) for k, v in re.findall(r"(L_[TPM]) = (\d+)", text))
    if set(got) != {"L_T", "L_P", "L_M"} or got["L_T"] or got["L_P"]:
        raise ValueError(f"edmg.trn_field#trn_r_packet is not 'L_T = 0, L_P = 0, L_M = <n>': {text!r}")
    return got["L_M"]


@lru_cache(maxsize=4)
def trn(units: int, blocks: tuple) -> Field:
    """The TRN field of an EDMG TRN-R packet, ``units`` units long, with the correlation windows on the ``blocks``
    (1-based) of every subfield.

    The basic subfield is six Golay blocks, [Ga, -Gb, Ga, Gb, Ga, -Gb] (``edmg.trn_field``).  For a path delayed
    by up to one block, a block's window holds the tail of the block before it; correlating blocks 2-5 with their
    own signs, the two Ga and the two Gb give R_a + R_b = 0 off the peak and the four neighbour terms come out
    -C_ab - C_ba + C_ab + C_ba = 0: no sidelobe within a block either way, blocks 1 and 6 standing guard.  With
    all six the neighbour terms leave -(C_ab + C_ba).  Built by mmw60_trx's ``compose`` and ``rotate`` (vendored)
    from the entry's block list."""
    n = int(S.num("edmg.trn_field#trn_bl@trn_seq_length=normal").value)
    ga, gb = tabulated(f"ga{n}"), tabulated(f"gb{n}")
    items = list(S.raw("edmg.trn_field#basic_subfield"))
    sub = P.compose({"Ga_N": ga, "Gb_N": gb}, items)
    count = int(units) * trn_r_subfields_per_unit()
    chips = P.rotate(np.tile(sub, count))                       # n = 0 at the field's first chip: the packet before it is a multiple of four chips
    chips.setflags(write=False)
    runs, start = [], None                                      # consecutive blocks become one window
    for b in sorted(int(x) for x in blocks) + [None]:
        if start is None:
            start = prev = b
        elif b is not None and b == prev + 1:
            prev = b
        else:
            runs.append((start, prev))
            start = prev = b
    wins = tuple((k * sub.size + (b0 - 1) * n, chips[k * sub.size + (b0 - 1) * n: k * sub.size + b1 * n])
                 for k in range(count) for b0, b1 in runs)
    if not all(1 <= b <= len(items) for b in blocks):
        raise ValueError(f"blocks must be 1 .. {len(items)}")
    return Field(chips=chips, n0=0, windows=wins, zero_sidelobe_chips=n)


@lru_cache(maxsize=4)
def designed(units: int) -> Field:
    """A waveform for the TRN field made for sensing (``jcas60.sensing_waveform``), as long as the TRN field of
    ``units`` units: Ga / Gb alternating, the signs periodic over the pattern's blocks.  With the pattern built on
    the length-4 perfect binary sequence the neighbours' cross-correlations cancel over one period, so every block is
    correlated -- all but the guard blocks at the two ends, chosen so that a whole number of periods is left (a
    partial period does not cancel).  The windows are one period long; the clean zone is ``clean_blocks`` blocks
    either way, which a test holds against the chips."""
    n = int(S.num("edmg.trn_field#trn_bl@trn_seq_length=normal").value)
    ga, gb = tabulated(f"ga{n}"), tabulated(f"gb{n}")
    items = list(S.raw("jcas60.sensing_waveform#pattern"))
    guard = int(S.num("jcas60.sensing_waveform#guard_blocks").value)
    period = P.compose({"Ga_N": ga, "Gb_N": gb}, items)
    blocks = int(units) * trn_r_subfields_per_unit() * len(S.raw("edmg.trn_field#basic_subfield"))
    if blocks % len(items):
        raise ValueError(f"the TRN field's {blocks} blocks are not a whole number of the pattern's {len(items)}")
    chips = P.rotate(np.tile(period, blocks // len(items)))
    chips.setflags(write=False)
    first, last = guard, blocks - guard
    starts = range(first, last, len(items))
    if (last - first) % len(items):
        raise ValueError("the guard blocks do not leave a whole number of periods to correlate")
    wins = tuple((b * n, chips[b * n: (b + len(items)) * n]) for b in starts)
    return Field(chips=chips, n0=0, windows=wins, zero_sidelobe_chips=int(S.num("jcas60.sensing_waveform#clean_blocks").value) * n)
