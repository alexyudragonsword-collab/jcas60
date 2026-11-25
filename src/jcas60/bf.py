"""The sensing frame of this table against spec_sources' 802.11bf entries, item by item (law 7).

spec_sources has entries for 802.11bf at 60 GHz since spec_sources#54 --
taken from the task group's public documents, not from the standard, which
this library does not read either.  They do not make this table an 802.11bf
implementation: no sensing session, no measurement exchange, no frame of the
procedure is built here.  What can be said is narrower and is said here: for
each choice of the sensing frame, what the entry holds, what the table does,
and whether the one fits the other.  Nothing in this module feeds a number of
the table; it reads the same :class:`~jcas60.config.Setup` and reports.

Field widths give the largest value a field can carry (2**(8 octets) - 1);
the entries' source gives no range of its own, so those bounds are derived.
"""
from __future__ import annotations

from . import sources as S
from .config import Setup

TRN = "edmg.trn_field"
TYPES = "ieee80211bf.dmg_sensing_types"
SCHED = "ieee80211bf.dmg_sensing_scheduling"
MONO = "ieee80211bf.dmg_monostatic_ppdu"
BW = "ieee80211bf.dmg_sensing_bandwidth"

FITS, OPEN, NO_ENTRY = "fits", "left open by the entry", "no entry"


def _of_type(field: str, kind: str) -> str:
    return str(S.raw(f"{TYPES}#{field}@type={kind}"))


def _field_max(name: str) -> int:
    """The largest value the scheduling field ``name`` can carry, from its width in octets."""
    return 2 ** (8 * int(S.num(f"{SCHED}#octets@field={name}").value)) - 1


def checks(s: Setup) -> list[dict]:
    """One dict per item: ``row`` (B, C or both), ``item``, ``key``, ``entry`` (what it holds), ``here`` (what the
    table does), ``verdict`` (``fits`` | ``left open by the entry`` | ``no entry``)."""
    unit_us = S.num(f"{SCHED}#time_unit_s").value * 1e6
    intra_us, inter_us = s.interval_s.value * 1e6 / unit_us, 1e6 / s.frames_per_s.value / unit_us
    n = int(s.n_packets.value)
    lim = {k: _field_max(k) for k in ("intraburst_interval", "interburst_interval", "number_of_instances_per_burst")}
    whole = abs(inter_us - round(inter_us)) < 1e-9
    out = [
        {"row": "B", "item": "who transmits and who receives", "key": TYPES,
         "entry": f"monostatic: transmitter {_of_type('transmitter', 'monostatic')}, receiver {_of_type('receiver', 'monostatic')}",
         "here": "one device transmits the packet and correlates its own echo", "verdict": FITS},
        {"row": "B", "item": "receiving while transmitting", "key": TYPES,
         "entry": f"full-duplex requirement: {S.raw(TYPES + '#full_duplex_requirement')}",
         "here": "the receiver is on during the packet (section 5 is what that costs)", "verdict": OPEN},
        {"row": "B", "item": "the PPDU", "key": MONO,
         "entry": f"{S.raw(MONO + '#ppdu')}",
         "here": f"a single-channel EDMG SC PPDU of {s.packet_chips.value:g} chips: preamble and the two headers, no data, no TRN field", "verdict": FITS},
        {"row": "B", "item": "the field the channel is read on", "key": MONO,
         "entry": f"sensing on the TRN field: {S.raw(MONO + '#trn_based_sensing')}; which field otherwise is not a value of the entry",
         "here": "the L-CEF (jcas60.sensing_field)", "verdict": OPEN},
        {"row": "B", "item": "sensing on the TRN field", "key": MONO,
         "entry": f"sensing on the TRN field: {S.raw(MONO + '#trn_based_sensing')}",
         "here": f"the TRN rows read the channel on the TRN field ({s.trn_field.n_correlated} chips of it a packet)", "verdict": FITS},
        {"row": "B, C", "item": "the TRN field's length", "key": TRN,
         "entry": f"EDMG_TRN_LEN at most {S.num(TRN + '#edmg_trn_len_max').value:g} units; a TRN-R packet's unit: {S.raw(TRN + '#trn_r_packet')}",
         "here": f"{s.trn_units.value:g} units of a TRN-R packet, {s.trn_field.chips.size} chips",
         "verdict": FITS if s.trn_units.value <= S.num(TRN + '#edmg_trn_len_max').value else "does not fit"},
        {"row": "B, C", "item": "the TRN field's Golay pair and rotation", "key": TRN,
         "entry": f"blocks {' '.join(S.raw(TRN + '#basic_subfield'))} of length {S.raw(TRN + '#golay_length')}; the first stream's pair and the field's rotation are not values of the entry",
         "here": "Ga128 / Gb128 of dmg.golay_sequences, turned like the preamble (jcas60.trn_sensing)", "verdict": OPEN},
        {"row": "C", "item": "who transmits and who receives", "key": TYPES,
         "entry": f"bistatic: transmitter {_of_type('transmitter', 'bistatic')}, receiver {_of_type('receiver', 'bistatic')}",
         "here": "one device transmits, the other correlates", "verdict": FITS},
        {"row": "C", "item": "the PPDU", "key": TYPES,
         "entry": "the entries give a PPDU for the monostatic types and for the multistatic type, none for the bistatic type",
         "here": "the same PPDUs as row B: the plain SC packet read on its L-CEF, or with the TRN field behind it", "verdict": NO_ENTRY},
        {"row": "B, C", "item": "packets of one frame, start to start", "key": SCHED,
         "entry": f"intraburst interval: {S.num(SCHED + '#octets@field=intraburst_interval').value:g} octet, in us: at most {lim['intraburst_interval']} us",
         "here": f"{intra_us:g} us", "verdict": FITS if intra_us <= lim["intraburst_interval"] else "does not fit"},
        {"row": "B, C", "item": "frames, start to start", "key": SCHED,
         "entry": f"interburst interval: {S.num(SCHED + '#octets@field=interburst_interval').value:g} octets, in us: at most {lim['interburst_interval']} us",
         "here": f"{inter_us:.2f} us" + ("" if whole else f" (not a whole number of us: {round(inter_us)} us is {1e6 / round(inter_us):.4f} frames a second)"),
         "verdict": FITS if inter_us <= lim["interburst_interval"] else "does not fit"},
        {"row": "B, C", "item": "packets per frame", "key": SCHED,
         "entry": f"measurement exchanges per burst: {S.num(SCHED + '#octets@field=number_of_instances_per_burst').value:g} octet: at most {lim['number_of_instances_per_burst']}",
         "here": f"{n}, taking one packet for one measurement exchange", "verdict": FITS if n <= lim["number_of_instances_per_burst"] else "does not fit"},
        {"row": "B, C", "item": "bandwidth", "key": BW,
         "entry": f"sensing by a non-EDMG DMG STA limited to one 2.16 GHz channel: {S.raw(BW + '#dmg_non_edmg')}; multistatic on {S.raw(BW + '#multistatic_channel_ghz')} GHz channels",
         "here": f"one {s.channel_width_hz.value / 1e9:g} GHz channel", "verdict": OPEN},
    ]
    return out
