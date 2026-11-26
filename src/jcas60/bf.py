"""The sensing frame of this table against spec_sources' 802.11bf entries, item by item (law 7).

spec_sources has entries for 802.11bf at 60 GHz since spec_sources#54, and
since spec_sources#57 the TRN field's first-stream Golay pair and rotation,
what the adopted text says of the monostatic PPDU's TRN type and data field,
and the bistatic type's PPDU -- taken from the task groups' public documents,
not from the standard, which this library does not read either.  They do not make this table an 802.11bf
implementation: no sensing session, no measurement exchange, no frame of the
procedure is built here.  What can be said is narrower and is said here: for
each choice of the sensing frame, what the entry holds, what the table does,
and whether the one fits the other.  Nothing in this module feeds a number of
the table; it reads the same :class:`~jcas60.config.Setup` and reports.

Field widths give the largest value a field can carry (2**(8 octets) - 1);
the entries' source gives no range of its own, so those bounds are derived.

Four verdicts, each written as found: ``fits``; ``does not fit`` (the entry
says one thing, the table does another -- row C's packets against the
bistatic type's PPDU and exchange); ``left open by the entry`` (the entry
says its source does not state it); ``no entry``.
"""
from __future__ import annotations

from . import sources as S
from .config import Setup

TRN = "edmg.trn_field"
TYPES = "ieee80211bf.dmg_sensing_types"
SCHED = "ieee80211bf.dmg_sensing_scheduling"
MONO = "ieee80211bf.dmg_monostatic_ppdu"
BW = "ieee80211bf.dmg_sensing_bandwidth"
BI = "ieee80211bf.dmg_bistatic_ppdu"
GOLAY = "dmg.golay_sequences"

FITS, DIFFERS, OPEN, NO_ENTRY = "fits", "does not fit", "left open by the entry", "no entry"


def _of_type(field: str, kind: str) -> str:
    return str(S.raw(f"{TYPES}#{field}@type={kind}"))


def _field_max(name: str) -> int:
    """The largest value the scheduling field ``name`` can carry, from its width in octets."""
    return 2 ** (8 * int(S.num(f"{SCHED}#octets@field={name}").value)) - 1


def checks(s: Setup) -> list[dict]:
    """One dict per item: ``row`` (B, C or both), ``item``, ``key``, ``entry`` (what it holds), ``here`` (what the
    table does), ``verdict`` (``fits`` | ``does not fit`` | ``left open by the entry`` | ``no entry``)."""
    unit_us = S.num(f"{SCHED}#time_unit_s").value * 1e6
    intra_us, inter_us = s.interval_s.value * 1e6 / unit_us, 1e6 / s.frames_per_s.value / unit_us
    n = int(s.n_packets.value)
    lim = {k: _field_max(k) for k in ("intraburst_interval", "interburst_interval", "number_of_instances_per_burst")}
    whole = abs(inter_us - round(inter_us)) < 1e-9
    before_trn = int(s.packet_chips.value)                       # chips of the packet ahead of the TRN field
    bi_ppdu = " / ".join(S.raw(BI + "#ppdu"))
    brpifs_us = S.num("dmg.phy_timing#brpifs_s").value * 1e6
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
         "verdict": FITS if s.trn_units.value <= S.num(TRN + '#edmg_trn_len_max').value else DIFFERS},
        {"row": "B, C", "item": "the TRN field's Golay pair", "key": TRN,
         "entry": f"blocks {' '.join(S.raw(TRN + '#basic_subfield'))} of length {S.raw(TRN + '#golay_length')}; first stream: {S.raw(TRN + '#stream1_golay')}",
         "here": f"Ga128 / Gb128 read from {GOLAY} by key (the pair that entry holds since spec_sources#57)",
         "verdict": FITS if GOLAY in str(S.raw(TRN + "#stream1_golay")) else OPEN},
        {"row": "B, C", "item": "the TRN field's rotation", "key": TRN,
         "entry": f"{S.raw(TRN + '#modulation')}; where the count of chips starts: {S.raw(TRN + '#rotation_origin')}",
         "here": f"a quarter turn a chip; the field starts at chip {before_trn} = 4 x {before_trn // 4} of the packet, so counting from the packet's "
                 "first chip or from the field's gives the same chips",
         "verdict": FITS if "pi/2" in str(S.raw(TRN + "#modulation")) and before_trn % 4 == 0 else OPEN},
        {"row": "B", "item": "which TRN packet", "key": MONO,
         "entry": f"TRN type of the monostatic PPDU: {S.raw(MONO + '#trn_type')}",
         "here": "TRN-R (jcas60.trn_sensing)", "verdict": FITS if "TRN-R" in str(S.raw(MONO + "#trn_type")) else OPEN},
        {"row": "B", "item": "the data field", "key": MONO,
         "entry": f"PSDU: {S.raw(MONO + '#psdu')}; for the monostatic type on its own the entry holds nothing",
         "here": "none: preamble, the two headers, then the TRN field -- so not the coordinated type's PPDU, which this table does not build", "verdict": OPEN},
        {"row": "C", "item": "who transmits and who receives", "key": TYPES,
         "entry": f"bistatic: transmitter {_of_type('transmitter', 'bistatic')}, receiver {_of_type('receiver', 'bistatic')}",
         "here": "one device transmits, the other correlates", "verdict": FITS},
        {"row": "C", "item": "the PPDU, read on the L-CEF", "key": BI,
         "entry": f"{bi_ppdu} PPDU with a TRN field; Sync field: {S.raw(BI + '#sync_field')}",
         "here": "the plain SC packet of row B, which has no TRN field", "verdict": DIFFERS},
        {"row": "C", "item": "the PPDU, with the TRN field", "key": BI,
         "entry": f"{bi_ppdu} PPDU with a TRN field; Sync field: {S.raw(BI + '#sync_field')}",
         "here": "a TRN-R field behind the packet and no Sync field, as a BRP-RX PPDU has -- but no data field, so no BRP frame is in it "
                 f"(a BRP frame's SC packet carries at least {S.num(BI + '#min_data_sc_blocks').value:g} SC blocks)", "verdict": DIFFERS},
        {"row": "C", "item": "the exchange", "key": BI,
         "entry": f"{S.raw(BI + '#tx_initiator')}; {S.raw(BI + '#spacing')}; feedback: {S.raw(BI + '#feedback')}",
         "here": f"{n} packets one way; no BRP frame comes back, no BRPIFS ({brpifs_us:g} us) and no feedback are in the air time of section 6 "
                 "(example 01 section 14 prices the exchange)", "verdict": DIFFERS},
        {"row": "C", "item": "the PHY mode", "key": BI,
         "entry": f"{' or '.join(S.raw(BI + '#phy_modes'))} mode",
         "here": "SC", "verdict": FITS if "SC" in S.raw(BI + "#phy_modes") else DIFFERS},
        {"row": "C", "item": "the TRN field's format", "key": BI,
         "entry": f"set by {S.raw(BI + '#trn_format_set_by')}: no value is fixed",
         "here": f"{s.trn_units.value:g} units of a TRN-R packet, Golay length 128 (jcas60.trn_sensing)", "verdict": FITS},
        {"row": "B, C", "item": "packets of one frame, start to start", "key": SCHED,
         "entry": f"intraburst interval: {S.num(SCHED + '#octets@field=intraburst_interval').value:g} octet, in us: at most {lim['intraburst_interval']} us",
         "here": f"{intra_us:g} us", "verdict": FITS if intra_us <= lim["intraburst_interval"] else DIFFERS},
        {"row": "B, C", "item": "frames, start to start", "key": SCHED,
         "entry": f"interburst interval: {S.num(SCHED + '#octets@field=interburst_interval').value:g} octets, in us: at most {lim['interburst_interval']} us",
         "here": f"{inter_us:.2f} us" + ("" if whole else f" (not a whole number of us: {round(inter_us)} us is {1e6 / round(inter_us):.4f} frames a second)"),
         "verdict": FITS if inter_us <= lim["interburst_interval"] else DIFFERS},
        {"row": "B, C", "item": "packets per frame", "key": SCHED,
         "entry": f"measurement exchanges per burst: {S.num(SCHED + '#octets@field=number_of_instances_per_burst').value:g} octet: at most {lim['number_of_instances_per_burst']}",
         "here": f"{n}, taking one packet for one measurement exchange", "verdict": FITS if n <= lim["number_of_instances_per_burst"] else DIFFERS},
        {"row": "B, C", "item": "bandwidth", "key": BW,
         "entry": f"sensing by a non-EDMG DMG STA limited to one 2.16 GHz channel: {S.raw(BW + '#dmg_non_edmg')}; multistatic on {S.raw(BW + '#multistatic_channel_ghz')} GHz channels",
         "here": f"one {s.channel_width_hz.value / 1e9:g} GHz channel", "verdict": OPEN},
    ]
    return out
