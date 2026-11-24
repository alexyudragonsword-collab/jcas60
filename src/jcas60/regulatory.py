"""Which EIRP limit each row is held to, and what else the rule asks (law 3).

47 CFR 15.255 has two sets of limits in 57-71 GHz: paragraph (c)(1) for
"devices other than field disturbance sensors" (the communication limit, an
average and a peak EIRP, no off-time), and (c)(2) for field disturbance
sensors / radars (tiers by band, some with an off-time).  Row A is a radar:
(c)(2), as fmcw_radar's own check reports.  Whether a communication device
that senses with its own waveform is one or the other is not settled here
(spec_sources#54), so each ay row is computed under both readings and both
are printed.  No row is ever compared at a limit that is not its own.

The radar reading runs fmcw_radar's own check (vendored
``link/compliance.py``, which holds the rule's structure and reads every
limit through ``jcas60.sources.fcc_number``) on a description of the ay
burst: the channel as the occupied band, one sensing packet where the radar
has one chirp.  The check reads a config object; it is given the few
attributes it reads and nothing else.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from . import sources as S
from .config import Setup
from .sources import Num, derive
from .vendor.fmcw_radar.link import compliance as _compliance


@dataclass(frozen=True)
class Reading:
    name: str                    # "communication limit" | "radar limit"
    rule: str
    cap_dbm: Num                 # the EIRP limit this reading holds the row to
    cap_kind: str                # "average EIRP" | "peak EIRP"
    eirp_dbm: Num                # what the shared front end radiates under it
    backoff_db: Num
    off_time: str                # what the rule asks besides
    off_time_required_s: Num | None
    off_time_window_s: Num | None
    off_time_have_s: float | None
    tier: str
    detail: object = None        # the vendored check's result (radar reading)


def communication(s: Setup) -> Reading:
    kind = str(S.raw("jcas60.comm_limit_kind"))
    cap = S.num(f"fcc15255.eirp_{'avg' if kind == 'average' else 'peak'}_dbm", "dBm")
    return Reading("communication limit", "47 CFR 15.255(c)(1)(i)", cap, f"{kind} EIRP", s.fe.eirp_dbm(cap), s.fe.backoff_db(cap),
                   "none in the rule", None, None, None, "devices other than field disturbance sensors")


def _view(band_hz, n_bursts: int, burst_s: float, pri_s: float, frame_period_s: float, tx_dbm: Num, backoff_db: float, gain: Num):
    """What ``compliance.check`` reads of a config, and nothing more."""
    lo, hi = band_hz
    values = {"tx_power_dbm": (tx_dbm.value, tx_dbm.tag), "tx_backoff_db": (backoff_db, "derived"), "tx_gain_dbi": (gain.value, gain.tag)}
    span = (n_bursts - 1) * pri_s + burst_s
    return SimpleNamespace(
        band=SimpleNamespace(center_hz=0.5 * (lo + hi), sweep_hz=hi - lo),
        chirp=SimpleNamespace(n_chirps=n_bursts, pri_s=pri_s, duration_s=burst_s, frame_on_s=span, frame_period_s=frame_period_s),
        mimo=SimpleNamespace(n_tx=1), regulatory=SimpleNamespace(scenario=None),
        frontend_value=lambda name: values[name])


def radar_check(band_hz, n_bursts: int, burst_s: float, pri_s: float, frame_period_s: float, tx_dbm: Num, gain: Num):
    """fmcw_radar's check, run twice: at full power to learn the governing tier's cap, then backed off to it."""
    first = _compliance.check(_view(band_hz, n_bursts, burst_s, pri_s, frame_period_s, tx_dbm, 0.0, gain))
    if first.eirp_cap_dbm is None:
        return first
    backoff = max(0.0, tx_dbm.value + gain.value - first.eirp_cap_dbm)
    return _compliance.check(_view(band_hz, n_bursts, burst_s, pri_s, frame_period_s, tx_dbm, backoff, gain))


def radar(s: Setup) -> Reading:
    res = radar_check(s.band_hz, int(s.n_packets.value), s.packet_s.value, s.interval_s.value, 1.0 / s.frames_per_s.value,
                      s.fe.tx_power_dbm, s.fe.tx_gain_dbi)
    g = res.tier
    if g is None or not g.compliant:
        raise ValueError(f"the ay sensing burst has no usable tier of {res.rule}: {res.verdict}")
    cap = S.num(f"fcc15255.{g.tier}#peak_eirp_dbm", "dBm")
    need = None if g.off_time_required_s is None else S.num(f"fcc15255.{g.tier}#off_time_sum_min_s", "s")
    win = None if g.off_time_required_s is None else S.num(f"fcc15255.{g.tier}#off_time_window_s", "s")
    return Reading("radar limit", res.rule, cap, "peak EIRP", s.fe.eirp_dbm(cap), s.fe.backoff_db(cap),
                   g.off_time, need, win, g.off_time_have_s, g.label, res)


def readings(s: Setup) -> tuple[Reading, Reading]:
    return communication(s), radar(s)


# ----------------------------------------------------------------------------- air time
@dataclass(frozen=True)
class AirTime:
    fraction: Num                # share of the time the front end is not available to communication
    fraction_upper: Num | None   # radar reading: if the off-time governs every emission of the device
    rate_mbps: Num               # the communication rate that share is taken from
    loss_mbps: Num
    loss_upper_mbps: Num | None


def comm_rate_mbps() -> Num:
    """mmw60_trx's comparison table, the ay_sc row on channel 2: the rate at the cap its calibrated chain allows."""
    import re
    from .vendor.mmw60 import COMPARISON_TABLE
    for line in COMPARISON_TABLE.read_text(encoding="utf-8").splitlines():
        m = re.match(r"ay_sc EDMG channel 2 @60\.48 GHz\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+\S+\s+\S+\s+\S+\s+\S+\s+(MCS\d+)\s+(\d+)\s", line)
        if m:
            return Num(float(m.group(2)), "derived", f"mmw60_trx examples/expected/04_comparison_table.txt: ay_sc EDMG channel 2, calibrated cap {m.group(1)}",
                       "Mb/s", "assumed")
    raise ValueError("mmw60_trx's comparison table no longer has the ay_sc row in the form read here")


def air_time(s: Setup, reading: Reading) -> AirTime:
    rate = comm_rate_mbps()
    per_packet = s.packet_s.value + S.num("jcas60.packet_overhead#sifs_per_packet").value * s.sifs_s.value
    frac = derive(per_packet * s.n_packets.value * s.frames_per_s.value, "(packet + SIFS) x packets x frames per second",
                  s.packet_s, s.sifs_s, s.n_packets, s.frames_per_s)
    upper = None
    if reading.off_time_required_s is not None:
        upper = derive(reading.off_time_required_s.value / reading.off_time_window_s.value, "off-time required / window",
                       reading.off_time_required_s, reading.off_time_window_s, S.num("jcas60.sensing_burst#frames_per_s"))
    loss = derive(frac.value * rate.value, "fraction x rate", frac, rate, unit="Mb/s")
    loss_up = None if upper is None else derive(upper.value * rate.value, "fraction x rate", upper, rate, unit="Mb/s")
    return AirTime(frac, upper, rate, loss, loss_up)
