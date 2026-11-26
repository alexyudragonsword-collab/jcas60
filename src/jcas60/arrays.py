"""Milestone 1, first step: a transmit array.

One antenna cannot reach the communication limit: 12 dBm into 10 dBi is
22 dBm of EIRP against 40 dBm.  N elements on a line, each behind its own PA
and steered in phase, add 20 log10(N) on boresight -- and with that a beam
that has to be pointed.  This module answers, for each N and under both
regulatory readings: the EIRP and the back-off, the SNR of the rows of the
comparison table, the beam's width, how many directions cover the field of
view, what a sweep of them costs in air time, and what the phase steering
does to a waveform 2.2 GHz wide (the squint) -- and what N power amplifiers
beside the receiver do to the monostatic row's leakage.

The array itself is arraytx's (``geometry``, ``beam``; vendored through
mmw60_trx): its geometry, its array factor, its steering and its far field.
Two roads where a conclusion is drawn (law 4): the textbook closed forms
against arraytx's numerical array factor for the width and the count of
beams; the closed-form link budget against the waveform simulation for the
SNR; a band integral against the sensing waveform driven element by element
and summed in the far field for the squint; the three leakage mechanisms in
closed form against the simulated receiver for the leakage.

The receiver stays one antenna here, and the table's rows are not replaced:
this is the comparison table's question asked again with N as a variable.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq

from . import budget as B
from . import config
from . import regulatory as R
from . import sources as S
from . import table as T
from . import wavesim as W
from .physconst import db
from .sources import derive
from .studies import _leak, _with_leakage_closed
from .vendor.mmw60.arraytx import beam as BM
from .vendor.mmw60.arraytx import geometry as G

KEY = "jcas60.array"


def half_power() -> float:
    """The half-power width of a uniform aperture's beam in units of lambda / aperture and in sin(theta): twice the x
    at which sin(pi x) / (pi x) has fallen to 1 / sqrt(2) -- the textbook's 0.886 (arraytx's width is the same level)."""
    return 2.0 * brentq(lambda x: np.sinc(x) - math.sqrt(0.5), 1e-9, 1.0 - 1e-9)


def geometry(s, n: int, element: str | None = None) -> G.ArrayGeometry:
    return G.LinearArray(int(n), S.num(KEY + "#spacing_lambda").value, s.carrier_hz.value, element or str(S.raw(KEY + "#element_pattern")))


def beamwidth_closed_deg(n: int, spacing_lambda: float, theta0_deg: float = 0.0) -> float:
    """Half-power width of a continuous aperture as long as ``n`` elements, steered to ``theta0``:
    0.886 lambda / (n d cos theta0), in degrees (small angles: the array's own width is wider at small n)."""
    return math.degrees(half_power() / (n * spacing_lambda * math.cos(math.radians(theta0_deg))))


def beams_closed(n: int, spacing_lambda: float, fov_deg: float) -> float:
    """How many beams, crossing at half power, cover +-fov: in sin(theta) every beam is 0.886 / (n d / lambda) wide."""
    return 2.0 * math.sin(math.radians(fov_deg)) * n * spacing_lambda / half_power()


def beam_half_width_u(geom: G.ArrayGeometry, mode: str) -> float:
    """Half the -3 dB width of the array's own beam in u = sin(theta), read on boresight from arraytx's array factor.
    With phase steering at the carrier the pattern depends on u - u0 alone, so this is the same wherever the beam points."""
    iso = G.LinearArray(geom.n, geom.spacing_lambda, geom.f0_hz, "isotropic")
    return math.sin(math.radians(0.5 * BM.beamwidth_3db_deg(iso, 0.0, mode)))


def beams_placed(geom: G.ArrayGeometry, fov_deg: float, mode: str) -> list[float]:
    """The beams' directions: evenly spaced in sin(theta) by one -3 dB width (so that neighbours cross at -3 dB),
    centred, as many as it takes to cover +-fov."""
    if geom.n < 2:
        return [0.0]
    du = beam_half_width_u(geom, mode)
    m = math.ceil(math.sin(math.radians(fov_deg)) / du - 1e-9)
    return [math.degrees(math.asin((2 * k - (m - 1)) * du)) for k in range(m)]


def squint_closed_db(s, geom: G.ArrayGeometry, theta0_deg: float, mode: str) -> float:
    """What the correlation peak loses in the beam's own direction to the steering's frequency response: the
    raised-cosine spectrum weighted by arraytx's ``beam_response`` (1 everywhere for true time delay)."""
    f = np.linspace(-1.0, 1.0, 4001)
    g = B.rc_spectrum(f, s.rolloff.value)
    h = BM.beam_response(geom, theta0_deg, f * s.chip_rate_hz.value, mode)
    return -2.0 * db(abs(np.trapezoid(g * h, f)) / np.trapezoid(g, f))


def squint_sim_db(s, sim: W.Simulator, geom: G.ArrayGeometry, theta0_deg: float, mode: str) -> float:
    """The same from the waveform: the sensing field driven element by element (arraytx's ``steer``), summed in the
    far field in the beam's direction (``far_field``), correlated; against the same array on boresight."""
    x = np.fft.ifft(sim.tx_f)
    j = -sim.lag_lo

    def peak(theta):
        e = BM.far_field(BM.steer(x, geom, theta, sim.fs, mode), geom, theta, sim.fs) / geom.element_gain(theta)
        return abs(sim.correlate(e)[j])
    return -2.0 * db(peak(theta0_deg) / peak(0.0))


def leakage(s, bench, rows, readings) -> dict:
    """The monostatic row's leakage with N power amplifiers beside the receiver.  Every element is taken to leak into
    the receive antenna at the isolation one antenna has; with phase steering the N leakages are one waveform N times,
    so they add as phasors whose phases this library has no model of: between nothing and all in phase (the table's
    own formula, EIRP - G_t - isolation), with the sum of their powers 10 log10(N) under that.  The criteria are the
    comparison table's (compression, ADC, LO phase noise), at the level each sum leaves; at one array and one reading
    both sums are simulated with all three mechanisms on."""
    fe, sim = s.fe, bench.sim
    env = sim.unit_envelope(s.leak_delay_s.value)
    label, iso = fe.isolation_db[-1]
    sums = (("sum of powers", 1.0), ("all in phase", 0.0))

    def effective(eirp, n, share):
        return derive(eirp.value - share * db(n), "EIRP - 10 log10(N) (the leakages' powers add)" if share else "EIRP (the leakages add in phase)", eirp, unit="dBm")
    out = []
    for r in rows:
        for p in r["readings"]:
            per, link = [], B.link(s, "monostatic", p["eirp_dbm"])
            for name, share in sums:
                cost = B.monostatic_cost(s, effective(p["eirp_dbm"], r["n"], share), env)
                leak = cost.conducted_dbm.value - iso.value
                per.append({"sum": name, "leak_dbm": leak, "need_db": cost.need_db.value, "binding": cost.binding, "shortfall_db": cost.need_db.value - iso.value,
                            "snr_closed_db": link.snr_grid_db.value + _with_leakage_closed(s, link, cost, env, leak)["snr_change_db"]})
            out.append({"n": r["n"], "limit": p["limit"], "eirp_dbm": p["eirp_dbm"].value, "snr_no_leak_db": link.snr_grid_db.value, "sums": per})
    nw, rd = int(S.num(KEY + "#simulated_elements").value), readings[0]
    eirp = next(p["eirp_dbm"] for r in rows if r["n"] == nw for p in r["readings"] if p["limit"] == rd.name)
    link = B.link(s, "monostatic", eirp)
    g, tgt, _ = bench.paths("monostatic", link, 0.0)
    near = int(round(g.delay_s * sim.fs))
    opt = bench.options("monostatic", adc_enob=s.adc_enob.value, lo="shared", rx_p1db_dbm=fe.p1db_inband_dbm.value)
    alone = bench.read([tgt], bench.options("monostatic"), near, strong=())
    roads = []
    for name, share in sums:
        cost = B.monostatic_cost(s, effective(eirp, nw, share), env)
        leak = cost.conducted_dbm.value - iso.value
        c = _with_leakage_closed(s, link, cost, env, leak)
        got = bench.read([tgt, _leak(s, leak)], opt, near)
        roads.append({"sum": name, "leak_dbm": leak, "snr_closed_db": link.snr_grid_db.value + c["snr_change_db"], "snr_sim_db": got["snr_db"],
                      "echo_db": c["echo_db"], "pn_over_thermal_db": c["pn_over_thermal_db"], "adc_over_thermal_db": c["adc_over_thermal_db"]})
    return {"rows": out, "roads": roads, "isolation_label": label, "isolation_db": iso.value, "n": nw, "limit": rd.name, "eirp_dbm": eirp.value,
            "snr_no_leak_closed_db": link.snr_grid_db.value, "snr_no_leak_sim_db": alone["snr_db"]}


def _above_rounding(loss_db: float) -> float:
    """A loss that is nothing but floating rounding (true time delay: the steering has no frequency response) is
    reported as zero, not as a reading: its last digits differ between numpy versions."""
    return loss_db if loss_db > 1e-9 else 0.0


@lru_cache(maxsize=1)
def build() -> dict:
    s = config.setup()
    fe = s.fe
    spacing, mode = S.num(KEY + "#spacing_lambda").value, str(S.raw(KEY + "#steering"))
    fov, readings, drop = S.num(KEY + "#field_of_view_deg").value, R.readings(s), S.num(KEY + "#crossover_db").value
    if drop != 3 or S.num(KEY + "#rx_elements").value != 1:
        raise ValueError("arraytx's beam width is the half-power one (3 dB), and the receiver here is one antenna")
    fields = (("L-CEF", s.field), ("TRN field", s.trn_field))
    rows = []
    for n in (int(x) for x in S.raw(KEY + "#tx_elements")):
        geom = geometry(s, n)
        full = derive(fe.eirp_max_dbm.value + 2.0 * db(n), "EIRP of one element + 20 log10(N)", fe.tx_power_dbm, fe.tx_gain_dbi, unit="dBm")
        beams = beams_placed(geom, fov, mode)
        per = []
        for rd in readings:
            eirp = derive(min(full.value, rd.cap_dbm.value), "min(EIRP of one element + 20 log10(N), cap)", full, rd.cap_dbm, unit="dBm")
            per.append({"limit": rd.name, "cap_dbm": rd.cap_dbm.value, "eirp_dbm": eirp, "backoff_db": full.value - eirp.value,
                        "fields": [{"field": name, "snr_cell_db": B.link(s, "monostatic", eirp, field=f).snr_cell_db.value,
                                    "air_one_beam": R.air_time(s, rd, f).fraction.value,
                                    "air_sweep": R.air_time(s, rd, f).fraction.value * len(beams)} for name, f in fields]})
        rows.append({"n": n, "eirp_full_dbm": full.value, "readings": per,
                     "width_closed_deg": beamwidth_closed_deg(n, spacing) if n > 1 else None,
                     "width_deg": BM.beamwidth_3db_deg(G.LinearArray(n, spacing, s.carrier_hz.value, "isotropic"), 0.0, mode) if n > 1 else None,
                     "beams_closed": beams_closed(n, spacing, fov) if n > 1 else 1.0,
                     "beams_array": math.sin(math.radians(fov)) / beam_half_width_u(geom, mode) if n > 1 else 1.0, "beams": len(beams), "beam_directions_deg": beams,
                     "sweeps_per_s_at_one_beams_air_time": s.frames_per_s.value / len(beams),
                     "scan_loss_edge_db": -2.0 * db(geom.element_gain(fov))})
    # the squint, on the widest array: phase steering against true time delay, both roads
    ns = int(S.num(KEY + "#squint_elements").value)
    gs = geometry(s, ns)
    bench = T._Bench(s)
    squint = [{"field": name, "n": ns, "theta_deg": th, "mode": m, "closed_db": _above_rounding(squint_closed_db(s, gs, th, m)),
               "sim_db": _above_rounding(squint_sim_db(s, bench.sim_for(f), gs, th, m))}
              for name, f in fields for th in (0.5 * fov, fov) for m in (mode, "ttd")]
    edge_hz = 0.5 * (1.0 + s.rolloff.value) * s.chip_rate_hz.value
    # the SNR on the second road, at the array that reaches the communication limit (no leakage: the link alone)
    nw = int(S.num(KEY + "#simulated_elements").value)
    rd = readings[0]
    eirp = next(p["eirp_dbm"] for r in rows if r["n"] == nw for p in r["readings"] if p["limit"] == rd.name)
    roads = []
    for name, f in fields:
        sim = bench.sim_for(f)
        link = B.link(s, "monostatic", eirp, field=f)
        g, tgt, _ = bench.paths("monostatic", link, 0.0)
        got = bench.read([tgt], bench.options("monostatic"), int(round(g.delay_s * sim.fs)), strong=(), sim=sim)
        roads.append({"field": name, "n": nw, "limit": rd.name, "eirp_dbm": eirp.value, "snr_grid_db": link.snr_grid_db.value, "snr_sim_db": got["snr_db"]})
    return {"setup": s, "rows": rows, "squint": squint, "squint_band_edge_deg": math.degrees(math.tan(math.radians(fov)) * edge_hz / s.carrier_hz.value),
            "roads": roads, "leakage": leakage(s, bench, rows, readings), "fov_deg": fov, "air_limit": readings[0].name, "mode": mode, "spacing_lambda": spacing, "element": str(S.raw(KEY + "#element_pattern")),
            "one_element_eirp_dbm": fe.eirp_max_dbm.value, "frames_per_s": s.frames_per_s.value}


def render(t: dict | None = None) -> str:
    from . import ENGINE_BANNER, WORDING
    t = t or build()
    L: list[str] = []
    w = L.append
    w(ENGINE_BANNER)
    w(WORDING)
    w("")
    w(f"== 1. a transmit array of N elements (line, {t['spacing_lambda']:g} lambda apart, {t['element']} element pattern, {t['mode']} steering; one receive antenna; "
      f"one element radiates {t['one_element_eirp_dbm']:g} dBm EIRP at full power) ==")
    lim = [p["limit"] for p in t["rows"][0]["readings"]]
    w(f"{'N':>4s} {'full EIRP':>11s} " + " ".join(f"{name + ': EIRP, back-off, SNR L-CEF / TRN':>58s}" for name in lim))
    for r in t["rows"]:
        cells = " ".join(f"{p['eirp_dbm'].value:22.2f} dBm {p['backoff_db']:6.2f} dB {p['fields'][0]['snr_cell_db']:8.2f} / {p['fields'][1]['snr_cell_db']:.2f} dB" for p in r["readings"])
        w(f"{r['n']:4d} {r['eirp_full_dbm']:7.2f} dBm {cells}")
    w("  SNR(cell) of row B at 5 m, closed form; the receiver is one antenna, so the array's gain is on the way out only")
    w("")
    w(f"== 2. the beam, and how many of them cover +-{t['fov_deg']:g} deg (adjacent beams crossing at -3 dB) ==")
    w(f"{'N':>4s} {'-3 dB width on boresight: closed':>34s} {'arraytx':>10s} {'beams: closed':>15s} {'arraytx':>9s} {'placed':>8s} {'loss at the edge of view':>26s} "
      f"{'air time of a sweep, ' + format(t['frames_per_s'], 'g') + ' sweeps/s: L-CEF / TRN':>52s} {'sweeps/s at one beam’s air time':>34s}")
    for r in t["rows"]:
        a = r["readings"][0]["fields"]
        w(f"{r['n']:4d} " + (f"{r['width_closed_deg']:30.2f} deg {r['width_deg']:6.2f} deg" if r["width_deg"] is not None else f"{'—':>34s} {'—':>10s}")
          + f" {r['beams_closed']:15.2f} {r['beams_array']:9.2f} {r['beams']:8d} {r['scan_loss_edge_db']:23.2f} dB {a[0]['air_sweep'] * 100:41.2f} / {a[1]['air_sweep'] * 100:.2f} % "
          f"{r['sweeps_per_s_at_one_beams_air_time']:28.2f}")
    w("  closed: 0.886 lambda / (N d) wide, 2 sin(fov) N d / (0.886 lambda) of them; arraytx: the same from its array factor's own width; placed: evenly in sin(theta), one width apart, centred")
    w(f"  the loss at the edge of view is the element's pattern, one way; a sweep sends one sensing frame per beam (air time under the {t['air_limit']})")
    w("  closed is the continuous aperture's width: the array's own is wider at small N (two elements: 60 deg, not 50.8)")
    w("")
    sq = t["squint"]
    w(f"== 3. phase steering over the occupied band: the squint (N = {sq[0]['n']}; the band's edge points {t['squint_band_edge_deg']:.2f} deg off at the edge of view) ==")
    w(f"{'field':>12s} {'steered to':>12s} {'steering':>10s} {'peak lost: closed':>20s} {'waveform':>10s}")
    for q in sq:
        w(f"{q['field']:>12s} {q['theta_deg']:8.1f} deg {q['mode']:>10s} {q['closed_db']:17.4f} dB {q['sim_db']:7.4f} dB")
    w("  closed: the raised-cosine spectrum weighted by the steering's response in the beam's direction; waveform: the field driven element by element and summed in the far field")
    w("")
    w("== 4. the link on the two roads at the array that reaches the communication limit (no leakage) ==")
    for q in t["roads"]:
        w(f"  {q['field']:10s} N = {q['n']}, {q['limit']}, EIRP {q['eirp_dbm']:.2f} dBm: closed on the grid {q['snr_grid_db']:.2f} dB, simulated {q['snr_sim_db']:.2f} dB, "
          f"difference {q['snr_sim_db'] - q['snr_grid_db']:+.2f} dB")
    w("")
    k = t["leakage"]
    w(f"== 5. the monostatic row: N power amplifiers beside the receiver (each element leaks at the isolation one antenna has; {k['isolation_label']}: {k['isolation_db']:g} dB) ==")
    w(f"{'N':>4s} {'limit':>22s} {'EIRP':>10s} {'SNR, no leakage':>16s} " + " ".join(f"{name + ': leakage, isolation needed, short by, SNR left':>62s}" for name in (q["sum"] for q in k["rows"][0]["sums"])))
    for r in k["rows"]:
        w(f"{r['n']:4d} {r['limit']:>22s} {r['eirp_dbm']:6.2f} dBm {r['snr_no_leak_db']:13.2f} dB " + " ".join(
            f"{q['leak_dbm']:22.2f} dBm {q['need_db']:8.2f} dB {q['shortfall_db']:+7.2f} dB {q['snr_closed_db']:8.2f} dB" for q in r["sums"]))
    w(f"  row B (L-CEF) at 5 m, SNR on the grid, closed form with the limiter, the LO's phase noise and the ADC at once; what binds the isolation everywhere: {', '.join(sorted({q['binding'] for r in k['rows'] for q in r['sums']}))}")
    w("  the N leakages are one waveform N times: they add as phasors, anywhere from nothing to all in phase, and the phases move with the beam; this library has no model of them")
    w(f"  both roads at N = {k['n']}, {k['limit']}, EIRP {k['eirp_dbm']:.2f} dBm, with the limiter, the LO's phase noise and the ADC on "
      f"(no leakage: closed {k['snr_no_leak_closed_db']:.2f} dB, simulated {k['snr_no_leak_sim_db']:.2f} dB):")
    for q in k["roads"]:
        w(f"    {q['sum']:14s} leakage {q['leak_dbm']:7.2f} dBm: closed {q['snr_closed_db']:.2f} dB, simulated {q['snr_sim_db']:.2f} dB, difference {q['snr_sim_db'] - q['snr_closed_db']:+.2f} dB "
          f"(closed: echo {q['echo_db']:+.2f} dB, phase noise {q['pn_over_thermal_db']:+.2f} dB and ADC {q['adc_over_thermal_db']:+.2f} dB against the thermal floor)")
    return "\n".join(L) + "\n"
