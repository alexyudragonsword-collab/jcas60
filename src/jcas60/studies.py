"""Studies beside the comparison table: what the table raises and does not answer (the ROADMAP's list, in its order).

Nothing here changes a number of the table.  Each study takes the same
:class:`~jcas60.config.Setup`, varies one thing, and answers by two roads
where it draws a conclusion (law 4):

1. :func:`limiter` -- the receiver's compression is an analogue thing and the
   table computes it at the ADC's two samples a chip.  The same criterion with
   everything ahead of the ADC at 4 and 8 samples a chip (closed form: the
   envelope sampled as densely; simulation: the limiter at that rate, an ideal
   low-pass, then the ADC), and with a softer and a harder knee.
2. :func:`gain_settings` -- the receiver's gain setting as a variable: each
   setting of the radar chip's figure has its own noise figure and in-band
   compression point, and what the leakage costs is read as the SNR that is
   left, not as three thresholds relative to a floor that moves with them.
3. :func:`beyond_window` -- a target farther than the clean window: where in
   the map it lands and how strong (chip-domain correlation of the field
   against a copy of itself; the waveform simulation at the same delays), and
   the two roads' difference across the window up to its edge.
4. :func:`fast_targets` -- what a fast target loses inside one packet, as a
   closed form next to the simulation.
5. :func:`undithered_adc` -- the ADC at the low isolation, where the thermal
   noise is a small fraction of an LSB: the ridge and the peak's wobble as a
   statistic over the target's delay rather than one cell's reading, and
   what a dither buys and costs.
6. :func:`leakage_paths` -- the leakage's delay as a variable, and a near
   reflector as a second path whose phase-noise residue adds in amplitude.
7. :func:`same_lo` -- every row on one LO: fmcw_radar's own report of its
   x10-LO config next to its default's, and the ay rows' phase-noise terms
   with that LO's profile in the place of the model.
8. :func:`bistatic_sync` -- what the two devices' leftover carrier offset
   costs inside a packet, and the direct path's phase noise on the TRN field.
9. :func:`range_window` -- the range response with the window row A's range
   FFT has: width, sidelobes, SNR, on both roads.
10. :func:`second_target` -- what a strong target (or a static reflector)
    leaves at another target's place: in range, with and without that
    window; in Doppler; and through the LO's phase noise.
11. :func:`designed_waveform` -- a waveform made for sensing in the TRN
    field's place (the entries allow one, on four conditions): every block
    correlated, the clean zone twice as long, against the standard subfield.
12. :func:`measured_knee` -- the limiter's knee had no source: two measured
    60 GHz compression curves fitted with the same Rapp model, and what each
    knee makes of the compression criterion and the ghost.
13. :func:`long_trn` -- the standard's long TRN subfield (TRN_BL = 256): the
    same field length and SNR, twice the clean zone.
14. :func:`bistatic_exchange` -- row C laid on the exchange the 802.11bf
    entry has (BRP frames, SIFS, BRPIFS, a reply): closed form only.
"""
from __future__ import annotations

import dataclasses
import math
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from . import budget as B
from . import config
from . import field as F
from . import fmcw_row as A
from . import frontend as FE
from . import regulatory as R
from . import sources as S
from . import table as T
from . import wavesim as W
from .physconst import C0, db, undb
from .sources import Num
from .vendor.fmcw_radar.link.compression import rapp

KEY = "jcas60.studies"
FIG = "rx60.radar_iwr6843_fig7_2"


def _bench(s) -> T._Bench:
    bench = T._Bench(s)
    bench.frames = int(S.num(KEY + "#frames").value)
    return bench


def _leak(s, level_dbm: float) -> W.Path:
    return W.Path(float(level_dbm), s.leak_delay_s.value, 0.0, "leakage")


# ----------------------------------------------------------------------------- 1. the limiter
def _limiter_roads(s, bench, analog: int, knee: float, levels=None) -> dict:
    """The compression criterion and the image at the best isolation found, for one analogue rate and one knee."""
    sp = dataclasses.replace(s, rx_rapp_p=Num(float(knee), "assumed", KEY + "#knees", "-"))
    sim = W.Simulator(sp, analog_osr=analog)
    eirp = R.readings(sp)[0].eirp_dbm
    link = B.link(sp, "monostatic", eirp)
    g, tgt, _ = bench.paths("monostatic", link, 0.0)
    near = int(round(g.delay_s * sim.fs))
    cond = eirp.value - sp.fe.tx_gain_dbi.value
    p1, loss = sp.fe.p1db_inband_dbm.value, sp.desense_allowed_db.value
    env = sim.unit_envelope(sp.leak_delay_s.value)
    level_closed = B.strong_level_for_snr_loss(env, sp, loss)
    x = np.array([level_closed + float(o) for o in S.raw(KEY + "#level_offsets_db")])
    opt = bench.options("monostatic", rx_p1db_dbm=p1)
    lin = bench.read([tgt, _leak(sp, x[0])], bench.options("monostatic"), near, sim=sim)
    y = np.array([bench.read([tgt, _leak(sp, v)], opt, near, sim=sim)["snr_db"] - lin["snr_db"] for v in x])
    level_sim = T._crossing(x, y, -loss)
    iso = sp.fe.isolation_db[-1][1].value
    c = B.compression(env, cond - iso, sp)

    def image(level):
        q = sim.read(sim.frame([tgt, _leak(sp, level)], bench.options("monostatic", rx_p1db_dbm=p1, noise=False), 0), near, (0,))
        return db(q["image_mw"] / q["peak_mw"]), db(q["peak_mw"])

    out = {"analog": analog, "knee": float(knee), "inverse_envelope_mean": float(np.mean(1.0 / env)),
           "need_closed_db": cond - level_closed, "need_sim_db": None if level_sim is None else cond - level_sim,
           "limit_closed_db": B.compression(env, p1 + 60.0, sp)["snr"], "isolation_db": iso,
           "image_bound_dbc": c["image"] - c["signal"], "image_sim_dbc": image(cond - iso)[0]}
    if levels is not None:                                   # the image against the level: where it meets the map's floor
        floor_dbc = -link.snr_grid_db.value
        at = B.strong_level_for_image(env, sp, floor_dbc)
        got = T._crossing(np.asarray(levels), np.array([image(v)[0] for v in levels]), floor_dbc)
        out.update(floor_dbc=floor_dbc, ghost_closed_db=None if at is None else cond - at, ghost_sim_db=None if got is None else cond - got)
    return out


def limiter(s, bench) -> dict:
    rates = [int(a) for a in S.raw(KEY + "#analog_samples_per_chip")]
    knees = [float(p) for p in S.raw(KEY + "#knees")]
    rate = int(S.num(KEY + "#knee_rate").value)
    p1 = s.fe.p1db_inband_dbm.value
    span, step = float(s.sim["compression_span_db"]), float(s.sim["compression_step_db"])
    levels = p1 + np.arange(-span, span + step / 2, step)
    return {"by_rate": [_limiter_roads(s, bench, a, s.rx_rapp_p.value) for a in rates],
            "by_knee": [_limiter_roads(s, bench, rate, p, levels) for p in knees], "knee_rate": rate}


# ----------------------------------------------------------------------------- 2. the gain setting
def _with_leakage_closed(sg, link: B.Link, cost: B.Monostatic, env: np.ndarray, leak_dbm: float) -> dict:
    """The SNR change with the leakage at ``leak_dbm`` and all three mechanisms at once.  The limiter comes first:
    it takes the echo and the thermal noise down by its small-signal gains; the phase-noise residue rides on the
    leakage, in quadrature to it, and leaves with the leakage's own gain; the ADC's range follows the leakage as it
    leaves the limiter.  Powers add at the correlator."""
    c = B.compression(env, leak_dbm, sg)
    n0 = link.noise_chip_dbm.value
    thermal = undb(c["noise"])
    pn = undb(leak_dbm + c["strong"] + cost.pn_residual_dbc.value - n0)
    adc = undb(B.adc_noise_inband_dbm(leak_dbm + c["strong"], sg.adc_enob.value, sg) - n0)
    return {"echo_db": c["signal"], "thermal_db": c["noise"], "pn_over_thermal_db": db(pn / thermal), "adc_over_thermal_db": db(adc / thermal),
            "snr_change_db": c["signal"] - db(thermal + pn + adc)}


def gain_settings(s, bench) -> dict:
    eirp = R.readings(s)[0].eirp_dbm
    rows = []
    for gdb in S.raw(FIG + "#gain_db"):
        at = f"@gain_db={float(gdb):g}"
        fe = dataclasses.replace(s.fe, nf_db=S.num(FIG + "#nf_db" + at, "dB"), rx_gain_db=S.num(FIG + "#gain_db" + at, "dB"),
                                 p1db_inband_dbm=S.num(FIG + "#p1db_inband_dbm" + at, "dBm"))
        sg = dataclasses.replace(s, fe=fe)
        sim = W.Simulator(sg)
        env = sim.unit_envelope(sg.leak_delay_s.value)
        link = B.link(sg, "monostatic", eirp)
        cost = B.monostatic_cost(sg, eirp, env)
        g, tgt, _ = bench.paths("monostatic", link, 0.0)
        near = int(round(g.delay_s * sim.fs))
        have = []
        for label, iso, _short, leak, _enob in cost.have:
            c = _with_leakage_closed(sg, link, cost, env, leak.value)
            opt = bench.options("monostatic", adc_enob=sg.adc_enob.value, lo="shared", rx_p1db_dbm=fe.p1db_inband_dbm.value)
            got = bench.read([tgt, _leak(sg, leak.value)], opt, near, sim=sim)
            have.append({"label": label, "isolation_db": iso.value, "leak_dbm": leak.value, **c,
                         "snr_closed_db": link.snr_grid_db.value + c["snr_change_db"], "snr_sim_db": got["snr_db"]})
        rows.append({"gain_db": float(gdb), "nf_db": fe.nf_db.value, "p1db_dbm": fe.p1db_inband_dbm.value, "snr_no_leak_db": link.snr_grid_db.value,
                     "need_compression_db": cost.need_compression_db.value, "need_pn_db": cost.need_pn_db.value, "need_adc_db": cost.need_adc_db.value,
                     "need_db": cost.need_db.value, "binding": cost.binding, "have": have})
    best = [max(range(len(rows)), key=lambda k: rows[k]["have"][j]["snr_closed_db"]) for j in range(len(rows[0]["have"]))]
    return {"rows": rows, "best": best, "table_nf_db": s.fe.nf_db.value, "table_gain_db": s.fe.rx_gain_db.value,
            "snr_no_leak_max_db": max(r["snr_no_leak_db"] for r in rows)}


# ----------------------------------------------------------------------------- 3. beyond the clean window
def _chip_profile(field, reach: int) -> np.ndarray:
    """|C(k)|^2 for k = -reach .. 0 (chips), C(k) = the summed window correlations of the field against a copy of
    itself that arrives -k chips later, nothing ahead of the field and nothing behind it, normalised to the peak:
    what a target k chips before the lag being read puts there.  Exact, in chips."""
    n = 1 << int(math.ceil(math.log2(field.chips.size + reach + field.zero_sidelobe_chips + 2)))
    rx = np.zeros(n, dtype=complex)
    rx[: field.chips.size] = field.chips
    tmpl = np.zeros(n, dtype=complex)
    for w0, ref in field.windows:
        tmpl[w0: w0 + len(ref)] = ref
    c = np.fft.ifft(np.fft.fft(rx) * np.conj(np.fft.fft(tmpl))) / field.n_correlated       # c[k mod n] = C(k)
    return np.abs(c[np.arange(-reach, 1) % n]) ** 2


def _ghost(profile: np.ndarray, delay: int, far: int) -> tuple[float, int]:
    """(level relative to an in-window target of the same strength [linear], the lag it shows at) for a target
    ``delay`` chips away, read on lags 0 .. ``far``."""
    reach = profile.size - 1
    lags = np.arange(0, far + 1)
    k = lags - delay                                             # C's argument at each lag
    v = np.where(k >= -reach, profile[np.clip(k + reach, 0, reach)], 0.0)
    j = int(np.argmax(v))
    return float(v[j]), int(lags[j])


def _beyond_rows(s, prof, sim, delays, far) -> list:
    """A target at each of ``delays`` chips: the strongest it shows on lags 0 .. ``far``, by the chips and simulated."""
    chip0 = (-sim.lag_lo) % sim.osr                                  # the lag index of the chip-spaced lags
    floor = undb(sim.rounding_dbc)                                   # under it a level is the arithmetic's rounding, not a reading
    rows = []
    for d in delays:
        closed, lag = _ghost(prof, d, far)
        cir = sim.frame([W.Path(0.0, d / s.chip_rate_hz.value, 0.0)], W.Options(noise=False, n_packets=1)).cir[0]
        p = np.abs(cir[chip0:: sim.osr]) ** 2
        lags = sim.lags[chip0:: sim.osr] // sim.osr
        keep = (lags >= 0) & (lags <= far)
        j = int(np.argmax(np.where(keep, p, 0.0)))
        rows.append({"delay_chips": d, "range_m": d * s.chip_m, "level_closed": max(closed, floor), "lag_closed": lag if closed > floor else None,
                     "level_sim": max(float(p[j]), floor), "lag_sim": int(lags[j]) if p[j] > floor else None})
    return rows


def beyond_window(s, bench) -> dict:
    delays = [int(d) for d in S.raw(KEY + "#beyond_delays_chips")]
    margin = int(S.num(KEY + "#floor_margin_chips").value)
    long_buf = dataclasses.replace(s, sim={**s.sim, "buffer_chips": int(S.num(KEY + "#beyond_buffer_chips").value)})
    eirp = R.readings(s)[0].eirp_dbm
    out = {"fields": [], "ranges": []}
    for name, field in (("L-CEF", s.field), ("TRN field", s.trn_field)):
        far = field.zero_sidelobe_chips - margin
        prof = _chip_profile(field, max(delays) + field.zero_sidelobe_chips)
        sim = W.Simulator(long_buf, field)
        rows = _beyond_rows(s, prof, sim, delays, far)
        # the worst a target anywhere past the window does, chip by chip (the chip road alone: it is exact)
        worst = max(((_ghost(prof, d, far)[0], d) for d in range(field.zero_sidelobe_chips + 1, max(delays) + 1)))
        out["fields"].append({"name": name, "window_chips": field.zero_sidelobe_chips, "window_m": field.zero_sidelobe_chips * s.chip_m,
                              "read_to_chips": far, "rows": rows, "worst_level": worst[0], "worst_delay_chips": worst[1]})
        # inside the window, out to its edge: the two roads at several ranges
        sim_w = bench.sim_for(field)
        for r in S.raw(KEY + "#window_ranges_m"):
            link = B.link(s, "monostatic", eirp, float(r), field)
            g = B.geometry(s, "monostatic", float(r))
            tgt = W.Path(link.p_rx_dbm.value, g.delay_s, g.path_rate(s.velocity_mps.value), "target")
            strong = T.strong_path_dbm(s, "monostatic", eirp, s.fe.isolation_db[-1][1])
            got = bench.read([tgt, _leak(s, strong.value)], bench.options("monostatic"), int(round(g.delay_s * sim_w.fs)), sim=sim_w)
            out["ranges"].append({"field": name, "range_m": float(r), "delay_chips": g.delay_s * s.chip_rate_hz.value,
                                  "snr_closed_db": link.snr_grid_db.value, "snr_sim_db": got["snr_db"]})
    return out


# ----------------------------------------------------------------------------- 4. fast targets
def fast_targets(s, bench) -> dict:
    v_max = s.wavelength_m / (4.0 * s.interval_s.value)
    vs = [s.velocity_mps.value, v_max] + [float(v) for v in S.raw(KEY + "#fast_velocities_mps")]
    g = B.geometry(s, "monostatic")
    delay = round(g.delay_s * s.chip_rate_hz.value) / s.chip_rate_hz.value       # on a chip: the peak is one lag's reading
    out = []
    for name, field in (("L-CEF", s.field), ("TRN field", s.trn_field)):
        sim = bench.sim_for(field)
        j = int(round(delay * sim.fs)) - sim.lag_lo

        def peak(v, sim=sim, j=j):
            return abs(sim.frame([W.Path(0.0, delay, 2.0 * v)], W.Options(noise=False, n_packets=1)).cir[0][j]) ** 2
        at_rest = peak(0.0)
        span = (field.windows[-1][0] + len(field.windows[-1][1]) - field.windows[0][0]) / s.chip_rate_hz.value
        out.append({"name": name, "span_s": span,
                    "rows": [{"velocity_mps": v, "loss_closed_db": -B.inpacket_loss_db(s, field, 2.0 * v), "loss_sim_db": -db(peak(v) / at_rest)} for v in vs]})
    return {"fields": out, "v_unambiguous_mps": v_max}


# ----------------------------------------------------------------------------- 5. the ADC that nothing dithers
def _spread(values) -> dict:
    a = np.asarray(values, dtype=float)
    return {"median": float(np.median(a)), "min": float(a.min()), "max": float(a.max()), "rms": float(np.sqrt(np.mean(a ** 2))), "mean": float(a.mean())}


def undithered_adc(s, bench) -> dict:
    eirp = R.readings(s)[0].eirp_dbm
    link = B.link(s, "monostatic", eirp)
    sim = bench.sim
    cost = B.monostatic_cost(s, eirp, sim.unit_envelope(s.leak_delay_s.value))
    _label, iso, _short, leak, _enob = cost.have[0]
    g = B.geometry(s, "monostatic")
    n, step, enob = int(S.num(KEY + "#adc_delays").value), S.num(KEY + "#adc_delay_step_chips").value, s.adc_enob.value
    n0 = link.noise_chip_dbm.value
    lin: dict = {}

    def at(k, dither=None):
        tau = g.delay_s + (k - n // 2) * step / s.chip_rate_hz.value
        paths = [W.Path(link.p_rx_dbm.value, tau, g.path_rate(s.velocity_mps.value), "target"), _leak(s, leak.value)]
        near = int(round(tau * sim.fs))
        if k not in lin:
            lin[k] = bench.read(paths, bench.options("monostatic"), near, sim=sim)
        adc = bench.read(paths, bench.options("monostatic", adc_enob=enob, adc_dither_lsb=dither), near, sim=sim)
        return {"ridge_under_target_db": adc["ridge_dbm"] - lin[k]["peak_dbm"], "peak_change_db": adc["peak_dbm"] - lin[k]["peak_dbm"],
                "rise_db": adc["floor_dbm"] - lin[k]["floor_dbm"], "snr_change_db": adc["snr_db"] - lin[k]["snr_db"]}

    def gather(ks, dither=None):
        rows = [at(k, dither) for k in ks]
        q = sim.frame([W.Path(link.p_rx_dbm.value, g.delay_s, g.path_rate(s.velocity_mps.value), "target"), _leak(s, leak.value)],
                      bench.options("monostatic", adc_enob=enob, adc_dither_lsb=dither), 0).quant
        d = 0.0 if dither is None else float(dither)
        # the quantiser's own noise is 1/12 LSB^2 a rail; a dither of d LSB rms adds d^2 to it
        closed = B.floor_rise_db(B.adc_noise_inband_dbm(leak.value, enob, s) + db(1.0 + 3.0 * 4.0 * d * d), n0)
        return {"dither_lsb": d, "delays": len(rows), "repeating_share": q["repeating_share"], "noise_lsb": q["noise_lsb"], "rise_closed_db": closed,
                **{k: _spread([r[k] for r in rows]) for k in rows[0]}}
    return {"isolation_db": iso.value, "leak_dbm": leak.value, "step_chips": step,
            "cases": [gather(range(n))] + [gather(range(0, n, 2), float(d)) for d in S.raw(KEY + "#adc_dither_lsb")]}


# ----------------------------------------------------------------------------- 6. the leakage's delay, and a second path
def leakage_paths(s, bench) -> dict:
    eirp = R.readings(s)[0].eirp_dbm
    link = B.link(s, "monostatic", eirp)
    sim = bench.sim
    g, tgt, _ = bench.paths("monostatic", link, 0.0)
    near = int(round(g.delay_s * sim.fs))
    iso = s.fe.isolation_db[-1][1].value
    level = eirp.value - s.fe.tx_gain_dbi.value - iso
    n0 = link.noise_chip_dbm.value

    def rise_sim(paths):
        lin = bench.read([tgt] + paths, bench.options("monostatic"), near, sim=sim)
        return bench.read([tgt] + paths, bench.options("monostatic", lo="shared"), near, sim=sim)["floor_dbm"] - lin["floor_dbm"]
    rows = []
    for tau in (float(x) for x in S.raw(KEY + "#leak_delays_s")):
        st = dataclasses.replace(s, leak_delay_s=Num(tau, "assumed", KEY + "#leak_delays_s", "s"))
        cost = B.monostatic_cost(st, eirp, sim.unit_envelope(tau))
        rows.append({"delay_s": tau, "residual_dbc": cost.pn_residual_dbc.value, "need_pn_db": cost.need_pn_db.value,
                     "rise_closed_db": B.floor_rise_db(level + cost.pn_residual_dbc.value, n0), "rise_sim_db": rise_sim([W.Path(level, tau, 0.0, "leakage")])})
    # a near reflector: a second path, weaker and later, off the same LO
    extra, rel = S.num(KEY + "#reflector_range_m").value, S.num(KEY + "#reflector_level_db").value
    tau1 = s.leak_delay_s.value
    tau2 = tau1 + 2.0 * extra / C0
    a2 = math.sqrt(undb(rel))
    phase = lambda tau: -2.0 * math.pi * s.carrier_hz.value * tau              # noqa: E731
    alone = B.pn_residual_shared(s, tau1)
    as_built = B.pn_residual_paths(s, ((1.0, phase(tau1), tau1), (a2, phase(tau2), tau2)))
    span = [B.pn_residual_paths(s, ((1.0, 0.0, tau1), (a2, th, tau2))) for th in (0.0, math.pi)]
    average = alone + a2 * a2 * B.pn_residual_shared(s, tau2)
    base = B.monostatic_cost(s, eirp, sim.unit_envelope(tau1)).need_pn_db.value
    both = [W.Path(level, tau1, 0.0, "leakage"), W.Path(level + rel, tau2, 0.0, "reflector")]
    return {"isolation_db": iso, "leak_dbm": level, "rows": rows,
            "reflector": {"range_m": extra, "level_db": rel, "delay_s": tau2, "turns": (tau2 - tau1) * s.carrier_hz.value,
                          "need_pn_alone_db": base, "need_pn_db": base + db(as_built / alone), "need_pn_average_db": base + db(average / alone),
                          "need_pn_min_db": base + db(min(span) / alone), "need_pn_max_db": base + db(max(span) / alone),
                          "rise_closed_db": B.floor_rise_db(level + db(as_built), n0), "rise_sim_db": rise_sim(both),
                          "rise_alone_closed_db": B.floor_rise_db(level + db(alone), n0)}}


# ----------------------------------------------------------------------------- 7. one LO for every row
def _bistatic_pn(sv, bench, field, sim, frames: int | None = None) -> dict:
    """The direct path's phase noise in the bistatic rows: two LOs, each packet turned back by the phase read on the
    direct path over the correlated field.  ``frames``: how many frames the floor is read on (the TRN field's phase
    noise is slow to synthesise, and a floor is the mean of thousands of cells: one frame reads it)."""
    keep, bench.frames = bench.frames, frames or bench.frames
    try:
        return _bistatic_pn_read(sv, bench, field, sim)
    finally:
        bench.frames = keep


def _bistatic_pn_read(sv, bench, field, sim) -> dict:
    eirp = R.readings(sv)[0].eirp_dbm
    link = B.link(sv, "bistatic", eirp, field=field)
    direct = T.strong_path_dbm(sv, "bistatic", eirp)
    g, tgt, st = bench.paths("bistatic", link, direct.value)
    near = int(round(g.delay_s * sim.fs))
    lin = bench.read([tgt, st], bench.options("bistatic"), near, sim=sim)
    pn = bench.read([tgt, st], bench.options("bistatic", lo="independent"), near, sim=sim)
    span = (field.windows[-1][0] + len(field.windows[-1][1]) - field.windows[0][0]) / sv.chip_rate_hz.value
    unit = len(field.windows[0][1]) / sv.chip_rate_hz.value           # one correlation window: a complementary unit
    res, whole = B.pn_residual_independent(sv, unit), B.pn_residual_independent(sv, span)
    return {"reference_s": unit, "span_s": span, "residual_dbc": db(res), "rise_whole_field_closed_db": B.floor_rise_db(direct.value + db(whole), link.noise_chip_dbm.value), "rise_closed_db": B.floor_rise_db(direct.value + db(res), link.noise_chip_dbm.value),
            "rise_sim_db": pn["floor_dbm"] - lin["floor_dbm"], "snr_change_sim_db": pn["snr_db"] - lin["snr_db"]}


def same_lo(s, bench) -> dict:
    sx = dataclasses.replace(s, fe=FE.with_x10_lo(s.fe))
    eirp = R.readings(s)[0].eirp_dbm
    rep = {"default": A.report(), "x10lo": A.variant("x10lo")["report"]}
    a = {k: {"lo_source": r["radar"]["lo_source"], "snr_cell_db": r["detection"]["snr_cell_db"], "lo_rms_rad": r["phase_noise"]["lo_rms_rad"],
             "pn_mean_cell_dbc": r["phase_noise"]["mean_cell_dbc"], "pn_local_cell_dbc": r["phase_noise"]["local_cell_dbc"]} for k, r in rep.items()}
    iso = s.fe.isolation_db[-1][1].value
    rows = []
    for name, sv in (("the model the table's rows use", s), ("mmw60_trx's x10 LO", sx)):
        sim = W.Simulator(sv)
        cost = B.monostatic_cost(sv, eirp, sim.unit_envelope(sv.leak_delay_s.value))
        link = B.link(sv, "monostatic", eirp)
        g, tgt, _ = bench.paths("monostatic", link, 0.0)
        near = int(round(g.delay_s * sim.fs))
        level = cost.conducted_dbm.value - iso
        lin = bench.read([tgt, _leak(sv, level)], bench.options("monostatic"), near, sim=sim)
        pn = bench.read([tgt, _leak(sv, level)], bench.options("monostatic", lo="shared"), near, sim=sim)
        rows.append({"lo": name, "at_1mhz_dbchz": float(sv.fe.pn_dbchz(1e6)), "at_100mhz_dbchz": float(sv.fe.pn_dbchz(1e8)),
                     "residual_dbc": cost.pn_residual_dbc.value, "need_pn_db": cost.need_pn_db.value, "need_db": cost.need_db.value, "binding": cost.binding,
                     "rise_closed_db": B.floor_rise_db(level + cost.pn_residual_dbc.value, link.noise_chip_dbm.value), "rise_sim_db": pn["floor_dbm"] - lin["floor_dbm"],
                     "bistatic": [{"field": fname, **_bistatic_pn(sv, bench, field, W.Simulator(sv, field), frames)}
                                  for fname, field, frames in (("L-CEF", sv.field, None), ("TRN field", sv.trn_field, 1))]})
    return {"row_a": a, "isolation_db": iso, "rows": rows}


# ----------------------------------------------------------------------------- 8. the bistatic rows' synchronisation
def bistatic_sync(s, bench) -> dict:
    eirp = R.readings(s)[0].eirp_dbm
    out = []
    for name, field in (("L-CEF", s.field), ("TRN field", s.trn_field)):
        sim = bench.sim_for(field)
        link = B.link(s, "bistatic", eirp, field=field)
        direct = T.strong_path_dbm(s, "bistatic", eirp)
        g, tgt, st = bench.paths("bistatic", link, direct.value)
        near = int(round(g.delay_s * sim.fs))

        def peak(cfo, sim=sim, tgt=tgt, st=st, near=near):
            return sim.read(sim.frame([tgt, st], bench.options("bistatic", noise=False, cfo_hz=cfo), 0), near, (0,))["peak_mw"]
        clean = peak(0.0)
        rate = tgt.path_rate_mps                                    # the target's own Doppler is in the same sum: the offset adds to it
        rows = [{"cfo_hz": float(c), "turns": float(c) * field.chips.size / s.chip_rate_hz.value,
                 "loss_closed_db": B.inpacket_loss_db(s, field, rate) - B.inpacket_loss_db(s, field, rate - float(c) * s.wavelength_m),
                 "loss_sim_db": -db(peak(float(c)) / clean)}
                for c in S.raw(KEY + "#cfo_hz")]
        out.append({"name": name, "rows": rows})
    return {"fields": out, "equivalent_mps_per_hz": s.wavelength_m / 2.0}


# ----------------------------------------------------------------------------- 9. a range window
def _response_sim(s, sim, offsets) -> np.ndarray:
    """The range response at ``offsets`` chips from a target, normalised to its peak: one lag's reading as a
    noise-free target is moved the other way (the lag grid is half a chip; the target is not bound to it)."""
    g = B.geometry(s, "monostatic")
    chips = round(g.delay_s * s.chip_rate_hz.value)
    j = chips * sim.osr - sim.lag_lo
    p = np.array([abs(sim.frame([W.Path(0.0, (chips - d) / s.chip_rate_hz.value, 0.0)], W.Options(noise=False, n_packets=1)).cir[0][j]) ** 2
                  for d in offsets])
    return p / p[0]


def _shape(offsets, resp) -> dict:
    """-3 dB width (both sides) and the highest sidelobe past the first null, from a response sampled from 0 outward."""
    half = float(np.interp(0.5, resp[::-1][-np.argmax(resp < 0.5) - 1:][:2], offsets[::-1][-np.argmax(resp < 0.5) - 1:][:2])) \
        if np.any(resp < 0.5) else float("nan")
    k = 1
    while k < len(resp) - 1 and resp[k + 1] < resp[k]:
        k += 1                                                   # the first minimum
    return {"width_chips": 2.0 * half, "sidelobe_db": db(float(resp[k:].max())), "first_null_chips": float(offsets[k])}


def range_window(s, bench) -> dict:
    name = A.get("radar.range_window")
    step = S.num(KEY + "#scan_step_chip").value
    marks = [float(x) for x in S.raw(KEY + "#mask_offsets_chips")]
    offsets = np.arange(0.0, max(marks) + 1.0 + step / 2, step)
    eirp = R.readings(s)[0].eirp_dbm
    link = B.link(s, "monostatic", eirp)
    g, tgt, _ = bench.paths("monostatic", link, 0.0)
    cases, snr0 = [], None
    off_grid = abs(g.delay_s * s.chip_rate_hz.value * s.adc_osr.value - round(g.delay_s * s.chip_rate_hz.value * s.adc_osr.value)) / s.adc_osr.value
    straddle = lambda w: -db(float(B.range_response(off_grid, s.rolloff.value, w)[0]))      # noqa: E731  what the lag grid takes off this target's peak
    for w in (None, name):
        sim = bench.sim if w is None else W.Simulator(s, range_window=w)
        closed, got = B.range_response(offsets, s.rolloff.value, w), _response_sim(s, sim, offsets)
        snr = bench.read([tgt], bench.options("monostatic"), int(round(g.delay_s * sim.fs)), strong=(), sim=sim)["snr_db"]
        snr0 = snr if snr0 is None else snr0
        masks = [{"offset_chips": m, "closed_db": db(float(closed[np.abs(offsets - m) <= 0.5].max())), "sim_db": db(float(got[np.abs(offsets - m) <= 0.5].max()))}
                 for m in marks]
        cases.append({"window": w or "none", "closed": _shape(offsets, closed), "sim": _shape(offsets, got),
                      "loss_closed_db": B.range_window_loss_db(s.rolloff.value, w), "straddle_db": straddle(w),
                      "loss_grid_closed_db": B.range_window_loss_db(s.rolloff.value, w) + straddle(w) - straddle(None), "loss_sim_db": snr0 - snr, "masks": masks})
    return {"window": name, "cases": cases, "chip_m": s.chip_m, "off_grid_chips": off_grid, "row_a_width_m": A.get("radar.range_res_m"), "row_a_sidelobe_db": A.get("radar.window_sidelobe_db")}


# ----------------------------------------------------------------------------- 10. a second target, a static reflector
def second_target(s, bench) -> dict:
    sim = bench.sim
    eirp = R.readings(s)[0].eirp_dbm
    link = B.link(s, "monostatic", eirp)
    g, tgt, _ = bench.paths("monostatic", link, 0.0)
    m = int(s.n_packets.value)
    chips = round(g.delay_s * s.chip_rate_hz.value)
    j = chips * sim.osr - sim.lag_lo
    bins = [int(b) for b in S.raw(KEY + "#doppler_offsets_bins")]
    floor = sim.rounding_dbc

    def doppler_row(v):
        """The map's row through a noise-free target on a chip at velocity ``v``, normalised to its strongest cell."""
        p = sim.frame([W.Path(0.0, chips / s.chip_rate_hz.value, 2.0 * v)], W.Options(noise=False)).power[:, j]
        return p / p.max(), int(np.argmax(p))
    out = {"doppler": []}
    dv = s.wavelength_m / (2.0 * m * s.interval_s.value)                      # one Doppler bin, in m/s
    for label, v in (("static", 0.0), ("moving", s.velocity_mps.value)):
        row, i0 = doppler_row(v)
        frac = v / dv - round(v / dv)                                         # how far off its bin's centre the target sits
        k = np.arange(m) - i0
        k = (k + m // 2) % m - m // 2                                         # circular offset from the peak's bin
        rows = []
        for b in bins:
            far = np.abs(k) >= b
            # the target's Doppler is -2v/lambda (receding): bin offsets from its true place are k + frac
            closed = float(B.doppler_response(k[far] + frac, s.doppler_window, m).max())
            rows.append({"beyond_bins": b, "closed_db": max(db(closed + 1e-40), floor), "sim_db": max(db(float(row[far].max()) + 1e-40), floor)})
        out["doppler"].append({"what": label, "velocity_mps": v, "off_bin": frac, "rows": rows})
    # a static reflector far stronger than the target: what its phase noise puts on the floor (shared LO)
    over, rng_m = S.num(KEY + "#clutter_over_target_db").value, S.num(KEY + "#clutter_range_m").value
    gc = B.geometry(s, "monostatic", rng_m)
    level = link.p_rx_dbm.value + over
    clutter = W.Path(level, gc.delay_s, 0.0, "reflector")
    near = int(round(g.delay_s * sim.fs))
    strong = (int(round(gc.delay_s * sim.fs)),)
    lin = bench.read([tgt, clutter], bench.options("monostatic"), near, strong=strong)
    pn = bench.read([tgt, clutter], bench.options("monostatic", lo="shared"), near, strong=strong)
    res = B.pn_residual_shared(s, gc.delay_s)
    out["clutter"] = {"over_target_db": over, "range_m": rng_m, "level_dbm": level, "residual_dbc": db(res),
                      "rise_closed_db": B.floor_rise_db(level + db(res), link.noise_chip_dbm.value), "rise_sim_db": pn["floor_dbm"] - lin["floor_dbm"],
                      "snr_linear_db": lin["snr_db"], "snr_alone_db": bench.read([tgt], bench.options("monostatic"), near, strong=())["snr_db"]}
    out["bin_mps"] = dv
    return out


# ----------------------------------------------------------------------------- 11. a waveform made for sensing
def _psd_10mhz(s, chips, width_hz) -> np.ndarray:
    """The chips' power spectrum averaged over ``width_hz`` (a running mean around each frequency), on a fine grid."""
    n = 1 << int(math.ceil(math.log2(chips.size))) + 1
    p = np.abs(np.fft.fft(chips, n)) ** 2
    m = max(1, int(round(width_hz / (s.chip_rate_hz.value / n))))
    c = np.cumsum(np.concatenate([p[-m:], p, p[:m]]))
    return (c[2 * m:] - c[:-2 * m])[: n] / (2 * m)


def designed_waveform(s, bench) -> dict:
    std, new = s.trn_field, F.designed(int(s.trn_units.value))
    n = int(S.num("edmg.trn_field#trn_bl@trn_seq_length=normal").value)
    eirp = R.readings(s)[0].eirp_dbm
    items = list(S.raw("jcas60.sensing_waveform#pattern"))
    kinds = [x.strip("+-") for x in items]
    run = max(sum(1 for _ in g) for g in (list(v) for _, v in __import__("itertools").groupby(kinds * 2)))
    width = S.num("ieee80211bf.dmg_monostatic_ppdu#trn_waveform_psd_average_hz").value
    ratio = db(1.0) + 10.0 * np.log10(_psd_10mhz(s, new.chips, width) / _psd_10mhz(s, std.chips, width))
    conditions = {"length_chips": (int(new.chips.size), int(std.chips.size)), "modulus_max_error": float(np.abs(np.abs(new.chips) - 1.0).max()),
                  "consecutive_same_sequence": (run, int(S.num("ieee80211bf.dmg_monostatic_ppdu#max_consecutive_ga128").value)),
                  "psd_average_hz": width, "psd_excess_max_db": float(ratio.max()), "psd_excess_min_db": float(ratio.min()),
                  "psd_total_db": db(float(np.sum(np.abs(new.chips) ** 2) / np.sum(np.abs(std.chips) ** 2)))}
    sim = W.Simulator(s, new)
    req = B.cfar_snr_required(s, A.rep("detection.pfa", "assumed"))
    rows = []
    for name, field, fsim in (("standard TRN subfield, blocks 2-5", std, bench.sim_trn), ("designed, every block", new, sim)):
        link = B.link(s, "monostatic", eirp, field=field)
        g, tgt, _ = bench.paths("monostatic", link, 0.0)
        strong = T.strong_path_dbm(s, "monostatic", eirp, s.fe.isolation_db[-1][1])
        got = bench.read([tgt, _leak(s, strong.value)], bench.options("monostatic"), int(round(g.delay_s * fsim.fs)), sim=fsim)
        golay = T.golay_under_doppler(s, bench, (0.0, s.velocity_mps.value), fsim, field.chips.size)
        rows.append({"name": name, "correlated_chips": field.n_correlated, "clean_chips": field.zero_sidelobe_chips, "clean_m": field.zero_sidelobe_chips * s.chip_m,
                     "snr_cell_db": link.snr_cell_db.value, "snr_grid_db": link.snr_grid_db.value, "snr_sim_db": got["snr_db"],
                     "range_snr_m": B.detection_range_m(s, "monostatic", eirp, req, field).value,
                     "sidelobe_rest_dbc": golay[0][1], "sidelobe_moving_dbc": golay[1][1],
                     "inpacket_loss_db": -B.inpacket_loss_db(s, field, 2.0 * s.wavelength_m / (4.0 * s.interval_s.value))})
    margin = int(S.num(KEY + "#floor_margin_chips").value)
    delays = [int(d) for d in S.raw(KEY + "#beyond_delays_chips")]
    far = new.zero_sidelobe_chips - margin
    prof = _chip_profile(new, max(delays) + new.zero_sidelobe_chips + n * len(items))
    beyond = _beyond_rows(s, prof, sim, delays, far)
    k = np.arange(-(prof.size - 1), 1)
    zones = [{"from_chips": a, "to_chips": b, "level": max(float(prof[(k <= -a) & (k >= -b)].max()), undb(sim.rounding_dbc))}
             for a, b in ((1, new.zero_sidelobe_chips), (new.zero_sidelobe_chips + 1, n * (len(items) - 2)), (n * len(items), n * len(items)))]
    return {"pattern": items, "guard_blocks": int(S.num("jcas60.sensing_waveform#guard_blocks").value), "period_chips": n * len(items),
            "conditions": conditions, "rows": rows, "beyond": beyond, "read_to_chips": far, "zones": zones, "gain_db": db(new.n_correlated / std.n_correlated)}


# ----------------------------------------------------------------------------- 12. the limiter's knee from measured curves
CURVES = (("rx60.amp_spasaro2024_compression", "gt_db", "an amplifier"), ("rx60.rx_moody2023_compression", "conversion_gain_db", "a receiver"))
CHIP2 = "rx60.radar_bgt60tr13c"
BI = "ieee80211bf.dmg_bistatic_ppdu"


def _rapp_gain_db(pin_dbm, ip1db_dbm: float, p: float) -> np.ndarray:
    """The Rapp limiter's gain [dB] at each input power, its input 1 dB point at ``ip1db_dbm`` (fmcw_radar's ``rapp``)."""
    a = np.sqrt(10.0 ** (np.atleast_1d(np.asarray(pin_dbm, dtype=float)) / 10.0))
    return 20.0 * np.log10(np.abs(rapp(a, math.sqrt(undb(ip1db_dbm)), p)) / a)


def _fit_knee(pin, gain, ip1db: float, free_p1db: bool) -> dict:
    """Least squares of the Rapp gain curve on every point read off the figure: the knee ``p`` and the small-signal
    gain, with the input 1 dB point the paper's own figure, or (``free_p1db``) fitted as well."""
    tight = {"xatol": 1e-10}

    def left(p, p1):
        r = gain - _rapp_gain_db(pin, p1, p)
        return r - float(np.mean(r)), float(np.mean(r))

    def cost(p, p1):
        return float(np.sum(left(p, p1)[0] ** 2))

    def best_p1(p):
        return minimize_scalar(lambda p1: cost(p, p1), bounds=(ip1db - 3.0, ip1db + 3.0), method="bounded", options=tight)
    if free_p1db:
        p = float(minimize_scalar(lambda q: best_p1(q).fun, bounds=(0.3, 8.0), method="bounded", options=tight).x)
        p1 = float(best_p1(p).x)
    else:
        p = float(minimize_scalar(lambda q: cost(q, ip1db), bounds=(0.3, 8.0), method="bounded", options=tight).x)
        p1 = ip1db
    r, g0 = left(p, p1)
    return {"p": p, "ip1db_dbm": p1, "small_signal_gain_db": g0, "rms_db": float(np.sqrt(np.mean(r ** 2))), "worst_db": float(np.abs(r).max())}


def measured_knee(s, bench) -> dict:
    """ROADMAP item 1's last part: the knee had no source.  Two measured 60 GHz compression curves (spec_sources#59;
    read off the papers' figures) fitted with the same Rapp model, and what each fitted knee makes of the comparison
    table's compression criterion and of the ghost.  Only the shape is taken from them: the 1 dB point stays the
    radar chip's, and the table's own knee stays what fmcw_radar's is (law 2)."""
    over = S.num(KEY + "#knee_point_over_p1db_db").value
    rate = int(S.num(KEY + "#knee_rate").value)
    p1 = s.fe.p1db_inband_dbm.value
    span, step = float(s.sim["compression_span_db"]), float(s.sim["compression_step_db"])
    levels = p1 + np.arange(-span, span + step / 2, step)
    eirp = R.readings(s)[0].eirp_dbm
    cost = B.monostatic_cost(s, eirp, bench.sim.unit_envelope(s.leak_delay_s.value))
    rows = []
    for key, column, what in CURVES:
        pin, gain = np.asarray(S.raw(f"{key}#pin_dbm"), dtype=float), np.asarray(S.raw(f"{key}#{column}"), dtype=float)
        ip1db = S.num(f"{key}#ip1db_dbm").value
        fixed, free = _fit_knee(pin, gain, ip1db, False), _fit_knee(pin, gain, ip1db, True)
        # the other road: one point, the compression read 3 dB above the 1 dB point (how spec_sources compared them)
        at = float(gain[0] - np.interp(ip1db + over, pin, gain))
        one = brentq(lambda p: -_rapp_gain_db(ip1db + over, ip1db, p)[0] - at, 0.3, 8.0)
        roads = _limiter_roads(s, bench, rate, fixed["p"], levels)
        needs = {"in-band compression": roads["need_closed_db"], "ADC quantisation": cost.need_adc_db.value, "LO phase noise on the leakage": cost.need_pn_db.value}
        rows.append({"key": key, "what": what, "points": int(pin.size), "ip1db_dbm": ip1db, "fixed": fixed, "free": free,
                     "compression_over_db": at, "p_one_point": float(one), "roads": roads,
                     "binding": max(needs, key=needs.get), "need_db": max(needs.values())})
    placeholder = {"p": s.rx_rapp_p.value, "need_compression_db": cost.need_compression_db.value, "need_pn_db": cost.need_pn_db.value,
                   "need_adc_db": cost.need_adc_db.value, "need_db": cost.need_db.value, "binding": cost.binding,
                   "compression_over_db": float(-_rapp_gain_db(over, 0.0, s.rx_rapp_p.value)[0])}
    iso = S.num(CHIP2 + "#tx_rx_isolation_db_typ")
    return {"rows": rows, "placeholder": placeholder, "over_db": over, "rate": rate, "isolation_best_db": s.fe.isolation_db[-1][1].value,
            "chip2": {"isolation_db": iso.value, "p1db_typ_dbm": S.num(CHIP2 + "#rx_p1db_dbm_typ").value, "spare_db": iso.value - cost.need_db.value}}


# ----------------------------------------------------------------------------- 13. the long TRN subfield
def long_trn(s, bench) -> dict:
    """ROADMAP item 2's other way to a longer clean window: TRN_SEQ_LENGTH = long, every block a Golay sequence of
    256 chips (spec_sources#59 brought the pair).  Half as many units keep the field's length, so the air time and the
    correlated chips stay what the table's TRN rows have."""
    blocks = tuple(int(b) for b in S.raw("jcas60.trn_sensing#correlated_blocks"))
    n_std = int(S.num("edmg.trn_field#trn_bl@trn_seq_length=normal").value)
    n_long = int(S.num("edmg.trn_field#trn_bl@trn_seq_length=long").value)
    units = s.trn_units.value * n_std / n_long
    if units != int(units):
        raise ValueError("the table's TRN units do not make a whole number of long units of the same length")
    std, new = s.trn_field, F.trn(int(units), blocks, "long")
    long_buf = dataclasses.replace(s, sim={**s.sim, "buffer_chips": int(S.num(KEY + "#beyond_buffer_chips").value)})
    sim = W.Simulator(long_buf, new)
    eirp = R.readings(s)[0].eirp_dbm
    req = B.cfar_snr_required(s, A.rep("detection.pfa", "assumed"))
    ga, gb = F.tabulated(f"ga{n_long}"), F.tabulated(f"gb{n_long}")
    ra, rb = F.from_recursion(n_long)
    acf = np.correlate(ga, ga, "full") + np.correlate(gb, gb, "full")
    rows = []
    for name, field, fsim in ((f"TRN_BL = {n_std}, {s.trn_units.value:g} units", std, bench.sim_trn), (f"TRN_BL = {n_long}, {units:g} units", new, sim)):
        link = B.link(s, "monostatic", eirp, field=field)
        g, tgt, _ = bench.paths("monostatic", link, 0.0)
        strong = T.strong_path_dbm(s, "monostatic", eirp, s.fe.isolation_db[-1][1])
        got = bench.read([tgt, _leak(s, strong.value)], bench.options("monostatic"), int(round(g.delay_s * fsim.fs)), sim=fsim)
        golay = T.golay_under_doppler(s, bench, (0.0, s.velocity_mps.value), fsim, field.chips.size)
        rows.append({"name": name, "field_chips": int(field.chips.size), "correlated_chips": field.n_correlated, "clean_chips": field.zero_sidelobe_chips,
                     "clean_m": field.zero_sidelobe_chips * s.chip_m, "snr_cell_db": link.snr_cell_db.value, "snr_grid_db": link.snr_grid_db.value,
                     "snr_sim_db": got["snr_db"], "range_snr_m": B.detection_range_m(s, "monostatic", eirp, req, field).value,
                     "sidelobe_rest_dbc": golay[0][1], "sidelobe_moving_dbc": golay[1][1]})
    margin = int(S.num(KEY + "#floor_margin_chips").value)
    delays = [int(d) for d in S.raw(KEY + "#long_delays_chips")]
    far = new.zero_sidelobe_chips - margin
    prof = _chip_profile(new, max(delays) + new.zero_sidelobe_chips)
    worst = max(((_ghost(prof, d, far)[0], d) for d in range(new.zero_sidelobe_chips + 1, max(delays) + 1)))
    return {"rows": rows, "units": int(units), "beyond": _beyond_rows(s, prof, sim, delays, far), "read_to_chips": far,
            "worst_level": worst[0], "worst_delay_chips": worst[1], "subfield_chips": len(S.raw("edmg.trn_field#basic_subfield")) * n_long,
            "pair": {"length": n_long, "recursion_equals_table": bool(np.array_equal(ra, ga) and np.array_equal(rb, gb)),
                     "acf_peak": float(acf[n_long - 1]), "acf_worst_sidelobe": float(np.abs(np.delete(acf, n_long - 1)).max())}}


# ----------------------------------------------------------------------------- 14. the bistatic exchange as the entry has it
def bistatic_exchange(s) -> dict:
    """ROADMAP item 13: row C sends its packets one way.  The entry's exchange (spec_sources#59): BRP frames, whose SC
    packet carries at least 18 SC blocks of data, SIFS between the frames of one station, a reply from the other after
    at most BRPIFS.  Closed form only -- it is arithmetic on the packet lengths -- for two ways of laying a frame of
    the table's 64 packets on it."""
    rc, n, frames = s.chip_rate_hz.value, int(s.n_packets.value), s.frames_per_s.value
    sifs, brpifs = s.sifs_s.value, S.num("dmg.phy_timing#brpifs_s").value
    if (S.raw(BI + "#same_sta_spacing"), S.raw(BI + "#reply_spacing")) != ("SIFS", "BRPIFS"):
        raise ValueError("the bistatic entry's spacings are not SIFS within a station and BRPIFS before the reply")
    blocks = int(S.num(BI + "#min_data_sc_blocks").value)
    data = blocks * S.num("dmg.phy_params#sc_block_chips").value + S.num("dmg.phy_params#sc_gi_chips").value
    if str(S.raw("dmg.sc_ppdu_chips#data_chips")).replace(" ", "") != "N_BLKS*512+64":
        raise ValueError("dmg.sc_ppdu_chips: the data field is no longer N_BLKS blocks and one guard interval")
    sound = (s.trn_packet_chips.value + data) / rc                # a BRP frame with the TRN field behind it
    reply = (s.packet_chips.value + data) / rc                    # the BRP frame that comes back: no TRN field
    here = R.air_time(s, R.readings(s)[0], s.trn_field).fraction.value
    lam = s.wavelength_m

    def doppler(pri):
        return {"interval_s": pri, "velocity_max_mps": lam / (4.0 * pri), "velocity_res_mps": lam / (2.0 * n * pri)}
    # (a) one exchange a frame: the packets back to back at SIFS, one reply after BRPIFS (at most) -- or after SIFS (at least)
    burst = n * sound + (n - 1) * sifs
    a = {"frame_s": burst + brpifs + reply, "air": (burst + brpifs + reply) * frames, "air_min": (burst + sifs + reply) * frames, **doppler(sound + sifs)}
    # (b) one exchange a packet, at the table's packet interval: sounding, BRPIFS, reply
    each = sound + brpifs + reply
    b = {"exchange_s": each, "fits_the_interval": bool(each + sifs <= s.interval_s.value), "air": each * n * frames,
         "air_emitted": (sound + reply) * n * frames, **doppler(s.interval_s.value)}
    body = S.raw(BI + "#body_octets")
    return {"data_blocks": blocks, "data_chips": data, "data_s": data / rc, "sounding_s": sound, "sounding_now_s": s.packet_seconds(s.trn_field).value,
            "reply_s": reply, "sifs_s": sifs, "brpifs_s": brpifs, "air_now": here, "one_exchange_a_frame": a, "one_exchange_a_packet": b,
            "table": doppler(s.interval_s.value), "body_octets": int(sum(body.values())), "frame_type": str(S.raw(BI + "#frame_type"))}


# ----------------------------------------------------------------------------- all of them
@lru_cache(maxsize=1)
def build() -> dict:
    s = config.setup()
    bench = _bench(s)
    return {"setup": s, "limiter": limiter(s, bench), "gains": gain_settings(s, bench), "beyond": beyond_window(s, bench),
            "fast": fast_targets(s, bench), "adc": undithered_adc(s, bench), "leak": leakage_paths(s, bench),
            "lo": same_lo(s, bench), "sync": bistatic_sync(s, bench), "window": range_window(s, bench), "second": second_target(s, bench),
            "designed": designed_waveform(s, bench), "knee": measured_knee(s, bench), "long_trn": long_trn(s, bench), "exchange": bistatic_exchange(s)}


def _f(v, fmt="{:.2f}") -> str:
    return "—" if v is None else fmt.format(v)


def render(t: dict | None = None) -> str:
    from . import ENGINE_BANNER, WORDING
    t = t or build()
    s = t["setup"]
    L: list[str] = []
    w = L.append
    w(ENGINE_BANNER)
    w(WORDING)
    w("")
    lim = t["limiter"]
    w("== 1. the limiter in continuous time, and its knee (row B, L-CEF, communication limit; ROADMAP item 1) ==")
    w(f"{'samples a chip ahead of the ADC':>32s} {'mean 1/envelope':>16s} {'isolation for 1 dB: closed':>28s} {'simulated':>10s} "
      f"{'hard-limit SNR loss':>20s} {'image at ' + format(lim['by_rate'][0]['isolation_db'], 'g') + ' dB: bound':>24s} {'simulated':>10s}")
    for r in lim["by_rate"]:
        w(f"{r['analog']:32d} {r['inverse_envelope_mean']:16.4f} {r['need_closed_db']:25.2f} dB {_f(r['need_sim_db']):>7s} dB "
          f"{-r['limit_closed_db']:17.2f} dB {r['image_bound_dbc']:20.1f} dBc {r['image_sim_dbc']:6.1f} dBc")
    w(f"  the first row is the table's own (the ADC's rate); the limiter (jcas60.rx_rapp_p = {s.rx_rapp_p.value:g}) acts on the leakage delayed by "
      f"{s.leak_delay_s.value * s.chip_rate_hz.value:.3f} chip")
    w(f"  the knee, at {lim['knee_rate']} samples a chip (the map's floor for this target is {lim['by_knee'][0]['floor_dbc']:.1f} dBc):")
    w(f"{'Rapp p':>32s} {'isolation for 1 dB: closed':>28s} {'simulated':>10s} {'image at ' + format(lim['by_knee'][0]['isolation_db'], 'g') + ' dB: bound':>24s} "
      f"{'simulated':>10s} {'isolation for the image under the floor: closed (bound)':>58s} {'simulated':>10s}")
    for r in lim["by_knee"]:
        w(f"{r['knee']:32g} {r['need_closed_db']:25.2f} dB {_f(r['need_sim_db']):>7s} dB {r['image_bound_dbc']:20.1f} dBc {r['image_sim_dbc']:6.1f} dBc "
          f"{_f(r['ghost_closed_db']):>55s} dB {_f(r['ghost_sim_db']):>7s} dB")
    w("")
    gs = t["gains"]
    w(f"== 2. the receiver's gain setting as a variable ({FIG}; row B, L-CEF, communication limit; ROADMAP item 2) ==")
    w(f"  the table's rows use the data sheet's typical noise figure ({gs['table_nf_db']:g} dB) with the {gs['table_gain_db']:g} dB setting's compression point; "
      "here each setting has the noise figure of the same figure")
    head = "".join(f" {'SNR with ' + format(h['isolation_db'], 'g') + ' dB isolation: closed / simulated':>50s}" for h in gs["rows"][0]["have"])
    w(f"{'gain':>8s} {'NF':>8s} {'P1dB in band':>14s} {'SNR, no leakage':>16s} {'isolation for 1 dB: compression / PN / ADC':>44s}{head}")
    for k, r in enumerate(gs["rows"]):
        have = "".join(f" {h['snr_closed_db']:36.2f} / {h['snr_sim_db']:6.2f} dB{'*' if gs['best'][j] == k else ' '}" for j, h in enumerate(r["have"]))
        w(f"{r['gain_db']:5g} dB {r['nf_db']:5.1f} dB {r['p1db_dbm']:10g} dBm {r['snr_no_leak_db']:13.2f} dB "
          f"{r['need_compression_db']:28.1f} / {r['need_pn_db']:.1f} / {r['need_adc_db']:.1f} dB{have}")
    for j, h0 in enumerate(gs["rows"][0]["have"]):
        b = gs["rows"][gs["best"][j]]
        h = b["have"][j]
        w(f"  {h0['isolation_db']:g} dB of isolation: the best setting (closed form, *) is {b['gain_db']:g} dB, SNR {h['snr_closed_db']:.2f} dB -- "
          f"{gs['snr_no_leak_max_db'] - h['snr_closed_db']:.2f} dB under the best setting without leakage ({gs['snr_no_leak_max_db']:.2f} dB); there the echo "
          f"{h['echo_db']:+.2f} dB, the thermal noise {h['thermal_db']:+.2f} dB, the phase-noise residue {h['pn_over_thermal_db']:+.1f} dB and the ADC "
          f"{h['adc_over_thermal_db']:+.1f} dB against it")
    w("  the three isolation figures are each 'at most 1 dB over the thermal floor': a noisier setting passes them more easily, which is why the SNR that is left is the figure to read")
    w("  what the closed form adds to example 00: the limiter takes the thermal noise down more than it takes the leakage down, so the phase-noise residue and the "
      "ADC's noise, which ride on the leakage, stand higher against it than in a linear receiver")
    low = gs["rows"][0]["have"][0]
    w(f"  at {low['isolation_db']:g} dB of isolation less than 10 dB of SNR is left and the roads part by {low['snr_sim_db'] - low['snr_closed_db']:.1f} dB: the closed form over-states the "
      "ADC's noise where nothing dithers it (example 00, section 5) and the simulation's peak search reads high at that SNR; neither figure is a design point")
    w("")
    bw = t["beyond"]
    w("== 3. a target beyond the clean window (one packet, no noise, a target of the strength that reads 0 dB inside the window; ROADMAP items 3 and 15) ==")
    for f in bw["fields"]:
        w(f"  {f['name']}: clean window {f['window_chips']} chips = {f['window_m']:.2f} m of range, read to {f['read_to_chips']} chips")
        w(f"{'target at [chips]':>20s} {'range [m]':>10s} {'strongest in the window: chips road':>38s} {'at lag':>8s} {'simulation':>12s} {'at lag':>8s}")
        for r in f["rows"]:
            w(f"{r['delay_chips']:20d} {r['range_m']:10.2f} {db(r['level_closed']):35.1f} dB {_f(r['lag_closed'], '{:d}'):>8s} "
              f"{db(r['level_sim']):9.1f} dB {_f(r['lag_sim'], '{:d}'):>8s}")
        w(f"    the worst of any delay past the window up to {f['rows'][-1]['delay_chips']} chips: {db(f['worst_level']):.1f} dB, a target at {f['worst_delay_chips']} chips "
          f"({f['worst_delay_chips'] * s.chip_m:.2f} m)")
    w("  a level at the arithmetic's rounding bound (example 00, section 7) is an exact cancellation: the far edge of the window is clean")
    w("  the L-CEF has nothing ahead of its cyclic prefix here; in a packet the STF is there (its layout: spec_sources#59), so its figures past 128 chips are this model's, not a packet's")
    w("  inside the window, the two roads (closed form on the grid against the waveform simulation):")
    for r in bw["ranges"]:
        w(f"    {r['field']:10s} {r['range_m']:5.1f} m ({r['delay_chips']:6.1f} chips): closed {r['snr_closed_db']:6.2f} dB, simulated {r['snr_sim_db']:6.2f} dB, "
          f"difference {r['snr_sim_db'] - r['snr_closed_db']:+.2f} dB")
    w("")
    ft = t["fast"]
    w(f"== 4. what a fast target loses inside one packet (one packet, no noise, the target on a chip; unambiguous velocity {ft['v_unambiguous_mps']:.1f} m/s; ROADMAP item 4) ==")
    for f in ft["fields"]:
        w(f"  {f['name']}: the correlated chips span {f['span_s'] * 1e6:.2f} us")
        w(f"{'radial velocity [m/s]':>26s} {'closed form':>14s} {'simulated':>12s}")
        for r in f["rows"]:
            w(f"{r['velocity_mps']:26.1f} {r['loss_closed_db']:11.4f} dB {r['loss_sim_db']:9.4f} dB")
    w("  closed form: |mean over the correlated chips of exp(j 2 pi f_d t)|^2 from where the correlation windows lie; inside the unambiguous velocity the table needs no such term")
    w("")
    ad = t["adc"]
    w(f"== 5. the ADC that nothing dithers ({ad['isolation_db']:g} dB of isolation: leakage {ad['leak_dbm']:g} dBm; linear receiver; the target moved over delays "
      f"{ad['step_chips']:g} chip apart; ROADMAP item 5) ==")
    w(f"{'dither [LSB rms]':>18s} {'delays':>7s} {'error the same in every packet':>32s} {'floor rise: closed':>20s} {'simulated':>10s} "
      f"{'zero-Doppler ridge under the target: median (min ... max)':>60s} {'the peak moves: rms (min ... max)':>36s} {'SNR change':>12s}")
    for c in ad["cases"]:
        rg, pk = c["ridge_under_target_db"], c["peak_change_db"]
        w(f"{c['dither_lsb']:18.2f} {c['delays']:7d} {c['repeating_share'] * 100:30.0f} % {c['rise_closed_db']:17.2f} dB {c['rise_db']['mean']:7.2f} dB "
          f"{rg['median']:38.1f} dB ({rg['min']:6.1f} ... {rg['max']:6.1f}) {pk['rms']:16.2f} dB ({pk['min']:+5.2f} ... {pk['max']:+5.2f}) {c['snr_change_db']['mean']:9.2f} dB")
    w(f"  thermal noise is {ad['cases'][0]['noise_lsb']:.2f} LSB rms; the dither is Gaussian, added ahead of the quantiser and not taken out again; the closed form counts the "
      "quantiser's 1/12 LSB^2 and the dither as noise")
    w("")
    lk = t["leak"]
    w(f"== 6. the leakage's delay, and a near reflector (shared LO; floor rise at {lk['isolation_db']:g} dB of isolation, leakage {lk['leak_dbm']:g} dBm; ROADMAP item 6) ==")
    w(f"{'leakage delay [ns]':>20s} {'phase-noise residue':>22s} {'isolation for 1 dB':>20s} {'floor rise: closed':>20s} {'simulated':>10s}")
    for r in lk["rows"]:
        w(f"{r['delay_s'] * 1e9:20.2f} {r['residual_dbc']:18.1f} dBc {r['need_pn_db']:17.1f} dB {r['rise_closed_db']:17.2f} dB {r['rise_sim_db']:7.2f} dB")
    rf = lk["reflector"]
    w(f"  a reflector {rf['range_m'] * 100:g} cm out, {rf['level_db']:g} dB under the direct leakage (delay {rf['delay_s'] * 1e9:.3f} ns, {rf['turns']:.2f} carrier turns behind it): "
      "the two residues come from one LO and add in amplitude")
    w(f"      isolation for 1 dB: {rf['need_pn_alone_db']:.1f} dB without it; {rf['need_pn_db']:.1f} dB with it at this phase; {rf['need_pn_min_db']:.1f} ... {rf['need_pn_max_db']:.1f} dB "
      f"over the phase between the paths, {rf['need_pn_average_db']:.1f} dB on average")
    w(f"      floor rise with it: closed {rf['rise_closed_db']:.2f} dB / simulated {rf['rise_sim_db']:.2f} dB (without it: {rf['rise_alone_closed_db']:.2f} dB)")
    w("")
    lo = t["lo"]
    w("== 7. every row on one LO: mmw60_trx's x10 LO in the place of each row's own (the table keeps fmcw_radar's default config; ROADMAP item 6) ==")
    for k, label in (("default", "A   fmcw_radar's default config"), ("x10lo", "A'  its radar_60g_x10lo.yaml")):
        r = lo["row_a"][k]
        w(f"  {label:32s} SNR(cell) {r['snr_cell_db']:.2f} dB; LO {r['lo_rms_rad']:.3f} rad rms; phase-noise floor of the map: mean {r['pn_mean_cell_dbc']:.1f} dBc, "
          f"at the target {r['pn_local_cell_dbc']:.1f} dBc   ({r['lo_source']})")
    w(f"  the ay rows (communication limit; floor rises at {lo['isolation_db']:g} dB of isolation for B, on the direct path for C):")
    w(f"{'LO':>44s} {'L(1 MHz)':>10s} {'L(100 MHz)':>11s} {'B: residue':>12s} {'isolation for 1 dB':>20s} {'of all three':>14s} {'floor rise closed / sim':>26s} "
      f"{'C, L-CEF: closed / sim':>24s} {'C, TRN field: closed / sim':>28s}")
    for r in lo["rows"]:
        c0, c1 = r["bistatic"]
        w(f"{r['lo']:>44s} {r['at_1mhz_dbchz']:10.1f} {r['at_100mhz_dbchz']:11.1f} {r['residual_dbc']:8.1f} dBc {r['need_pn_db']:17.1f} dB {r['need_db']:11.1f} dB "
          f"{r['rise_closed_db']:15.2f} / {r['rise_sim_db']:.2f} dB {c0['rise_closed_db']:13.2f} / {c0['rise_sim_db']:.2f} dB {c1['rise_closed_db']:17.2f} / {c1['rise_sim_db']:.2f} dB")
    for r in lo["rows"]:
        c0, c1 = r["bistatic"]
        w(f"      {r['lo']}: B is bound by {r['binding']}; C's residue {c0['residual_dbc']:.1f} dBc, the phase noise less its mean over one correlation window "
          f"({c0['reference_s'] * 1e6:.2f} us); taken over the whole correlated span instead the closed form would say {c0['rise_whole_field_closed_db']:.2f} dB "
          f"(L-CEF, {c0['span_s'] * 1e6:.2f} us) and {c1['rise_whole_field_closed_db']:.2f} dB (TRN field, {c1['span_s'] * 1e6:.2f} us)")
    w("")
    sy = t["sync"]
    w("== 8. the bistatic rows: the carrier offset the two devices are left with (one frame, no noise; every packet's phase taken from the direct path; ROADMAP item 7) ==")
    for f in sy["fields"]:
        w(f"  {f['name']}:")
        w(f"{'offset':>16s} {'as a velocity':>16s} {'turns across the field':>24s} {'the target loses: closed':>26s} {'simulated':>12s}")
        for r in f["rows"]:
            w(f"{r['cfo_hz'] / 1e3:12.0f} kHz {r['cfo_hz'] * sy['equivalent_mps_per_hz']:12.1f} m/s {r['turns']:24.3f} {r['loss_closed_db']:23.4f} dB {r['loss_sim_db']:9.4f} dB")
    w("  between packets nothing is left: the offset is common to the target and to the direct path the phase is read on; inside a packet it is the fast target's loss of section 4, "
      "the offset added to the target's own Doppler")
    w("")
    rw = t["window"]
    w(f"== 9. a range window: row A's ({rw['window']}) across the occupied band, after the matched filter (ROADMAP item 10) ==")
    w(f"{'window':>10s} {'-3 dB width: closed':>22s} {'simulated':>12s} {'in range':>10s} {'highest sidelobe: closed':>26s} {'simulated':>10s} {'SNR lost at the peak':>22s} {'on the lag grid: closed':>25s} {'simulated':>10s}")
    for c in rw["cases"]:
        w(f"{c['window']:>10s} {c['closed']['width_chips']:16.3f} chips {c['sim']['width_chips']:6.3f} chips {c['sim']['width_chips'] * rw['chip_m'] * 100:7.2f} cm "
          f"{c['closed']['sidelobe_db']:23.1f} dB {c['sim']['sidelobe_db']:7.1f} dB {c['loss_closed_db']:19.2f} dB {c['loss_grid_closed_db']:22.2f} dB {c['loss_sim_db']:7.2f} dB")
    wide = rw["cases"][-1]["sim"]["width_chips"] * rw["chip_m"]
    w(f"  row A: {rw['row_a_width_m'] * 100:.2f} cm with the same window, sidelobes {rw['row_a_sidelobe_db']:.1f} dB; like for like the ay rows are "
      f"{wide / rw['row_a_width_m']:.2f} times as wide (without a window on either side the bandwidths give 2.27)")
    w(f"  on the lag grid the window costs less than at the peak: this target sits {rw['off_grid_chips']:.2f} chip off a lag, and a wider main lobe loses less there "
      f"({rw['cases'][0]['straddle_db']:.2f} dB without the window, {rw['cases'][-1]['straddle_db']:.2f} dB with it)")
    w("  what a target leaves for a second one at the same velocity (the highest level within half a chip of each offset):")
    w(f"{'offset [chips]':>16s} " + " ".join(f"{c['window'] + ': closed / simulated':>30s}" for c in rw["cases"]))
    for k, mk in enumerate(rw["cases"][0]["masks"]):
        w(f"{mk['offset_chips']:16.1f} " + " ".join(f"{c['masks'][k]['closed_db']:19.1f} / {c['masks'][k]['sim_db']:6.1f} dB" for c in rw["cases"]))
    w("")
    st2 = t["second"]
    w(f"== 10. what a strong target leaves in Doppler, and a static reflector (one frame, no noise; a Doppler bin is {st2['bin_mps']:.2f} m/s; ROADMAP item 9) ==")
    w(f"{'the strong one':>16s} {'off its bin':>12s} " + " ".join(f"{'beyond ' + str(r['beyond_bins']) + ' bins: closed / sim':>34s}" for r in st2["doppler"][0]["rows"]))
    for d in st2["doppler"]:
        w(f"{d['what'] + ' ' + format(d['velocity_mps'], 'g') + ' m/s':>16s} {d['off_bin']:+12.2f} " + " ".join(f"{r['closed_db']:23.1f} / {r['sim_db']:6.1f} dB" for r in d["rows"]))
    w("  a static path sits on a bin and the window's nulls fall on every other bin: past the main lobe it leaves the arithmetic's rounding; a mover is off its bin and leaves the window's sidelobes")
    cl = st2["clutter"]
    w(f"  a static reflector {cl['over_target_db']:g} dB over the target, at {cl['range_m']:g} m ({cl['level_dbm']:.1f} dBm): linear receiver and ideal LO, the target's SNR is "
      f"{cl['snr_linear_db']:.2f} dB (alone: {cl['snr_alone_db']:.2f} dB); with the shared LO's phase noise (residue {cl['residual_dbc']:.1f} dBc at that delay) the floor rises "
      f"{cl['rise_closed_db']:.2f} dB closed / {cl['rise_sim_db']:.2f} dB simulated")
    w("")
    dw = t["designed"]
    cd = dw["conditions"]
    w("== 11. a waveform made for sensing in the TRN field's place (jcas60.sensing_waveform; the entry ieee80211bf.dmg_monostatic_ppdu allows one on four conditions; ROADMAP item 12) ==")
    w(f"  one period: {' '.join(dw['pattern'])} ({dw['period_chips']} chips); {dw['guard_blocks']} guard blocks at each end, every block between them correlated")
    w(f"{'':>36s} {'correlated chips':>18s} {'clean zone':>22s} {'SNR(cell)':>11s} {'on the grid: closed':>21s} {'simulated':>10s} {'range at Pd 0.9 by SNR':>24s} "
      f"{'sidelobes at rest':>19s} {'at 2 m/s':>10s} {'in-packet loss at v_max':>25s}")
    for r in dw["rows"]:
        w(f"{r['name']:>36s} {r['correlated_chips']:18d} {r['clean_chips']:6d} chips = {r['clean_m']:5.2f} m {r['snr_cell_db']:8.2f} dB {r['snr_grid_db']:18.2f} dB "
          f"{r['snr_sim_db']:7.2f} dB {r['range_snr_m']:22.2f} m {r['sidelobe_rest_dbc']:15.1f} dBc {r['sidelobe_moving_dbc']:6.1f} dBc {r['inpacket_loss_db']:22.3f} dB")
    w(f"  the entry's conditions: length {cd['length_chips'][0]} chips against the TRN field's {cd['length_chips'][1]}: fits; power over every 128 Tc: the same (unit modulus, "
      f"worst error {cd['modulus_max_error']:.0e}): fits; at most {cd['consecutive_same_sequence'][1]} of one sequence in a row: {cd['consecutive_same_sequence'][0]}: fits")
    w(f"      spectral density averaged over {cd['psd_average_hz'] / 1e6:g} MHz, no higher than the TRN field's: {cd['psd_excess_max_db']:+.2f} dB at the worst frequency "
      f"({cd['psd_excess_min_db']:+.2f} dB at the best, {cd['psd_total_db']:+.2f} dB in all): does not fit as built -- the excess is where the TRN field's own spectrum dips")
    w(f"      backed off by that {cd['psd_excess_max_db']:.2f} dB it fits, and the {dw['gain_db']:.2f} dB of more chips become {dw['gain_db'] - cd['psd_excess_max_db']:+.2f} dB against the standard subfield; "
      "the longer clean zone and what follows do not depend on the level")
    w("  what a target leaves for a lag this many chips from it, either side (the chips alone; rounding bound = exact cancellation): "
      + "; ".join(f"{z['from_chips']}" + ("" if z['from_chips'] == z['to_chips'] else f"-{z['to_chips']}") + f": {db(z['level']):.1f} dB" for z in dw["zones"]))
    w(f"{'target at [chips]':>20s} {'range [m]':>10s} {'strongest on lags 0-' + str(dw['read_to_chips']) + ': chips road':>38s} {'simulation':>12s}")
    for r in dw["beyond"]:
        w(f"{r['delay_chips']:20d} {r['range_m']:10.2f} {db(r['level_closed']):35.1f} dB {db(r['level_sim']):9.1f} dB")
    w("  against section 3: the standard subfield shows a target 256 chips past its place 6.1 dB down and folds at 768 chips; this one is clean to 256 chips, "
      "noise-like about 21 dB down from there to 768, and folds at 1024 chips")
    w("")
    kn = t["knee"]
    ph = kn["placeholder"]
    w("== 12. the limiter's knee from measured curves (ROADMAP item 1; spec_sources#59: two 60 GHz compression curves read off the papers' figures, fitted with the same Rapp model) ==")
    w(f"{'curve':>36s} {'points':>7s} {'its P1dB':>10s} {'p, P1dB as published':>22s} {'residual rms / worst':>22s} {'p, P1dB fitted too':>20s} {'P1dB':>10s} {'rms':>9s} "
      f"{'compression ' + format(kn['over_db'], 'g') + ' dB over P1dB':>30s} {'p from that point':>19s}")
    for r in kn["rows"]:
        w(f"{r['key']:>36s} {r['points']:7d} {r['ip1db_dbm']:6.1f} dBm {r['fixed']['p']:22.3f} {r['fixed']['rms_db']:11.3f} / {r['fixed']['worst_db']:.3f} dB {r['free']['p']:20.3f} "
          f"{r['free']['ip1db_dbm']:6.2f} dBm {r['free']['rms_db']:6.3f} dB {r['compression_over_db']:27.2f} dB {r['p_one_point']:19.3f}")
    w(f"  the table's knee is p = {ph['p']:g}, fmcw_radar's placeholder: {ph['compression_over_db']:.2f} dB of compression at that point.  Both measured curves bend more gently")
    w(f"  what each knee makes of row B (communication limit; the radar chip's own P1dB kept, only the shape exchanged; {kn['rate']} samples a chip):")
    w(f"{'knee':>36s} {'compression needs: closed':>27s} {'simulated':>10s} {'phase noise':>12s} {'ADC':>9s} {'isolation needed':>18s} {'what binds':>30s} "
      f"{'against ' + format(kn['isolation_best_db'], 'g') + ' dB':>15s} {'ghost under the floor: bound':>30s} {'simulated':>10s}")
    w(f"{'p = ' + format(ph['p'], 'g') + ' (the table)':>36s} {ph['need_compression_db']:24.2f} dB {'':>10s} {ph['need_pn_db']:9.2f} dB {ph['need_adc_db']:6.2f} dB {ph['need_db']:15.2f} dB "
      f"{ph['binding']:>30s} {ph['need_db'] - kn['isolation_best_db']:+12.2f} dB")
    for r in kn["rows"]:
        q = r["roads"]
        w(f"{'p = ' + format(r['fixed']['p'], '.3f') + ' (' + r['what'] + ')':>36s} {q['need_closed_db']:24.2f} dB {_f(q['need_sim_db']):>7s} dB {ph['need_pn_db']:9.2f} dB {ph['need_adc_db']:6.2f} dB "
          f"{r['need_db']:15.2f} dB {r['binding']:>30s} {r['need_db'] - kn['isolation_best_db']:+12.2f} dB {_f(q['ghost_closed_db']):>27s} dB {_f(q['ghost_sim_db']):>7s} dB")
    w("  the ghost's closed form is the small-signal image, an upper bound once the leakage is deep in the limiter (the table's section 5); the simulated figure is the reading")
    w("  a gentler knee lets the compression criterion go and leaves the LO's phase noise on the leakage as what binds: the isolation needed hardly moves, the reason does")
    w("  not known from these curves: AM/PM (neither paper gives it), memory, and whether the radar chip's receiver bends like either device")
    c2 = kn["chip2"]
    w(f"  another radar chip's datasheet ({CHIP2}): TX-to-RX isolation {c2['isolation_db']:g} dB typ at the die, package excluded -- {c2['spare_db']:+.2f} dB against the table's need; "
      f"its receiver's P1dB is {c2['p1db_typ_dbm']:g} dBm typ, one number, no curve")
    w("")
    lt = t["long_trn"]
    w(f"== 13. the long TRN subfield: TRN_BL = {lt['pair']['length']} (ROADMAP item 2; the pair is spec_sources#59's edmg.golay_sequences_long) ==")
    w(f"  Ga{lt['pair']['length']} / Gb{lt['pair']['length']}: the entry's D and W through the recursion give the tabulated sequences: {'yes' if lt['pair']['recursion_equals_table'] else 'NO'}; "
      f"aperiodic autocorrelations summed: {lt['pair']['acf_peak']:g} at zero shift, {lt['pair']['acf_worst_sidelobe']:g} at the worst other")
    w(f"{'':>28s} {'field':>12s} {'correlated chips':>18s} {'clean zone':>22s} {'SNR(cell)':>11s} {'on the grid: closed':>21s} {'simulated':>10s} {'range at Pd 0.9 by SNR':>24s} "
      f"{'sidelobes at rest':>19s} {'at 2 m/s':>10s}")
    for r in lt["rows"]:
        w(f"{r['name']:>28s} {r['field_chips']:6d} chips {r['correlated_chips']:18d} {r['clean_chips']:6d} chips = {r['clean_m']:5.2f} m {r['snr_cell_db']:8.2f} dB {r['snr_grid_db']:18.2f} dB "
          f"{r['snr_sim_db']:7.2f} dB {r['range_snr_m']:22.2f} m {r['sidelobe_rest_dbc']:15.1f} dBc {r['sidelobe_moving_dbc']:6.1f} dBc")
    w(f"{'target at [chips]':>20s} {'range [m]':>10s} {'strongest on lags 0-' + str(lt['read_to_chips']) + ': chips road':>38s} {'simulation':>12s} {'shows at [chips]':>18s}")
    for r in lt["beyond"]:
        w(f"{r['delay_chips']:20d} {r['range_m']:10.2f} {db(r['level_closed']):35.1f} dB {db(r['level_sim']):9.1f} dB {_f(r['lag_closed'], '{:d}'):>18s}")
    w(f"  the worst a target anywhere past the window does: {db(lt['worst_level']):.1f} dB, from {lt['worst_delay_chips']} chips; the subfield is {lt['subfield_chips']} chips, "
      "where a target folds back whole.  The same field length, air time and SNR as the table's TRN rows, twice the clean zone, the ghost twice as far")
    w("")
    ex = t["exchange"]
    xa, xb, xt = ex["one_exchange_a_frame"], ex["one_exchange_a_packet"], ex["table"]
    w("== 14. the bistatic exchange as the entry has it (ROADMAP item 13; ieee80211bf.dmg_bistatic_ppdu after spec_sources#59; closed form only: arithmetic on packet lengths) ==")
    w(f"  a BRP frame ({ex['frame_type']}; {ex['body_octets']} octets of body known) in an SC packet carries at least {ex['data_blocks']} SC blocks: {ex['data_chips']:g} chips, {ex['data_s'] * 1e6:.2f} us of data field")
    w(f"  the sounding packet {ex['sounding_now_s'] * 1e6:.2f} -> {ex['sounding_s'] * 1e6:.2f} us; the BRP frame that comes back {ex['reply_s'] * 1e6:.2f} us; "
      f"SIFS {ex['sifs_s'] * 1e6:g} us between one station's frames, at most BRPIFS {ex['brpifs_s'] * 1e6:g} us before the other's")
    w(f"{'':>52s} {'air time':>12s} {'packet interval':>17s} {'unambiguous to +-':>22s} {'velocity resolution':>21s}")
    w(f"{'row C as the table has it (one way, no data field)':>52s} {ex['air_now'] * 100:10.2f} % {xt['interval_s'] * 1e6:14.2f} us {xt['velocity_max_mps']:18.2f} m/s {xt['velocity_res_mps']:17.2f} m/s")
    w(f"{'one exchange a frame: 64 frames at SIFS, one reply':>52s} {xa['air'] * 100:10.2f} % {xa['interval_s'] * 1e6:14.2f} us {xa['velocity_max_mps']:18.2f} m/s {xa['velocity_res_mps']:17.2f} m/s")
    w(f"{'one exchange a packet, at the table interval':>52s} {xb['air'] * 100:10.2f} % {xb['interval_s'] * 1e6:14.2f} us {xb['velocity_max_mps']:18.2f} m/s {xb['velocity_res_mps']:17.2f} m/s")
    w(f"  one exchange a frame: {xa['air_min'] * 100:.2f} % if the reply comes after SIFS; the packets cannot be spread out, so the Doppler axis is the packet's own length")
    w(f"  one exchange a packet: sounding, BRPIFS and reply take {xb['exchange_s'] * 1e6:.2f} us of the {xt['interval_s'] * 1e6:g} us ({'fits' if xb['fits_the_interval'] else 'DOES NOT FIT'}); "
      f"{xb['air_emitted'] * 100:.2f} % is emitted, the rest of the {xb['air'] * 100:.2f} % is the medium held through BRPIFS")
    w("  not found in public material (spec_sources#59): the MAC header's length, a reply once per several exchanges, an exchange with no reply")
    return "\n".join(L) + "\n"
