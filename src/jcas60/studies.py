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
"""
from __future__ import annotations

import dataclasses
import math
from functools import lru_cache

import numpy as np

from . import budget as B
from . import config
from . import regulatory as R
from . import sources as S
from . import table as T
from . import wavesim as W
from .physconst import db, undb
from .sources import Num

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
        chip0 = (-sim.lag_lo) % sim.osr                              # the lag index of the chip-spaced lags
        rows = []
        floor = undb(sim.rounding_dbc)                               # under it a level is the arithmetic's rounding, not a reading
        for d in delays:
            closed, lag = _ghost(prof, d, far)
            cir = sim.frame([W.Path(0.0, d / s.chip_rate_hz.value, 0.0)], W.Options(noise=False, n_packets=1)).cir[0]
            p = np.abs(cir[chip0:: sim.osr]) ** 2
            lags = sim.lags[chip0:: sim.osr] // sim.osr
            keep = (lags >= 0) & (lags <= far)
            j = int(np.argmax(np.where(keep, p, 0.0)))
            rows.append({"delay_chips": d, "range_m": d * s.chip_m, "level_closed": max(closed, floor), "lag_closed": lag if closed > floor else None,
                         "level_sim": max(float(p[j]), floor), "lag_sim": int(lags[j]) if p[j] > floor else None})
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


# ----------------------------------------------------------------------------- all of them
@lru_cache(maxsize=1)
def build() -> dict:
    s = config.setup()
    bench = _bench(s)
    return {"setup": s, "limiter": limiter(s, bench), "gains": gain_settings(s, bench), "beyond": beyond_window(s, bench)}


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
    return "\n".join(L) + "\n"
