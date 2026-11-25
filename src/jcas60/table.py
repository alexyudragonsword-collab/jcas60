"""Milestone 0: the comparison table, built by code and printed by ``examples/00_comparison_table.py``.

Five rows on one target and one set of front-end numbers:

    A        fmcw_radar's default radar, as it reports itself (nothing recomputed);
    B        the ay single-carrier Golay training field, monostatic -- under the communication limit and
             under the radar limit (the rule does not say which applies: both are rows);
    C        the same, bistatic -- under both limits.

``build()`` returns the rows as data (every number a :class:`~jcas60.sources.Num`); ``render()`` turns
them into the text.  The closed form (``budget``) and the waveform simulation (``wavesim``) are both
run and their difference is a column.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from . import ENGINE_BANNER, WORDING, __version__
from . import bf
from . import budget as B
from . import fmcw_row as A
from . import regulatory as R
from . import sources as S
from . import wavesim as W
from .config import Setup, setup
from .physconst import db
from .sources import Num, derive
from .vendor.fmcw_radar.processing import rd as _rd

N = Num


# ----------------------------------------------------------------------------- row A
def row_a() -> dict:
    g = A.governing_tier()
    pfa = A.rep("detection.pfa", "assumed")
    snr = A.rep("detection.snr_cell_db", unit="dB")
    if (A.get("detection.range_at_pd.pd"), A.get("detection.range_at_pd.pfa")) != (setup().pd.value, pfa.value):
        raise ValueError("fmcw_radar's report gives its detection range at another Pd / Pfa than this table's")
    fs, tc, nch = A.rep("radar.adc_fs_msps", "assumed"), A.rep("radar.chirp_us", "assumed"), A.rep("radar.n_chirps", "assumed")
    on, period = A.rep("compliance.frame_on_s", "assumed", "s"), A.rep("compliance.frame_period_s", "assumed", "s")
    rate = R.comm_rate_mbps()
    frac = derive(on.value / period.value, "frame_on_s / frame_period_s", on, period)
    return {
        "name": "A  FMCW (fmcw_radar default)", "limit": "radar limit", "kind": "monostatic",
        "bandwidth_hz": N(A.get("radar.sweep_ghz") * 1e9, "assumed", "fmcw_radar report: radar.sweep_ghz", "Hz"),
        "range_res_m": A.rep("radar.range_res_rect_m", unit="m"),
        "range_width_m": A.rep("radar.range_res_m", unit="m"), "range_window": A.get("radar.range_window"),
        "range_sidelobe_db": A.rep("radar.window_sidelobe_db", unit="dB"),
        "range_max_m": A.rep("radar.max_range_m", unit="m"), "range_max_why": "the IF band",
        "velocity_res_mps": A.rep("radar.velocity_res_mps", unit="m/s"),
        "velocity_max_mps": A.rep("radar.max_velocity_mps", unit="m/s"),
        "eirp_dbm": A.budget_term("EIRP"),
        "processing_gain_db": derive(db(fs.value * tc.value * nch.value), "10 log10(samples per chirp x chirps)", fs, tc, nch, unit="dB"),
        "processing_what": f"{fs.value * tc.value:.0f} samples x {nch.value:.0f} chirps",
        "window_loss_db": A.budget_term("window loss"),
        "snr_db": snr, "pd": A.rep("detection.pd"), "pfa": pfa,
        "snr_required_db": A.rep("detection.range_at_pd.snr_required_db", unit="dB"),
        "range_pd_m": A.rep("detection.range_at_pd.snr_limited_range_m", unit="m"),
        "range_actual_m": A.rep("detection.range_at_pd.range_m", unit="m"), "range_limited_by": A.get("detection.range_at_pd.limited_by"),
        "two_paths": None, "paths_delta_db": A.engines_difference_db(),
        "rule": A.get("compliance.rule"), "tier": g["label"],
        "cap_dbm": N(g["peak_eirp_limit_dbm"], g["sources"]["peak_eirp_dbm"], "fmcw_radar report: compliance.tiers[governing].peak_eirp_limit_dbm", "dBm"),
        "cap_kind": "peak EIRP", "backoff_db": A.rep_input("tx_backoff_db", "dB"),
        "off_time_required_s": N(g["off_time_required_s"], g["sources"]["off_time"], "fmcw_radar report: compliance.tiers[governing].off_time_required_s", "s"),
        "off_time_window_s": N(float(S.num("fcc15255.radar_57_64#off_time_window_s").value), "spec", "fcc15255.radar_57_64#off_time_window_s", "s"),
        "off_time_have_s": A.rep("compliance.tiers.5.off_time_have_s", unit="s") if A.get("compliance.tiers.5.tier") == g["tier"] else None,
        "adc_fs_hz": N(fs.value * 1e6, "assumed", "fmcw_radar report: radar.adc_fs_msps", "Hz"), "adc_channels": "1 (real IF)",
        "adc_bits": A.rep("radar.adc_bits", "assumed", "bit"), "adc_bits_what": "as configured",
        "leakage": A.get("leakage"),
        "air_fraction": frac, "air_fraction_upper": None, "rate_mbps": rate,
        "loss_mbps": derive(frac.value * rate.value, "fraction x rate", frac, rate, unit="Mb/s"), "loss_upper_mbps": None,
    }


# ----------------------------------------------------------------------------- the ay rows
class _Bench:
    """One simulator, the paths of each geometry, and the averaged readings."""

    def __init__(self, s: Setup):
        self.s, self.sim = s, W.Simulator(s)
        self.sim_trn = W.Simulator(s, s.trn_field)
        self.frames = int(s.sim["frames_averaged"])
        self.seed = int(s.sim["seed"])

    def sim_for(self, field) -> W.Simulator:
        return self.sim_trn if field is self.s.trn_field else self.sim

    def paths(self, kind: str, link: B.Link, strong_dbm: float, velocity: float | None = None):
        g = B.geometry(self.s, kind)
        v = self.s.velocity_mps.value if velocity is None else velocity
        tgt = W.Path(link.p_rx_dbm.value, g.delay_s, g.path_rate(v), "target")
        strong = W.Path(strong_dbm, self.s.leak_delay_s.value if kind == "monostatic" else 0.0, 0.0, "leakage" if kind == "monostatic" else "direct path")
        return g, tgt, strong

    def options(self, kind: str, **kw) -> W.Options:
        if kind == "bistatic":
            kw = {"random_packet_phase": True, "reference_lag": 0, **kw}
        return W.Options(seed=self.seed, **kw)

    def read(self, paths, opt, near, strong=(0,), sim=None) -> dict:
        sim = sim or self.sim
        rs = [sim.read(sim.frame(paths, opt, k), near, strong) for k in range(self.frames)]
        return {"snr_db": float(np.mean([r["snr_db"] for r in rs])), "peak_dbm": db(float(np.mean([r["peak_mw"] for r in rs]))),
                "floor_dbm": db(float(np.mean([r["floor_mw"] for r in rs]))), "ridge_dbm": db(float(np.mean([r["ridge_mw"] for r in rs]))),
                "first": rs[0]}


def strong_path_dbm(s: Setup, kind: str, eirp: Num, isolation: Num | None = None) -> Num:
    """The static path the receiver cannot avoid: the transmit leakage (monostatic), the direct path (bistatic)."""
    if kind == "monostatic":
        return derive(eirp.value - s.fe.tx_gain_dbi.value - isolation.value, "EIRP - G_t - isolation", eirp, s.fe.tx_gain_dbi, isolation, unit="dBm")
    loss = B.one_way_loss_db(s.wavelength_m, s.bi_baseline_m.value) + s.oxygen_db_per_km.value * s.bi_baseline_m.value / 1e3
    return derive(eirp.value + s.fe.rx_gain_dbi.value - loss, "EIRP + G_r - (4 pi L / lambda)^2 - oxygen over L", eirp, s.fe.rx_gain_dbi, s.bi_baseline_m, s.carrier_hz, unit="dBm")


def row_ay(s: Setup, bench: _Bench, kind: str, reading: R.Reading, pfa: Num, snr_req: Num, field=None) -> dict:
    field = s.field if field is None else field
    on_trn = field is s.trn_field
    fe, sim = s.fe, bench.sim_for(field)
    eirp = reading.eirp_dbm
    link = B.link(s, kind, eirp, field=field)
    g = B.geometry(s, kind)
    lam = s.wavelength_m
    m, t_int = s.n_packets.value, s.interval_s.value
    width_chips = B.width_at_half_power_chips(s.rolloff.value)
    zone_m = field.zero_sidelobe_chips * s.chip_m * 2.0            # of path length
    r_max = zone_m / 2.0 if kind == "monostatic" else (zone_m + g.baseline_m) / 2.0
    iso_best = fe.isolation_db[-1][1]
    strong = strong_path_dbm(s, kind, eirp, iso_best if kind == "monostatic" else None)
    _, tgt, st = bench.paths(kind, link, strong.value)
    near = int(round(g.delay_s * sim.fs))
    meas = bench.read([tgt, st], bench.options(kind), near, sim=sim)
    fr = sim.frame([tgt, st], bench.options(kind), 0)
    i, lag = meas["first"]["doppler_index"], meas["first"]["lag"]
    j = lag - sim.lag_lo
    p_db = _rd.power_db(fr.rd)
    r_read = _rd.read_peak(sim.range_axis(g)[j - 3: j + 4], p_db[i, j - 3: j + 4])["peak_pos"]
    v_read = _rd.read_peak(sim.velocity_axis(g)[i - 3: i + 4], p_db[i - 3: i + 4, j])["peak_pos"]
    det = sim.cfar(fr, pfa.value)
    jc = (j - (-sim.lag_lo) % sim.osr) // sim.osr
    r_snr = B.detection_range_m(s, kind, eirp, snr_req, field)
    air = R.air_time(s, reading, field)
    sim_snr = derive(meas["snr_db"], f"waveform simulation, mean of {bench.frames} frames: peak / floor of the range-Doppler map", eirp, s.n_packets, unit="dB")
    row = {
        "name": f"{'B' if kind == 'monostatic' else 'C'}  ay {'TRN field' if on_trn else 'L-CEF'}, {kind}", "limit": reading.name, "kind": kind,
        "field": "TRN field" if on_trn else "L-CEF", "packet_s": s.packet_seconds(field),
        "bandwidth_hz": s.chip_rate_hz,
        "range_res_m": derive(s.chip_m / g.cos_half_beta, "c / (2 Rc cos(beta/2))", s.chip_rate_hz, s.bi_baseline_m if kind == "bistatic" else s.chip_rate_hz, unit="m"),
        "range_width_m": derive(width_chips * s.chip_m / g.cos_half_beta, "-3 dB width of the raised-cosine pulse", s.chip_rate_hz, s.rolloff, unit="m"),
        "range_window": "none (the pulse)",
        "range_max_m": derive(r_max, "cyclic prefix x c / Rc" + (" / 2" if kind == "monostatic" else ", as distance from both devices"),
                              *((s.chip_rate_hz,) if kind == "monostatic" else (s.chip_rate_hz, s.bi_baseline_m)), unit="m"),
        "range_max_why": f"the {field.zero_sidelobe_chips}-chip " + ("guard block either side of the correlated blocks" if on_trn else "cyclic prefix"),
        "velocity_res_mps": derive(lam / (2.0 * m * t_int * g.cos_half_beta), "lambda / (2 M T cos(beta/2))", s.n_packets, s.interval_s, s.carrier_hz, unit="m/s"),
        "velocity_max_mps": derive(lam / (4.0 * t_int * g.cos_half_beta), "lambda / (4 T cos(beta/2))", s.interval_s, s.carrier_hz, unit="m/s"),
        "eirp_dbm": eirp,
        "processing_gain_db": link.processing_gain_db, "processing_what": f"{field.n_correlated} chips x {m:.0f} packets",
        "window_loss_db": link.window_loss_db,
        "snr_db": link.snr_cell_db, "pd": derive(B.pd_at(link.snr_cell_db.value, pfa.value), "Pd(SNR, Pfa), Swerling 0", link.snr_cell_db, pfa), "pfa": pfa,
        "range_pd_m": r_snr,
        "two_paths": {"closed_grid_db": link.snr_grid_db, "sim_db": sim_snr, "straddle_db": link.straddle_db,
                      "range_read_m": r_read, "velocity_read_mps": v_read, "cfar_target": bool(det[i, jc]),
                      "strong_dbm": strong, "strong_what": st.name},
        "paths_delta_db": derive(meas["snr_db"] - link.snr_grid_db.value, "waveform simulation - closed form on the grid", sim_snr, link.snr_grid_db, unit="dB"),
        "rule": reading.rule, "tier": reading.tier, "cap_dbm": reading.cap_dbm, "cap_kind": reading.cap_kind, "backoff_db": reading.backoff_db,
        "off_time_required_s": reading.off_time_required_s, "off_time_window_s": reading.off_time_window_s,
        "off_time_have_s": None if reading.off_time_have_s is None else derive(reading.off_time_have_s, "least off-time in any window (fmcw_radar's check)", s.n_packets, s.interval_s, s.frames_per_s, unit="s"),
        "adc_fs_hz": s.adc_fs_hz, "adc_channels": "2 (I and Q)", "adc_bits": s.adc_enob, "adc_bits_what": "ENOB of mmw60_trx's converter model at this rate",
        "air_fraction": air.fraction, "air_fraction_upper": air.fraction_upper, "rate_mbps": air.rate_mbps,
        "loss_mbps": air.loss_mbps, "loss_upper_mbps": air.loss_upper_mbps,
        "link": link,
    }
    return row


# ----------------------------------------------------------------------------- the receiver's cost, both roads
def _crossing(x: np.ndarray, y: np.ndarray, level: float) -> float | None:
    """Where the rising or falling curve ``y(x)`` first crosses ``level`` (linear between the points)."""
    for k in range(len(x) - 1):
        if (y[k] - level) * (y[k + 1] - level) <= 0.0 and y[k] != y[k + 1]:
            return float(x[k] + (level - y[k]) * (x[k + 1] - x[k]) / (y[k + 1] - y[k]))
    return None


def compression_sweep(s: Setup, bench: _Bench, link: B.Link) -> dict:
    """The receiver's compression against the level of the leakage, on both roads: the echo's SNR (closed form:
    fmcw_radar's small-signal theory; simulation: noisy frames, against a linear receiver) and the image it throws
    at the mirrored Doppler (closed form: an upper bound; simulation: one noise-free frame per level)."""
    sim, env = bench.sim, bench.sim.unit_envelope()
    p1 = s.fe.p1db_inband_dbm.value
    span, step = float(s.sim["compression_span_db"]), float(s.sim["compression_step_db"])
    levels = p1 + np.arange(-span, span + step / 2, step)
    g, tgt, _ = bench.paths("monostatic", link, 0.0)
    near = int(round(g.delay_s * sim.fs))
    leak = lambda lvl: W.Path(float(lvl), s.leak_delay_s.value, 0.0, "leakage")      # noqa: E731
    lin = bench.read([tgt, leak(levels[0])], bench.options("monostatic"), near)
    rows = []
    for lvl in levels:
        c = B.compression(env, float(lvl), s)
        nl = bench.read([tgt, leak(lvl)], bench.options("monostatic", rx_p1db_dbm=p1), near)
        quiet = sim.read(sim.frame([tgt, leak(lvl)], bench.options("monostatic", rx_p1db_dbm=p1, noise=False), 0), near, (0,))
        rows.append({"level_dbm": float(lvl), "snr_closed_db": c["snr"], "snr_sim_db": nl["snr_db"] - lin["snr_db"],
                     "signal_closed_db": c["signal"], "signal_sim_db": nl["peak_dbm"] - lin["peak_dbm"],
                     "image_closed_dbc": c["image"] - c["signal"], "image_sim_dbc": db(quiet["image_mw"] / quiet["peak_mw"]),
                     "row_sidelobe_sim_dbc": db(quiet["row_sidelobe_mw"] / quiet["peak_mw"])})
    x = np.array([r["level_dbm"] for r in rows])
    loss = s.desense_allowed_db.value
    return {"rows": rows, "env": env,
            "level_snr_closed_dbm": B.strong_level_for_snr_loss(env, s, loss),
            "level_snr_sim_dbm": _crossing(x, np.array([r["snr_sim_db"] for r in rows]), -loss),
            "snr_limit_closed_db": B.compression(env, p1 + 60.0, s)["snr"],
            "level_snr_closed_by_p": [(p, B.strong_level_for_snr_loss(env, s, loss, p)) for p in (1.0, 3.0)]}


def monostatic_block(s: Setup, bench: _Bench, reading: R.Reading, sweep: dict) -> dict:
    eirp = reading.eirp_dbm
    link = B.link(s, "monostatic", eirp)
    cost = B.monostatic_cost(s, eirp, sweep["env"])
    near = int(round(B.geometry(s, "monostatic").delay_s * bench.sim.fs))
    p1 = s.fe.p1db_inband_dbm.value
    checks = []
    for label, iso, short, leak, enob_need in cost.have:
        _, tgt, st = bench.paths("monostatic", link, leak.value)
        lin = bench.read([tgt, st], bench.options("monostatic"), near)
        adc = bench.read([tgt, st], bench.options("monostatic", adc_enob=s.adc_enob.value), near)
        pn = bench.read([tgt, st], bench.options("monostatic", lo="shared"), near)
        q = bench.sim.frame([tgt, st], bench.options("monostatic", adc_enob=s.adc_enob.value), 0).quant
        comp = bench.read([tgt, st], bench.options("monostatic", rx_p1db_dbm=p1), near)
        quiet = bench.sim.read(bench.sim.frame([tgt, st], bench.options("monostatic", rx_p1db_dbm=p1, noise=False), 0), near, (0,))
        c = B.compression(sweep["env"], leak.value, s)
        checks.append({
            "noise_lsb": q["noise_lsb"], "repeating_share": q["repeating_share"],
            "adc_rise_varying_db": B.floor_rise_db(db(q["varying_mw"]) - db(s.adc_osr.value), link.noise_chip_dbm.value),
            "label": label, "isolation_db": iso, "shortfall_db": short, "leak_dbm": leak, "enob_needed": enob_need,
            "snr_linear_db": lin["snr_db"],
            "adc_rise_closed_db": B.floor_rise_db(B.adc_noise_inband_dbm(leak.value, s.adc_enob.value, s), link.noise_chip_dbm.value),
            "adc_rise_sim_db": adc["floor_dbm"] - lin["floor_dbm"], "adc_ridge_dbm": adc["ridge_dbm"], "adc_peak_change_db": adc["peak_dbm"] - lin["peak_dbm"],
            "pn_rise_closed_db": B.floor_rise_db(leak.value + cost.pn_residual_dbc.value, link.noise_chip_dbm.value),
            "pn_rise_sim_db": pn["floor_dbm"] - lin["floor_dbm"],
            "ridge_linear_dbm": lin["ridge_dbm"], "peak_dbm": lin["peak_dbm"], "floor_dbm": lin["floor_dbm"],
            "comp_snr_closed_db": c["snr"], "comp_snr_sim_db": comp["snr_db"] - lin["snr_db"],
            "comp_signal_closed_db": c["signal"], "comp_signal_sim_db": comp["peak_dbm"] - lin["peak_dbm"],
            "comp_image_closed_dbc": c["image"] - c["signal"], "comp_image_sim_dbc": db(quiet["image_mw"] / quiet["peak_mw"]),
            "comp_row_sidelobe_sim_dbc": db(quiet["row_sidelobe_mw"] / quiet["peak_mw"]), "comp_ridge_dbm": comp["ridge_dbm"]})
    _, tgt, _ = bench.paths("monostatic", link, 0.0)
    alone = bench.read([tgt], bench.options("monostatic"), near, strong=())
    above = cost.compression_level_dbm.value - p1                 # the criterion's level over the 1 dB point: the same at every gain setting
    gains = []
    for gdb, nf, p1g in zip(S.raw("rx60.radar_iwr6843_fig7_2#gain_db"), S.raw("rx60.radar_iwr6843_fig7_2#nf_db"), S.raw("rx60.radar_iwr6843_fig7_2#p1db_inband_dbm")):
        gains.append((float(gdb), float(nf), float(p1g), cost.conducted_dbm.value - (float(p1g) + above)))
    cond = cost.conducted_dbm.value
    # the limiter's image of this row's target at the level of the map's floor: where the two roads put the leakage
    floor_dbc = -link.snr_grid_db.value
    x = np.array([r["level_dbm"] for r in sweep["rows"]])
    lvl_sim = _crossing(x, np.array([r["image_sim_dbc"] for r in sweep["rows"]]), floor_dbc)
    ghost = {"floor_dbc": floor_dbc, "need_closed_db": cond - B.strong_level_for_image(sweep["env"], s, floor_dbc),
             "need_sim_db": None if lvl_sim is None else cond - lvl_sim}
    return {"reading": reading.name, "cost": cost, "checks": checks, "snr_no_leak_db": alone["snr_db"], "gains": gains,
            "need_compression_sim_db": None if sweep["level_snr_sim_dbm"] is None else cond - sweep["level_snr_sim_dbm"],
            "need_by_p": [(p, None if lvl is None else cond - lvl) for p, lvl in sweep["level_snr_closed_by_p"]], "ghost": ghost}


def bistatic_block(s: Setup, bench: _Bench, reading: R.Reading) -> dict:
    eirp = reading.eirp_dbm
    link = B.link(s, "bistatic", eirp)
    direct = strong_path_dbm(s, "bistatic", eirp)
    g, tgt, st = bench.paths("bistatic", link, direct.value)
    near = int(round(g.delay_s * bench.sim.fs))
    lin = bench.read([tgt, st], bench.options("bistatic"), near)
    pn = bench.read([tgt, st], bench.options("bistatic", lo="independent"), near)
    res = B.pn_residual_independent(s, s.field.n_correlated / s.chip_rate_hz.value)
    comp = B.compression(bench.sim.unit_envelope(), direct.value, s)
    return {"reading": reading.name, "direct_dbm": direct, "comp_snr_closed_db": comp["snr"], "comp_image_closed_dbc": comp["image"] - comp["signal"],
            "enob_needed": derive(B.enob_required(direct.value, link.noise_chip_dbm.value, s), "ENOB for the quantisation noise to cost the allowed desensitisation",
                                  direct, link.noise_chip_dbm, s.adc_backoff_db, s.adc_osr, s.desense_allowed_db, unit="bit"),
            "under_p1db_db": derive(s.fe.p1db_inband_dbm.value - direct.value, "P1dB(in band) - direct path", s.fe.p1db_inband_dbm, direct, unit="dB"),
            "pn_residual_dbc": db(res),
            "pn_rise_closed_db": B.floor_rise_db(direct.value + db(res), link.noise_chip_dbm.value),
            "pn_rise_sim_db": pn["floor_dbm"] - lin["floor_dbm"], "pn_snr_change_db": pn["snr_db"] - lin["snr_db"]}


def golay_under_doppler(s: Setup, bench: _Bench, velocities, sim=None, span_chips: int | None = None) -> list:
    """Peak sidelobe of the summed correlations at the chip-spaced lags of the zero-sidelobe zone, one packet, no
    noise, the target on a chip: what the pair's complementarity is worth when the target moves.  The last column
    is the carrier's turn over ``span_chips`` (default: one correlation window of the simulator's field).  A level
    under the arithmetic's rounding (``Simulator.rounding_dbc``) is reported as that bound: at rest the pair cancels
    exactly and what is left is not a reading."""
    sim, out = sim or bench.sim, []
    span = len(sim.field.windows[0][1]) if span_chips is None else span_chips
    delay_chips = round(B.geometry(s, "monostatic").delay_s * s.chip_rate_hz.value)
    for v in velocities:
        fr = sim.frame([W.Path(0.0, delay_chips / s.chip_rate_hz.value, 2.0 * v)], W.Options(noise=False, n_packets=1))
        h = np.abs(fr.cir[0]) ** 2
        chip = h[fr.lag0 % sim.osr:: sim.osr]
        lags = sim.lags[fr.lag0 % sim.osr:: sim.osr] // sim.osr
        zone = (lags >= 0) & (lags <= sim.field.zero_sidelobe_chips) & (lags != delay_chips)
        peak = float(chip[lags == delay_chips][0])
        turn = 2.0 * math.pi * (2.0 * v / s.wavelength_m) * span / s.chip_rate_hz.value
        out.append((float(v), max(db(float(chip[zone].max()) / peak + 1e-40), sim.rounding_dbc), db(peak), turn))
    return out


def pulse_sidelobes(s: Setup, bench: _Bench) -> list:
    """The range response between the chips: the raised-cosine pulse's own sidelobes at 1.5, 2.5, 3.5 chips."""
    sim = bench.sim
    delay_chips = round(B.geometry(s, "monostatic").delay_s * s.chip_rate_hz.value)
    fr = sim.frame([W.Path(0.0, delay_chips / s.chip_rate_hz.value, 0.0)], W.Options(noise=False, n_packets=1))
    h = np.abs(fr.cir[0]) ** 2
    j = delay_chips * sim.osr - sim.lag_lo
    return [(k / sim.osr, db(float(h[j + k]) / float(h[j]))) for k in (3, 5, 7)]


def measured_width_chips(s: Setup, bench: _Bench) -> float:
    """-3 dB width of the range response, measured: one lag's reading as a noise-free target is moved past it."""
    sim = bench.sim
    delay_chips = round(B.geometry(s, "monostatic").delay_s * s.chip_rate_hz.value)
    j = delay_chips * sim.osr - sim.lag_lo
    offs = np.linspace(-0.75, 0.75, 97)
    resp = np.array([abs(sim.frame([W.Path(0.0, (delay_chips + d) / s.chip_rate_hz.value, 0.0)], W.Options(noise=False, n_packets=1)).cir[0][j]) ** 2
                     for d in offs])
    y = 10.0 * np.log10(resp / resp.max())
    left = np.interp(-3.0103, y[: offs.size // 2 + 1], offs[: offs.size // 2 + 1])
    right = np.interp(-3.0103, y[offs.size // 2:][::-1], offs[offs.size // 2:][::-1])
    return float(right - left)


def _all_blocks(s: Setup):
    """The TRN field with every block of the subfield correlated, guards included: what the choice of blocks 2-5 avoids."""
    from . import field as F
    return F.trn(int(s.trn_units.value), tuple(range(1, len(S.raw("edmg.trn_field#basic_subfield")) + 1)))


def trn_sweep(s: Setup, reading: R.Reading, snr_a_db: float) -> list:
    """Closed form, monostatic under the communication limit: (TRN units, correlated chips, SNR, SNR over row A's,
    packet duration, share of the air time) for 1 .. 8 units, the blocks correlated as in the table."""
    from . import field as F
    blocks = tuple(int(b) for b in S.raw("jcas60.trn_sensing#correlated_blocks"))
    base = B.link(s, "monostatic", reading.eirp_dbm)
    per_packet = S.num("jcas60.packet_overhead#sifs_per_packet").value * s.sifs_s.value
    out = []
    for units in range(1, 9):
        f = F.trn(units, blocks)
        snr = base.snr_cell_db.value + db(f.n_correlated / s.field.n_correlated)
        t_packet = (s.packet_chips.value + f.chips.size) / s.chip_rate_hz.value
        out.append((units, f.n_correlated, snr, snr - snr_a_db, t_packet, (t_packet + per_packet) * s.n_packets.value * s.frames_per_s.value))
    return out


def sensitivity(s: Setup, reading: R.Reading, pfa: float) -> list:
    """Closed form, monostatic under the communication limit, one choice changed at a time: (what, d SNR, note)."""
    base = B.link(s, "monostatic", reading.eirp_dbm).snr_cell_db.value
    stf = S.num("edmg.sc_siso_preamble#l_stf_chips").value
    comm = {k: S.num(f"rx60.comm_asada2011#{k}").value for k in ("antenna_gain_dbi", "nf_db", "tx_p1db_dbm")}
    d_ant = 2.0 * (comm["antenna_gain_dbi"] - s.fe.tx_gain_dbi.value)
    d_nf = s.fe.nf_db.value - comm["nf_db"]
    d_tx = comm["tx_p1db_dbm"] - s.fe.tx_power_dbm.value
    rcs_hi = S.num("target_rcs.deep2019_pedestrian_77g#co_pol_max_dbsm").value - s.rcs_dbsm.value
    rows = [
        ("packets per frame x 2 (jcas60.sensing_burst)", db(2.0), "velocity resolution twice as fine, air time x 2"),
        ("rectangular Doppler window (jcas60.doppler_window)", B.link(s, "monostatic", reading.eirp_dbm).window_loss_db.value, "Doppler sidelobes -13 dB instead of -31 dB"),
        ("the L-STF correlated as well (jcas60.sensing_field)", db((stf + s.field.n_correlated) / s.field.n_correlated), "its composition has no public source"),
        ("RCS at the top of the pedestrian range (target_rcs.deep2019_pedestrian_77g#co_pol_max_dbsm)", rcs_hi, "all three rows move together"),
        ("the communication transceiver's numbers (rx60.comm_asada2011: antenna, NF, TX P1dB)", d_ant + d_nf + d_tx,
         f"antennas {d_ant:+.1f} dB, NF {d_nf:+.1f} dB, TX power {d_tx:+.1f} dB; row A is not recomputed on them"),
    ]
    return [(what, d, base + d, note) for what, d, note in rows]


@lru_cache(maxsize=1)
def build() -> dict:
    S.TRACE.clear()
    for cached in (setup,):
        cached.cache_clear()
    from . import field as _field
    from . import frontend as _frontend
    _frontend.shared.cache_clear()
    _field.l_cef.cache_clear()
    s = setup()
    bench = _Bench(s)
    a = row_a()
    comm, radar = R.readings(s)
    comm_t, radar_t = R.readings(s, s.trn_field)                 # the same limits; the burst is longer by the TRN field
    rows = [a] + [row_ay(s, bench, kind, rd_, a["pfa"], a["snr_required_db"]) for kind in ("monostatic", "bistatic") for rd_ in (comm, radar)]
    rows += [row_ay(s, bench, kind, rd_, a["pfa"], a["snr_required_db"], s.trn_field) for kind in ("monostatic", "bistatic") for rd_ in (comm_t, radar_t)]
    v_max = rows[1]["velocity_max_mps"].value
    return {
        "setup": s, "rows": rows, "readings": (comm, radar),
        "compression": (sweep := compression_sweep(s, bench, rows[1]["link"])),
        "monostatic": [monostatic_block(s, bench, r, sweep) for r in (comm, radar)],
        "bistatic": [bistatic_block(s, bench, r) for r in (comm, radar)],
        "golay": golay_under_doppler(s, bench, (0.0, s.velocity_mps.value, round(v_max, 1), 50.0, 200.0)),
        "golay_trn": golay_under_doppler(s, bench, (0.0, s.velocity_mps.value, round(v_max, 1), 50.0, 200.0), bench.sim_trn, s.trn_field.chips.size),
        "rounding_dbc": {"l_cef": bench.sim.rounding_dbc, "trn": bench.sim_trn.rounding_dbc},
        "golay_trn_all_blocks": golay_under_doppler(s, bench, (0.0,), W.Simulator(s, _all_blocks(s)))[0][1],
        "trn_sweep": trn_sweep(s, comm, a["snr_db"].value),
        "pulse_sidelobes": pulse_sidelobes(s, bench), "width_measured_chips": measured_width_chips(s, bench),
        "sensitivity": sensitivity(s, comm, a["pfa"].value),
        "bf": bf.checks(s),
        "snr_required_db": a["snr_required_db"].value,
        "snr_required_known_noise_db": B.snr_required_db(s.pd.value, a["pfa"].value),
        "range_known_noise_m": B.detection_range_m(s, "monostatic", comm.eirp_dbm, derive(B.snr_required_db(s.pd.value, a["pfa"].value),
                                                   "SNR for Pd at Pfa, noise level known", s.pd, a["pfa"], unit="dB")).value,
        "trace": dict(S.TRACE),
    }


# ----------------------------------------------------------------------------- text
def _n(x: Num | None, fmt: str, scale: float = 1.0) -> str:
    return "—" if x is None else f"{format(x.value * scale, fmt)} {x.label}"


def _short(row: dict) -> str:
    return f"{row['name'][:1]} {row['kind'][:4]} / {'comm' if row['limit'].startswith('comm') else 'radar'} limit" if row["name"][0] != "A" else "A FMCW / radar limit"


def render(t: dict | None = None) -> str:
    t = t or build()
    s: Setup = t["setup"]
    rows = t["rows"]
    a = rows[0]
    fe = s.fe
    L = []
    w = L.append
    w(ENGINE_BANNER)
    w(f"jcas60 {__version__}; milestone 0: {WORDING}")
    w(f"one target for every row: a person at {s.range_m.show('g')} m receding at {s.velocity_mps.show('g')} m/s, RCS {s.rcs_dbsm.show('g')} dBsm; "
      f"detection: Pd {s.pd.show('g')} at Pfa {a['pfa'].show('g')}, Swerling 0, at the cell SNR fmcw_radar's report says its map's CA-CFAR needs, "
      f"{a['snr_required_db'].show('.2f')} dB, for every row (with the noise level known: {t['snr_required_known_noise_db']:.2f} dB)")
    w("tags: [spec] the source's own text, [survey] not verified, [assumed] this library's or a sibling's choice, [derived] computed from "
      "spec numbers, [derived/assumed] computed from at least one assumed number")
    w("")
    w("== 0. the shared front-end numbers (one instance; row A's own inputs, from its report, beside them) ==")
    pairs = [("TX power, conducted [dBm]", fe.tx_power_dbm, A.rep_input("tx_power_dbm")), ("TX antenna gain [dBi]", fe.tx_gain_dbi, A.rep_input("tx_gain_dbi")),
             ("RX antenna gain [dBi]", fe.rx_gain_dbi, A.rep_input("rx_gain_dbi")), ("RX noise figure [dB]", fe.nf_db, A.rep_input("nf_db")),
             ("RX gain setting [dB]", fe.rx_gain_db, A.rep_input("rx_gain_db")), ("RX in-band P1dB at that gain [dBm]", fe.p1db_inband_dbm, A.rep_input("p1db_inband_dbm"))]
    for label, mine, theirs in pairs:
        w(f"  {label:38s} {mine.show('g'):>16s}   row A's report: {theirs.show('g'):<16s} {'equal' if mine.value == theirs.value else 'DIFFERENT'}   ({mine.ref})")
    w(f"  {'PA back-off':38s} one rule for every row: max(0, TX power + antenna gain - the row's EIRP limit); row A's report has {a['backoff_db'].show('g')} dB, the rule gives "
      f"{fe.backoff_db(a['cap_dbm']).value:g} dB")
    for label, iso in fe.isolation_db:
        w(f"  {'TX-RX isolation [dB]':38s} {iso.show('g'):>16s}   {label} ({iso.ref})")
    w(f"  {'LO phase noise':38s} ay rows: PSD0 {fe.pn_psd0_dbchz.show('g')} dBc/Hz, pole {fe.pn_f_pole_hz.value / 1e6:g} MHz, zero {fe.pn_f_zero_hz.value / 1e6:g} MHz "
      f"({fe.pn_dbchz(1e6):.1f} dBc/Hz at 1 MHz; {fe.pn_psd0_dbchz.ref.partition('#')[0]}). NOT the same in row A: its LO is fmcw_radar's chirping ADPLL "
      f"({A.get('radar.lo_dco_ghz'):g} GHz x {A.get('radar.lo_mult')}), the one front-end number the rows do not share")
    w("")
    w("== 1. resolution and ambiguity ==")
    w(f"{'row':34s} {'bandwidth':>22s} {'c/2B':>22s} {'-3 dB width as processed':>42s} {'clean range':>24s} {'velocity res.':>24s} {'unambiguous v':>26s}")
    for r in (rows[0], rows[1], rows[3]):
        w(f"{r['name']:34s} {_n(r['bandwidth_hz'], '.2f', 1e-9) + ' GHz':>22s} {_n(r['range_res_m'], '.2f', 100) + ' cm':>22s} "
          f"{_n(r['range_width_m'], '.2f', 100) + ' cm, ' + r['range_window']:>42s} {_n(r['range_max_m'], '.2f') + ' m':>24s} "
          f"{_n(r['velocity_res_mps'], '.3f') + ' m/s':>24s} {'+-' + _n(r['velocity_max_mps'], '.2f') + ' m/s':>26s}")
    wm = t["width_measured_chips"]
    w(f"  clean range: A by {a['range_max_why']}; the ay rows by {rows[1]['range_max_why']} (beyond it a window no longer holds a cyclic shift of its sequence; "
      f"the sequences themselves repeat after {len(s.field.windows[0][1])} chips = {len(s.field.windows[0][1]) * s.chip_m:.1f} m)")
    w(f"  ay -3 dB width: closed form {B.width_at_half_power_chips(s.rolloff.value):.4f} chips, waveform simulation {wm:.4f} chips ({wm * s.chip_m * 100:.2f} cm monostatic); "
      f"roll-off {s.rolloff.show('g')} (assumed.ay_sc_pulse, mmw60_trx)")
    w(f"  range sidelobes: A {a['range_sidelobe_db'].value:.1f} dB ({a['range_window']} window); ay, between the chips, the pulse's own: "
      + ", ".join(f"{lvl:.1f} dB at {off:g} chips" for off, lvl in t["pulse_sidelobes"]) + " (waveform simulation; no range window is applied)")
    w(f"  the TRN rows (B and C sensed on an EDMG TRN-R field of {s.trn_units.show('g')} units = {s.trn_field.chips.size} chips behind the same packet) have these same figures: the field changes the gain, "
      f"not the pulse or the burst; their clean range is by {rows[5]['range_max_why']} ({rows[5]['range_max_m'].value:.2f} m)")
    w(f"  bistatic (C): legs {s.bi_r_tx_m.value:g} m and {s.bi_r_rx_m.value:g} m, baseline {s.bi_baseline_m.value:g} m, bistatic angle "
      f"{2 * math.degrees(math.acos(B.geometry(s, 'bistatic').cos_half_beta)):.0f} deg; resolutions are along the bisector (jcas60.bistatic_geometry)")
    w("")
    w("== 2. link: EIRP, processing gain, SNR of the person at 5 m, range at Pd 0.9 ==")
    w(f"{'row':34s} {'limit':>20s} {'EIRP [dBm]':>22s} {'processing gain':>40s} {'window loss':>20s} {'SNR(cell) [dB]':>24s} {'Pd':>8s} {'range at Pd 0.9 (SNR)':>30s} {'within the clean range':>24s}")
    for r in rows:
        inside = r["range_pd_m"].value <= r["range_max_m"].value
        w(f"{r['name']:34s} {r['limit']:>20s} {_n(r['eirp_dbm'], '.2f'):>22s} {_n(r['processing_gain_db'], '.2f') + ' dB = ' + r['processing_what']:>40s} "
          f"{_n(r['window_loss_db'], '.2f') + ' dB':>20s} {_n(r['snr_db'], '.2f'):>24s} {r['pd'].value:8.4f} {_n(r['range_pd_m'], '.2f') + ' m':>30s} "
          f"{'yes' if inside else 'no: ' + format(r['range_max_m'].value, '.2f') + ' m':>24s}")
    w(f"  SNR(cell): thermal noise only, at the peak (no straddle).  Row A's range is its report's detection.range_at_pd: {a['range_pd_m'].value:.2f} m by SNR, "
      f"{a['range_actual_m'].value:.2f} m as the report concludes (limited by: {a['range_limited_by']})")
    w(f"  the ay rows are held to the same required SNR; what a CA-CFAR on their own map would need is not computed (jcas60.detection): with the noise level known the ranges are "
      f"{(t['range_known_noise_m'] / rows[1]['range_pd_m'].value - 1.0) * 100:.1f} % longer (row B; {t['range_known_noise_m']:.2f} m)")
    for r in (rows[1], rows[2]):
        gap = a["snr_db"].value - r["snr_db"].value
        chips = 10.0 ** (gap / 10.0) * s.field.n_correlated * s.n_packets.value
        w(f"  row A over row B, {r['limit']}: {gap:+.2f} dB with {r['eirp_dbm'].value - a['eirp_dbm'].value:g} dB less EIRP: its frame is on the target for {A.get('compliance.frame_on_s') * 1e3:g} ms, "
          f"the ay frame's correlated chips for {s.field.n_correlated * s.n_packets.value / s.chip_rate_hz.value * 1e6:.1f} us; equal SNR takes x{10.0 ** (gap / 10.0):.1f} the chips, "
          f"{chips / s.chip_rate_hz.value * 1e3:.2f} ms a frame ({chips / s.chip_rate_hz.value * s.frames_per_s.value * 100:.1f} % of the time, training alone)")
    for r in (rows[5], rows[6]):
        w(f"  row B on the TRN field against row A, {r['limit']}: {r['snr_db'].value - a['snr_db'].value:+.2f} dB with {r['eirp_dbm'].value - a['eirp_dbm'].value:g} dB more EIRP: "
          f"{s.trn_field.n_correlated} chips a packet, {s.trn_field.n_correlated * s.n_packets.value / s.chip_rate_hz.value * 1e3:.3f} ms of correlated chips a frame")
    b1 = rows[1]["link"]
    w(f"  the ay budget, monostatic at {rows[1]['eirp_dbm'].value:g} dBm: spreading {b1.spreading_db.value:.2f} dB, oxygen {b1.oxygen_db.value:.2f} dB, echo {b1.p_rx_dbm.value:.2f} dBm, "
      f"noise in {s.chip_rate_hz.value / 1e9:g} GHz {b1.noise_chip_dbm.value:.2f} dBm, SNR per chip {b1.snr_chip_db.value:.2f} dB")
    w("")
    w("== 3. the two roads (law 4): closed form against waveform simulation ==")
    w(f"{'row':34s} {'limit':>20s} {'closed, on the grid':>22s} {'simulation':>22s} {'difference':>12s} {'straddle':>10s} {'read at':>26s} {'CA-CFAR':>9s}  static path in the simulation")
    for r in rows:
        tp = r["two_paths"]
        if tp is None:
            w(f"{r['name']:34s} {r['limit']:>20s} {'—':>22s} {'—':>22s} {r['paths_delta_db'].value:+10.2f} dB {'':>10s} {'':>26s} {'':>9s}  "
              f"fmcw_radar's own two engines (time - statistical), from its example 00 as printed")
            continue
        w(f"{r['name']:34s} {r['limit']:>20s} {tp['closed_grid_db'].value:19.2f} dB {tp['sim_db'].value:19.2f} dB {r['paths_delta_db'].value:+10.2f} dB "
          f"{tp['straddle_db'].value:7.2f} dB {format(tp['range_read_m'], '.3f') + ' m, ' + format(tp['velocity_read_mps'], '+.2f') + ' m/s':>26s} "
          f"{'detected' if tp['cfar_target'] else 'missed':>9s}  {tp['strong_what']} at {tp['strong_dbm'].value:.1f} dBm")
    w(f"  simulation: mean of {s.sim['frames_averaged']} frames; the floor is read off the zero-Doppler ridge of the static path; the ADC and the LO phase noise are off here (section 5 turns them on)")
    w("")
    w("== 4. the limit each row is held to (law 3: each at its own; the ay rows under both readings of the rule) ==")
    w(f"{'row':34s} {'limit':>20s} {'rule':>24s} {'EIRP limit [dBm]':>26s} {'back-off [dB]':>22s} {'off-time required':>30s} {'off-time the frame leaves':>26s}  tier")
    for r in rows:
        need = "none in the rule" if r["off_time_required_s"] is None else f"{r['off_time_required_s'].value * 1e3:g} ms in any {r['off_time_window_s'].value * 1e3:g} ms {r['off_time_required_s'].label}"
        have = "—" if r["off_time_have_s"] is None else f"{r['off_time_have_s'].value * 1e3:.2f} ms"
        w(f"{r['name']:34s} {r['limit']:>20s} {r['rule']:>24s} {_n(r['cap_dbm'], 'g') + ' ' + r['cap_kind'].split()[0]:>26s} {_n(r['backoff_db'], 'g'):>22s} {need:>30s} {have:>26s}  {r['tier']}")
    w(f"  the front end's full power is {fe.eirp_max_dbm.show('g')} dBm EIRP: one antenna does not reach the communication limit; the two readings differ by "
      f"{rows[1]['eirp_dbm'].value - rows[2]['eirp_dbm'].value:g} dB of EIRP and by the off-time")
    w(f"  radar reading of the ay rows: fmcw_radar's own check (vendored) on the {s.channel_width_hz.value / 1e9:g} GHz channel "
      f"{s.band_hz[0] / 1e9:g}-{s.band_hz[1] / 1e9:g} GHz, a burst of {s.n_packets.value:g} packets of {s.packet_s.value * 1e6:.3f} us every {s.interval_s.value * 1e6:g} us, {s.frames_per_s.value:g} frames a second")
    w("")
    w("== 5. what the receiver pays ==")
    w(f"{'row':34s} {'ADC rate':>26s} {'channels':>14s} {'bits':>26s}   what the bits are")
    for r in (rows[0], rows[1], rows[3]):
        w(f"{r['name']:34s} {_n(r['adc_fs_hz'], 'g', 1e-6) + ' MS/s':>26s} {r['adc_channels']:>14s} {_n(r['adc_bits'], '.2f'):>26s}   {r['adc_bits_what']}")
    w(f"  A: its leakage budget is {'not in its default report (the model is off in that config; its example 08 turns it on)' if a['leakage'] is None else 'in its report'}; "
      "an FMCW radar's leakage beats to near zero frequency and is filtered ahead of the ADC")
    for blk in t["monostatic"]:
        c = blk["cost"]
        w(f"  B, {blk['reading']} (conducted {c.conducted_dbm.value:g} dBm): isolation needed for the leakage to cost at most {s.desense_allowed_db.value:g} dB of sensitivity")
        gh = blk["ghost"]
        w(f"      in-band compression  {c.need_compression_db.value:6.1f} dB {c.need_compression_db.label}  (closed form; simulated {blk['need_compression_sim_db']:.1f} dB: the leakage at {c.compression_level_dbm.value:.1f} dBm, "
          f"{c.compression_level_dbm.value - fe.p1db_inband_dbm.value:+.1f} dB from the in-band P1dB; {c.need_compression_at_p1db_db.value:.1f} dB puts it at P1dB itself)")
        w(f"      ADC quantisation     {c.need_adc_db.value:6.1f} dB {c.need_adc_db.label}  (ENOB {s.adc_enob.value:.2f} at {s.adc_fs_hz.value / 1e9:g} GS/s, {s.adc_backoff_db.value:g} dB under full scale)")
        w(f"      LO phase noise       {c.need_pn_db.value:6.1f} dB {c.need_pn_db.label}  (shared LO, leakage delayed {s.leak_delay_s.value * 1e9:g} ns: residual {c.pn_residual_dbc.value:.1f} dBc after range correlation)")
        w(f"      -> {c.need_db.value:.1f} dB, set by the {c.binding}")
        w(f"      besides sensitivity, the compression mirrors the target to minus its Doppler: for that image to stay under the map's floor for this target ({gh['floor_dbc']:.1f} dBc) the isolation is "
          f"{gh['need_closed_db']:.1f} dB by the closed form (an upper bound on the image) / {gh['need_sim_db']:.1f} dB simulated")
        for ck in blk["checks"]:
            w(f"      have {ck['isolation_db'].show('g')} dB, {ck['label']}: short by {ck['shortfall_db'].value:.1f} dB; leakage {ck['leak_dbm'].value:g} dBm; ENOB needed {ck['enob_needed'].value:.2f}")
            w(f"          both roads, floor rise: ADC closed {ck['adc_rise_closed_db']:.2f} dB / simulated {ck['adc_rise_sim_db']:.2f} dB; LO phase noise closed {ck['pn_rise_closed_db']:.2f} dB / simulated {ck['pn_rise_sim_db']:.2f} dB")
            w(f"          the ADC's error: thermal noise is {ck['noise_lsb']:.2f} LSB rms, {ck['repeating_share'] * 100:.0f} % of the error is the same in every packet (the closed form takes all of it for noise); "
              f"the part that is not gives a floor rise of {ck['adc_rise_varying_db']:.2f} dB, the rest is the zero-Doppler ridge")
            w(f"          compression (limiter, p = {s.rx_rapp_p.value:g}): SNR closed {ck['comp_snr_closed_db']:+.2f} dB / simulated {ck['comp_snr_sim_db']:+.2f} dB (the echo itself {ck['comp_signal_closed_db']:+.2f} / {ck['comp_signal_sim_db']:+.2f} dB); "
              f"image at the mirrored Doppler {ck['comp_image_sim_dbc']:.1f} dBc (closed-form bound {ck['comp_image_closed_dbc']:.1f} dBc); the target's range sidelobes {ck['comp_row_sidelobe_sim_dbc']:.1f} dBc; ridge {ck['comp_ridge_dbm']:.1f} dBm")
            w(f"          simulation, linear receiver and ideal ADC: SNR {ck['snr_linear_db']:.2f} dB with the leakage, {blk['snr_no_leak_db']:.2f} dB without; "
              f"with the ADC the zero-Doppler ridge at the target's range is {ck['adc_ridge_dbm']:.1f} dBm (target {ck['peak_dbm']:.1f} dBm, floor {ck['floor_dbm']:.1f} dBm) and the target's peak moves {ck['adc_peak_change_db']:+.2f} dB")
    g0 = t["monostatic"][0]
    w(f"  B, compression at each gain setting of the receiver (rx60.radar_iwr6843_fig7_2), {t['monostatic'][0]['reading']}: "
      + "; ".join(f"{gd:g} dB gain: NF {nf:g} dB, P1dB {p1:g} dBm, needs {need:.1f} dB" for gd, nf, p1, need in g0["gains"]))
    cs = t["compression"]
    w(f"  B, the limiter against the leakage's level (in-band P1dB {fe.p1db_inband_dbm.value:g} dBm, p = {s.rx_rapp_p.show('g')}; the echo's SNR never loses more than {-cs['snr_limit_closed_db']:.2f} dB: the limiter takes the noise down with it):")
    w(f"      {'leakage [dBm]':>14s} {'SNR closed':>12s} {'SNR simulated':>14s} {'echo closed':>12s} {'echo simulated':>15s} {'image bound [dBc]':>18s} {'image simulated':>16s} {'range sidelobes':>16s}")
    for r in cs["rows"]:
        w(f"      {r['level_dbm']:14.0f} {r['snr_closed_db']:12.2f} {r['snr_sim_db']:14.2f} {r['signal_closed_db']:12.2f} {r['signal_sim_db']:15.2f} {r['image_closed_dbc']:18.1f} {r['image_sim_dbc']:16.1f} {r['row_sidelobe_sim_dbc']:16.1f}")
    w(f"      with another knee the {s.desense_allowed_db.value:g} dB criterion moves ({t['monostatic'][0]['reading']}): "
      + "; ".join(f"p = {p:g}: {'never reached' if need is None else format(need, '.1f') + ' dB'}" for p, need in t["monostatic"][0]["need_by_p"]))
    w(f"  the TRN rows pay the same three: the leakage (or the direct path) is at the same level whatever field is sensed on; the receiver is on for "
      f"{s.trn_field.chips.size / s.chip_rate_hz.value * 1e6:.1f} us of TRN field a packet instead of {s.field.chips.size / s.chip_rate_hz.value * 1e6:.2f} us of L-CEF")
    for blk in t["bistatic"]:
        w(f"  C, {blk['reading']}: direct path {blk['direct_dbm'].value:.1f} dBm {blk['direct_dbm'].label}, {blk['under_p1db_db'].value:.1f} dB under the in-band P1dB; ENOB needed {blk['enob_needed'].value:.2f}; "
          f"independent LOs after referencing to the direct path: residual {blk['pn_residual_dbc']:.1f} dBc, floor rise closed {blk['pn_rise_closed_db']:.2f} dB / simulated {blk['pn_rise_sim_db']:.2f} dB (SNR {blk['pn_snr_change_db']:+.2f} dB); "
          f"compression at that level (closed form): SNR {blk['comp_snr_closed_db']:+.3f} dB, image bound {blk['comp_image_closed_dbc']:.1f} dBc")
    w("     C pays elsewhere: the two devices share no clock, so every packet's phase is taken from the direct path (simulated) and its timing and carrier offset from the preamble (not simulated)")
    w("")
    w(f"== 6. air time: {s.frames_per_s.value:g} sensing frames a second, against a link at {a['rate_mbps'].show('g')} Mb/s (mmw60_trx's comparison table, ay_sc on channel 2) ==")
    w(f"{'row':34s} {'limit':>20s} {'share of the time':>40s} {'throughput given up [Mb/s]':>40s}")
    for r in rows:
        share = f"{r['air_fraction'].value * 100:.2f} %" + ("" if r["air_fraction_upper"] is None else f" ... {r['air_fraction_upper'].value * 100:.0f} %") + f" {r['air_fraction'].label}"
        loss = f"{r['loss_mbps'].value:.0f}" + ("" if r["loss_upper_mbps"] is None else f" ... {r['loss_upper_mbps'].value:.0f}") + f" {r['loss_mbps'].label}"
        w(f"{r['name']:34s} {r['limit']:>20s} {share:>40s} {loss:>40s}")
    w(f"  A: the frame's {A.get('compliance.frame_on_s') * 1e3:g} ms of chirps, during which this front end cannot communicate (jcas60.fmcw_airtime)")
    w(f"  ay: {s.n_packets.value:g} packets of {s.packet_chips.value:g} chips ({s.packet_s.value * 1e6:.3f} us; with the TRN field {s.trn_packet_chips.value:g} chips, {rows[5]['packet_s'].value * 1e6:.3f} us) "
      f"plus one SIFS ({s.sifs_s.value * 1e6:g} us) each; under the radar reading the upper figure is the off-time "
      "taken to bind every emission of the device (jcas60.radar_reading_airtime); C's figure is the transmitting device's")
    w("")
    w("== 7. a Golay pair's sidelobes when the target moves (waveform simulation: one packet, no noise, target on a chip) ==")
    w(f"{'radial velocity [m/s]':>24s} {'peak sidelobe in the zero zone [dBc]':>40s} {'peak [dB]':>12s} {'carrier turn over one sequence [rad]':>40s}")
    for v, side, peak, turn in t["golay"]:
        w(f"{v:24.1f} {side:40.1f} {peak:12.4f} {turn:40.2e}")
    w(f"  0 m/s: the pair cancels exactly; what the arithmetic leaves is rounding and is printed as its bound (one unit in the last place on every "
      f"sample of the buffer: {t['rounding_dbc']['l_cef']:.1f} dBc, {t['rounding_dbc']['trn']:.1f} dBc for the TRN field's); "
      f"{t['golay'][2][0]:g} m/s is the unambiguous velocity of section 1")
    w(f"  the TRN field, blocks {'-'.join(str(b) for b in (S.raw('jcas60.trn_sensing#correlated_blocks')[0], S.raw('jcas60.trn_sensing#correlated_blocks')[-1]))} of every subfield correlated "
      f"(all six blocks at rest: {t['golay_trn_all_blocks']:.1f} dBc):")
    w(f"{'radial velocity [m/s]':>24s} {'peak sidelobe in the zero zone [dBc]':>40s} {'peak [dB]':>12s} {'carrier turn over the TRN field [rad]':>40s}")
    for v, side, peak, turn in t["golay_trn"]:
        w(f"{v:24.1f} {side:40.1f} {peak:12.4f} {turn:40.2e}")
    w("  the peak column is the coherent loss inside one packet: the TRN field is long enough for a fast target to turn the carrier across it")
    w("")
    w("== 8. sensitivity: row B under the communication limit, closed form, one choice changed at a time ==")
    for what, d, snr, note in t["sensitivity"]:
        w(f"  {what:92s} {d:+6.2f} dB -> SNR(cell) {snr:6.2f} dB   {note}")
    w(f"  the TRN field's length, closed form, communication limit (the table uses {s.trn_units.value:g} units: the fewest that reach row A's SNR):")
    w(f"      {'TRN units':>10s} {'correlated chips':>17s} {'SNR(cell) [dB]':>15s} {'over row A [dB]':>16s} {'packet [us]':>12s} {'air time [%]':>13s}")
    for units, chips, snr, over, t_pkt, share in t["trn_sweep"]:
        w(f"      {units:10d} {chips:17d} {snr:15.2f} {over:+16.2f} {t_pkt * 1e6:12.2f} {share * 100:13.2f}")
    w("")
    w("== 9. this library's own choices (every entry of data/assumed.yaml; what would move is in each entry's note) ==")
    for k, e in S.assumed_entries().items():
        w(f"  {k:32s} {e.value}  [{e.source}]")
    w("")
    w("== 10. every reference read to build this table ==")
    for ref_, n in sorted(t["trace"].items()):
        val = "" if math.isnan(n.value) else f"{n.value:g}"
        w(f"  {ref_:72s} {val:>12s} [{n.tag}]")
    w(f"  row A: fmcw_radar@{A.provenance()['sha'][:7]} `{A.provenance()['command']}` (src/jcas60/data/upstream/fmcw_radar_budget.json)")
    w("")
    w("== 11. the ay sensing frame against spec_sources' 802.11bf entries, item by item (law 7; the entries are from the task group's public documents) ==")
    for c in t["bf"]:
        w(f"  {c['row']:5s} {c['item']:38s} {c['verdict']:24s} here: {c['here']}")
        w(f"  {'':5s} {'':38s} {'':24s} entry ({c['key']}): {c['entry']}")
    w("  this is not a statement of conformance: no sensing session, measurement exchange or frame of the procedure is built here, and the standard itself was not read")
    return "\n".join(L) + "\n"
