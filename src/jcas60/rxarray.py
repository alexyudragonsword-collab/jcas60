"""Milestone 1, second step: a receive array.

Example 02 put N power amplifiers on the transmit side: 20 log10(N) of EIRP,
capped by the rule, paid for with a beam that has to be swept.  This module
puts N receive chains behind N elements instead -- each with its own ADC,
all on one LO, the beams formed in arithmetic afterwards -- and asks:

1. what each way of spending N elements buys, under both regulatory readings
   and at the same air time: one transmit antenna with N receivers, N
   transmitters swept with one receiver, and both;
2. whether the simulated chains, combined, give the SNR the closed form says,
   and how well the combination tells the target's direction (against the
   Cramer-Rao bound);
3. what the transmit leakage does when it enters every chain: each chain
   compresses on its own, the phase-noise residue riding on the leakage is
   common to all of them (one LO, one waveform) and adds as the leakage's
   phases across the elements make it;
4. what a combiner that puts a null on the leakage gets back: the leakage is
   one waveform, so across the chains it is one complex vector whatever its
   phases are, and the data give it -- by projection (the vector read off
   the zero-Doppler ridge) and adaptively (the covariance of the map's floor).

The array is arraytx's geometry (vendored through mmw60_trx), the element and
the spacing are example 02's (``jcas60.array``), and the receive chains are
the table's receiver N times.  Two roads where a conclusion is drawn (law 4);
the table's rows are not replaced.
"""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from . import arrays as AR
from . import budget as B
from . import config
from . import regulatory as R
from . import sources as S
from . import table as T
from . import wavesim as W
from .physconst import db, undb
from .studies import _with_leakage_closed
from .vendor.fmcw_radar.processing import rd as _rd
from .vendor.mmw60.arraytx import beam as BM

KEY = "jcas60.rx_array"


# ----------------------------------------------------------------------------- the chains and their combination
def weights(s, geom, theta_deg: float) -> np.ndarray:
    """The combiner's weights for a look at ``theta``: each chain turned back by the carrier phase its element's
    position gives a wave from there (arraytx's ``delays_s``: element k hears it that much earlier)."""
    return np.exp(-2j * np.pi * s.carrier_hz.value * geom.delays_s(theta_deg))


def chains(s, sim: W.Simulator, geom, sources, index: int, **options) -> np.ndarray:
    """The correlator outputs of every chain for one frame, (elements, packets, lags).  ``sources`` is a list of
    (path, direction): each chain receives every path advanced by its element's position towards that direction.
    The chains share the frame's LO (one ``seed``); each has its own thermal noise (``noise_seed``)."""
    seed = int(s.sim["seed"])
    out = []
    for k in range(geom.n):
        paths = [W.Path(p.power_dbm, p.delay_s - float(geom.delays_s(th)[k]), p.path_rate_mps, p.name) for p, th in sources]
        out.append(sim.frame(paths, W.Options(seed=seed, noise_seed=seed + 1 + k, **options), index).cir)
    return np.array(out)


def combined(s, sim: W.Simulator, geom, cirs: np.ndarray, theta_deg: float) -> W.Frame:
    """The chains summed with the look's weights: one frame, as one receiver would have delivered it."""
    cir = np.tensordot(weights(s, geom, theta_deg), cirs, axes=1)
    return W.Frame(cir=cir, lags_per_chip=sim.osr, lag0=-sim.lag_lo, rd=_rd.doppler_fft(cir, s.doppler_window, 1))


def angle_bound_deg(s, geom, theta_deg: float, snr_element_db: float) -> float:
    """The Cramer-Rao bound on the direction of one target from one frame of a uniform line array: the phase step
    between neighbours is a frequency along the array, var >= 6 / (SNR N (N^2 - 1)), and d(step)/d(theta) is
    2 pi (d / lambda) cos(theta).  ``snr_element_db``: the target's SNR in its cell on one chain."""
    n = geom.n
    step = math.sqrt(6.0 / (undb(snr_element_db) * n * (n * n - 1)))
    return math.degrees(step / (2.0 * math.pi * geom.spacing_lambda * math.cos(math.radians(theta_deg))))


def angle_estimate_deg(s, geom, x: np.ndarray, near_deg: float, half_width_deg: float) -> float:
    """The direction at which the combined cell is largest (maximum likelihood for one target), searched inside the
    main lobe around ``near_deg``."""
    def miss(th):
        return -abs(np.dot(weights(s, geom, th), x))
    return float(minimize_scalar(miss, bounds=(near_deg - half_width_deg, near_deg + half_width_deg), method="bounded", options={"xatol": 1e-10}).x)


# ----------------------------------------------------------------------------- 1. three ways to spend N elements
def ways(s) -> list:
    mode, fov = str(S.raw(AR.KEY + "#steering")), S.num(AR.KEY + "#field_of_view_deg").value
    fields = (("L-CEF", s.field), ("TRN field", s.trn_field))
    lam, m, t_int = s.wavelength_m, s.n_packets.value, s.interval_s.value
    res = lam / (2.0 * m * t_int)
    rows = []
    for n in (int(x) for x in S.raw(KEY + "#rx_elements")):
        beams = len(AR.beams_placed(AR.geometry(s, n), fov, mode))
        per = []
        for rd in R.readings(s):
            one = rd.eirp_dbm.value                                   # one antenna under this reading: the table's row
            swept = min(s.fe.eirp_max_dbm.value + 2.0 * db(n), rd.cap_dbm.value)
            cells = []
            for name, f in fields:
                snr1 = B.link(s, "monostatic", rd.eirp_dbm, field=f).snr_cell_db.value
                cells.append({"field": name, "one_antenna_db": snr1,
                              "rx_only_db": snr1 + db(n),
                              "tx_swept_db": snr1 + swept - one,
                              "tx_swept_equal_air_db": snr1 + swept - one - db(beams),
                              "both_db": snr1 + swept - one + db(n)})
            per.append({"limit": rd.name, "eirp_one_dbm": one, "eirp_swept_dbm": swept, "fields": cells})
        rows.append({"n": n, "beams": beams, "readings": per, "velocity_res_mps": res, "velocity_res_equal_air_mps": res * beams,
                     "chains_rx_only": n, "chains_tx_swept": 1})
    return rows


# ----------------------------------------------------------------------------- 2. the combination on the second road
def combination(s, bench) -> dict:
    n, theta = int(S.num(KEY + "#simulated_elements").value), S.num(KEY + "#target_angle_deg").value
    trials = int(S.num(KEY + "#angle_trials").value)
    geom = AR.geometry(s, n)
    g = geom.element_gain(theta)
    eirp = R.readings(s)[0].eirp_dbm
    sim = bench.sim
    link = B.link(s, "monostatic", eirp)
    geo, tgt, _ = bench.paths("monostatic", link, 0.0)
    near = int(round(geo.delay_s * sim.fs))
    pattern_db = 4.0 * db(g)                                          # the element's pattern both ways, in power
    target = W.Path(tgt.power_dbm + pattern_db, tgt.delay_s, tgt.path_rate_mps, "target")
    half = 0.5 * BM.beamwidth_3db_deg(AR.geometry(s, n, "isotropic"), theta, "phase")
    snr_all, snr_one, got = [], [], []
    for k in range(trials):
        cirs = chains(s, sim, geom, [(target, theta)], k)
        fr = combined(s, sim, geom, cirs, theta)
        r = sim.read(fr, near)
        snr_all.append((r["peak_mw"], r["floor_mw"]))
        one = sim.read(W.Frame(cir=cirs[0], lags_per_chip=sim.osr, lag0=-sim.lag_lo, rd=_rd.doppler_fft(cirs[0], s.doppler_window, 1)), near)
        snr_one.append((one["peak_mw"], one["floor_mw"]))
        x = np.array([_rd.doppler_fft(c, s.doppler_window, 1)[r["doppler_index"], r["lag"] - sim.lag_lo] for c in cirs])
        got.append(angle_estimate_deg(s, geom, x, theta, half))

    def ratio(pairs):
        p, f = np.mean(pairs, axis=0)
        return db(float(p / f))
    err = np.asarray(got) - theta
    snr_el = link.snr_grid_db.value + pattern_db
    return {"n": n, "theta_deg": theta, "trials": trials, "pattern_db": pattern_db,
            "snr_one_closed_db": snr_el, "snr_one_sim_db": ratio(snr_one),
            "snr_closed_db": snr_el + db(n), "snr_sim_db": ratio(snr_all), "gain_closed_db": db(n), "gain_sim_db": ratio(snr_all) - ratio(snr_one),
            "angle_bound_deg": angle_bound_deg(s, geom, theta, snr_el), "angle_rms_deg": float(np.sqrt(np.mean(err ** 2))), "angle_mean_deg": float(np.mean(err)),
            "beamwidth_deg": 2.0 * half}


# ----------------------------------------------------------------------------- 3. the leakage in every chain
def _direction_at(geom, look_deg: float, level: float, mode: str) -> float:
    """The direction on the main lobe's skirt, beyond the look, from which a wave comes through the combiner at
    ``level`` (power, 1 at the look): where the leakage's phases across the elements leave that much of it."""
    def af2(th):
        a = BM.array_factor_cw(geom, look_deg, 0.0, np.array([th, look_deg]), mode)
        return float((a[0] / a[1]) ** 2)
    null = math.degrees(math.asin(math.sin(math.radians(look_deg)) + 1.0 / (geom.n * geom.spacing_lambda)))
    if level <= 0.0:
        return null
    return brentq(lambda th: af2(th) - level, look_deg, null, xtol=1e-12)


def leakage(s, bench) -> dict:
    """One transmit antenna, N receive chains, the table's leakage in each (the best isolation found).  Combined
    for a look at the target: the echo adds in amplitude, each chain's thermal and ADC noise in power, and the
    leakage with what rides on it as its phases across the elements make it -- written here as the direction a wave
    with those phases would come from: the look itself, the skirt where 1 / N of it is left, and the first null."""
    n, theta = int(S.num(KEY + "#simulated_elements").value), S.num(KEY + "#target_angle_deg").value
    frames = int(S.num(KEY + "#leak_frames").value)
    geom, iso = AR.geometry(s, n), AR.geometry(s, n, "isotropic")
    eirp = R.readings(s)[0].eirp_dbm
    sim = bench.sim
    link = B.link(s, "monostatic", eirp)
    geo, tgt, _ = bench.paths("monostatic", link, 0.0)
    near = int(round(geo.delay_s * sim.fs))
    pattern_db = 4.0 * db(geom.element_gain(theta))
    target = W.Path(tgt.power_dbm + pattern_db, tgt.delay_s, tgt.path_rate_mps, "target")
    env = sim.unit_envelope(s.leak_delay_s.value)
    cost = B.monostatic_cost(s, eirp, env)
    label, isolation = s.fe.isolation_db[-1]
    level = cost.conducted_dbm.value - isolation.value
    c = _with_leakage_closed(s, link, cost, env, level)
    thermal = undb(c["thermal_db"])
    pn, adc = thermal * undb(c["pn_over_thermal_db"]), thermal * undb(c["adc_over_thermal_db"])
    no_leak = link.snr_grid_db.value + pattern_db
    on = {"adc_enob": s.adc_enob.value, "lo": "shared", "rx_p1db_dbm": s.fe.p1db_inband_dbm.value}
    leak = W.Path(level, s.leak_delay_s.value, 0.0, "leakage")

    def read(cirs, strong):
        return sim.read(combined(s, sim, geom, cirs, theta), near, strong)

    def mean_snr(sources, strong, **options):
        rs = [read(chains(s, sim, geom, sources, k, **options), strong) for k in range(frames)]
        return db(float(np.mean([r["peak_mw"] for r in rs]) / np.mean([r["floor_mw"] for r in rs])))
    cases = []
    for name, share in (("at the look", 1.0), ("on the skirt, 1 / N left", 1.0 / n), ("in the first null", 0.0)):
        th = _direction_at(iso, theta, share, "phase")
        closed = no_leak + db(n) + c["echo_db"] - db(thermal + n * share * pn + adc)
        cases.append({"name": name, "direction_deg": th, "share": share, "snr_closed_db": closed,
                      "snr_sim_db": mean_snr([(target, theta), (leak, th)], (0,), **on)})
    one = np.mean([[r[k] for k in ("peak_mw", "floor_mw")] for r in
                   (sim.read(sim.frame([target, leak], W.Options(seed=int(s.sim["seed"]), noise_seed=int(s.sim["seed"]) + 1, **on), k), near, (0,)) for k in range(frames))], axis=0)
    return {"n": n, "theta_deg": theta, "frames": frames, "isolation_label": label, "isolation_db": isolation.value, "leak_dbm": level,
            "need_db": cost.need_db.value, "binding": cost.binding, "need_compression_db": cost.need_compression_db.value,
            "one_chain_no_leak_db": no_leak, "one_chain_closed_db": no_leak + c["snr_change_db"], "one_chain_sim_db": db(float(one[0] / one[1])),
            "no_leak_closed_db": no_leak + db(n), "no_leak_sim_db": mean_snr([(target, theta)], ()),
            "echo_db": c["echo_db"], "pn_over_thermal_db": c["pn_over_thermal_db"], "adc_over_thermal_db": c["adc_over_thermal_db"], "cases": cases}


# ----------------------------------------------------------------------------- 4. a null on the leakage
METHODS = ("conventional", "projected", "adaptive")


def leak_signature(sim: W.Simulator, cirs: np.ndarray, leak_delay_s: float) -> np.ndarray:
    """The leakage across the chains, from the data: it is static, so each chain's packets averaged at the leakage's
    lag hold its complex gain there (the echo and the noise average down, the phase-noise residue with them)."""
    return cirs[:, :, int(round(leak_delay_s * sim.fs)) - sim.lag_lo].mean(axis=1)


def floor_cells(s, sim: W.Simulator, shape: tuple, cell: tuple) -> np.ndarray:
    """The cells of a map that hold neither the zero-Doppler ridge, nor the leakage's own lags, nor the target and
    its neighbours, nor the lags past the clean zone: what the adaptive combiner learns the floor's covariance on."""
    m, n_lags = shape
    ex_d, ex_r = int(s.sim["exclude_doppler_bins"]), int(s.sim["exclude_range_lags"]) * sim.osr
    i, j = cell
    rows = (np.abs(np.arange(m) - m // 2) > ex_d) & (np.abs(np.arange(m) - i) > ex_d)
    cols = (sim.lags > ex_r) & (sim.lags <= (sim.field.zero_sidelobe_chips - int(s.sim["exclude_range_lags"])) * sim.osr) & (np.abs(np.arange(n_lags) - j) > ex_r)
    return rows[:, None] & cols[None, :]


def null_weights(s, sim: W.Simulator, geom, cirs: np.ndarray, theta_deg: float, how: str, cell: tuple) -> np.ndarray:
    """The combiner u (the output is u^H x) for a look at ``theta``.  ``conventional``: the look's own signature a.
    ``projected``: a with its part along the leakage's signature taken out -- a null on the leakage, whatever its
    phases.  ``adaptive``: R^-1 a with R the chains' covariance over the map's floor cells (thermal and ADC noise
    on the diagonal, the phase-noise residue along the leakage's signature): the most SNR any weights give."""
    a = np.conj(weights(s, geom, theta_deg))
    if how == "conventional":
        return a
    if how == "projected":
        v = leak_signature(sim, cirs, s.leak_delay_s.value)
        return a - v * (np.vdot(v, a) / np.vdot(v, v))
    if how == "adaptive":
        rds = np.array([_rd.doppler_fft(c, s.doppler_window, 1) for c in cirs])
        x = rds[:, floor_cells(s, sim, rds.shape[1:], cell)]
        return np.linalg.solve(x @ x.conj().T / x.shape[1], a)
    raise ValueError(how)


def beam_shares(s, geom, v: np.ndarray) -> np.ndarray:
    """What share of a signature ``v`` each of the N orthogonal beams of the array receives.  They sum to one
    (Parseval), whatever v is: a null on the leakage costs one beam's worth of gain in all, spread as v spreads it."""
    n = geom.n
    looks = [math.degrees(math.asin(k / (n * geom.spacing_lambda))) for k in range(-(n // 2), n - n // 2)]
    return np.array([abs(np.vdot(np.conj(weights(s, geom, th)), v)) ** 2 for th in looks]) / (n * float(np.vdot(v, v).real))


def _level_db(sim: W.Simulator, share: float) -> float:
    """A share of power in dB as a reading: one that is the arithmetic's rounding is recorded as what it stands for --
    zero dB for a share of one, the simulator's rounding bound for a share of nothing."""
    level = db(max(share, undb(sim.rounding_dbc)))
    return 0.0 if abs(level) < 1e-9 else level


def nulling(s, bench) -> dict:
    """Section 3's set-up -- one transmit antenna, N chains, the table's leakage in each -- with the leakage coming
    through the conventional combiner at several shares of its power, and three combiners on the same frames."""
    n, theta = int(S.num(KEY + "#simulated_elements").value), S.num(KEY + "#target_angle_deg").value
    frames = int(S.num(KEY + "#leak_frames").value)
    geom, iso = AR.geometry(s, n), AR.geometry(s, n, "isotropic")
    eirp = R.readings(s)[0].eirp_dbm
    sim = bench.sim
    link = B.link(s, "monostatic", eirp)
    geo, tgt, _ = bench.paths("monostatic", link, 0.0)
    near = int(round(geo.delay_s * sim.fs))
    pattern_db = 4.0 * db(geom.element_gain(theta))
    target = W.Path(tgt.power_dbm + pattern_db, tgt.delay_s, tgt.path_rate_mps, "target")
    env = sim.unit_envelope(s.leak_delay_s.value)
    cost = B.monostatic_cost(s, eirp, env)
    level = cost.conducted_dbm.value - s.fe.isolation_db[-1][1].value
    c = _with_leakage_closed(s, link, cost, env, level)
    thermal = undb(c["thermal_db"])
    pn, own = thermal * undb(c["pn_over_thermal_db"]), thermal * (1.0 + undb(c["adc_over_thermal_db"]))     # the residue; each chain's own noise
    top = link.snr_grid_db.value + pattern_db + db(n) + c["echo_db"]                # the echo combined, over one chain's undisturbed floor
    on = {"adc_enob": s.adc_enob.value, "lo": "shared", "rx_p1db_dbm": s.fe.p1db_inband_dbm.value}
    leak = W.Path(level, s.leak_delay_s.value, 0.0, "leakage")
    a = np.conj(weights(s, geom, theta))
    rows = []
    for share in (float(x) for x in S.raw(KEY + "#null_shares")):
        th = theta if share >= 1.0 else _direction_at(iso, theta, share, "phase")
        v_true = np.conj(weights(s, geom, th))
        closed = {"conventional": top - db(own + n * share * pn),
                  "projected": None if share >= 1.0 else top + db(1.0 - share) - db(own),
                  "adaptive": top - db(own) + db(1.0 - share * n * pn / (own + n * pn))}
        got = {how: {"peak": [], "floor": [], "gain": [], "left": []} for how in METHODS}
        for k in range(frames):
            cirs = chains(s, sim, geom, [(target, theta), (leak, th)], k, **on)
            first = sim.read(combined(s, sim, geom, cirs, theta), near, (0,))
            cell = (first["doppler_index"], first["lag"] - sim.lag_lo)
            for how in METHODS:
                u = null_weights(s, sim, geom, cirs, theta, how, cell)
                cir = np.tensordot(np.conj(u), cirs, axes=1)
                r = sim.read(W.Frame(cir=cir, lags_per_chip=sim.osr, lag0=-sim.lag_lo, rd=_rd.doppler_fft(cir, s.doppler_window, 1)), near, (0,))
                norm = float(np.vdot(u, u).real) * n
                got[how]["peak"].append(r["peak_mw"])
                got[how]["floor"].append(r["floor_mw"])
                got[how]["gain"].append(abs(np.vdot(u, a)) ** 2 / norm)
                got[how]["left"].append(abs(np.vdot(u, v_true)) ** 2 / norm)
        methods = {}
        for how in METHODS:
            g = got[how]
            blind = closed[how] is None
            methods[how] = {"snr_closed_db": closed[how], "snr_sim_db": None if blind else db(float(np.mean(g["peak"]) / np.mean(g["floor"]))),
                            "look_gain_db": _level_db(sim, float(np.mean(g["gain"]))), "leak_left_db": _level_db(sim, float(np.mean(g["left"])))}
        rows.append({"share": share, "direction_deg": th, "methods": methods,
                     "look_gain_projected_closed_db": None if share >= 1.0 else db(1.0 - share)})
    rng = np.random.default_rng(int(s.sim["seed"]))
    random_v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    plane = np.conj(weights(s, geom, _direction_at(iso, theta, 1.0 / n, "phase")))
    sums = {"a plane wave's phases": beam_shares(s, geom, plane), "random gains and phases": beam_shares(s, geom, random_v)}
    return {"n": n, "theta_deg": theta, "frames": frames, "rows": rows, "no_leak_closed_db": link.snr_grid_db.value + pattern_db + db(n), "null_limit_closed_db": top - db(own),
            "floor_cells": int(floor_cells(s, sim, (int(s.n_packets.value), sim.lags.size), (0, 0)).sum()),
            "beam_sums": {k: {"sum": float(v.sum()), "largest": float(v.max()), "mean_gain_db": db(float(np.mean(1.0 - v)))} for k, v in sums.items()},
            "mean_gain_closed_db": db(1.0 - 1.0 / n)}


# ----------------------------------------------------------------------------- all of it
@lru_cache(maxsize=1)
def build() -> dict:
    s = config.setup()
    bench = T._Bench(s)
    if S.num(AR.KEY + "#crossover_db").value != 3:
        raise ValueError("arraytx's beam width is the half-power one (3 dB)")
    return {"setup": s, "ways": ways(s), "combination": combination(s, bench), "leakage": leakage(s, bench), "nulling": nulling(s, bench),
            "adc_fs_hz": s.adc_fs_hz.value, "adc_enob": s.adc_enob.value, "fov_deg": S.num(AR.KEY + "#field_of_view_deg").value,
            "spacing_lambda": S.num(AR.KEY + "#spacing_lambda").value, "element": str(S.raw(AR.KEY + "#element_pattern"))}


def _d(v, fmt: str = "{:.2f}") -> str:
    return "blind" if v is None else fmt.format(v)


def render(t: dict | None = None) -> str:
    from . import ENGINE_BANNER, WORDING
    t = t or build()
    L: list[str] = []
    w = L.append
    w(ENGINE_BANNER)
    w(WORDING)
    w("")
    w(f"== 1. three ways to spend N elements (line, {t['spacing_lambda']:g} lambda apart, {t['element']} element pattern; row B at 5 m, SNR(cell), closed form; L-CEF / TRN field) ==")
    w("  rx only: one transmit antenna, N receive chains, beams formed afterwards -- every direction from one frame")
    w("  tx swept: N transmit elements, one receiver, a frame per beam (example 02); at the same air time each beam gets 1 / beams of the packets")
    w("  both: N transmit elements swept, N receive chains")
    for i, name in enumerate(p["limit"] for p in t["ways"][0]["readings"]):
        w(f"  {name}:")
        w(f"{'N':>6s} {'beams':>6s} {'one antenna':>18s} {'rx only':>18s} {'tx swept':>18s} {'tx swept, same air time':>26s} {'both':>18s}")
        for r in t["ways"]:
            a, b = r["readings"][i]["fields"]
            w(f"{r['n']:6d} {r['beams']:6d} " + " ".join(f"{a[k]:{wd}.2f} / {b[k]:.2f} dB" for k, wd in
                                                       (("one_antenna_db", 6), ("rx_only_db", 6), ("tx_swept_db", 6), ("tx_swept_equal_air_db", 14), ("both_db", 6))))
    r8 = next(r for r in t["ways"] if r["n"] == t["combination"]["n"])
    w(f"  air time: rx only x1, tx swept and both x beams.  Velocity resolution {r8['velocity_res_mps']:.2f} m/s; tx swept at the same air time, N = {r8['n']}: {r8['velocity_res_equal_air_mps']:.2f} m/s")
    w(f"  receive chains: rx only and both N, tx swept 1; each is the table's receiver and ADC ({t['adc_fs_hz'] / 1e9:g} GS/s complex, ENOB {t['adc_enob']:.2f})")
    w("")
    c = t["combination"]
    w(f"== 2. the chains combined, on the second road (N = {c['n']}, one transmit antenna, the target at {c['theta_deg']:g} deg; L-CEF, communication limit; {c['trials']} frames) ==")
    w(f"  one chain: closed on the grid {c['snr_one_closed_db']:.2f} dB (the element's pattern both ways: {c['pattern_db']:+.2f} dB), simulated {c['snr_one_sim_db']:.2f} dB")
    w(f"  combined:  closed {c['snr_closed_db']:.2f} dB, simulated {c['snr_sim_db']:.2f} dB; the gain: closed {c['gain_closed_db']:.2f} dB, simulated {c['gain_sim_db']:.2f} dB")
    w(f"  the direction from one frame (the combined cell at its largest, searched inside the {c['beamwidth_deg']:.2f} deg beam): rms error {c['angle_rms_deg']:.4f} deg, "
      f"mean {c['angle_mean_deg']:+.4f} deg; Cramer-Rao bound {c['angle_bound_deg']:.4f} deg")
    w("")
    k = t["leakage"]
    w(f"== 3. the leakage in every chain (one transmit antenna, N = {k['n']}; each chain at {k['isolation_db']:g} dB of isolation, {k['isolation_label']}: leakage {k['leak_dbm']:g} dBm; "
      f"limiter, LO phase noise and ADC on; {k['frames']} frames) ==")
    w(f"  each chain compresses on its own: the isolation needed is one antenna's, {k['need_db']:.2f} dB ({k['binding']}); the array does not relax it")
    w(f"  one chain: without leakage {k['one_chain_no_leak_db']:.2f} dB; with it, closed {k['one_chain_closed_db']:.2f} dB, simulated {k['one_chain_sim_db']:.2f} dB "
      f"(echo {k['echo_db']:+.2f} dB, phase noise {k['pn_over_thermal_db']:+.2f} dB and ADC {k['adc_over_thermal_db']:+.2f} dB against the thermal floor)")
    w(f"  combined, without leakage: closed {k['no_leak_closed_db']:.2f} dB, simulated {k['no_leak_sim_db']:.2f} dB")
    w(f"{'the leakage comes through the combiner':>42s} {'as if from':>12s} {'left of it':>12s} {'SNR closed':>12s} {'simulated':>11s} {'difference':>12s}")
    for q in k["cases"]:
        w(f"{q['name']:>42s} {q['direction_deg']:8.2f} deg {q['share']:12.3f} {q['snr_closed_db']:9.2f} dB {q['snr_sim_db']:8.2f} dB {q['snr_sim_db'] - q['snr_closed_db']:+9.2f} dB")
    w("  combined, the echo adds in amplitude; each chain's thermal and ADC noise in power; the phase-noise residue is the same in every chain (one LO, one waveform) and follows the leakage")
    w("  through the combiner.  The leakage's phases across the elements are not modelled: the three rows span what they can do")
    w("")
    z = t["nulling"]
    w("== 4. a null on the leakage (the same set-up and frames for three combiners; the leakage's vector across the chains is taken from the data, not from a model of its phases) ==")
    w("  conventional: the look's own weights.  projected: the look's weights with their part along the leakage taken out (the leakage read off each chain's zero-Doppler ridge)")
    w(f"  adaptive: R^-1 a, R the chains' covariance over the map's {z['floor_cells']} floor cells (no ridge, no leakage lags, no target)")
    w(f"{'leakage left by the':>22s} {'as if':>10s} {'conventional':>22s} {'projected':>22s} {'adaptive':>22s} {'gain at the look, projected':>30s} {'adaptive':>10s} {'leakage left: projected':>25s} {'adaptive':>10s}")
    w(f"{'conventional combiner':>22s} {'from':>10s} {'closed / simulated':>22s} {'closed / simulated':>22s} {'closed / simulated':>22s} {'closed / simulated':>30s}")
    for r in z["rows"]:
        m = r["methods"]
        cells = " ".join(f"{_d(m[h]['snr_closed_db']):>9s} / {_d(m[h]['snr_sim_db']):>6s} dB" for h in METHODS)
        w(f"{r['share']:22.3f} {r['direction_deg']:6.2f} deg {cells} {_d(r['look_gain_projected_closed_db']):>17s} / {m['projected']['look_gain_db']:6.2f} dB {m['adaptive']['look_gain_db']:7.2f} dB "
          f"{m['projected']['leak_left_db']:22.1f} dB {m['adaptive']['leak_left_db']:7.1f} dB")
    w(f"  without leakage {z['no_leak_closed_db']:.2f} dB; the most a null can give back {z['null_limit_closed_db']:.2f} dB (each chain's limiter has already taken its part: that is not undone)")
    w("  the projection is blind where the leakage sits at the look, and loses to the conventional combiner near it; the adaptive weights are never under either")
    b = z["beam_sums"]
    w("  whatever the leakage's vector, the shares its N orthogonal beams receive sum to one: " + "; ".join(f"{k}: sum {v['sum']:.6f}, largest {v['largest']:.3f}" for k, v in b.items()))
    w(f"  so a projected null costs one beam's worth of gain in all: {z['mean_gain_closed_db']:.2f} dB on the average beam (" + ", ".join(f"{v['mean_gain_db']:.2f}" for v in b.values()) + " dB counted)")
    return "\n".join(L) + "\n"
