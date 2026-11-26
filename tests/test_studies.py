"""The studies beside the table (ROADMAP items 1-3; example 01): what each concluded, held."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from jcas60 import ENGINE_BANNER, studies
from jcas60.physconst import db

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def st():
    return studies.build()


def test_example_01_prints_the_expected_output(st):
    want = (REPO / "examples" / "expected" / "01_receiver_and_window.txt").read_text(encoding="utf-8")
    assert want.splitlines()[0] == ENGINE_BANNER and studies.render(st) == want


# ----------------------------------------------------------------------------- 1. the limiter
def test_the_limiter_in_continuous_time_is_where_the_table_put_it(st):
    """Law 4 at every rate; and the answer to 'what is it in continuous time': converged by 4 samples a chip, 0.04 dB
    from the table's two."""
    by = {r["analog"]: r for r in st["limiter"]["by_rate"]}
    assert set(by) == {2, 4, 8}
    assert abs(by[4]["need_closed_db"] - by[8]["need_closed_db"]) < 0.01 and abs(by[4]["inverse_envelope_mean"] - by[8]["inverse_envelope_mean"]) < 1e-3
    assert 0.0 < by[8]["need_closed_db"] - by[2]["need_closed_db"] < 0.1
    assert by[2]["inverse_envelope_mean"] < by[8]["inverse_envelope_mean"]        # the ADC-rate samples sit where the envelope is high
    for r in by.values():
        assert abs(r["need_sim_db"] - r["need_closed_db"]) < 0.3 and r["image_sim_dbc"] < r["image_bound_dbc"]
    assert by[8]["need_closed_db"] == pytest.approx(43.28, abs=0.006)


def test_the_knee_moves_the_sensitivity_criterion_and_hardly_the_ghost(st):
    knees = st["limiter"]["by_knee"]
    assert [r["knee"] for r in knees] == [1.0, 2.0, 3.0]
    need = [r["need_closed_db"] for r in knees]
    assert need == sorted(need) and need[2] - need[0] > 7.0                        # a soft knee is 8 dB kinder to the sensitivity ...
    ghost = [r["ghost_sim_db"] for r in knees]
    assert max(ghost) - min(ghost) < 2.0                                           # ... and 1.4 dB kinder to the ghost
    for r in knees:
        assert abs(r["need_sim_db"] - r["need_closed_db"]) < 0.3
        assert r["image_sim_dbc"] < r["image_bound_dbc"] and r["ghost_sim_db"] < r["ghost_closed_db"]     # the closed form bounds the image
        assert r["image_sim_dbc"] > r["floor_dbc"]                                 # at 40 dB of isolation the ghost shows, whatever the knee
    assert [round(x, 1) for x in need] == [37.0, 43.3, 45.0] and [round(x, 1) for x in ghost] == [47.2, 48.6, 48.6]


# ----------------------------------------------------------------------------- 2. the gain setting
def test_the_gain_setting_with_the_leakage_there_both_roads(st):
    gs = st["gains"]
    rows = gs["rows"]
    assert [r["gain_db"] for r in rows] == [30.0, 36.0, 40.0, 44.0, 48.0]
    for r in rows:
        lo, hi = r["have"]
        assert (lo["isolation_db"], hi["isolation_db"]) == (20.0, 40.0)
        assert abs(hi["snr_sim_db"] - hi["snr_closed_db"]) < 0.5                    # law 4 where a design could sit
        # where less than 10 dB is left the closed form over-states the undithered ADC and the peak search reads high
        assert lo["snr_closed_db"] < 10.0 and 0.0 < lo["snr_sim_db"] - lo["snr_closed_db"] < 5.0
        # riding on the leakage, the phase-noise residue and the ADC's noise lose less to the limiter than the thermal noise does
        assert r["need_db"] == max(r["need_compression_db"], r["need_pn_db"], r["need_adc_db"])
    best = rows[gs["best"][1]]
    assert best["gain_db"] == 36.0 and best["have"][1]["snr_closed_db"] == pytest.approx(25.57, abs=0.006)
    assert gs["snr_no_leak_max_db"] - best["have"][1]["snr_closed_db"] == pytest.approx(6.00, abs=0.006)
    # the thresholds flatter a noisy setting: 30 dB of gain passes all three at 40 dB of isolation and still is not the best
    assert rows[0]["need_db"] < 40.0 < rows[1]["need_db"] and rows[0]["have"][1]["snr_closed_db"] < best["have"][1]["snr_closed_db"]
    assert rows[2]["binding"] == "LO phase noise on the leakage"                   # at its own noise figure the 40 dB setting is bound by the phase noise


# ----------------------------------------------------------------------------- 3. beyond the window
def test_a_target_beyond_the_window_lands_where_the_chips_say(st):
    lcef, trn = st["beyond"]["fields"]
    assert (lcef["name"], trn["name"]) == ("L-CEF", "TRN field")
    for f in (lcef, trn):
        rows = {r["delay_chips"]: r for r in f["rows"]}
        for r in rows.values():                                                     # two roads: the chips alone and the waveform simulation
            assert db(r["level_sim"]) == pytest.approx(db(r["level_closed"]), abs=0.01) and r["lag_sim"] == r["lag_closed"]
        assert rows[58]["lag_sim"] == 58 and db(rows[58]["level_sim"]) == pytest.approx(0.0, abs=0.01)
        assert rows[128]["lag_sim"] is None and db(rows[128]["level_sim"]) < -200  # the window's far edge is an exact cancellation
    t = {r["delay_chips"]: r for r in trn["rows"]}
    # two blocks farther: three of the four correlated blocks meet their like, one its opposite -- half the amplitude
    # (79 / 160 over 40 subfields: the first one has no subfield ahead of it to supply its first block)
    assert db(t[256]["level_sim"]) == pytest.approx(20 * math.log10(79 / 160), abs=0.01) and t[256]["lag_sim"] == 0
    assert t[314]["lag_sim"] == 58 and db(t[314]["level_sim"]) == pytest.approx(db(t[256]["level_sim"]), abs=0.01)
    assert db(t[826]["level_sim"]) > -0.5 and t[826]["lag_sim"] == 58               # one subfield farther: the whole target, at the wrong range
    assert db(trn["worst_level"]) > -0.5
    c = {r["delay_chips"]: r for r in lcef["rows"]}
    for d in (256, 512, 768):                                                       # the L-CEF: one unpaired 128-chip block of 1024
        assert db(c[d]["level_sim"]) == pytest.approx(20 * math.log10(128 / 1024), abs=0.1)
    assert db(lcef["worst_level"]) == pytest.approx(-18.1, abs=0.1)


def test_inside_the_window_the_two_roads_agree_out_to_its_edge(st):
    for r in st["beyond"]["ranges"]:
        assert abs(r["snr_sim_db"] - r["snr_closed_db"]) < 0.4, r
        assert r["delay_chips"] < 120


# ----------------------------------------------------------------------------- 4. fast targets
def test_the_in_packet_loss_is_a_closed_form(st):
    lcef, trn = st["fast"]["fields"]
    for f in (lcef, trn):
        for r in f["rows"]:
            assert r["loss_sim_db"] == pytest.approx(r["loss_closed_db"], abs=1e-3), (f["name"], r)
    assert st["fast"]["v_unambiguous_mps"] == pytest.approx(12.4, abs=0.05)
    by = {round(r["velocity_mps"], 1): r["loss_closed_db"] for r in trn["rows"]}
    assert by[2.0] < 0.01 and by[12.4] == pytest.approx(0.109, abs=0.001) and by[50.0] == pytest.approx(1.849, abs=0.001)
    assert max(r["loss_closed_db"] for r in lcef["rows"]) < 0.05                    # the L-CEF is 30 times shorter
    assert trn["span_s"] > 25 * lcef["span_s"]


# ----------------------------------------------------------------------------- 5. the ADC that nothing dithers
def test_the_undithered_adc_as_a_statistic_and_what_a_dither_costs(st):
    plain, quarter, half = st["adc"]["cases"]
    assert st["adc"]["isolation_db"] == 20.0 and plain["noise_lsb"] < 0.1 and (plain["dither_lsb"], quarter["dither_lsb"], half["dither_lsb"]) == (0.0, 0.25, 0.5)
    # one cell's reading was never the quantity: over the target's delay the ridge under it spans 18 dB
    rg = plain["ridge_under_target_db"]
    assert rg["max"] - rg["min"] > 10.0 and rg["max"] > 0.0 and -6.0 < rg["median"] < 0.0
    # and the quantiser's gain on the echo is not 1 for this leakage: the peak reads high at every delay
    assert 0.3 < plain["peak_change_db"]["min"] and plain["peak_change_db"]["max"] < 2.0
    assert plain["repeating_share"] > 0.5 and plain["rise_db"]["mean"] < plain["rise_closed_db"] - 3.0
    # a dither makes the error noise -- the closed form holds again -- and the noise is worse than what it replaced
    for c in (quarter, half):
        assert c["repeating_share"] < 0.1 and abs(c["rise_db"]["mean"] - c["rise_closed_db"]) < 0.4
        assert c["snr_change_db"]["mean"] < plain["snr_change_db"]["mean"] - 5.0
    assert half["repeating_share"] < 0.01


# ----------------------------------------------------------------------------- 6. the leakage's delay
def test_the_phase_noise_residue_against_the_leakage_delay_and_a_second_path(st):
    lk = st["leak"]
    rows = lk["rows"]
    assert [round(r["delay_s"] * 1e9, 2) for r in rows] == [0.1, 0.2, 0.5, 1.0]
    for r in rows:
        assert abs(r["rise_sim_db"] - r["rise_closed_db"]) < 0.2                     # law 4
    need = [r["need_pn_db"] for r in rows]
    assert need == sorted(need) and need[1] - need[0] == pytest.approx(5.5, abs=0.1)
    assert need[3] - need[2] < 0.5 and need[3] == pytest.approx(53.2, abs=0.06)      # it saturates: beyond 0.5 ns the residue is the LO's own noise twice
    rf = lk["reflector"]
    assert rf["need_pn_min_db"] < rf["need_pn_alone_db"] < rf["need_pn_average_db"] < rf["need_pn_db"] <= rf["need_pn_max_db"]
    assert rf["need_pn_db"] == pytest.approx(47.0, abs=0.06) and rf["need_pn_average_db"] == pytest.approx(44.6, abs=0.06)
    assert abs(rf["rise_sim_db"] - rf["rise_closed_db"]) < 0.3 and rf["rise_closed_db"] > rf["rise_alone_closed_db"] + 1.0


def test_two_paths_off_one_lo_reduce_to_one():
    """The several-path residue with one path is the one-path residue; with two equal paths at the same delay and
    opposite phase it is nothing."""
    from jcas60 import budget, config
    s = config.setup()
    tau = s.leak_delay_s.value
    one = budget.pn_residual_shared(s, tau)
    assert budget.pn_residual_paths(s, ((1.0, 0.3, tau),)) == pytest.approx(one, rel=1e-12)
    assert budget.pn_residual_paths(s, ((1.0, 0.0, tau), (1.0, math.pi, tau))) < 1e-20 * one + 1e-30
    assert budget.pn_residual_paths(s, ((1.0, 0.0, tau), (1.0, 0.0, tau))) == pytest.approx(4.0 * one, rel=1e-12)


# ----------------------------------------------------------------------------- 7. one LO for every row
def test_every_row_on_the_x10_lo(st):
    lo = st["lo"]
    a, ax = lo["row_a"]["default"], lo["row_a"]["x10lo"]
    assert ax["snr_cell_db"] == pytest.approx(a["snr_cell_db"], abs=0.01)          # row A is thermal-noise limited at 5 m: its LO does not show
    assert ax["lo_rms_rad"] < a["lo_rms_rad"] and "x10" in ax["lo_source"]
    model, x10 = lo["rows"]
    assert model["need_pn_db"] == pytest.approx(43.0, abs=0.06) and x10["need_pn_db"] == pytest.approx(39.0, abs=0.06)
    assert model["binding"] == x10["binding"] == "in-band compression" and x10["need_db"] == model["need_db"]   # the shortfall does not move
    for r in (model, x10):
        assert abs(r["rise_sim_db"] - r["rise_closed_db"]) < 0.2                   # law 4, monostatic
        for c in r["bistatic"]:
            assert abs(c["rise_sim_db"] - c["rise_closed_db"]) < 0.3, c            # ... and bistatic, on both fields
    # what the closed form said while it referenced over the whole correlated span: more than 1 dB from the simulation on the TRN field
    trn = model["bistatic"][1]
    assert trn["rise_whole_field_closed_db"] - trn["rise_sim_db"] > 1.0 and trn["span_s"] > 30 * trn["reference_s"]


def test_the_x10_lo_is_the_vendored_profile_and_the_table_does_not_use_it():
    import numpy as np

    from jcas60 import config, frontend
    fe = config.setup().fe
    assert fe.pn_profile is None and fe is frontend.shared()                       # law 2: the table's one instance is untouched
    x = frontend.with_x10_lo(fe)
    assert float(x.pn_dbchz(1e6)) == pytest.approx(-98.5, abs=0.1) and float(x.pn_dbchz(1e8)) == pytest.approx(-133.5, abs=0.1)
    assert float(x.pn_dbchz(1.0)) == float(x.pn_dbchz(100.0)) and float(fe.pn_dbchz(1e6)) == pytest.approx(-93.0, abs=0.05)
    assert np.all(np.diff(x.pn_profile[0]) > 0) and dataclass_fields_equal(fe, x)


def dataclass_fields_equal(a, b) -> bool:
    import dataclasses
    return all(getattr(a, f.name) == getattr(b, f.name) for f in dataclasses.fields(a) if f.name != "pn_profile")


# ----------------------------------------------------------------------------- 8. the bistatic rows' synchronisation
def test_the_leftover_carrier_offset_is_the_fast_targets_loss(st):
    lcef, trn = st["sync"]["fields"]
    for f in (lcef, trn):
        rows = {round(r["cfo_hz"]): r for r in f["rows"]}
        assert set(rows) == {5000, 50000, 2420000}
        assert abs(rows[5000]["loss_sim_db"] - rows[5000]["loss_closed_db"]) < 0.002
        assert abs(rows[50000]["loss_sim_db"] - rows[50000]["loss_closed_db"]) < 0.01
        assert abs(rows[2420000]["loss_sim_db"] - rows[2420000]["loss_closed_db"]) < 0.1 and rows[2420000]["loss_closed_db"] > 10.0
    t5 = {round(r["cfo_hz"]): r for r in trn["rows"]}[5000]
    assert 0.05 < t5["loss_closed_db"] < 0.11                                       # 5 kHz across the TRN field: under a tenth of a dB
    assert st["sync"]["equivalent_mps_per_hz"] * 5000 == pytest.approx(12.4, abs=0.05)


# ----------------------------------------------------------------------------- 9. a range window
def test_the_range_window_trades_sidelobes_for_width_on_both_roads(st):
    rw = st["window"]
    plain, hann = rw["cases"]
    assert (rw["window"], plain["window"], hann["window"]) == ("hann", "none", "hann")          # row A's own window, by name
    for c in (plain, hann):
        assert c["sim"]["width_chips"] == pytest.approx(c["closed"]["width_chips"], rel=1e-3)
        assert c["sim"]["sidelobe_db"] == pytest.approx(c["closed"]["sidelobe_db"], abs=0.05)
        assert abs(c["loss_sim_db"] - c["loss_grid_closed_db"]) < 0.15
        for m in c["masks"]:
            assert m["sim_db"] == pytest.approx(m["closed_db"], abs=0.05)
    assert plain["closed"]["width_chips"] == pytest.approx(0.8726, abs=1e-3) and plain["closed"]["sidelobe_db"] == pytest.approx(-14.3, abs=0.06)
    assert hann["closed"]["width_chips"] == pytest.approx(1.197, abs=1e-3) and hann["closed"]["sidelobe_db"] == pytest.approx(-27.3, abs=0.06)
    assert hann["loss_closed_db"] == pytest.approx(0.95, abs=0.006) and plain["loss_closed_db"] == pytest.approx(0.0, abs=1e-9)
    assert hann["straddle_db"] < plain["straddle_db"]                               # the wider lobe straddles the lag grid better
    wide = hann["sim"]["width_chips"] * rw["chip_m"]
    assert wide / rw["row_a_width_m"] == pytest.approx(1.91, abs=0.01)             # like for like, not the 2.27 of the bandwidths
    # past the main lobe the window buys 6-11 dB against a second target; inside 1.5 chips it costs
    by = lambda c: {m["offset_chips"]: m["closed_db"] for m in c["masks"]}          # noqa: E731
    assert by(hann)[1.5] > by(plain)[1.5] and all(by(hann)[k] < by(plain)[k] - 6.0 for k in (2.5, 3.5, 5.0, 8.0))


def test_the_range_window_closed_forms():
    import numpy as np

    from jcas60 import budget
    a = 0.25
    assert budget.range_response(0.0, a)[0] == pytest.approx(1.0) and budget.range_response(1.0, a)[0] < 1e-6      # the raised cosine: zero at a chip
    assert budget.range_response(np.array([0.3, 0.7]), a) == pytest.approx(budget.rc(np.array([0.3, 0.7]), a) ** 2, abs=1e-5)
    assert budget.range_window_loss_db(a, None) == pytest.approx(0.0, abs=1e-9) and budget.range_window_loss_db(a, "hann") > 0.5
    # the slow-time window: Hann's nulls on every bin past the first, its first sidelobe between bins 2 and 3
    d = budget.doppler_response(np.array([0.0, 1.0, 2.0, 5.0, 2.5]), "hann", 64)
    assert d[0] == pytest.approx(1.0) and d[1] == pytest.approx(0.25, abs=1e-9) and d[2] < 1e-25 and d[3] < 1e-25
    assert 10 * np.log10(d[4]) == pytest.approx(-32.3, abs=0.3)


# ----------------------------------------------------------------------------- 10. a second target
def test_what_a_strong_target_leaves_in_doppler_and_a_static_reflector(st):
    sec = st["second"]
    static, moving = sec["doppler"]
    assert static["off_bin"] == 0.0 and abs(moving["off_bin"]) > 0.1
    for r in static["rows"]:                                                        # a static path: the rounding bound, on both roads
        assert r["sim_db"] == r["closed_db"] < -200.0
    for r in moving["rows"]:
        assert abs(r["sim_db"] - r["closed_db"]) < 0.5
    lv = [r["closed_db"] for r in moving["rows"]]
    assert lv == sorted(lv, reverse=True) and lv[0] > -32.0 and lv[1] < -48.0     # the window's sidelobes, falling away
    cl = sec["clutter"]
    assert abs(cl["rise_sim_db"] - cl["rise_closed_db"]) < 0.2 and cl["rise_closed_db"] == pytest.approx(1.08, abs=0.006)
    assert abs(cl["snr_linear_db"] - cl["snr_alone_db"]) < 0.1                     # an ideal LO: the reflector costs the mover nothing


# ----------------------------------------------------------------------------- 11. a waveform made for sensing
def test_the_designed_waveform_against_the_standard_subfield(st):
    dw = st["designed"]
    std, new = dw["rows"]
    # every block between the guards, a whole number of periods; the clean zone twice the standard subfield's
    assert (std["correlated_chips"], new["correlated_chips"]) == (20480, 29 * 1024) and dw["gain_db"] == pytest.approx(1.61, abs=0.006)
    assert (std["clean_chips"], new["clean_chips"]) == (128, 256) and new["clean_m"] == pytest.approx(21.80, abs=0.006)
    for r in (std, new):
        assert abs(r["snr_sim_db"] - r["snr_grid_db"]) < 0.15                       # law 4
        assert r["sidelobe_rest_dbc"] < -200.0 and r["sidelobe_moving_dbc"] < -75.0
    assert new["snr_cell_db"] - std["snr_cell_db"] == pytest.approx(dw["gain_db"], abs=1e-9)
    # the entry's four conditions: three fit, the spectral one does not -- and backing off to meet it takes the SNR gain back
    cd = dw["conditions"]
    assert cd["length_chips"][0] == cd["length_chips"][1] and cd["modulus_max_error"] < 1e-12 and cd["consecutive_same_sequence"][0] <= cd["consecutive_same_sequence"][1]
    assert cd["psd_average_hz"] == 10e6 and cd["psd_excess_max_db"] == pytest.approx(1.86, abs=0.006) and abs(cd["psd_total_db"]) < 1e-9
    assert dw["gain_db"] - cd["psd_excess_max_db"] < 0.0
    # either side of a target: nothing out to two blocks, noise-like from there to six, the whole target again at eight
    clean, noisy, fold = dw["zones"]
    assert (clean["from_chips"], clean["to_chips"]) == (1, 256) and db(clean["level"]) < -200.0
    assert -22.0 < db(noisy["level"]) < -20.0 and db(fold["level"]) > -0.5 and fold["from_chips"] == 1024
    rows = {r["delay_chips"]: r for r in dw["beyond"]}
    for r in rows.values():                                                         # the chips alone and the waveform simulation
        assert db(r["level_sim"]) == pytest.approx(db(r["level_closed"]), abs=0.01)
    for d in (58, 128, 160, 192):                                                   # targets the standard subfield's window ends before
        assert rows[d]["lag_sim"] == d and db(rows[d]["level_sim"]) == pytest.approx(0.0, abs=0.01)
    assert rows[256]["lag_sim"] is None and all(db(rows[d]["level_sim"]) < -20.0 for d in (314, 384, 512, 640, 768, 826))


def test_the_designed_pattern_is_built_on_the_perfect_sequence():
    """The sign pattern's two halves: each a length-4 sequence whose periodic autocorrelation is zero off the peak,
    and whose cross-correlation is zero at the two shifts that neighbouring blocks meet at."""
    import numpy as np

    from jcas60 import sources as S_
    items = list(S_.raw("jcas60.sensing_waveform#pattern"))
    assert [x.strip("+-") for x in items] == ["Ga_N", "Gb_N"] * 4
    sign = np.array([1 if x[0] == "+" else -1 for x in items])
    e, o = sign[0::2], sign[1::2]
    acf = lambda x, m: int(np.sum(x * np.roll(x, m)))                              # noqa: E731
    assert [acf(e, m) for m in range(4)] == [4, 0, 0, 0] and [acf(o, m) for m in range(4)] == [4, 0, 0, 0]
    assert int(np.sum(e * o)) == 0 and int(np.sum(e * np.roll(o, 1))) == 0
