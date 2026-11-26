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
