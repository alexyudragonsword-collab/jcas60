"""Milestone 1's first step, a transmit array (example 02): what it concluded, held."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from jcas60 import ENGINE_BANNER, arrays, config, table
from jcas60 import sources as S
from jcas60.physconst import db
from jcas60.vendor.mmw60.arraytx import beam as BM

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def ar():
    return arrays.build()


def _row(ar, n):
    return next(r for r in ar["rows"] if r["n"] == n)


def _cell(ar, n, limit, field):
    return next(f for p in _row(ar, n)["readings"] if p["limit"] == limit for f in p["fields"] if f["field"] == field)


def test_example_02_prints_the_expected_output(ar):
    want = (REPO / "examples" / "expected" / "02_transmit_array.txt").read_text(encoding="utf-8")
    assert want.splitlines()[0] == ENGINE_BANNER and arrays.render(ar) == want


def test_one_element_is_the_comparison_tables_row(ar):
    """Law 2: N = 1 is the same front end at the same EIRP, so the same closed-form SNR as the table's rows."""
    t = table.build()
    assert _row(ar, 1)["eirp_full_dbm"] == config.setup().fe.eirp_max_dbm.value
    for p in _row(ar, 1)["readings"]:
        rows = [r for r in t["rows"] if r.get("limit") == p["limit"] and "field" in r and r["kind"] == "monostatic"]
        assert len(rows) == 2                                                            # L-CEF and TRN field
        assert sorted(f["snr_cell_db"] for f in p["fields"]) == pytest.approx(sorted(r["link"].snr_cell_db.value for r in rows), abs=1e-9)


def test_law3_each_reading_caps_the_array_at_its_own_limit(ar):
    """The array reaches the communication limit at 8 elements and can go no further; under the radar reading it never
    buys EIRP, only back-off: 20 log10(N) more of it on every PA."""
    for r in ar["rows"]:
        assert r["eirp_full_dbm"] == pytest.approx(ar["one_element_eirp_dbm"] + 2.0 * db(r["n"]))
        for p in r["readings"]:
            assert p["eirp_dbm"].value == min(r["eirp_full_dbm"], p["cap_dbm"]) and p["eirp_dbm"].tag == "derived"
            assert p["backoff_db"] == pytest.approx(r["eirp_full_dbm"] - p["eirp_dbm"].value) and p["backoff_db"] >= 0.0
    comm, radar = (p["limit"] for p in ar["rows"][0]["readings"])
    assert [r["n"] for r in ar["rows"] if any(p["limit"] == comm and p["backoff_db"] > 0 for p in r["readings"])] == [8, 16]
    assert _cell(ar, 8, comm, "L-CEF")["snr_cell_db"] == pytest.approx(48.74, abs=0.006)
    assert _cell(ar, 8, comm, "TRN field")["snr_cell_db"] == pytest.approx(61.75, abs=0.006)
    assert _cell(ar, 16, comm, "L-CEF")["snr_cell_db"] == _cell(ar, 8, comm, "L-CEF")["snr_cell_db"]
    assert len({_cell(ar, r["n"], radar, "L-CEF")["snr_cell_db"] for r in ar["rows"] if r["n"] > 1}) == 1
    assert _cell(ar, 8, comm, "L-CEF")["snr_cell_db"] - _cell(ar, 8, radar, "L-CEF")["snr_cell_db"] == pytest.approx(20.0)


def test_law4_the_beam_and_the_count_of_beams_on_two_roads(ar):
    """The continuous aperture's closed form against arraytx's array factor: within 1 % from 8 elements up, and the
    array's own beam is the wider one at small N (two elements: 60 deg)."""
    for r in ar["rows"]:
        if r["n"] == 1:
            assert r["beams"] == 1 and r["width_deg"] is None
            continue
        assert r["width_deg"] >= r["width_closed_deg"]
        if r["n"] >= 8:
            assert r["width_deg"] / r["width_closed_deg"] - 1.0 < 0.01
        assert r["beams"] == math.ceil(r["beams_closed"]) == math.ceil(r["beams_array"]) == r["n"]
        assert r["beams_array"] <= r["beams_closed"]
        d = r["beam_directions_deg"]
        assert d[0] > -ar["fov_deg"] and d == sorted(d) and np.allclose(d, -np.asarray(d[::-1]), atol=1e-9)
    assert _row(ar, 2)["width_deg"] == pytest.approx(60.0, abs=0.1) and _row(ar, 8)["width_deg"] == pytest.approx(12.79, abs=0.006)
    assert arrays.half_power() == pytest.approx(0.886, abs=5e-4)


def test_adjacent_beams_cross_where_the_entry_says(ar):
    """Each pair of neighbours meets at -3 dB of both (arraytx's array factor, evaluated between them), and the outermost
    beams' -3 dB edges are outside the field of view."""
    s = ar["setup"]
    drop = S.num("jcas60.array#crossover_db").value
    for n in (2, 8, 16):
        g = arrays.geometry(s, n, "isotropic")
        d = _row(ar, n)["beam_directions_deg"]
        u = np.sin(np.radians(d))
        for a, b, ua, ub in zip(d, d[1:], u, u[1:]):
            mid = np.array([math.degrees(math.asin(0.5 * (ua + ub)))])
            for th in (a, b):
                af = abs(BM.array_factor_cw(g, th, 0.0, mid, ar["mode"]))[0] / abs(BM.array_factor_cw(g, th, 0.0, np.array([th]), ar["mode"]))[0]
                assert 2.0 * db(af) == pytest.approx(-drop, abs=0.05)
        assert u[-1] + arrays.beam_half_width_u(g, ar["mode"]) >= math.sin(math.radians(ar["fov_deg"]))


def test_a_sweep_costs_one_frame_per_beam(ar):
    for r in ar["rows"]:
        for p in r["readings"]:
            for f in p["fields"]:
                assert f["air_sweep"] == pytest.approx(f["air_one_beam"] * r["beams"])
        assert r["sweeps_per_s_at_one_beams_air_time"] == pytest.approx(ar["frames_per_s"] / r["beams"])
        assert r["scan_loss_edge_db"] == pytest.approx(6.02, abs=0.006)
    comm = ar["rows"][0]["readings"][0]["limit"]
    assert _cell(ar, 8, comm, "TRN field")["air_sweep"] == pytest.approx(0.3611, abs=5e-5)
    assert _cell(ar, 8, comm, "L-CEF")["air_sweep"] == pytest.approx(0.0930, abs=5e-5)


def test_law4_the_squint_on_two_roads(ar):
    """The band integral against the waveform driven element by element: the same to 1e-3 dB; under 0.1 dB at 16
    elements and 60 deg, nothing with true time delay, and four times less at half the angle's tangent squared."""
    for q in ar["squint"]:
        assert abs(q["closed_db"] - q["sim_db"]) < 1e-3
        assert (q["closed_db"] < 1e-6) if q["mode"] == "ttd" else (0.0 < q["closed_db"] < 0.1)
    by = {(q["field"], q["theta_deg"]): q["closed_db"] for q in ar["squint"] if q["mode"] == "phase"}
    assert by[("L-CEF", 60.0)] == pytest.approx(0.0499, abs=5e-5) and by[("TRN field", 60.0)] == pytest.approx(by[("L-CEF", 60.0)], abs=1e-4)
    assert by[("L-CEF", 60.0)] / by[("L-CEF", 30.0)] == pytest.approx((math.tan(math.radians(60.0)) / math.tan(math.radians(30.0))) ** 2 / 3.0, rel=0.02)


def test_law4_the_link_on_two_roads_at_the_communication_limit(ar):
    for q in ar["roads"]:
        assert q["n"] == 8 and q["eirp_dbm"] == 40.0 and abs(q["snr_sim_db"] - q["snr_grid_db"]) < 0.3


def test_the_leakage_takes_most_of_the_arrays_gain_back(ar):
    """What was not expected: with N PAs beside the receiver the monostatic row keeps +1.3 dB (leakages in phase) to
    +10.0 dB (their powers add) of the array's 18 dB at the communication limit, and both roads say so."""
    k = ar["leakage"]
    row = {(r["n"], r["limit"]): r for r in k["rows"]}
    comm, radar = k["limit"], next(r["limit"] for r in k["rows"] if r["limit"] != k["limit"])
    one, eight = row[(1, comm)], row[(8, comm)]
    assert one["sums"][0] == {**one["sums"][1], "sum": "sum of powers"}                  # one element: the two sums are one
    assert one["sums"][0]["need_db"] == pytest.approx(43.24, abs=0.006)                  # the comparison table's own figure
    assert eight["snr_no_leak_db"] - one["snr_no_leak_db"] == pytest.approx(18.0, abs=0.006)
    gain = [q["snr_closed_db"] - one["sums"][0]["snr_closed_db"] for q in eight["sums"]]
    assert gain[0] == pytest.approx(10.04, abs=0.02) and gain[1] == pytest.approx(1.31, abs=0.02)
    # the needed isolation follows the conducted power: +10 log10(N) (powers) or +20 log10(N) (in phase) at full power
    assert eight["sums"][0]["need_db"] - one["sums"][0]["need_db"] == pytest.approx(eight["eirp_dbm"] - 22.0 - db(8), abs=0.01)
    assert eight["sums"][1]["need_db"] - one["sums"][0]["need_db"] == pytest.approx(eight["eirp_dbm"] - 22.0, abs=0.01)
    # under the radar reading the array only backs the PAs off: the leakage falls by up to 10 log10(N) and the 40 dB found are enough
    assert row[(8, radar)]["sums"][0]["shortfall_db"] < 0.0 < row[(1, radar)]["sums"][0]["shortfall_db"]
    assert row[(8, radar)]["sums"][1]["need_db"] == row[(1, radar)]["sums"][0]["need_db"]
    assert all(q["binding"] == "in-band compression" for r in k["rows"] for q in r["sums"])
    for q, r in zip(k["roads"], eight["sums"]):                                          # law 4
        assert q["snr_closed_db"] == r["snr_closed_db"] and abs(q["snr_sim_db"] - q["snr_closed_db"]) < 0.5
    assert abs(k["snr_no_leak_sim_db"] - k["snr_no_leak_closed_db"]) < 0.3


def test_the_entry_is_read_and_its_limits_are_refused():
    assert S.raw("jcas60.array#steering") == "phase" and S.num("jcas60.array#rx_elements").value == 1
    assert arrays.beams_placed(arrays.geometry(config.setup(), 1), 60.0, "phase") == [0.0]
