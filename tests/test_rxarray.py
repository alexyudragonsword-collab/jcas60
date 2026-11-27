"""Milestone 1's second step, a receive array (example 03): what it concluded, held."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from jcas60 import ENGINE_BANNER, arrays, config, rxarray, table
from jcas60 import wavesim as W
from jcas60.physconst import db

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def rx():
    return rxarray.build()


def _cell(rx, n, limit, field):
    return next(f for r in rx["ways"] if r["n"] == n for p in r["readings"] if p["limit"] == limit for f in p["fields"] if f["field"] == field)


def test_example_03_prints_the_expected_output(rx):
    want = (REPO / "examples" / "expected" / "03_receive_array.txt").read_text(encoding="utf-8")
    assert want.splitlines()[0] == ENGINE_BANNER and rxarray.render(rx) == want


def test_one_element_is_the_table_and_the_swept_column_is_example_02(rx):
    """Law 2: N = 1 is the comparison table's row in every column; the transmit-swept column is example 02's SNR."""
    t, ar = table.build(), arrays.build()
    for p in rx["ways"][0]["readings"]:
        rows = {r["field"]: r for r in t["rows"] if r.get("limit") == p["limit"] and "field" in r and r["kind"] == "monostatic"}
        for f in p["fields"]:
            for k in ("rx_only_db", "tx_swept_db", "tx_swept_equal_air_db", "both_db"):
                assert f[k] == pytest.approx(f["one_antenna_db"], abs=1e-9), k
            assert f["one_antenna_db"] == pytest.approx(rows[f["field"]]["link"].snr_cell_db.value, abs=1e-9)
    for r in rx["ways"]:
        a = next(x for x in ar["rows"] if x["n"] == r["n"])
        assert r["beams"] == a["beams"]
        for p, q in zip(r["readings"], a["readings"]):
            for f, g in zip(p["fields"], q["fields"]):
                assert f["tx_swept_db"] == pytest.approx(g["snr_cell_db"], abs=1e-9)


def test_law3_receive_gain_is_not_capped_and_transmit_gain_is(rx):
    """The receive array adds 10 log10(N) under both readings; the transmit array nothing under the radar reading and
    nothing past the communication limit."""
    comm, radar = (p["limit"] for p in rx["ways"][0]["readings"])
    for r in rx["ways"]:
        for limit in (comm, radar):
            for field in ("L-CEF", "TRN field"):
                f = _cell(rx, r["n"], limit, field)
                assert f["rx_only_db"] - f["one_antenna_db"] == pytest.approx(db(r["n"])) == pytest.approx(f["both_db"] - f["tx_swept_db"])
                assert f["tx_swept_equal_air_db"] == pytest.approx(f["tx_swept_db"] - db(r["beams"]))
        assert _cell(rx, r["n"], radar, "L-CEF")["tx_swept_db"] == _cell(rx, 1, radar, "L-CEF")["one_antenna_db"]
    assert _cell(rx, 8, comm, "L-CEF")["rx_only_db"] == pytest.approx(39.77, abs=0.006) and _cell(rx, 8, radar, "L-CEF")["rx_only_db"] == pytest.approx(37.77, abs=0.006)
    assert _cell(rx, 8, comm, "TRN field")["both_db"] == pytest.approx(70.78, abs=0.006) and _cell(rx, 16, comm, "L-CEF")["tx_swept_db"] == _cell(rx, 8, comm, "L-CEF")["tx_swept_db"]


def test_at_the_same_air_time_sweeping_buys_what_receiving_buys_and_costs_the_doppler_axis(rx):
    """Until the cap: N transmitters swept, each beam given 1 / N of the packets, come to the receive array's SNR
    (the cap's 0.06 dB apart at N = 8) with N times coarser velocity resolution; past the cap and under the radar
    reading the sweep falls behind."""
    comm, radar = (p["limit"] for p in rx["ways"][0]["readings"])
    for n in (2, 4):
        f = _cell(rx, n, comm, "L-CEF")
        assert f["tx_swept_equal_air_db"] == pytest.approx(f["rx_only_db"], abs=1e-9)
    f8 = _cell(rx, 8, comm, "L-CEF")
    assert f8["rx_only_db"] - f8["tx_swept_equal_air_db"] == pytest.approx(0.06, abs=0.006)
    f16 = _cell(rx, 16, comm, "L-CEF")
    assert f16["rx_only_db"] - f16["tx_swept_equal_air_db"] == pytest.approx(6.08, abs=0.006)
    r8 = _cell(rx, 8, radar, "L-CEF")
    assert r8["rx_only_db"] - r8["tx_swept_equal_air_db"] == pytest.approx(2.0 * db(8), abs=1e-9)          # 18 dB, not the 9 expected
    row = next(r for r in rx["ways"] if r["n"] == 8)
    assert row["velocity_res_equal_air_mps"] == pytest.approx(8 * row["velocity_res_mps"]) == pytest.approx(3.10, abs=0.006)


def test_law4_the_combined_chains_give_the_closed_forms_snr(rx):
    c = rx["combination"]
    assert c["n"] == 8 and c["theta_deg"] == 20.0 and c["pattern_db"] == pytest.approx(4.0 * db(math.cos(math.radians(20.0))))
    assert abs(c["snr_one_sim_db"] - c["snr_one_closed_db"]) < 0.3 and abs(c["snr_sim_db"] - c["snr_closed_db"]) < 0.3
    assert c["gain_closed_db"] == pytest.approx(9.03, abs=0.006) and abs(c["gain_sim_db"] - c["gain_closed_db"]) < 0.3


def test_law4_the_direction_from_one_frame_meets_the_cramer_rao_bound(rx):
    """Maximum likelihood on the combined cell against the bound: the rms error within the scatter of 32 frames
    (about 13 %) above it, no bias, and a small fraction of the beam."""
    c = rx["combination"]
    assert c["angle_bound_deg"] == pytest.approx(0.0770, abs=5e-4)
    assert 0.8 < c["angle_rms_deg"] / c["angle_bound_deg"] < 1.35
    assert abs(c["angle_mean_deg"]) < 3.0 * c["angle_bound_deg"] / math.sqrt(c["trials"]) and c["angle_rms_deg"] < c["beamwidth_deg"] / 100.0
    s = rx["setup"]
    g = arrays.geometry(s, 8)
    assert rxarray.angle_bound_deg(s, g, 20.0, 30.0) / rxarray.angle_bound_deg(s, g, 20.0, 36.0) == pytest.approx(10 ** 0.3, rel=1e-9)


def test_law4_the_leakage_in_every_chain(rx):
    """Each chain compresses on its own, so the isolation needed is one antenna's; combined, the phase-noise residue
    follows the leakage through the combiner.  Closed form against the simulated chains in three cases."""
    k = rx["leakage"]
    t = table.build()
    assert k["need_db"] == t["monostatic"][0]["cost"].need_db.value == pytest.approx(43.24, abs=0.006) and k["binding"] == "in-band compression"
    assert abs(k["one_chain_sim_db"] - k["one_chain_closed_db"]) < 0.5 and abs(k["no_leak_sim_db"] - k["no_leak_closed_db"]) < 0.3
    look, skirt, null = k["cases"]
    for q in k["cases"]:
        assert abs(q["snr_sim_db"] - q["snr_closed_db"]) < 0.5, q["name"]
    assert (look["share"], skirt["share"], null["share"]) == (1.0, 1.0 / 8, 0.0) and look["direction_deg"] == 20.0 < skirt["direction_deg"] < null["direction_deg"]
    assert look["snr_closed_db"] < skirt["snr_closed_db"] < null["snr_closed_db"] < k["no_leak_closed_db"]
    gain = [q["snr_closed_db"] - k["one_chain_closed_db"] for q in k["cases"]]
    assert gain == pytest.approx([3.29, 9.03, 11.20], abs=0.02)             # the array's 9.03 dB: less, all of it, more
    assert skirt["snr_closed_db"] - k["one_chain_closed_db"] == pytest.approx(db(8), abs=1e-6)
    # the skirt's direction is where arraytx's array factor has fallen to 1 / N
    from jcas60.vendor.mmw60.arraytx import beam as BM
    iso = arrays.geometry(rx["setup"], 8, "isotropic")
    a = BM.array_factor_cw(iso, 20.0, 0.0, np.array([skirt["direction_deg"], null["direction_deg"], 20.0]), "phase")
    assert (a[0] / a[2]) ** 2 == pytest.approx(1.0 / 8, rel=1e-6) and (a[1] / a[2]) ** 2 < 1e-20


def test_the_chains_share_the_lo_and_nothing_else():
    """The simulator's second noise stream: with it unset a frame is what it was; two chains with the same seed and
    different noise seeds have the same phase noise and different thermal noise."""
    s = config.setup()
    sim = W.Simulator(s)
    p = [W.Path(-60.0, 3.0e-8, 0.0)]
    quiet = [sim.frame(p, W.Options(seed=1, noise_seed=k, lo="shared", noise=False, n_packets=2)).cir for k in (2, 3)]
    assert np.array_equal(quiet[0], quiet[1])
    noisy = [sim.frame(p, W.Options(seed=1, noise_seed=k, n_packets=2)).cir for k in (2, 3)]
    assert not np.allclose(noisy[0], noisy[1])
    assert np.array_equal(sim.frame(p, W.Options(seed=1, n_packets=2)).cir, sim.frame(p, W.Options(seed=1, noise_seed=None, n_packets=2)).cir)
