"""Acceptance of milestone 0's table: row A is fmcw_radar's report item by item; the example prints what is kept."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from jcas60 import ENGINE_BANNER, budget, config, regulatory, table
from jcas60 import fmcw_row as A

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def t():
    return table.build()


#: the table's cell -> where fmcw_radar's report has it (a path, or a function of the report)
ROW_A = {
    "bandwidth_hz": lambda g: g("radar.sweep_ghz") * 1e9,
    "range_res_m": "radar.range_res_rect_m", "range_width_m": "radar.range_res_m", "range_sidelobe_db": "radar.window_sidelobe_db",
    "range_max_m": "radar.max_range_m", "velocity_res_mps": "radar.velocity_res_mps", "velocity_max_mps": "radar.max_velocity_mps",
    "eirp_dbm": "compliance.eirp_dbm", "snr_db": "detection.snr_cell_db", "pd": "detection.pd", "pfa": "detection.pfa",
    "snr_required_db": "detection.range_at_pd.snr_required_db", "range_pd_m": "detection.range_at_pd.snr_limited_range_m",
    "range_actual_m": "detection.range_at_pd.range_m",
    "cap_dbm": "compliance.limit_dbm", "backoff_db": "compliance.tx_backoff_db", "off_time_required_s": "compliance.off_time_required_s",
    "adc_fs_hz": lambda g: g("radar.adc_fs_msps") * 1e6, "adc_bits": "radar.adc_bits",
    "air_fraction": lambda g: g("compliance.frame_on_s") / g("compliance.frame_period_s"),
}


def test_row_a_equals_fmcw_radars_report_item_by_item(t):
    a = t["rows"][0]
    for cell, where in ROW_A.items():
        want = where(A.get) if callable(where) else A.get(where)
        assert a[cell].value == want, cell
    assert a["window_loss_db"].value == next(r["value"] for r in A.get("budget") if r["term"] == "window loss")
    assert a["rule"] == A.get("compliance.rule") and a["range_window"] == A.get("radar.range_window")
    assert a["off_time_have_s"].value == A.governing_tier()["off_time_have_s"] and a["tier"] == A.governing_tier()["label"]
    assert a["leakage"] is None and A.get("leakage") is None
    assert a["paths_delta_db"].value == 0.03                     # fmcw_radar's example 00 as printed


def test_row_a_cells_not_in_the_report_are_derived_from_it_and_say_so(t):
    a = t["rows"][0]
    for cell in ("processing_gain_db", "loss_mbps"):
        assert a[cell].tag == "derived" and "fmcw_radar report" not in a[cell].ref
    # the range at Pd 0.9 is read from the report since fmcw_radar#42; this library's old extrapolation -- the report's
    # SNR carried out by R^4 and the oxygen term -- lands on the same figure (the default config has no IF high-pass)
    r = a["range_pd_m"].value
    snr_there = a["snr_db"].value - 40 * __import__("math").log10(r / 5.0) - A.get("inputs.oxygen_db_per_km.value") * 2 * (r - 5.0) / 1e3
    assert snr_there == pytest.approx(a["snr_required_db"].value, abs=1e-6)
    assert a["range_limited_by"] == "if_band" and a["range_actual_m"].value == a["range_max_m"].value


def test_every_row_is_held_to_the_same_required_snr(t):
    """One criterion for the three ways: the report's CA-CFAR figure; each ay row's range is where its closed-form
    cell SNR falls to it."""
    req = t["rows"][0]["snr_required_db"]
    assert req.value == t["snr_required_db"] == pytest.approx(12.24, abs=0.005)
    assert t["snr_required_db"] - t["snr_required_known_noise_db"] == pytest.approx(0.49, abs=0.005)
    s = t["setup"]
    for r in t["rows"][1:]:
        field = s.trn_field if r["field"] == "TRN field" else s.field
        at = budget.link(s, r["kind"], r["eirp_dbm"], r["range_pd_m"].value, field).snr_cell_db.value
        assert at == pytest.approx(req.value, abs=1e-6), (r["name"], r["limit"])


def test_the_recorded_report_is_what_the_pinned_fmcw_radar_writes():
    """Regenerates the report from the pinned commit when the sibling checkout is there (CI checks it out)."""
    import upstream_report as U
    if U.sibling() is None:
        pytest.skip("fmcw_radar checkout not present")
    assert U.main(["check", "--strict"]) == 0
    assert A.provenance() == U.UPSTREAM


def test_the_numbers_the_conclusions_quote(t):
    """cairn/comparison-table.md quotes these; a change that moves one has to revisit the conclusions."""
    rows = t["rows"]
    got = {"A snr": rows[0]["snr_db"].value, "B comm snr": rows[1]["snr_db"].value, "B radar snr": rows[2]["snr_db"].value,
           "B comm range": rows[1]["range_pd_m"].value, "A range": rows[0]["range_pd_m"].value,
           "B res cm": rows[1]["range_res_m"].value * 100, "A res cm": rows[0]["range_res_m"].value * 100,
           "B width cm": rows[1]["range_width_m"].value * 100, "B clean": rows[1]["range_max_m"].value,
           "need": t["monostatic"][0]["cost"].need_db.value, "ghost": t["monostatic"][0]["ghost"]["need_sim_db"],
           "B trn snr": rows[5]["snr_db"].value, "B trn radar snr": rows[6]["snr_db"].value, "air trn": rows[5]["air_fraction"].value * 100,
           "comp 40": t["monostatic"][0]["checks"][1]["comp_snr_closed_db"], "image 40": t["monostatic"][0]["checks"][1]["comp_image_sim_dbc"], "need adc": t["monostatic"][0]["cost"].need_adc_db.value,
           "need pn": t["monostatic"][0]["cost"].need_pn_db.value, "air": rows[1]["air_fraction"].value * 100, "air A": rows[0]["air_fraction"].value * 100}
    want = {"A snr": 43.34, "B comm snr": 30.74, "B radar snr": 28.74, "B comm range": 14.28, "A range": 28.78, "B res cm": 8.52, "A res cm": 3.75,
            "B width cm": 7.43, "B clean": 10.90, "need": 43.2, "ghost": 48.7, "B trn snr": 43.75, "B trn radar snr": 41.75, "air trn": 4.51, "comp 40": -2.17, "image 40": -11.8, "need adc": 39.9, "need pn": 43.0, "air": 1.16, "air A": 19.20}
    for k, v in want.items():
        assert got[k] == pytest.approx(v, abs=0.006 if (abs(v) < 20 and k not in ("ghost", "image 40")) or "snr" in k or "range" in k else 0.06), k
    assert regulatory.comm_rate_mbps().value == 4620.0
    assert budget.snr_required_db(config.setup().pd.value, 1e-4) == pytest.approx(11.749, abs=1e-3)


def test_golay_sidelobes_rise_with_velocity_from_exact_cancellation(t):
    g = t["golay"]
    # at rest the pair cancels exactly: what is left is under the arithmetic's rounding, and is reported as that bound
    assert g[0][0] == 0.0 and g[0][1] == t["rounding_dbc"]["l_cef"] < -200.0 and g[1][1] > g[0][1] + 100.0
    assert g[1][0] == 2.0 and -90.0 < g[1][1] < -85.0
    levels = [x[1] for x in g]
    assert levels == sorted(levels)
    # 20 dB a decade: the sidelobe goes as the carrier's turn over the field
    assert (g[4][1] - g[3][1]) == pytest.approx(20 * __import__("math").log10(g[4][0] / g[3][0]), abs=0.3)


def test_example_prints_the_expected_output_with_the_banner_first(t):
    want = (REPO / "examples" / "expected" / "00_comparison_table.txt").read_text(encoding="utf-8")
    assert want.splitlines()[0] == ENGINE_BANNER
    assert table.render(t) == want                                # in process
    got = subprocess.run([sys.executable, str(REPO / "examples" / "00_comparison_table.py")], capture_output=True, text=True, check=True).stdout
    assert got == want                                            # and as a script (tools/example_drift.py compare does the same)


def test_every_number_cell_of_the_table_carries_a_tag(t):
    from jcas60.sources import TAGS, Num
    n = 0
    for r in t["rows"]:
        for k, v in r.items():
            if isinstance(v, Num):
                n += 1
                assert v.tag in TAGS and v.ref, (r["name"], k)
    assert n > 100


def test_the_trn_field_is_built_from_the_entry_and_is_clean_on_blocks_2_to_5(t):
    """The TRN rows' field: a TRN-R packet's units from edmg.trn_field, correlated on blocks 2-5 of every subfield."""
    import numpy as np

    from jcas60 import field as F
    from jcas60 import sources as S
    s = t["setup"]
    f = s.trn_field
    per_unit, n = F.trn_r_subfields_per_unit(), int(S.num("edmg.trn_field#trn_bl@trn_seq_length=normal").value)
    assert (per_unit, n) == (10, 128) and f.chips.size == int(s.trn_units.value) * per_unit * 6 * n == 30720
    assert f.n_correlated == int(s.trn_units.value) * per_unit * 4 * n == 20480 and len(f.windows) == 40
    assert s.trn_packet_chips.value == s.packet_chips.value + f.chips.size and s.trn_packet_chips.value % 4 == 0   # the rotation carries on
    # a delayed copy, any delay up to one block either way: the correlation over blocks 2-5 is the delay's and nothing else
    for d in (0, 1, 57, 127, -1, -100):
        got = sum(np.vdot(ref, f.chips[w0 + n * 6 - d: w0 + n * 6 - d + len(ref)]) for w0, ref in f.windows[:-1])
        assert abs(got) < 1e-9 if d else abs(got - (len(f.windows) - 1) * 4 * n) < 1e-9
    assert t["golay_trn"][0][1] == t["rounding_dbc"]["trn"] < -200 and -38 < t["golay_trn_all_blocks"] < -34   # all six blocks leave sidelobes
    assert t["golay_trn"][1][0] == 2.0 and t["golay_trn"][1][1] < -80 and t["golay_trn"][1][2] > -0.01
    # the number of units is the fewest that reach row A under the communication limit
    a = t["rows"][0]["snr_db"].value
    sweep = {u: snr for u, _, snr, *_ in t["trn_sweep"]}
    assert sweep[int(s.trn_units.value)] >= a > sweep[int(s.trn_units.value) - 1]
    assert sweep[int(s.trn_units.value)] == pytest.approx(t["rows"][5]["snr_db"].value, abs=1e-9)


def test_the_two_correlators_agree():
    """A long field is correlated through one FFT against a template, the L-CEF by gathering its two windows: the
    same sums, checked on the L-CEF where both can run."""
    import numpy as np

    from jcas60 import wavesim as W
    s = config.setup()
    a, b = W.Simulator(s), W.Simulator(s)
    assert not a._by_fft and W.Simulator(s, s.trn_field)._by_fft
    tmpl = np.zeros(b.n, dtype=complex)
    for w0, ref in s.field.windows:
        tmpl[(b.start + w0) * b.osr: (b.start + w0 + len(ref)) * b.osr: b.osr] = ref
    b._by_fft, b._tmpl_f = True, np.conj(np.fft.fft(tmpl))
    paths = [W.Path(-85.0, 58.7 / s.chip_rate_hz.value, 4.0), W.Path(-28.0, s.leak_delay_s.value)]
    x, y = a.frame(paths, W.Options(n_packets=2)).cir, b.frame(paths, W.Options(n_packets=2)).cir
    assert np.abs(x - y).max() < 1e-12 * np.abs(x).max()


def test_the_reader_finds_a_slow_target_next_to_the_ridge(t):
    """The bistatic target sits 4.5 Doppler bins from the static direct path; the peak is read on its own bin, not
    on the next one out (which reads 0.32 dB low)."""
    for r in (t["rows"][3], t["rows"][4], t["rows"][7], t["rows"][8]):
        assert abs(r["paths_delta_db"].value) < 0.2, (r["name"], r["limit"], r["paths_delta_db"].value)
