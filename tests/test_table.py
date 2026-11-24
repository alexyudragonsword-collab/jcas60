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
    for cell in ("range_pd_m", "processing_gain_db", "loss_mbps"):
        assert a[cell].tag == "derived"
    # the range at Pd 0.9: the report's SNR carried out by R^4 and the oxygen term to the required SNR
    r = a["range_pd_m"].value
    snr_there = a["snr_db"].value - 40 * __import__("math").log10(r / 5.0) - A.get("inputs.oxygen_db_per_km.value") * 2 * (r - 5.0) / 1e3
    assert snr_there == pytest.approx(t["snr_required_db"], abs=1e-6)


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
           "need": t["monostatic"][0]["cost"].need_db.value, "need adc": t["monostatic"][0]["cost"].need_adc_db.value,
           "need pn": t["monostatic"][0]["cost"].need_pn_db.value, "air": rows[1]["air_fraction"].value * 100, "air A": rows[0]["air_fraction"].value * 100}
    want = {"A snr": 43.34, "B comm snr": 30.75, "B radar snr": 28.75, "B comm range": 14.68, "A range": 29.56, "B res cm": 8.52, "A res cm": 3.75,
            "B width cm": 7.43, "B clean": 10.90, "need": 56.0, "need adc": 39.9, "need pn": 43.0, "air": 1.16, "air A": 19.20}
    for k, v in want.items():
        assert got[k] == pytest.approx(v, abs=0.006 if v < 20 or "snr" in k or "range" in k else 0.06), k
    assert regulatory.comm_rate_mbps().value == 4620.0
    assert budget.snr_required_db(config.setup().pd.value, 1e-4) == pytest.approx(11.749, abs=1e-3)


def test_golay_sidelobes_rise_with_velocity_from_exact_cancellation(t):
    g = t["golay"]
    assert g[0][0] == 0.0 and g[0][1] < -250.0                     # the pair's cancellation at double precision
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
