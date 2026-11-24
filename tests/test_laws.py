"""The seven laws of AGENTS.md, one group of tests each (law 6's vendor check is tests/test_vendor.py,
law 5's fingerprint tests/test_fingerprint.py; the laws' own statements are checked here)."""
from __future__ import annotations

import ast
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

from jcas60 import WORDING, budget, config, frontend, regulatory, table
from jcas60 import fmcw_row as A
from jcas60 import sources as S
from jcas60.vendor.spec_sources import audit_tags

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "src" / "jcas60"
OWN = sorted(p for p in SRC.glob("*.py"))


@pytest.fixture(scope="module")
def t():
    return table.build()


# ----------------------------------------------------------------------------- law 1: numbers come in through spec_sources, tagged
def test_law1_every_reference_read_resolves_and_carries_one_of_the_four_tags(t):
    assert len(t["trace"]) > 40
    for ref, n in t["trace"].items():
        key = ref.partition("#")[0].partition("[")[0]
        assert key in S.table(), ref
        assert n.tag == S.table()[key].source and n.tag in S.TAGS


def test_law1_the_three_tables_are_kept_apart_by_prefix():
    assert all(k.startswith("jcas60.") for k in S.local())
    assert not any(k.startswith("jcas60.") for k in list(S.shared()) + list(S.mmw60()))
    assert all(e.source != "assumed" for e in S.shared().values())          # the shared tables hold the outside world only


def test_law1_own_entries_are_all_assumed_and_say_what_would_move():
    assert S.local() and set(S.assumed_entries()) == set(S.local())
    for k, e in S.local().items():
        assert "换了会动什么" in e.note, k


def test_law1_grep_count_of_each_tag_equals_the_entries():
    assert audit_tags(S.LOCAL) == []
    for tag in S.TAGS:
        n = int(subprocess.run(["grep", "-c", f"source: {tag}", str(S.LOCAL)], capture_output=True, text=True).stdout or 0)
        assert n == sum(1 for e in S.local().values() if e.source == tag), tag


def _literals(path: Path):
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            yield node.lineno, node.value


def _one_digit(v: float) -> bool:
    """0, or a single significant digit (2, 0.5, 1e9, 40): arithmetic, not data."""
    if v == 0:
        return True
    return len(f"{abs(v):.6e}".split("e")[0].rstrip("0").rstrip(".").replace(".", "")) == 1


#: literals with more than one significant digit that are not physical numbers, by file
ALLOWED = {
    "budget.py": {1.5, 4001, 2001, 6001, 2000.0, 80, 60},        # 1.5 = SQNR of a full-scale complex tone; grid sizes; bisection bounds
    "wavesim.py": {16, 64},                                      # the phase-noise grid's length in buffers; pulse-tail guard in chips
    "table.py": {0.75, 97, 3.0103, 1e-40, 100, 180, 34, 92},    # scan grids, -3 dB, formatting widths
    "config.py": set(), "frontend.py": set(), "field.py": set(), "regulatory.py": set(), "fmcw_row.py": set(), "sources.py": set(),
    "__init__.py": set(),
}


def test_law1_no_stray_numbers_in_the_code():
    """A physical number written into the code would bypass the tables.  Outside ``physconst.py`` a literal has one
    significant digit (arithmetic: 2, 0.5, 1e9), is a formatting width inside an f-string, or is listed above."""
    stray = []
    for p in OWN:
        if p.name == "physconst.py":
            continue
        allowed = ALLOWED[p.name]
        src = p.read_text(encoding="utf-8").splitlines()
        for line, v in _literals(p):
            if _one_digit(v) or v in allowed or abs(v) in (1e-9, 1e-300):
                continue
            text = src[line - 1]
            if re.search(r"['\"]", text) and re.search(rf"[:>\s]{re.escape(str(v))}s?\}}", text):
                continue                                           # a width in a format spec
            stray.append((p.name, line, v, text.strip()[:80]))
    assert not stray, stray


# ----------------------------------------------------------------------------- law 2: one set of front-end numbers
def test_law2_one_front_end_instance_serves_every_row(t):
    fe = frontend.shared()
    assert t["setup"].fe is fe and config.setup().fe is fe
    for r in t["rows"][1:]:
        assert r["eirp_dbm"].value == fe.eirp_dbm(r["cap_dbm"]).value
        assert r["backoff_db"].value == fe.backoff_db(r["cap_dbm"]).value


def test_law2_row_a_ran_on_the_same_numbers():
    fe, s = frontend.shared(), config.setup()
    for mine, name in ((fe.tx_power_dbm, "tx_power_dbm"), (fe.tx_gain_dbi, "tx_gain_dbi"), (fe.rx_gain_dbi, "rx_gain_dbi"), (fe.nf_db, "nf_db"),
                       (fe.rx_gain_db, "rx_gain_db"), (fe.p1db_inband_dbm, "p1db_inband_dbm"), (s.rcs_dbsm, "rcs_dbsm"),
                       (s.oxygen_db_per_km, "oxygen_db_per_km"), (s.range_m, "range_m")):
        assert mine.value == A.get(f"inputs.{name}.value"), name
    assert s.velocity_mps.value == A.get("target.velocity_mps")
    # the back-off is one rule, and it gives row A's own figure
    cap = table.row_a()["cap_dbm"]
    assert fe.backoff_db(cap).value == A.get("inputs.tx_backoff_db.value") == 8.0
    assert fe.eirp_dbm(cap).value == A.budget_term("EIRP").value


def test_law2_front_end_numbers_are_read_in_one_module_only():
    """The front-end references appear in frontend.py and nowhere else in the library's own code."""
    keys = ("rx60.radar_iwr6843#tx_power_dbm", "rx60.radar_iwr6843#nf_db", "jcas60.antenna_gain_dbi", "jcas60.rx_gain_db",
            "isolation.tx_rx", "tgad_em.pn_params", "rx60.radar_iwr6843_fig7_2#p1db_inband_dbm@")
    for p in OWN:
        if p.name == "frontend.py":
            continue
        text = p.read_text(encoding="utf-8")
        assert not [k for k in keys if k in text], p.name


def test_law2_the_lo_is_named_as_the_one_number_not_shared(t):
    assert "the one front-end number the rows do not share" in table.render(t)


# ----------------------------------------------------------------------------- law 3: each row at its own limit, both in the table
def test_law3_every_row_states_its_rule_and_limit_and_both_readings_are_rows(t):
    rows = t["rows"]
    assert [r["limit"] for r in rows] == ["radar limit", "communication limit", "radar limit", "communication limit", "radar limit"]
    for r in rows:
        assert r["cap_dbm"].tag == "spec" and r["rule"].startswith("47 CFR 15.255")
        assert r["eirp_dbm"].value <= r["cap_dbm"].value + 1e-12            # never above its own limit
    text = table.render(t)
    assert "47 CFR 15.255(c)(1)(i)" in text and "47 CFR 15.255(c)(2)" in text


def test_law3_no_row_runs_at_another_rows_limit(t):
    """The radar rows are backed off to the radar tier; the communication rows are not held to it."""
    comm, radar = t["readings"]
    assert comm.cap_dbm.value == S.num("fcc15255.eirp_avg_dbm").value and radar.cap_dbm.value == S.num("fcc15255.radar_57_61p56_off_time#peak_eirp_dbm").value
    rows = t["rows"]
    assert rows[2]["eirp_dbm"].value == radar.cap_dbm.value < rows[1]["eirp_dbm"].value == frontend.shared().eirp_max_dbm.value < comm.cap_dbm.value
    assert rows[0]["eirp_dbm"].value == rows[0]["cap_dbm"].value == 14.0


def test_law3_the_vendored_check_reproduces_fmcw_radars_own_report():
    """fmcw_radar's check, as vendored and fed from spec_sources by this library's shim, on fmcw_radar's own default
    frame: every tier and the verdict equal its report."""
    c, rad = A.get("compliance"), A.get("radar")
    lo, hi = c["sweep_ghz"]
    res = regulatory.radar_check((lo * 1e9, hi * 1e9), int(rad["n_chirps"]), rad["chirp_us"] * 1e-6, (rad["chirp_us"] + rad["idle_us"]) * 1e-6,
                                 c["frame_period_s"], frontend.shared().tx_power_dbm, frontend.shared().tx_gain_dbi).as_dict()
    from upstream_report import differences                     # equal to a relative 1e-9: the frame's length is summed differently
    keys = ("governing", "eirp_cap_dbm", "eirp_dbm", "conducted_dbm", "tx_backoff_db", "min_frame_period_s", "frame_on_s", "verdict", "tiers")
    assert differences({k: res[k] for k in keys}, {k: c[k] for k in keys}) == []


# ----------------------------------------------------------------------------- law 4: two roads, and their difference is reported
def test_law4_the_two_roads_agree_within_1_db_on_every_ay_row(t):
    for r in t["rows"][1:]:
        d = r["paths_delta_db"].value
        assert d == pytest.approx(r["two_paths"]["sim_db"].value - r["two_paths"]["closed_grid_db"].value)
        assert abs(d) < 1.0, (r["name"], r["limit"], d)
        assert abs(d) < 0.3                                                   # where they are today
        assert r["two_paths"]["cfar_target"]
        assert abs(r["two_paths"]["range_read_m"] - 5.0) < 0.02 and abs(r["two_paths"]["velocity_read_mps"] - 2.0) < 0.05


def test_law4_the_difference_is_a_column_of_the_table(t):
    text = table.render(t)
    assert "== 3. the two roads" in text
    for r in t["rows"][1:]:
        assert f"{r['paths_delta_db'].value:+10.2f} dB" in text
    assert f"{t['rows'][0]['paths_delta_db'].value:+10.2f} dB" in text            # row A: fmcw_radar's own two engines


def test_law4_the_receivers_cost_has_both_roads_where_the_closed_form_claims_to_hold(t):
    """Quantisation and phase noise: closed form and simulation within 0.2 dB wherever the thermal noise dithers the
    ADC (at least half an LSB) -- and where it does not, the table says so and the split explains the gap."""
    for blk in t["monostatic"]:
        for ck in blk["checks"]:
            assert abs(ck["pn_rise_sim_db"] - ck["pn_rise_closed_db"]) < 0.2
            if ck["noise_lsb"] >= 0.5:
                assert abs(ck["adc_rise_sim_db"] - ck["adc_rise_closed_db"]) < 0.2
            else:
                assert ck["repeating_share"] > 0.5 and ck["adc_rise_sim_db"] < ck["adc_rise_closed_db"] - 1.0
            assert abs(ck["adc_rise_sim_db"] - ck["adc_rise_varying_db"]) < 0.2    # the part that is not the same in every packet is the floor
            assert abs(ck["snr_linear_db"] - blk["snr_no_leak_db"]) < 0.2          # a linear receiver does not lose sensitivity to the leakage
    for blk in t["bistatic"]:
        assert abs(blk["pn_rise_sim_db"] - blk["pn_rise_closed_db"]) < 1.0


# ----------------------------------------------------------------------------- law 7: wording
def test_law7_no_802_11bf_entry_no_conformance_claim():
    assert not [k for k in S.shared() if "80211bf" in k or k.startswith("bf.")], "spec_sources has an 802.11bf entry now: law 7's wording can be revisited"
    assert "no claim of 802.11bf conformance" in WORDING and WORDING in table.render()
    claim = re.compile(r"(符合|遵循|compliant with|conform\w* to|conformant)[^。.\n]{0,24}802\.11bf", re.I)
    files = [p for p in REPO.rglob("*") if p.is_file() and p.suffix in (".py", ".md", ".txt", ".yaml", ".yml", ".toml")
             and ".git" not in p.parts and "vendor" not in p.parts]
    hits = []
    for p in files:
        for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
            if claim.search(line) and not re.search(r"不说|不能说|no claim|not claim|does not claim|nothing here claims", line):
                hits.append((str(p.relative_to(REPO)), i))
    assert not hits, hits


# ----------------------------------------------------------------------------- the field itself
def test_the_golay_pair_is_complementary_within_the_cyclic_prefix_and_the_recursion_gives_the_tables():
    from jcas60 import field as F
    for n in (128, 64, 32):
        a, b = F.from_recursion(n)
        assert np.array_equal(a, F.tabulated(f"ga{n}")) and np.array_equal(b, F.tabulated(f"gb{n}"))
    f = F.l_cef()
    assert f.chips.size == 1280 and f.n_correlated == 1024 and np.allclose(np.abs(f.chips), 1.0)
    for d in (0, 1, 57, 128):                                    # a delayed copy: each window holds a cyclic shift
        got = sum(np.vdot(ref, f.chips[w0 - d: w0 - d + len(ref)]) for w0, ref in f.windows)
        assert abs(got) < 1e-9 if d else abs(got - 1024) < 1e-9
    assert budget.geometry(config.setup(), "monostatic").delay_s * config.setup().chip_rate_hz.value < f.zero_sidelobe_chips
    # Gu512 and Gv512 as a pair: their cyclic autocorrelations sum to zero at EVERY lag but 0.  (Inside the 128-chip
    # prefix each of them is already clean alone -- it is made of Ga128 / Gb128 -- so only the lags beyond tell a
    # true pair from one sequence sent twice; a break-check found that the test above did not.)
    (_, u), (_, v) = f.windows
    cyc = lambda x: np.fft.ifft(np.abs(np.fft.fft(x)) ** 2)      # noqa: E731
    alone, pair = np.abs(cyc(u))[1:], np.abs(cyc(u) + cyc(v))[1:]
    assert pair.max() < 1e-9 and alone[:128].max() < 1e-9 and alone.max() > 10
