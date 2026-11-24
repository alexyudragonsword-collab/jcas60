"""Law 6: the vendor tree is a verifiable copy -- full SHAs, every file classified, nothing but import lines changed."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

import vendor_check as V

REPO = Path(__file__).resolve().parent.parent
PINS = {"spec_sources": "f9ce3f2b035d4d83269aa4a4eb33982c53838286", "mmw60_trx": "cfb401aae3d4271f1ed1146f8c1f182065b26fe7",
        "fmcw_radar": "fb0fc165457b4b1ba40e3e777a987f27d842b796"}


def test_every_vendored_file_names_its_upstream_at_a_full_sha():
    files = V.vendored_files()
    assert len(files) == 15
    for f in files:
        m = V.HEADER.match(f.read_text(encoding="utf-8").splitlines()[0])
        assert m, f
        assert re.fullmatch(r"[0-9a-f]{40}", m["sha"]) and PINS[m["repo"]] == m["sha"], f


def test_the_pins_are_the_same_everywhere_they_are_written():
    init = (REPO / "src/jcas60/vendor/__init__.py").read_text(encoding="utf-8")
    ci = (REPO / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    import upstream_report as U
    for repo, sha in PINS.items():
        assert sha in init and sha in ci, repo
    assert U.UPSTREAM["sha"] == PINS["fmcw_radar"]


def test_vendor_check_strict_passes_against_the_siblings():
    sib = V.Siblings()
    if any(sib.path(r) is None for r in PINS):
        pytest.skip("sibling checkouts not present (CI's vendor job checks them out at the pins)")
    results = V.check(sib, {})
    assert V.report(results, strict=True, fail_on_skip=True) == 0
    kinds = {r.rel: r.status for r in results}
    assert sorted(k for k, v in kinds.items() if v == "adapted") == ["fmcw_radar/link/compliance.py", "fmcw_radar/processing/rd.py"]
    assert all(v in ("verbatim", "adapted") for v in kinds.values())


def test_only_what_is_read_is_vendored():
    """No mmw60 training.py / Frontend60, no fmcw_radar config layer, no second spec_sources."""
    names = {str(p.relative_to(V.VENDOR)) for p in V.VENDOR.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    assert not [n for n in names if "training" in n or "frontend" in n.lower().replace("frontend60.yaml", "") or "schema" in n]
    assert len([n for n in names if n.endswith("reader.py")]) == 1
    used = (REPO / "src/jcas60").glob("*.py")
    text = "\n".join(p.read_text(encoding="utf-8") for p in used)
    for mod in ("processing.detection", "processing.windows", "processing import cfar", "processing import rd", "link import compliance", "ay.golay"):
        assert mod in text, mod
