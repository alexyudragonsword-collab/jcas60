"""Law 6: the vendor tree is a verifiable copy -- full SHAs, every file classified, nothing but import lines changed."""
from __future__ import annotations

import re
from pathlib import Path

import pytest

import vendor_check as V

REPO = Path(__file__).resolve().parent.parent
PINS = {"spec_sources": "da5472183fc4b9447fbbd036d2e2d22ff30f7b0a", "mmw60_trx": "442fabd99c16b0e6cfc4432d35fe64dc053dc515",
        "fmcw_radar": "d1c48bd972fabc693517a30d976ac85feda3739c"}


def test_every_vendored_file_names_its_upstream_at_a_full_sha():
    files = V.vendored_files()
    assert len(files) == 21
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
    assert all(v == "verbatim" for v in kinds.values())                  # rd.py too, since fmcw_radar#44


def test_only_what_is_read_is_vendored():
    """No mmw60 training.py / Frontend60, no fmcw_radar config layer, no second spec_sources."""
    names = {str(p.relative_to(V.VENDOR)) for p in V.VENDOR.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    assert not [n for n in names if "training" in n or "frontend" in n.lower().replace("frontend60.yaml", "") or "schema" in n]
    assert len([n for n in names if n.endswith("reader.py")]) == 1
    # the regulatory check is copied as it is and never reached through a config object (fmcw_radar#42)
    own = (REPO / "src/jcas60/regulatory.py").read_text(encoding="utf-8")
    assert "check_emission(" in own and "_compliance.check(" not in own and "SimpleNamespace" not in own
    # the vendored preamble builder imports numpy and the Golay module beside it, nothing else (mmw60_trx#89)
    pre = (V.VENDOR / "mmw60/ay/preamble.py").read_text(encoding="utf-8")
    imports = [ln for ln in pre.splitlines() if ln.startswith(("import ", "from "))]
    assert imports == ["from __future__ import annotations", "import numpy as np", "from .golay import golay_pair"]
    # arraytx's two files import numpy and each other, nothing of mmw60_trx's own array link (which stands on Frontend60)
    for name, want in (("geometry.py", ["from __future__ import annotations", "from dataclasses import dataclass, replace", "import numpy as np"]),
                       ("beam.py", ["from __future__ import annotations", "import numpy as np", "from .geometry import C0, ArrayGeometry"])):
        text = (V.VENDOR / "mmw60/arraytx" / name).read_text(encoding="utf-8")
        assert [ln for ln in text.splitlines() if ln.startswith(("import ", "from "))] == want, name
    assert not [n for n in names if "arraylink" in n or "array60" in n]
    used = (REPO / "src/jcas60").glob("*.py")
    text = "\n".join(p.read_text(encoding="utf-8") for p in used)
    for mod in ("processing.detection", "processing.windows", "processing import cfar", "processing import rd", "link import compliance", "link.compression import", "ay import preamble", "arraytx import beam", "arraytx import geometry", "processing import cfarstats"):
        assert mod in text, mod
