"""Law 5: the fingerprint is recorded, kept, and exact."""
from __future__ import annotations

import json

import fingerprint as F


def test_the_kept_record_is_this_trees_fingerprint():
    kept = json.loads(F.KEPT.read_text(encoding="utf-8"))
    now = F.record()
    assert len(kept) > 400
    assert F.compare(kept, now) == [], "numbers moved: explain them in the commit message with (numbers-move: why) and re-record tools/fingerprint.json"


def test_compare_sees_one_moved_value_and_a_missing_one():
    a = {"x": "1.0", "y": "sha256:aa"}
    assert F.compare(a, dict(a)) == []
    assert F.compare(a, {"x": "1.0000000000000002", "y": "sha256:aa"}) == ["x: 1.0 -> 1.0000000000000002"]
    assert F.compare(a, {"x": "1.0"}) == ["y: only in the first"]
    # another environment: scalars to a relative 1e-6, hashes left out; a real move is still seen
    b = {"_environment": "elsewhere", "x": "1.0000000000000002 [derived]", "y": "sha256:bb"}
    assert F.compare({"_environment": "here", "x": "1.0 [derived]", "y": "sha256:aa"}, b) == []
    assert len(F.compare({"_environment": "here", "x": "1.01 [derived]", "y": "sha256:aa"}, b)) == 1
    assert len(F.compare({"_environment": "here", "x": "1.0 [spec]", "y": "sha256:aa"}, b)) == 1
