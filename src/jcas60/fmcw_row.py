"""Row A: fmcw_radar's default radar, as fmcw_radar reports it.

Nothing is recomputed: every number is read out of the budget report that
``tools/upstream_report.py`` regenerated from the pinned commit
(``data/upstream/fmcw_radar_budget.json``), and the two engines' difference
out of its example 00 as printed (vendored).  A test holds every cell of the
table's row A equal to the report, item by item.

Tags: the report's ``inputs`` carry their own (its ``config`` reads as this
table's ``assumed``: fmcw_radar's design choice); the limits of its
regulatory check are ``spec``; a design number of its config is ``assumed``;
anything its engines computed is ``derived`` on that basis.
"""
from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path

from .sources import Num
from .vendor.fmcw_radar import EXAMPLE_00

RECORD = Path(__file__).resolve().parent / "data" / "upstream" / "fmcw_radar_budget.json"


@lru_cache(maxsize=1)
def _loaded() -> dict:
    return json.loads(RECORD.read_text(encoding="utf-8"))


def provenance() -> dict:
    return dict(_loaded()["_provenance"])


def report() -> dict:
    return _loaded()["report"]


def get(path: str):
    """A value of the report by dotted path (``radar.max_range_m``)."""
    v = report()
    for part in path.split("."):
        v = v[int(part)] if isinstance(v, list) else v[part]
    return v


def rep(path: str, tag: str = "derived", unit: str = "") -> Num:
    """A scalar of the report as a Num.  ``derived`` numbers rest on fmcw_radar's own design choices."""
    return Num(float(get(path)), tag, f"fmcw_radar report: {path}", unit, "assumed" if tag == "derived" else "")


def rep_input(name: str, unit: str = "") -> Num:
    e = get(f"inputs.{name}")
    tag = "assumed" if e["source"] == "config" else e["source"]
    return Num(float(e["value"]), tag, f"fmcw_radar report: inputs.{name}", unit, "spec" if tag == "derived" else "")


def budget_term(term: str) -> Num:
    """A row of the report's radar-equation table, by its label."""
    for row in get("budget"):
        if row["term"] == term:
            return Num(float(row["value"]), "derived", f"fmcw_radar report: budget[{term}]", row["unit"], "assumed")
    raise KeyError(term)


def governing_tier() -> dict:
    return next(t for t in get("compliance.tiers") if t["tier"] == get("compliance.governing"))


def engines_difference_db() -> Num:
    """The time engine's cell SNR minus the statistical engine's, as fmcw_radar's example 00 prints it."""
    text = EXAMPLE_00.read_text(encoding="utf-8")
    m = re.search(r"measured cell SNR ([-\d.]+) dB vs statistical ([-\d.]+) dB expected on the grid.*?difference ([-+\d.]+) dB", text)
    if not m:
        raise ValueError("fmcw_radar's example 00 no longer prints the two engines' difference in the form read here")
    return Num(float(m.group(3)), "derived", "fmcw_radar examples/expected/00_single_target_rd_map.txt: time engine - statistical engine",
               "dB", "assumed")
