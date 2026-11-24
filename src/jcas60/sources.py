"""Where every number comes from (law 1).

Three tables, read by key through spec_sources' own reader:

* the vendored spec_sources tables -- everything about the outside world;
* ``data/assumed.yaml`` -- this library's own choices (``jcas60.*``), each
  ``assumed`` with what would move if it were chosen otherwise;
* mmw60_trx's own table, vendored -- its preset choices, read by key.

A number travels as a :class:`Num`: the value, its tag and the reference it
was read by.  A number computed here is ``derive(...)``: tagged ``derived``
and carrying the weakest tag among its inputs, so a result that leans on an
assumption prints as ``[derived/assumed]``.  Every reference read is kept in
``TRACE`` (the example prints it; a test resolves each one again).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .vendor import spec_sources as _shared
from .vendor.mmw60 import SURVEY as _MMW60_TABLE

LOCAL = Path(__file__).resolve().parent / "data" / "assumed.yaml"
TAGS = _shared.SOURCES                               # spec, survey, assumed, derived
_WEAKNESS = {"spec": 0, "derived": 0, "survey": 1, "assumed": 2}

#: every reference read since the last ``TRACE.clear()``: reference -> Num
TRACE: dict[str, "Num"] = {}


@dataclass(frozen=True)
class Num:
    value: float
    tag: str                 # spec | survey | assumed | derived
    ref: str                 # the reference it was read by; for a derived number, its formula
    unit: str = ""
    basis: str = ""          # derived only: the weakest tag among the inputs (spec | survey | assumed)

    def __post_init__(self):
        if self.tag not in TAGS:
            raise ValueError(f"{self.ref}: tag must be one of {TAGS}, got {self.tag!r}")

    def __float__(self) -> float:
        return float(self.value)

    @property
    def label(self) -> str:
        return f"[{self.tag}/{self.basis}]" if self.tag == "derived" and self.basis not in ("", "spec") else f"[{self.tag}]"

    def show(self, fmt: str = "g") -> str:
        return f"{format(self.value, fmt)} {self.label}"


@lru_cache(maxsize=1)
def shared() -> dict:
    return _shared.all_entries()


@lru_cache(maxsize=1)
def local() -> dict:
    return _shared.load_table(LOCAL)


@lru_cache(maxsize=1)
def mmw60() -> dict:
    return _shared.load_table(_MMW60_TABLE)


@lru_cache(maxsize=1)
def table() -> dict:
    """The three tables as one; a key in two of them is an error (the reader's rule)."""
    return _shared.merge_tables(shared(), local(), mmw60())


def entry(key: str):
    return table()[key]


def raw(reference: str):
    """The value behind a reference, whatever its type (a list, a text); recorded like ``num``."""
    v = _shared.ref(reference, table())
    e = entry(reference.partition("#")[0])
    TRACE[reference] = Num(float("nan"), e.source, reference, e.unit)
    return v


def num(reference: str, unit: str = "") -> Num:
    """One scalar by reference: ``key``, ``key#field`` or ``key#field@axis=value``."""
    v = _shared.ref(reference, table())
    e = entry(reference.partition("#")[0])
    n = Num(float(v), e.source, reference, unit or e.unit)
    TRACE[reference] = n
    return n


def derive(value: float, formula: str, *inputs: Num, unit: str = "") -> Num:
    worst = max((i.basis or i.tag for i in inputs), key=lambda t: _WEAKNESS[t], default="spec")
    return Num(float(value), "derived", formula, unit, "spec" if worst == "derived" else worst)


def fcc_number(key: str):
    """``fcc15255.<tier>.<field>`` (fmcw_radar's local key for a field of
    spec_sources' ``fcc15255.<tier>`` entry) or a plain ``fcc15255.<name>``:
    what the vendored ``link/compliance.py`` asks for, answered from the
    vendored spec_sources tables.  Returns ``.value`` and ``.source``."""
    if key in shared():
        n = num(key)
    else:
        head, _, field = key.rpartition(".")
        n = num(f"{head}#{field}")
    return _FccNumber(n.value, n.tag)


@dataclass(frozen=True)
class _FccNumber:
    value: float
    source: str


def assumed_entries() -> dict:
    """This library's own choices, for the list the example prints."""
    return {k: e for k, e in local().items() if e.source == "assumed"}
