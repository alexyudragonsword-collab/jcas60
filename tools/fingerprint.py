#!/usr/bin/env python3
"""Numerical fingerprint: prove a change did not move a number (law 5).

Tests assert relations and tolerances; a result that moves a little and still
satisfies them passes.  A refactor is therefore checked bit for bit (the
siblings' ``tools/fingerprint.py`` method): every number of the comparison
table, the closed-form budgets behind it, the transmitted field and the
waveform simulation's maps are recorded -- scalars by ``repr``, arrays by
SHA-256 -- and two records are compared exactly.

    python tools/fingerprint.py record out.json      # record this tree
    python tools/fingerprint.py compare a.json b.json  # list every differing value; exit 1 if there is one
    python tools/fingerprint.py check                # this tree against the record kept in the repository

The record kept in the repository is ``tools/fingerprint.json``: what the
numbers were when the last change that meant to move them was committed.
Bit-exact means the same library versions, so CI compares the base commit
and the head commit recorded side by side in one job.  The kept record names
the environment that wrote it (``_environment``); ``check`` is exact in that
environment -- which is how a refactor is proved locally -- and in another
one it compares the scalars to a relative 1e-6 and leaves the array hashes
out, saying so.  A change meant to move numbers lists them in its commit
message with the ``(numbers-move: why)`` marker CI looks for, and re-records
the kept file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
KEPT = Path(__file__).resolve().parent / "fingerprint.json"


def digest(a) -> str:
    a = np.ascontiguousarray(a)
    return f"sha256:{hashlib.sha256(a.tobytes()).hexdigest()} {a.dtype} {list(a.shape)}"


def _flat(prefix: str, v, out: dict) -> None:
    from jcas60.sources import Num
    if isinstance(v, Num):
        out[prefix] = f"{v.value!r} [{v.tag}{'/' + v.basis if v.basis else ''}]"
    elif isinstance(v, np.ndarray):
        out[prefix] = digest(v)
    elif isinstance(v, dict):
        for k in sorted(v, key=str):
            _flat(f"{prefix}.{k}", v[k], out)
    elif isinstance(v, (list, tuple)):
        for i, x in enumerate(v):
            _flat(f"{prefix}[{i}]", x, out)
    elif hasattr(v, "__dataclass_fields__"):
        for k in v.__dataclass_fields__:
            _flat(f"{prefix}.{k}", getattr(v, k), out)
    elif isinstance(v, (float, np.floating)):
        out[prefix] = repr(float(v))
    elif isinstance(v, (int, np.integer, bool, str)) or v is None:
        out[prefix] = repr(v)
    else:
        out[prefix] = f"<{type(v).__name__}>"


def environment() -> str:
    import platform

    import scipy
    return f"python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.machine()}"


def record() -> dict:
    from jcas60 import budget, config, table
    from jcas60 import wavesim as W
    out: dict = {"_environment": environment()}
    t = table.build()
    s = t["setup"]
    _flat("table.rows", [{k: v for k, v in r.items() if k not in ("link",)} for r in t["rows"]], out)
    _flat("table.links", [r["link"] for r in t["rows"][1:]], out)
    _flat("table.monostatic", [{k: v for k, v in b.items() if k != "reading"} for b in t["monostatic"]], out)
    _flat("table.compression", {k: v for k, v in t["compression"].items() if k != "env"}, out)
    _flat("sim.unit_envelope", t["compression"]["env"], out)
    _flat("table.bistatic", t["bistatic"], out)
    _flat("table.golay", t["golay"], out)
    _flat("table.golay_trn", {"blocks_2_to_5": t["golay_trn"], "all_blocks_at_rest": t["golay_trn_all_blocks"]}, out)
    _flat("table.rounding_dbc", t["rounding_dbc"], out)
    _flat("table.trn_sweep", t["trn_sweep"], out)
    _flat("field.trn_chips", np.asarray(s.trn_field.chips), out)
    _flat("table.pulse_sidelobes", t["pulse_sidelobes"], out)
    _flat("table.width_measured_chips", t["width_measured_chips"], out)
    _flat("table.sensitivity", t["sensitivity"], out)
    _flat("table.text", hashlib.sha256(table.render(t).encode()).hexdigest(), out)
    # the pieces under the table, on their own
    _flat("field.chips", np.asarray(s.field.chips), out)
    sim = W.Simulator(config.setup())
    _flat("sim.tx_f", sim.tx_f, out)
    g = budget.geometry(s, "monostatic")
    link = t["rows"][1]["link"]
    paths = [W.Path(link.p_rx_dbm.value, g.delay_s, g.path_rate(s.velocity_mps.value)), W.Path(-28.0, s.leak_delay_s.value)]
    for name, opt in (("linear", W.Options()), ("adc", W.Options(adc_enob=s.adc_enob.value)), ("pn_shared", W.Options(lo="shared")),
                      ("noiseless", W.Options(noise=False))):
        fr = sim.frame(paths, opt, 0)
        _flat(f"sim.frame.{name}.cir", fr.cir, out)
        _flat(f"sim.frame.{name}.rd", fr.rd, out)
    _flat("budget.pn_residual_shared", budget.pn_residual_shared(s, s.leak_delay_s.value), out)
    _flat("budget.pn_residual_independent", budget.pn_residual_independent(s, s.field.n_correlated / s.chip_rate_hz.value), out)
    _flat("budget.snr_required_db", t["snr_required_db"], out)
    return out


def _close(x: str, y: str) -> bool:
    """Two recorded scalars equal to a relative 1e-6 (and the same tag); array hashes are not comparable this way."""
    if x.startswith("sha256:") or y.startswith("sha256:"):
        return True
    (xv, _, xt), (yv, _, yt) = x.partition(" "), y.partition(" ")
    try:
        return xt == yt and math.isclose(float(xv), float(yv), rel_tol=1e-6, abs_tol=1e-9)
    except ValueError:
        return x == y


def compare(a: dict, b: dict, exact: bool | None = None) -> list[str]:
    """Every value that differs.  ``exact`` None: exact when both records come from one environment."""
    if exact is None:
        exact = a.get("_environment") == b.get("_environment")
    a = {k: v for k, v in a.items() if k != "_environment"}
    b = {k: v for k, v in b.items() if k != "_environment"}
    lines = [f"{k}: only in the {'first' if k in a else 'second'}" for k in sorted(set(a) ^ set(b))]
    lines += [f"{k}: {a[k]} -> {b[k]}" for k in sorted(set(a) & set(b)) if (a[k] != b[k] if exact else not _close(a[k], b[k]))]
    return lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("record").add_argument("out", type=Path)
    c = sub.add_parser("compare")
    c.add_argument("a", type=Path)
    c.add_argument("b", type=Path)
    sub.add_parser("check")
    args = ap.parse_args(argv)
    if args.cmd == "record":
        fp = record()
        args.out.write_text(json.dumps(fp, indent=0, sort_keys=True) + "\n", encoding="utf-8")
        print(f"recorded {len(fp)} values in {args.out}")
        return 0
    a, b = (json.loads(KEPT.read_text()), record()) if args.cmd == "check" else (json.loads(args.a.read_text()), json.loads(args.b.read_text()))
    diff = compare(a, b)
    if a.get("_environment") != b.get("_environment"):
        print(f"different environments ({a.get('_environment')} / {b.get('_environment')}): scalars to 1e-6, array hashes not compared")
    for line in diff:
        print(line[:300])
    print(f"{len(diff)} values differ of {len(set(a) | set(b))}" if diff else f"identical: {len(a)} values")
    return 1 if diff else 0


if __name__ == "__main__":
    raise SystemExit(main())
