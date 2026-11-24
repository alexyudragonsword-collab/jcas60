"""Verifiable copies of what jcas60 takes from its three upstreams.

A copy, not a fork: every file here carries an attribution header (repo, full
SHA, upstream path) and ``tools/vendor_check.py --strict`` re-derives its
relation to that commit on every CI run.  Only the files this library reads,
and what they import, are copied; code outside this tree imports them through
``jcas60.vendor...`` and edits nothing in it.

Pins
----
* ``spec_sources/`` -- spec_sources @ ``f9ce3f2b035d4d83269aa4a4eb33982c53838286``:
  the reader and five tables (``ieee80211`` -- the Golay sequences, the CEF
  layout, the EDMG SC numbers and PPDU lengths, the DMG timing; ``regulatory``
  -- 47 CFR 15.255, the communication limits and the radar tiers;
  ``frontend60`` -- the 60 GHz radar chip's figures and the TX-RX isolation;
  ``radar_targets`` -- the pedestrian RCS; ``atmosphere`` -- the oxygen
  attenuation).  Every number about the outside world enters here (law 1).

* ``mmw60/`` -- mmw60_trx @ ``cfb401aae3d4271f1ed1146f8c1f182065b26fe7`` (0.3.0):

  - ``ay/golay.py`` (the Golay recursion; numpy only),
  - ``data/spec/survey.yaml`` (its own table: the SC pulse and the converter
    model are read from it by key),
  - ``examples/expected/04_comparison_table.txt`` (its comparison table as
    printed: the ay_sc row's rate is read from it).

  ``ay/training.py`` is NOT copied: it builds the OFDM-rate EDMG-STF / EDMG-CEF,
  and it imports mmw60's ``spec`` shim (its local table plus a second, older
  copy of spec_sources).  The single-carrier Golay training field is built in
  ``jcas60.field`` from the two spec_sources entries until mmw60_trx#89 lands.
  ``Frontend60`` is not copied either (milestone 0's boundary: the same
  numbers, not the same object).

* ``fmcw_radar/`` -- fmcw_radar @ ``fb0fc165457b4b1ba40e3e777a987f27d842b796`` (0.2.32):

  - ``processing/detection.py``, ``processing/cfar.py``,
    ``processing/windows.py`` (verbatim),
  - ``processing/rd.py`` (adapted: one import line -- ``C0`` comes from this
    package's ``__init__`` instead of fmcw_radar's config schema, which is not
    copied; the ``RadarConfig`` name it also imported is an annotation only),
  - ``link/compliance.py`` (adapted: the two import lines -- the
    ``RadarConfig`` annotation import is dropped and ``number`` is
    ``jcas60.sources.fcc_number``, which answers the same
    ``fcc15255.<tier>.<field>`` keys from the vendored spec_sources tables
    instead of fmcw_radar's ``sources.yaml``),
  - ``examples/expected/00_single_target_rd_map.txt`` (its example 00 as
    printed: the two engines' difference is read from it).

  fmcw_radar's budget report of its default config is not a file of its
  repository; it is regenerated from the pinned commit by
  ``tools/upstream_report.py`` and kept in ``jcas60/data/upstream/``.

Local changes
-------------
None in any body: the three import lines above are the adapted copies' whole
difference, which ``tools/vendor_check.py`` classifies as such.  A change to
a vendored body goes upstream first (an issue there), then is re-vendored at
a new pin.
"""
