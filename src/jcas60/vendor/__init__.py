"""Verifiable copies of what jcas60 takes from its three upstreams.

A copy, not a fork: every file here carries an attribution header (repo, full
SHA, upstream path) and ``tools/vendor_check.py --strict`` re-derives its
relation to that commit on every CI run.  Only the files this library reads,
and what they import, are copied; code outside this tree imports them through
``jcas60.vendor...`` and edits nothing in it.

Pins
----
* ``spec_sources/`` -- spec_sources @ ``3a6e1905ef7c48c01338cf9cd808e3fd4bdad8a7``:
  the reader and five tables (``ieee80211`` -- the Golay sequences, the CEF
  layout, the EDMG SC numbers and PPDU lengths, the DMG timing; ``regulatory``
  -- 47 CFR 15.255, the communication limits and the radar tiers;
  ``frontend60`` -- the 60 GHz radar chip's figures and the TX-RX isolation;
  ``radar_targets`` -- the pedestrian RCS; ``atmosphere`` -- the oxygen
  attenuation, with a point at 60.48 GHz).  ``ieee80211`` also holds the
  802.11bf entries of spec_sources#54 and #57 (sensing types, scheduling,
  the monostatic, bistatic and multistatic PPDUs, bandwidth, the TRN
  fields), which ``jcas60.bf`` checks the sensing frame against.  Since this
  pin ``dmg.golay_sequences`` holds the final standard's 128 / 64 pairs (the
  earlier ones were a 2010 draft's; spec_sources#57) and the length-32 pair,
  unverified, is a separate entry this library does not read.  Every number
  about the outside world enters here (law 1).

* ``mmw60/`` -- mmw60_trx @ ``442fabd99c16b0e6cfc4432d35fe64dc053dc515`` (0.3.0 plus its unreleased
  ``ay/preamble.py``, mmw60_trx#89):

  - ``ay/golay.py`` (the Golay recursion; numpy only),
  - ``ay/preamble.py`` (the standard single-carrier training field: Gu512 /
    Gv512 from Ga128 / Gb128, the CEF, its cyclic prefix, the pi/2 rotation;
    numpy and ``golay.py`` only, the spec_sources entries are its arguments),
  - ``data/spec/survey.yaml`` (its own table: the SC pulse and the converter
    model are read from it by key),
  - ``examples/expected/04_comparison_table.txt`` (its comparison table as
    printed: the ay_sc row's rate is read from it).

  ``ay/training.py`` is NOT copied: it builds the OFDM-rate EDMG-STF / EDMG-CEF,
  and it imports mmw60's ``spec`` shim (its local table plus a second, older
  copy of spec_sources).  ``Frontend60`` is not copied either (milestone 0's
  boundary: the same numbers, not the same object).

* ``fmcw_radar/`` -- fmcw_radar @ ``c303566180c2ef22bd2bf6ffaa1f08e59174e110`` (0.2.33):

  - ``processing/detection.py``, ``processing/cfar.py``,
    ``processing/windows.py`` (verbatim),
  - ``processing/rd.py`` (adapted: one import line -- ``C0`` comes from this
    package's ``__init__`` instead of fmcw_radar's config schema, which is not
    copied; the ``RadarConfig`` name it also imported is an annotation only),
  - ``link/compliance.py`` (verbatim since fmcw_radar#42: nothing of
    fmcw_radar is imported when it loads; this library calls its
    ``check_emission`` with ``jcas60.sources.fcc_number`` as the function
    that looks a limit up, and never its ``check(cfg)``),
  - ``link/compression.py`` (verbatim; numpy only: the Rapp limiter the
    simulation puts in the receiver and the small-signal theory the closed
    form uses for the same thing),
  - ``examples/expected/00_single_target_rd_map.txt`` (its example 00 as
    printed: the two engines' difference is read from it).

  fmcw_radar's budget report of its default config is not a file of its
  repository; it is regenerated from the pinned commit by
  ``tools/upstream_report.py`` and kept in ``jcas60/data/upstream/``.

Local changes
-------------
None in any body: the one import line of ``processing/rd.py`` is the adapted
copy's whole difference, which ``tools/vendor_check.py`` classifies as such.
A change to a vendored body goes upstream first (an issue there), then is
re-vendored at a new pin.
"""
