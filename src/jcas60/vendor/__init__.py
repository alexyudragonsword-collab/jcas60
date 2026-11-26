"""Verifiable copies of what jcas60 takes from its three upstreams.

A copy, not a fork: every file here carries an attribution header (repo, full
SHA, upstream path) and ``tools/vendor_check.py --strict`` re-derives its
relation to that commit on every CI run.  Only the files this library reads,
and what they import, are copied; code outside this tree imports them through
``jcas60.vendor...`` and edits nothing in it.

Pins
----
* ``spec_sources/`` -- spec_sources @ ``da5472183fc4b9447fbbd036d2e2d22ff30f7b0a``:
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
  unverified, is a separate entry this library does not read.  With
  spec_sources#59 the pin also holds two measured 60 GHz compression curves and
  a second radar chip's datasheet figures (``rx60.*``), the length-256 / 512
  Golay pairs (``edmg.golay_sequences_long``), the bistatic BRP exchange's
  fields, and BRPIFS corrected from 20 to 44 us.  Every number about the
  outside world enters here (law 1).

* ``mmw60/`` -- mmw60_trx @ ``442fabd99c16b0e6cfc4432d35fe64dc053dc515`` (0.3.0 plus its unreleased
  ``ay/preamble.py``, mmw60_trx#89):

  - ``ay/golay.py`` (the Golay recursion; numpy only),
  - ``ay/preamble.py`` (the standard single-carrier training field: Gu512 /
    Gv512 from Ga128 / Gb128, the CEF, its cyclic prefix, the pi/2 rotation;
    numpy and ``golay.py`` only, the spec_sources entries are its arguments),
  - ``data/spec/survey.yaml`` (its own table: the SC pulse and the converter
    model are read from it by key),
  - ``examples/expected/04_comparison_table.txt`` (its comparison table as
    printed: the ay_sc row's rate is read from it),
  - ``data/clock_profiles/pllsim_wu19_x10_60p48g_floor150.yaml`` (its x10 LO's
    phase-noise profile: the one fmcw_radar's ``radar_60g_x10lo.yaml`` reads,
    for the study that gives every row that LO; the table's rows do not use it),
  - ``arraytx/geometry.py``, ``arraytx/beam.py`` (array_tx's array geometry and
    beam functions as mmw60_trx vendors them, at array_tx @ ``286609e``: their
    headers keep both attributions; numpy only, and each other.  The transmit
    array of ``jcas60.arrays`` is built on them: the geometry, the array factor,
    the element-by-element steering and the far-field sum).

  mmw60_trx's own array link (``arraylink.py``, ``frontend/array60.py``) is NOT
  copied: it stands on ``Frontend60`` and the whole of wifitrx.

  ``ay/training.py`` is NOT copied: it builds the OFDM-rate EDMG-STF / EDMG-CEF,
  and it imports mmw60's ``spec`` shim (its local table plus a second, older
  copy of spec_sources).  ``Frontend60`` is not copied either (milestone 0's
  boundary: the same numbers, not the same object).

* ``fmcw_radar/`` -- fmcw_radar @ ``d1c48bd972fabc693517a30d976ac85feda3739c`` (0.2.34):

  - ``processing/detection.py``, ``processing/cfar.py``,
    ``processing/windows.py`` (verbatim),
  - ``processing/rd.py`` (verbatim since fmcw_radar#44: nothing of its config
    layer is imported when it loads; this library calls ``doppler_fft``,
    ``power_db``, ``read_peak`` and ``noise_floor``, none of which takes a
    config),
  - ``processing/cfarstats.py`` (verbatim, new with fmcw_radar#44: what a
    map's windows and zero-padding do to a CA-CFAR -- the SNR a steady target
    needs for a given Pd on that map.  Imports ``detection``, ``cfar`` and
    ``windows`` beside it, numpy and scipy),
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
  ``tools/upstream_report.py`` and kept in ``jcas60/data/upstream/``, and so
  are the reports of two other configs of the same commit: ``x10lo`` (every
  row on one LO, a study) and ``leakage`` (the default radar with its leakage
  model on, fmcw_radar#44: row A's leakage in the table's section 5).

Local changes
-------------
None: every file is verbatim (``tools/vendor_check.py --strict``: 21 files,
21 verbatim).  A change to a vendored body goes upstream first (an issue there), then is
re-vendored at a new pin.
"""
