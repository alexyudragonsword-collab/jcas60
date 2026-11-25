"""jcas60: sensing with a 60 GHz communication waveform against a dedicated FMCW radar, on one set of front-end numbers.

Milestone 0 is one comparison table (``jcas60.table``): fmcw_radar's default
radar as it reports itself, and the 802.11ay single-carrier Golay training
field used for sensing, monostatic and bistatic.  The wording is "sensing
with the ay waveform": nothing here claims conformance to 802.11bf; what the
sensing frame was checked against, entry by entry, is ``jcas60.bf`` (law 7).
"""

__version__ = "0.1.0"

#: first line of every example (the siblings' convention): what kind of engine produced the numbers
ENGINE_BANNER = ("engine: closed-form radar budget + chip-level waveform simulation (linear channel, point target; "
                 "no PA / LNA nonlinearity) -- compression conclusions are closed-form only")
#: the wording law 7 allows: the frame is checked against spec_sources' 802.11bf entries item by item, and that is all
WORDING = ("sensing with the 802.11ay waveform (no claim of 802.11bf conformance: section 11 holds what was checked against "
           "spec_sources' 802.11bf entries, item by item)")
