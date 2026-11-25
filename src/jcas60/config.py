"""The one description of milestone 0's comparison (laws 1 and 2).

``setup()`` reads every number by reference -- the outside world from
spec_sources, mmw60_trx's preset choices from its table, this library's own
from ``data/assumed.yaml`` -- into one frozen :class:`Setup` holding the
shared :class:`~jcas60.frontend.FrontEnd` instance.  The budget, the
simulation, the regulatory check and the table read this object and nothing
else; the example builds it once.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from . import field as F
from . import sources as S
from .frontend import FrontEnd, shared
from .physconst import C0
from .sources import Num, derive


@dataclass(frozen=True)
class Setup:
    fe: FrontEnd
    # the target (the same one fmcw_radar's default config has)
    range_m: Num
    velocity_mps: Num
    rcs_dbsm: Num
    oxygen_db_per_km: Num
    # the ay single-carrier waveform
    carrier_hz: Num
    chip_rate_hz: Num
    channel_width_hz: Num
    rolloff: Num
    field: F.Field
    packet_chips: Num
    trn_field: F.Field             # the same packet with an EDMG TRN-R field behind it, sensed on that field
    trn_units: Num
    trn_packet_chips: Num
    sifs_s: Num
    # the sensing burst
    n_packets: Num
    interval_s: Num
    frames_per_s: Num
    doppler_window: str
    # the receiver
    adc_osr: Num
    desense_allowed_db: Num
    adc_enob_ref: Num
    adc_enob_per_octave: Num
    adc_enob_ref_hz: float
    adc_backoff_db: Num
    rx_rapp_p: Num
    leak_delay_s: Num
    # detection
    pd: Num
    # bistatic geometry
    bi_r_tx_m: Num
    bi_r_rx_m: Num
    bi_baseline_m: Num
    # simulation
    sim: dict

    @property
    def wavelength_m(self) -> float:
        return C0 / self.carrier_hz.value

    @property
    def band_hz(self) -> tuple[float, float]:
        half = self.channel_width_hz.value / 2
        return self.carrier_hz.value - half, self.carrier_hz.value + half

    @property
    def chip_m(self) -> float:
        """One chip of two-way delay as range: c / (2 chip rate)."""
        return C0 / (2 * self.chip_rate_hz.value)

    @property
    def packet_s(self) -> Num:
        return derive(self.packet_chips.value / self.chip_rate_hz.value, "packet chips / chip rate",
                      self.packet_chips, self.chip_rate_hz, unit="s")

    def packet_seconds(self, field: F.Field) -> Num:
        """The sensing packet's duration when ``field`` is the one sensed on (the TRN field rides behind the packet)."""
        chips = self.trn_packet_chips if field is self.trn_field else self.packet_chips
        return derive(chips.value / self.chip_rate_hz.value, "packet chips / chip rate", chips, self.chip_rate_hz, unit="s")

    @property
    def adc_fs_hz(self) -> Num:
        return derive(self.adc_osr.value * self.chip_rate_hz.value, "samples per chip x chip rate",
                      self.adc_osr, self.chip_rate_hz, unit="Hz")

    @property
    def adc_enob(self) -> Num:
        """mmw60_trx's converter model at this ADC's rate: enob_ref - bits_per_octave log2(fs / fs_ref)."""
        import math
        v = self.adc_enob_ref.value - self.adc_enob_per_octave.value * math.log2(self.adc_fs_hz.value / self.adc_enob_ref_hz)
        return derive(v, "enob_ref - bits_per_octave log2(fs / 1 GS/s)", self.adc_enob_ref, self.adc_enob_per_octave,
                      self.adc_fs_hz, unit="bit")


@lru_cache(maxsize=1)
def setup() -> Setup:
    channel = int(S.num("jcas60.channel#number").value)
    centres = S.raw("ieee80211ay.channel_centres_ghz")
    carrier = Num(float(centres[channel - 1]) * 1e9, S.entry("ieee80211ay.channel_centres_ghz").source,
                  f"ieee80211ay.channel_centres_ghz[{channel}]", "Hz")
    S.TRACE[carrier.ref] = carrier
    fields = S.raw("jcas60.sensing_packet#fields")
    chips = [S.num(f"edmg.sc_siso_preamble#{f}", "chips") for f in fields]
    packet = derive(sum(c.value for c in chips), " + ".join(fields), *chips, S.num("jcas60.sensing_packet#data_blocks"), unit="chips")
    oxygen = S.num(f"itu_p676.gamma_o_dry_sea_level#db_per_km@f_ghz={carrier.value / 1e9:g}", "dB/km")     # at the channel's centre
    if str(S.raw("assumed.ay_sc_pulse#shape")) != "root-raised-cosine":
        raise ValueError("mmw60_trx's SC pulse is no longer a root-raised cosine: jcas60's pulse model does not fit it")
    sim = {k: S.raw(f"jcas60.sim#{k}") for k in ("seed", "frames_averaged", "buffer_chips", "exclude_doppler_bins", "ridge_doppler_bins",
                                                "exclude_range_lags", "cfar_guard", "cfar_train", "compression_span_db", "compression_step_db")}
    units = S.num("jcas60.trn_sensing#units")
    trn = F.trn(int(units.value), tuple(int(b) for b in S.raw("jcas60.trn_sensing#correlated_blocks")))
    trn_packet = derive(packet.value + trn.chips.size, "packet + TRN units x subfields per unit x 6 TRN_BL", packet, units,
                        S.num("edmg.trn_field#trn_bl@trn_seq_length=normal"), unit="chips")
    return Setup(
        fe=shared(),
        range_m=S.num("jcas60.target#range_m", "m"), velocity_mps=S.num("jcas60.target#velocity_mps", "m/s"),
        rcs_dbsm=S.num("target_rcs.deep2019_pedestrian_77g#co_pol_min_dbsm", "dBsm"), oxygen_db_per_km=oxygen,
        carrier_hz=carrier, chip_rate_hz=S.num("edmg.sc_params#chip_rate_hz", "Hz"),
        channel_width_hz=S.num("jcas60.channel#width_hz", "Hz"), rolloff=S.num("assumed.ay_sc_pulse#rolloff"),
        field=F.l_cef(), packet_chips=packet, trn_field=trn, trn_units=units, trn_packet_chips=trn_packet, sifs_s=S.num("dmg.phy_timing#sifs_s", "s"),
        n_packets=S.num("jcas60.sensing_burst#n_packets"), interval_s=S.num("jcas60.sensing_burst#interval_s", "s"),
        frames_per_s=S.num("jcas60.sensing_burst#frames_per_s", "1/s"), doppler_window=str(S.raw("jcas60.doppler_window")),
        adc_osr=S.num("jcas60.adc#samples_per_chip"), desense_allowed_db=S.num("jcas60.adc#desense_allowed_db", "dB"),
        adc_enob_ref=S.num("assumed.converter_enob_ref", "bit"), adc_enob_per_octave=S.num("assumed.converter_bits_per_octave"),
        adc_enob_ref_hz=1e9, adc_backoff_db=S.num("assumed.converter_backoff_db", "dB"),
        rx_rapp_p=S.num("jcas60.rx_rapp_p"), leak_delay_s=S.num("jcas60.leakage_delay_s", "s"),
        pd=S.num("jcas60.detection#pd"),
        bi_r_tx_m=S.num("jcas60.bistatic_geometry#r_tx_m", "m"), bi_r_rx_m=S.num("jcas60.bistatic_geometry#r_rx_m", "m"),
        bi_baseline_m=S.num("jcas60.bistatic_geometry#baseline_m", "m"),
        sim=sim,
    )
