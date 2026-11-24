"""The one set of front-end numbers every row of the table reads (law 2).

Transmit power, the back-off rule, receive noise figure, antenna gains, the
in-band compression point, the TX-RX isolation and the LO phase noise are
configured here and nowhere else: ``shared()`` builds one frozen instance and
hands the same object to every budget and simulation.  Row A is fmcw_radar's
own report; a test holds its inputs equal to this instance, field by field.

Milestone 0's boundary: the same *numbers*, not mmw60_trx's ``Frontend60``
object (which stands on its whole wifitrx chain and is not vendored).  The
one number that cannot be the same is the LO: an FMCW radar's LO is the
chirping synthesiser itself, so row A keeps fmcw_radar's ADPLL and the
phase-noise model here serves the ay rows only (``jcas60.lo_phase_noise``).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

from . import sources as S
from .physconst import KT0_DBM_HZ
from .sources import Num, derive


@dataclass(frozen=True)
class FrontEnd:
    tx_power_dbm: Num          # what the transmitter can put out, conducted
    tx_gain_dbi: Num
    rx_gain_dbi: Num
    nf_db: Num
    rx_gain_db: Num            # the receiver's gain setting ...
    p1db_inband_dbm: Num       # ... and the in-band input compression point at that setting
    isolation_db: tuple        # ((label, Num), ...): what the literature offers for TX-RX isolation
    pn_psd0_dbchz: Num         # the ay rows' LO: PSD0 (1 + (f/fz)^2) / (1 + (f/fp)^2)
    pn_f_pole_hz: Num
    pn_f_zero_hz: Num

    @property
    def eirp_max_dbm(self) -> Num:
        return derive(self.tx_power_dbm.value + self.tx_gain_dbi.value, "tx_power_dbm + tx_gain_dbi",
                      self.tx_power_dbm, self.tx_gain_dbi, unit="dBm")

    def eirp_dbm(self, cap: Num) -> Num:
        """The EIRP under a regulatory cap: full power, or backed off to the cap -- the one rule for every row."""
        return derive(min(self.eirp_max_dbm.value, cap.value), "min(tx_power_dbm + tx_gain_dbi, cap)",
                      self.tx_power_dbm, self.tx_gain_dbi, cap, unit="dBm")

    def backoff_db(self, cap: Num) -> Num:
        return derive(max(0.0, self.eirp_max_dbm.value - cap.value), "max(0, tx_power_dbm + tx_gain_dbi - cap)",
                      self.tx_power_dbm, self.tx_gain_dbi, cap, unit="dB")

    @property
    def noise_density_dbm_hz(self) -> Num:
        return derive(KT0_DBM_HZ + self.nf_db.value, "k T0 F", self.nf_db, unit="dBm/Hz")

    def pn_dbchz(self, f_hz):
        """The LO's phase-noise PSD at offset ``f_hz`` (array or scalar), dBc/Hz: tgad_em.pn_model_shape."""
        import numpy as np
        f = np.abs(np.asarray(f_hz, dtype=float))
        return self.pn_psd0_dbchz.value + 10.0 * np.log10((1.0 + (f / self.pn_f_zero_hz.value) ** 2)
                                                          / (1.0 + (f / self.pn_f_pole_hz.value) ** 2))


@lru_cache(maxsize=1)
def shared() -> FrontEnd:
    gain = S.num("jcas60.rx_gain_db")
    return FrontEnd(
        tx_power_dbm=S.num("rx60.radar_iwr6843#tx_power_dbm", "dBm"),
        tx_gain_dbi=S.num("jcas60.antenna_gain_dbi#tx", "dBi"),
        rx_gain_dbi=S.num("jcas60.antenna_gain_dbi#rx", "dBi"),
        nf_db=S.num("rx60.radar_iwr6843#nf_db", "dB"),
        rx_gain_db=gain,
        p1db_inband_dbm=S.num(f"rx60.radar_iwr6843_fig7_2#p1db_inband_dbm@gain_db={gain.value:g}", "dBm"),
        isolation_db=(("antenna arrays without EBG (low end)", S.num("isolation.tx_rx_tue2018_no_ebg_min_db", "dB")),
                      ("antenna arrays with EBG", S.num("isolation.tx_rx_tue2018#with_ebg_db", "dB"))),
        pn_psd0_dbchz=S.num("tgad_em.pn_params#psd0_dbchz", "dBc/Hz"),
        pn_f_pole_hz=S.num("tgad_em.pn_params#f_pole_hz", "Hz"),
        pn_f_zero_hz=S.num("tgad_em.pn_params#f_zero_hz", "Hz"),
    )
