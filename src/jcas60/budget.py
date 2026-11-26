"""The closed-form path: the radar equation plus processing gain (law 4's first road).

One ay sensing frame is ``n_packets`` packets; in each the receiver
correlates the L-CEF's two 512-chip Golay sequences (1024 chips) against a
matched-filtered copy of what it received, and a windowed FFT across the
packets gives the Doppler axis.  With a root-raised-cosine pulse and its
matched filter the SNR of one chip at the peak is ``P_r / (N0 Rc)`` whatever
the roll-off, so

    SNR(cell) = P_r / (N0 Rc) x 1024 x n_packets / ENBW(window)

and what a sampled map shows is that less the straddling of the lag grid and
of the Doppler bins (``on_grid``: the window-weighted coherent sum of the
raised-cosine pulse over the packets, the target's range walk included).

Also here: what a monostatic receiver needs for the transmit leakage not
to cost sensitivity -- its in-band compression (fmcw_radar's small-signal
theory of a limiter, on the leakage's envelope), the ADC's quantisation
noise, and the LO phase noise the leakage carries after range correlation
-- each against the same allowed desensitisation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .config import Setup
from .physconst import C0, db, undb
from .sources import Num, derive
from .vendor.fmcw_radar.link.compression import small_signal
from .vendor.fmcw_radar.processing.detection import pd_swerling0, snr_for_pd
from .vendor.fmcw_radar.processing.windows import enbw_bins, window


# ----------------------------------------------------------------------------- the pulse
def rc(t, alpha: float):
    """Raised-cosine pulse at ``t`` chips (the root-raised-cosine pulse through its matched filter), 1 at t = 0."""
    t = np.asarray(t, dtype=float)
    den = 1.0 - (2.0 * alpha * t) ** 2
    safe = np.where(np.abs(den) < 1e-9, 1.0, den)
    out = np.sinc(t) * np.cos(np.pi * alpha * t) / safe
    return np.where(np.abs(den) < 1e-9, np.pi / 4.0 * np.sinc(1.0 / (2.0 * alpha)), out)


def rc_spectrum(f, alpha: float):
    """The raised-cosine spectrum at ``f`` in units of the chip rate: 1 in the flat part, 0 beyond (1 + alpha) / 2.
    The matched filter's power response, and the transmitted PSD's shape."""
    f = np.abs(np.asarray(f, dtype=float))
    lo, hi = (1.0 - alpha) / 2.0, (1.0 + alpha) / 2.0
    roll = 0.5 * (1.0 + np.cos(np.pi * (f - lo) / alpha))
    return np.where(f <= lo, 1.0, np.where(f >= hi, 0.0, roll))


def width_at_half_power_chips(alpha: float) -> float:
    """Full width of the pulse's main lobe at -3 dB, in chips."""
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if float(rc(mid, alpha)) ** 2 > 0.5 else (lo, mid)
    return 2.0 * lo


# ----------------------------------------------------------------------------- geometry
@dataclass(frozen=True)
class Geometry:
    kind: str                # monostatic | bistatic
    r_tx_m: float
    r_rx_m: float
    baseline_m: float        # 0 for monostatic
    cos_half_beta: float     # cos of half the bistatic angle (1 for monostatic)

    @property
    def path_m(self) -> float:
        return self.r_tx_m + self.r_rx_m

    @property
    def delay_s(self) -> float:
        """The target's delay on the receiver's lag axis: the round trip, or (bistatic) the excess over the direct path."""
        return (self.path_m - self.baseline_m) / C0

    def path_rate(self, v_mps: float) -> float:
        """d(path)/dt for a target receding at ``v_mps`` (bistatic: along the bisector)."""
        return 2.0 * v_mps * self.cos_half_beta


def geometry(s: Setup, kind: str, range_m: float | None = None) -> Geometry:
    if kind == "monostatic":
        r = s.range_m.value if range_m is None else range_m
        return Geometry(kind, r, r, 0.0, 1.0)
    if kind != "bistatic":
        raise ValueError(kind)
    base = s.bi_baseline_m.value
    if range_m is None:
        rt, rr = s.bi_r_tx_m.value, s.bi_r_rx_m.value
    else:                                   # the target moved out along the bisector, both legs equal
        rt = rr = range_m
    cos_beta = (rt * rt + rr * rr - base * base) / (2.0 * rt * rr)
    return Geometry(kind, rt, rr, base, math.sqrt(0.5 * (1.0 + cos_beta)))


def spreading_loss_db(wavelength_m: float, g: Geometry) -> float:
    """(4 pi)^3 R_t^2 R_r^2 / lambda^2 in dB: the radar equation's geometric term."""
    return 3.0 * db(4.0 * math.pi) + 2.0 * db(g.r_tx_m) + 2.0 * db(g.r_rx_m) - 2.0 * db(wavelength_m)


def one_way_loss_db(wavelength_m: float, r_m: float) -> float:
    return 2.0 * db(4.0 * math.pi * r_m / wavelength_m)


# ----------------------------------------------------------------------------- the link
@dataclass(frozen=True)
class Link:
    kind: str
    eirp_dbm: Num
    spreading_db: Num
    oxygen_db: Num
    p_rx_dbm: Num            # the echo at the receiver input (after the receive antenna)
    noise_chip_dbm: Num      # thermal noise in the chip-rate bandwidth
    snr_chip_db: Num
    processing_gain_db: Num  # 10 log10(correlated chips x packets)
    window_loss_db: Num
    snr_matched_db: Num      # chips x packets, no window (the matched filter in both dimensions)
    snr_cell_db: Num         # with the Doppler window, at the peak
    straddle_db: Num         # what the lag grid and the Doppler bins take off the peak (range walk included)
    snr_grid_db: Num         # what a sampled map is expected to show


def _window(s: Setup) -> np.ndarray:
    return window(s.doppler_window, int(s.n_packets.value))


def on_grid(s: Setup, g: Geometry, lags_per_chip: int) -> tuple[float, int, int]:
    """(amplitude, lag index, Doppler bin) of the strongest cell a sampled map shows for a unit-amplitude
    target: the window-weighted sum over the packets of the pulse at each packet's delay, turning at the
    Doppler.  Bins are ``fftfreq`` order (unshifted), the lag index is in 1 / lags_per_chip chips."""
    m = int(s.n_packets.value)
    t = np.arange(m) * s.interval_s.value
    rate = g.path_rate(s.velocity_mps.value)
    d = (g.delay_s + rate * t / C0) * s.chip_rate_hz.value                      # delay of each packet, chips
    phase = np.exp(-2j * np.pi * rate * t / s.wavelength_m)                   # carrier phase of each packet
    w = _window(s)
    centre = int(round(float(np.mean(d)) * lags_per_chip))
    best = (0.0, 0, 0)
    for lag in range(centre - 2, centre + 3):
        col = np.fft.fft(w * rc(d - lag / lags_per_chip, s.rolloff.value) * phase) / w.sum()
        k = int(np.argmax(np.abs(col)))
        if abs(col[k]) > best[0]:
            best = (float(abs(col[k])), lag, k)
    return best


def link(s: Setup, kind: str, eirp: Num, range_m: float | None = None, field=None) -> Link:
    g = geometry(s, kind, range_m)
    field = s.field if field is None else field
    fe = s.fe
    spread = derive(spreading_loss_db(s.wavelength_m, g), "(4 pi)^3 R_t^2 R_r^2 / lambda^2", s.carrier_hz, s.range_m, unit="dB")
    oxy = derive(s.oxygen_db_per_km.value * g.path_m / 1e3, "gamma_o (R_t + R_r)", s.oxygen_db_per_km, s.range_m, unit="dB")
    p_rx = derive(eirp.value + fe.rx_gain_dbi.value - spread.value - oxy.value + s.rcs_dbsm.value,
                  "EIRP + G_r - spreading - oxygen + RCS", eirp, fe.rx_gain_dbi, spread, oxy, s.rcs_dbsm, unit="dBm")
    n0 = fe.noise_density_dbm_hz
    noise = derive(n0.value + db(s.chip_rate_hz.value), "k T0 F Rc", n0, s.chip_rate_hz, unit="dBm")
    snr_chip = derive(p_rx.value - noise.value, "P_r / (N0 Rc)", p_rx, noise, unit="dB")
    pg = derive(db(field.n_correlated * s.n_packets.value), "10 log10(correlated chips x packets)", s.n_packets, s.chip_rate_hz,
                *((s.trn_units,) if field is s.trn_field else ()), unit="dB")
    wl = derive(db(enbw_bins(_window(s))), "10 log10(ENBW of the Doppler window)", s.n_packets, unit="dB")
    matched = derive(snr_chip.value + pg.value, "SNR(chip) + processing gain", snr_chip, pg, unit="dB")
    cell = derive(matched.value - wl.value, "SNR(chip) + processing gain - window loss", matched, wl, unit="dB")
    amp, _, _ = on_grid(s, g, int(s.adc_osr.value))
    straddle = derive(-2.0 * db(amp), "-20 log10(sampled peak of a unit target)", s.rolloff, s.adc_osr, s.velocity_mps, s.interval_s, unit="dB")
    grid = derive(cell.value - straddle.value, "SNR(cell) - straddle", cell, straddle, unit="dB")
    return Link(kind, eirp, spread, oxy, p_rx, noise, snr_chip, pg, wl, matched, cell, straddle, grid)


# ----------------------------------------------------------------------------- detection
def snr_required_db(pd: float, pfa: float) -> float:
    """The cell SNR a steady target needs for ``pd`` at ``pfa`` with the noise level known (fmcw_radar's
    ``snr_for_pd``).  The table's criterion is the one fmcw_radar's report uses instead -- what its map's CA-CFAR
    needs, 0.49 dB more -- and this figure is printed beside it."""
    return snr_for_pd(pd, pfa)


def pd_at(snr_db: float, pfa: float) -> float:
    return float(pd_swerling0(snr_db, pfa))


def range_for_snr(snr_of_range, snr_req_db: float, lo: float = 0.5, hi: float = 2000.0) -> float:
    """The range at which ``snr_of_range(R)`` falls to ``snr_req_db`` (it falls monotonically)."""
    if snr_of_range(hi) > snr_req_db:
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if snr_of_range(mid) > snr_req_db else (lo, mid)
    return 0.5 * (lo + hi)


def detection_range_m(s: Setup, kind: str, eirp: Num, snr_req_db: Num, field=None) -> Num:
    """The range at which the cell SNR falls to ``snr_req_db`` (the table hands every row the same one)."""
    r = range_for_snr(lambda x: link(s, kind, eirp, x, field).snr_cell_db.value, snr_req_db.value)
    return derive(r, "SNR(cell)(R) = the required SNR", eirp, snr_req_db, s.pd, s.fe.nf_db, s.fe.rx_gain_dbi, s.rcs_dbsm, s.n_packets, unit="m")


# ----------------------------------------------------------------------------- phase noise through the correlator
def overlap_weight(f_hz, s: Setup) -> np.ndarray:
    """W(f): the share of a strong signal's power that, moved by ``f_hz``, still passes the matched filter --
    (1 / Rc) Int RC(v - f) RC(v) dv.  Phase noise at offset f on that signal reaches the correlator with this weight."""
    f = np.atleast_1d(np.asarray(f_hz, dtype=float)) / s.chip_rate_hz.value
    v = np.linspace(-1.0, 1.0, 4001)
    a = s.rolloff.value
    return np.trapezoid(rc_spectrum(v[None, :] - f[:, None], a) * rc_spectrum(v[None, :], a), v, axis=1)


def _pn_grid(s: Setup) -> np.ndarray:
    top = (1.0 + s.rolloff.value) * s.chip_rate_hz.value
    return np.concatenate([np.linspace(0.0, 1e6, 2001)[1:], np.geomspace(1e6, top, 6001)[1:]])


def pn_residual_shared(s: Setup, tau_s: float) -> float:
    """Phase-noise power (rad^2, relative to the carrier) a path of delay ``tau_s`` carries into the correlator when
    transmitter and receiver share the LO: Int 4 sin^2(pi f tau) S_phi(f) W(f) df, both sidebands (range correlation)."""
    f = _pn_grid(s)
    y = 4.0 * np.sin(np.pi * f * tau_s) ** 2 * undb_arr(s.fe.pn_dbchz(f)) * overlap_weight(f, s)
    return 2.0 * float(np.trapezoid(y, f))


def pn_residual_paths(s: Setup, paths) -> float:
    """The same for several static paths that share the LO, ``paths`` = ((amplitude, carrier phase [rad], delay [s]),
    ...): their residues come from one phase process and add in amplitude --
    Int S_phi(f) W(f) |sum_i a_i e^{j theta_i} (e^{-j 2 pi f tau_i} - 1)|^2 df over both sidebands, in the units of
    the amplitudes squared (give the first path amplitude 1 to have it relative to that path)."""
    f = _pn_grid(s)
    w = undb_arr(s.fe.pn_dbchz(f)) * overlap_weight(f, s)
    total = 0.0
    for sign in (1.0, -1.0):
        z = sum(a * np.exp(1j * th) * (np.exp(-2j * np.pi * sign * f * tau) - 1.0) for a, th, tau in paths)
        total += float(np.trapezoid(np.abs(z) ** 2 * w, f))
    return total


def inpacket_loss_db(s: Setup, field, path_rate_mps: float) -> float:
    """What a path whose length changes at ``path_rate_mps`` loses inside one packet: the carrier turns across the
    correlated chips, and their sum is no longer in phase -- |mean over the correlated chips of e^{j 2 pi f_d t}|^2,
    f_d = rate / lambda, summed over the field's correlation windows as they lie in time."""
    fd = path_rate_mps / s.wavelength_m
    acc = sum(np.exp(2j * np.pi * fd * (w0 + np.arange(len(ref))) / s.chip_rate_hz.value).sum() for w0, ref in field.windows)
    return 2.0 * db(abs(acc) / field.n_correlated)


def pn_residual_independent(s: Setup, t_ref_s: float) -> float:
    """The same for two independent LOs (bistatic): Int 2 S_phi(f) (1 - sinc^2(f T)) W(f) df, both sidebands -- the
    phase noise less its mean over ``t_ref_s``.  T is one correlation window, not the whole field: a window's
    correlations cancel on their own (each is a complementary unit), so a phase that is constant over a window puts
    nothing off the peak whatever it is from window to window; only what moves inside a window does.  (The receiver's
    own referencing, the phase it reads on the direct path over the whole field, removes less than that and does not
    set the floor.)"""
    f = _pn_grid(s)
    y = 2.0 * undb_arr(s.fe.pn_dbchz(f)) * (1.0 - np.sinc(f * t_ref_s) ** 2) * overlap_weight(f, s)
    return 2.0 * float(np.trapezoid(y, f))


def undb_arr(x):
    return 10.0 ** (np.asarray(x, dtype=float) / 10.0)


# ----------------------------------------------------------------------------- the receiver's cost
def adc_noise_inband_dbm(p_in_dbm: float, enob: float, s: Setup) -> float:
    """Quantisation noise that reaches the correlator, for an input of ``p_in_dbm`` held ``adc_backoff_db`` under full
    scale on each rail: SQNR = 1.76 + 6.02 ENOB - back-off over the whole Nyquist band, of which the matched filter
    keeps 1 / (samples per chip)."""
    return p_in_dbm - (db(1.5) + 2.0 * db(2.0) * enob - s.adc_backoff_db.value) - db(s.adc_osr.value)


def enob_required(p_in_dbm: float, noise_chip_dbm: float, s: Setup) -> float:
    """The ENOB at which that noise sits ``desense_allowed_db`` of floor rise under the thermal noise."""
    allowed = noise_chip_dbm + db(undb(s.desense_allowed_db.value) - 1.0)
    return (p_in_dbm - allowed - db(1.5) + s.adc_backoff_db.value - db(s.adc_osr.value)) / (2.0 * db(2.0))


def floor_rise_db(extra_dbm: float, noise_chip_dbm: float) -> float:
    return db(1.0 + undb(extra_dbm - noise_chip_dbm))


def compression(env_unit: np.ndarray, strong_dbm: float, s: Setup, p: float | None = None) -> dict:
    """What the receiver's compression does to a weak echo and to the noise riding with a strong static path of
    ``strong_dbm`` whose envelope at unit mean power is ``env_unit``: fmcw_radar's small-signal result for a Rapp
    limiter whose input 1 dB point is the in-band P1dB (dB: ``signal``, ``noise``, ``snr``, ``image``, ``strong``)."""
    a1 = math.sqrt(undb(s.fe.p1db_inband_dbm.value))
    return small_signal(math.sqrt(undb(strong_dbm)) * np.asarray(env_unit), a1, s.rx_rapp_p.value if p is None else p)


def strong_level_for_snr_loss(env_unit: np.ndarray, s: Setup, loss_db: float, p: float | None = None) -> float | None:
    """The level [dBm] of the strong path at which the echo's SNR is down ``loss_db`` (it falls monotonically and
    stops at 3 dB, a hard limiter's figure); None when the limiter never costs that much."""
    lo, hi = s.fe.p1db_inband_dbm.value - 60.0, s.fe.p1db_inband_dbm.value + 60.0
    if compression(env_unit, hi, s, p)["snr"] > -loss_db:
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if compression(env_unit, mid, s, p)["snr"] > -loss_db else (lo, mid)
    return 0.5 * (lo + hi)


def strong_level_for_image(env_unit: np.ndarray, s: Setup, image_dbc: float) -> float | None:
    """The level [dBm] of the strong path at which the limiter's mirror image of the echo (at minus its Doppler)
    reaches ``image_dbc`` under the echo itself.  The small-signal image is an upper bound when the strong path is
    not a single tone, so this level is a lower bound; None when the image never gets that high."""
    rel = lambda lvl: (lambda c: c["image"] - c["signal"])(compression(env_unit, lvl, s))      # noqa: E731
    lo, hi = s.fe.p1db_inband_dbm.value - 60.0, s.fe.p1db_inband_dbm.value + 60.0
    if rel(hi) < image_dbc:
        return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if rel(mid) < image_dbc else (lo, mid)
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class Monostatic:
    conducted_dbm: Num
    noise_chip_dbm: Num
    allowed_extra_dbm: Num           # what one mechanism may add under the thermal floor (the allowed desensitisation)
    need_compression_db: Num         # isolation at which the receiver's compression costs the echo the allowed SNR
    need_compression_at_p1db_db: Num  # ... at which the leakage sits exactly at the in-band P1dB (for reference)
    compression_level_dbm: Num       # the leakage level of that criterion
    need_adc_db: Num                 # ... for the ADC's quantisation noise (the converter model's ENOB)
    need_pn_db: Num                  # ... for the leakage's residual LO phase noise
    pn_residual_dbc: Num
    need_db: Num                     # the largest of the three
    binding: str
    have: tuple                      # ((label, isolation Num, shortfall Num, leak dBm Num, ENOB needed Num), ...)


def monostatic_cost(s: Setup, eirp: Num, env_unit: np.ndarray) -> Monostatic:
    fe = s.fe
    cond = derive(eirp.value - fe.tx_gain_dbi.value, "EIRP - G_t", eirp, fe.tx_gain_dbi, unit="dBm")
    noise = derive(fe.noise_density_dbm_hz.value + db(s.chip_rate_hz.value), "k T0 F Rc", fe.nf_db, s.chip_rate_hz, unit="dBm")
    allowed = derive(noise.value + db(undb(s.desense_allowed_db.value) - 1.0), "N + 10 log10(10^(D/10) - 1)", noise, s.desense_allowed_db, unit="dBm")
    comp0 = derive(cond.value - fe.p1db_inband_dbm.value, "P_tx - P1dB(in band)", cond, fe.p1db_inband_dbm, unit="dB")
    level = strong_level_for_snr_loss(env_unit, s, s.desense_allowed_db.value)
    if level is None:
        raise ValueError("the limiter never costs the allowed desensitisation: the compression criterion has no level")
    comp_level = derive(level, "the leakage level at which the limiter's small-signal SNR is down by the allowed desensitisation",
                        fe.p1db_inband_dbm, s.rx_rapp_p, s.desense_allowed_db, s.rolloff, unit="dBm")
    comp = derive(cond.value - level, "P_tx - that level", cond, comp_level, unit="dB")
    enob = s.adc_enob
    # the leakage level at which the model ADC's in-band noise equals the allowance
    leak_max = allowed.value + (db(1.5) + 2.0 * db(2.0) * enob.value - s.adc_backoff_db.value) + db(s.adc_osr.value)
    adc = derive(cond.value - leak_max, "P_tx - (N_allowed + SQNR + 10 log10(samples per chip))", cond, allowed, enob, s.adc_backoff_db, s.adc_osr, unit="dB")
    res = pn_residual_shared(s, s.leak_delay_s.value)
    res_dbc = derive(db(res), "Int 4 sin^2(pi f tau) S_phi(f) W(f) df", fe.pn_psd0_dbchz, s.leak_delay_s, s.rolloff, unit="dBc")
    pn = derive(cond.value + res_dbc.value - allowed.value, "P_tx + residual PN - N_allowed", cond, res_dbc, allowed, unit="dB")
    needs = {"in-band compression": comp, "ADC quantisation": adc, "LO phase noise on the leakage": pn}
    binding = max(needs, key=lambda k: needs[k].value)
    have = []
    for label, iso in fe.isolation_db:
        leak = derive(cond.value - iso.value, "P_tx - isolation", cond, iso, unit="dBm")
        short = derive(needs[binding].value - iso.value, "isolation needed - isolation available", needs[binding], iso, unit="dB")
        b = derive(enob_required(leak.value, noise.value, s), "ENOB for the quantisation noise to cost the allowed desensitisation",
                   leak, noise, s.adc_backoff_db, s.adc_osr, s.desense_allowed_db, unit="bit")
        have.append((label, iso, short, leak, b))
    return Monostatic(cond, noise, allowed, comp, comp0, comp_level, adc, pn, res_dbc, needs[binding], binding, tuple(have))
