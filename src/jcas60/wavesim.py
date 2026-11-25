"""The waveform-level path (law 4's second road).

The L-CEF as transmitted -- the vendored Golay sequences, the standard's
layout, a root-raised-cosine pulse -- through a point-target channel, into a
correlation receiver, onto a range-Doppler map whose peak and floor are read
with fmcw_radar's own readers (vendored ``processing/rd.py``).

One packet is one buffer of ``buffer_chips`` chips at ``samples_per_chip``;
the field sits in its middle, the rest is silence (the STF before it and the
headers after it are not sent: for delays up to the cyclic prefix they never
enter a correlation window, and the pulse's tails from there are below
-100 dB).  Per packet:

1. chips -> impulses -> root-raised-cosine pulse (frequency domain, no
   truncation);
2. the LO's phase noise on the transmitter, if on;
3. each path: its delay (a phase ramp in frequency: any fraction of a
   sample), its carrier phase at that packet's time -- which is the Doppler
   across packets -- and the Doppler's turn inside the packet, which is what
   spoils a Golay pair's complementarity;
4. the receiver's LO (the same phase noise for a shared LO, another
   realisation for an independent one), thermal noise, the receiver's
   compression if on (fmcw_radar's memoryless Rapp limiter on the complex
   envelope, vendored), the ADC if on;
5. matched filter, then for each lag the two 512-chip cyclic correlations
   summed and divided by 1024: a path of power P reads P at its lag.

Across packets: the vendored Doppler FFT (window normalised to unit sum).
Bistatic: every packet is first turned back by the phase it reads at the
direct path's lag -- the two devices share no LO.

Amplitudes are in sqrt(mW): ``abs(x)**2`` is a power in mW.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .budget import Geometry, rc_spectrum
from .config import Setup
from .physconst import C0, db, undb
from .vendor.fmcw_radar.link.compression import rapp as _rapp
from .vendor.fmcw_radar.processing import cfar as _cfar
from .vendor.fmcw_radar.processing import rd as _rd


@dataclass(frozen=True)
class Path:
    power_dbm: float          # at the receiver input
    delay_s: float            # at the frame's first packet
    path_rate_mps: float = 0.0   # d(path length)/dt: 2 v for a monostatic target receding at v; 0 for leakage / direct path
    name: str = ""


@dataclass
class Options:
    noise: bool = True
    adc_enob: float | None = None      # None: no quantiser
    rx_p1db_dbm: float | None = None   # None: a linear receiver; else the input 1 dB point of a Rapp limiter ahead of the ADC
    lo: str | None = None              # None: no phase noise | "shared" | "independent"
    reference_lag: int | None = None   # bistatic: the lag (in samples) whose phase every packet is turned back by
    random_packet_phase: bool = False  # bistatic: an unknown carrier phase per packet (no common LO)
    seed: int = 1
    n_packets: int | None = None       # default: the setup's


@dataclass
class Frame:
    cir: np.ndarray            # (packets, lags) complex: the correlator's output, sqrt(mW)
    lags_per_chip: int
    lag0: int                  # index of lag 0 in ``cir``'s second axis
    rd: np.ndarray = field(default=None)   # (Doppler bins centred, lags) complex
    quant: dict = field(default=None)      # with the ADC on: its error over the field, split (see ``Simulator.frame``)

    @property
    def power(self) -> np.ndarray:
        return np.abs(self.rd) ** 2


class Simulator:
    def __init__(self, s: Setup, field=None):
        self.s = s
        self.field = f = s.field if field is None else field
        self.osr = int(s.adc_osr.value)
        self.n_chips = int(s.sim["buffer_chips"])
        while self.n_chips < f.chips.size + f.zero_sidelobe_chips + 2 * 64:        # a longer field gets the next power of two
            self.n_chips *= 2
        self.n = self.n_chips * self.osr
        self.fs = self.osr * s.chip_rate_hz.value
        self.start = (self.n_chips - f.chips.size) // 2          # chip index of the field's first chip in the buffer
        need = f.chips.size + f.zero_sidelobe_chips + 2 * 64
        if self.n_chips < need:
            raise ValueError(f"buffer of {self.n_chips} chips is too short for the field, the longest delay and the pulse tails ({need})")
        freq = np.fft.fftfreq(self.n, d=1.0 / self.fs)
        self.freq = freq
        self.rrc = np.sqrt(rc_spectrum(freq / s.chip_rate_hz.value, s.rolloff.value))    # unit gain in the flat part
        x = np.zeros(self.n, dtype=complex)
        x[self.start * self.osr: (self.start + f.chips.size) * self.osr: self.osr] = f.chips
        # a chip of unit amplitude -> a pulse whose matched-filter output peaks at 1: the impulse train carries
        # ``osr`` in its spectrum per unit chip amplitude, the two root filters multiply to the raised cosine
        self.tx_f = np.fft.fft(x) * self.rrc * math.sqrt(self.osr)
        self.tx_power = float(np.mean(np.abs(np.fft.ifft(self.tx_f)[self.start * self.osr:(self.start + f.chips.size) * self.osr]) ** 2))
        self.t = np.arange(self.n) / self.fs
        # correlation: lags -8 chips .. zero-sidelobe zone + 8 chips, in samples
        self.lag_lo = -8 * self.osr
        self.lags = np.arange(self.lag_lo, (f.zero_sidelobe_chips + 8) * self.osr + 1)
        self._norm = f.n_correlated
        # two ways to the same sums: gathering every window at every lag (the L-CEF: two windows), or one FFT against
        # a template that holds the references at their places (a long field: its gathers would be millions of samples)
        self._by_fft = self.lags.size * f.n_correlated > 1_000_000
        if self._by_fft:
            tmpl = np.zeros(self.n, dtype=complex)
            for w0, ref in f.windows:
                tmpl[(self.start + w0) * self.osr: (self.start + w0 + len(ref)) * self.osr: self.osr] = ref
            self._tmpl_f = np.conj(np.fft.fft(tmpl))
        else:
            self._idx = [((self.start + w0) * self.osr + np.arange(len(ref))[None, :] * self.osr + self.lags[:, None], np.conj(ref))
                         for w0, ref in f.windows]

    # ------------------------------------------------------------------ pieces
    def _phase_noise(self, rng: np.random.Generator) -> np.ndarray:
        """One realisation of the LO phase over the buffer [rad]: Gaussian with the front end's PSD, synthesised on a
        grid 16 buffers long so that offsets down to fs / (16 N) are in it, of which one buffer is cut."""
        m = 16 * self.n
        f = np.fft.rfftfreq(m, d=1.0 / self.fs)
        psd = 10.0 ** (self.s.fe.pn_dbchz(f) / 10.0)                # rad^2/Hz, two-sided: S_phi(f) = L(f)
        psd[0] = 0.0
        # a real phi whose two-sided PSD is S(f): E|rfft bin|^2 = S fs m, so that var(phi) = Int S df over both sides
        amp = np.sqrt(psd * self.fs * m)
        spec = amp * (rng.standard_normal(f.size) + 1j * rng.standard_normal(f.size)) / math.sqrt(2.0)
        phi = np.fft.irfft(spec, m)
        return phi[: self.n]

    def _quantise(self, x: np.ndarray, full_scale: float, enob: float) -> np.ndarray:
        step = 2.0 * full_scale / (2.0 ** enob)

        def q(v):
            return np.clip((np.floor(v / step) + 0.5) * step, -full_scale + step / 2, full_scale - step / 2)
        return q(x.real) + 1j * q(x.imag)

    @property
    def rounding_dbc(self) -> float:
        """Relative to a peak, the level under which a result is the arithmetic's own rounding and not a reading:
        every sample of the buffer the transforms run on off by one unit in the last place, all of them the same
        way.  What is left of an exact cancellation sits far under it and differs between FFT implementations, so
        it is reported as this bound."""
        return 2.0 * db(float(np.finfo(float).eps) * self.n)

    def unit_envelope(self) -> np.ndarray:
        """The transmitted field's envelope at unit mean power, one value per sample of the field: what a strong
        static path (the leakage) looks like to the receiver's compression, up to its level."""
        seg = slice(self.start * self.osr, (self.start + self.field.chips.size) * self.osr)
        return np.abs(np.fft.ifft(self.tx_f)[seg]) / math.sqrt(self.tx_power)

    def correlate(self, y: np.ndarray) -> np.ndarray:
        """Matched filter, then the summed cyclic correlations at every lag."""
        zf = np.fft.fft(y) * self.rrc
        if self._by_fft:                                           # sum_n z[n + lag] conj(t[n]) at every lag, the wanted ones picked out
            return np.fft.ifft(zf * self._tmpl_f)[self.lags % self.n] / self._norm
        z = np.fft.ifft(zf)
        return sum(z[idx] @ ref for idx, ref in self._idx) / self._norm

    # ------------------------------------------------------------------ one frame
    def frame(self, paths: list[Path], opt: Options | None = None, frame_index: int = 0) -> Frame:
        s, opt = self.s, opt or Options()
        m = int(opt.n_packets or s.n_packets.value)
        rng = np.random.default_rng([int(opt.seed), int(frame_index)])
        lam = s.wavelength_m
        sigma = math.sqrt(undb(s.fe.noise_density_dbm_hz.value) * self.fs / 2.0)        # per quadrature, per sample
        scale = 1.0 / math.sqrt(self.tx_power)                                           # the field at unit mean power
        cir = np.empty((m, self.lags.size), dtype=complex)
        full_scale = None
        errors = []
        seg = slice(self.start * self.osr, (self.start + self.field.chips.size) * self.osr)
        for k in range(m):
            t_k = k * s.interval_s.value
            tx_f = self.tx_f
            phi_tx = phi_rx = None
            if opt.lo is not None:
                phi_tx = self._phase_noise(rng)
                phi_rx = phi_tx if opt.lo == "shared" else self._phase_noise(rng)
                tx_f = np.fft.fft(np.fft.ifft(self.tx_f) * np.exp(1j * phi_tx))
            y = np.zeros(self.n, dtype=complex)
            for p in paths:
                tau = p.delay_s + p.path_rate_mps * t_k / C0
                fd = -p.path_rate_mps / lam                                             # receding: negative Doppler
                wave = np.fft.ifft(tx_f * np.exp(-2j * np.pi * self.freq * tau))
                y += (math.sqrt(undb(p.power_dbm)) * scale * np.exp(-2j * np.pi * s.carrier_hz.value * tau)) * wave * np.exp(2j * np.pi * fd * self.t)
            if opt.random_packet_phase:
                y *= np.exp(2j * np.pi * rng.uniform())
            if phi_rx is not None:
                y *= np.exp(-1j * phi_rx)
            if opt.noise:
                y += sigma * (rng.standard_normal(self.n) + 1j * rng.standard_normal(self.n))
            if opt.rx_p1db_dbm is not None:                                            # the receiver's compression, on signal and noise alike
                y = _rapp(y, math.sqrt(undb(opt.rx_p1db_dbm)), s.rx_rapp_p.value)
            if opt.adc_enob is not None:
                if full_scale is None:                                                  # the AGC sets the range once, on the first packet's field
                    full_scale = math.sqrt(float(np.mean(np.abs(y[seg]) ** 2))) * 10.0 ** (s.adc_backoff_db.value / 20.0)
                yq = self._quantise(y, full_scale, opt.adc_enob)
                errors.append(yq[seg] - y[seg])
                y = yq
            cir[k] = self.correlate(y)
        if opt.reference_lag is not None:
            ref = cir[:, opt.reference_lag - self.lag_lo]
            cir = cir * np.exp(-1j * np.angle(ref))[:, None]
        out = Frame(cir=cir, lags_per_chip=self.osr, lag0=-self.lag_lo)
        if errors:
            # the quantiser's error over the field: the part that is the same in every packet (it follows the static
            # path's waveform when the thermal noise is too small to dither it) and the part that is not
            err = np.array(errors)
            total = float(np.mean(np.abs(err) ** 2))
            vary = float(np.mean(np.abs(err - err.mean(axis=0)) ** 2))
            out.quant = {"total_mw": total, "varying_mw": vary, "repeating_share": 1.0 - vary / total,
                         "noise_lsb": sigma / (2.0 * full_scale / 2.0 ** opt.adc_enob)}
        out.rd = _rd.doppler_fft(cir, s.doppler_window, 1)
        return out

    # ------------------------------------------------------------------ reading a map
    def velocity_axis(self, g: Geometry, m: int | None = None) -> np.ndarray:
        """Velocity of each (centred) Doppler bin, receding positive (bistatic: along the bisector)."""
        m = int(m or self.s.n_packets.value)
        fd = np.fft.fftshift(np.fft.fftfreq(m, d=self.s.interval_s.value))
        return -fd * self.s.wavelength_m / (2.0 * g.cos_half_beta)

    def range_axis(self, g: Geometry) -> np.ndarray:
        """Range of each lag: monostatic c tau / 2; bistatic the target's distance from both devices (equal legs) at
        that excess delay."""
        tau = self.lags / self.fs
        return C0 * tau / 2.0 if g.kind == "monostatic" else (C0 * tau + g.baseline_m) / 2.0

    def read(self, fr: Frame, near_lag: int, strong_lags: tuple = ()) -> dict:
        """The peak near ``near_lag`` (samples) and the noise floor away from it.  The floor is read like fmcw_radar's
        time engine reads its own: the mean cell power away from the peak (``rd.noise_floor``).  With static paths at
        ``strong_lags`` (leakage, direct path) the map has a ridge along zero Doppler -- the pulse's skirt, which no
        sequence cancels -- and a spot at their lags; those rows and lags are left out of the floor, and the ridge's
        level at the peak's lag is returned on its own (``ridge_mw``)."""
        s = self.s
        p = fr.power
        m = p.shape[0]
        j0 = near_lag - self.lag_lo
        ex_d, ex_r = int(s.sim["exclude_doppler_bins"]), int(s.sim["exclude_range_lags"]) * self.osr
        rows, cols, search = np.ones(m, dtype=bool), np.ones(p.shape[1], dtype=bool), np.ones(m, dtype=bool)
        if strong_lags:
            off = np.abs(np.arange(m) - m // 2)                            # Doppler bins from zero (the map is centred)
            rows[off <= ex_d] = False                                      # the floor is read well away from the ridge ...
            search[off <= int(s.sim["ridge_doppler_bins"])] = False        # ... the peak only off the window's main lobe around it:
            # a slow target a few bins out is still a peak (the bistatic row's sits at 4.5 bins)
        win = np.where(search[:, None], p[:, j0 - 2 * self.osr: j0 + 2 * self.osr + 1], 0.0)
        i, dj = np.unravel_index(int(np.argmax(win)), win.shape)
        j = j0 - 2 * self.osr + int(dj)
        for lag in strong_lags:
            cols[max(0, lag - self.lag_lo - ex_r): lag - self.lag_lo + ex_r + 1] = False
        r_idx, c_idx = np.flatnonzero(rows), np.flatnonzero(cols)
        sub = p[np.ix_(r_idx, c_idx)]
        floor = _rd.noise_floor(sub, (int(np.argmin(np.abs(r_idx - i))), int(np.searchsorted(c_idx, j))), (ex_d, ex_r))
        # the target's own Doppler row away from its peak (range sidelobes: what a gain that follows the strong path's
        # envelope scatters the echo into) and the mirrored Doppler at the target's lag (a limiter's image)
        far = cols & (np.abs(np.arange(p.shape[1]) - j) > ex_r)
        mirror = (m - i) % m if m % 2 == 0 else m - 1 - i
        return {"peak_mw": float(p[i, j]), "floor_mw": floor, "snr_db": 10.0 * math.log10(p[i, j] / floor),
                "doppler_index": int(i), "lag": int(j + self.lag_lo), "ridge_mw": float(p[m // 2, j]),
                "row_sidelobe_mw": float(p[i, far].max()), "image_mw": float(p[mirror, j])}

    def cfar(self, fr: Frame, pfa: float) -> np.ndarray:
        """fmcw_radar's CA-CFAR on the map's chip-spaced lags (independent cells)."""
        p = fr.power[:, (-self.lag_lo) % self.osr:: self.osr]
        _, det = _cfar.ca_cfar_2d(p, int(self.s.sim["cfar_guard"]), int(self.s.sim["cfar_train"]), pfa)
        return det
