"""
idiophone_density.model
=======================

Theoretical (bibliography-based) spectral-density model for unpitched
percussion: cymbals, gongs / tam-tams (flexural plates) and bass drums
(membranes).

The model takes the instrument's physical parameters (diameter, thickness
or membrane tension, material) and generates an ERB-band spectral density
profile: modal count per ERB band, energy weight per band, and per-phase
composite density indices, schema-aligned with note-level density metadata.

Provenance classes (per value):
  primary_source     -- taken directly from a cited source
  derived            -- computed from primary-source values via stated theory
  literature_derived -- fitted/estimated from figures in cited sources
  internal_default   -- tunable engineering choice, documented

Principal sources
-----------------
[R]  Rossing, T. D. (2000). Science of Percussion Instruments. World
     Scientific. Ch. 9 (cymbals): Chladni-law fits f = c(m+2n)^p
     (Table 9.1), decay times (Fig. 9.5), spectral evolution (Fig. 9.6).
[FR] Fletcher, N. H., & Rossing, T. D. (1998). The Physics of Musical
     Instruments (2nd ed.). Springer. Plates (Ch. 3, 20), membranes and
     bass drum (Ch. 18).
[M]  Meyer, J. (2009). Acoustics and the Performance of Music (5th ed.).
     Springer. Comparative SPL / dynamic ranges / spectral extent.
[C]  Cremer, L., Heckl, M., & Petersson, B. A. T. (2005). Structure-Borne
     Sound (3rd ed.). Springer. Asymptotic modal density of plates.
[GM] Glasberg, B. R., & Moore, B. C. J. (1990). Derivation of auditory
     filter shapes from notched-noise data. Hearing Research, 47, 103-138.
     ERB scale.

Author: prototype prepared for L. Raimundo's doctoral toolchain, 2026.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.special import jn_zeros

# Reference acoustic impedance of air (internal_default, ISO-round).
_RHO_C = 413.0  # Pa·s/m
_P_REF = 20e-6  # Pa

# ----------------------------------------------------------------------
# 1. Materials  (provenance: primary_source, standard handbook values)
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class Material:
    name: str
    E: float        # Young's modulus [Pa]
    rho: float      # density [kg/m^3]
    nu: float       # Poisson ratio [-]

    @property
    def c_L(self) -> float:
        """Longitudinal (quasi-longitudinal) plate wave speed [m/s]."""
        return float(np.sqrt(self.E / (self.rho * (1.0 - self.nu ** 2))))


MATERIALS: Dict[str, Material] = {
    # B20 bell bronze (CuSn20) -- cymbals, tam-tams
    "bronze_B20": Material("bronze_B20", E=110e9, rho=8700.0, nu=0.33),
    # Brass (CuZn) -- some gongs / cheaper cymbals
    "brass":      Material("brass",      E=100e9, rho=8500.0, nu=0.35),
    # Steel -- some gongs
    "steel":      Material("steel",      E=200e9, rho=7850.0, nu=0.30),
}

# Mylar (PET) membrane material, expressed as areal density per thickness.
MYLAR_RHO = 1390.0  # kg/m^3 (primary_source, PET handbook value)


# ----------------------------------------------------------------------
# 2. ERB auditory scale  [GM]  (provenance: primary_source)
# ----------------------------------------------------------------------

def erb_bandwidth(f: np.ndarray) -> np.ndarray:
    """ERB bandwidth in Hz at centre frequency f [GM, eq. 3]."""
    return 24.7 * (4.37 * f / 1000.0 + 1.0)


def erb_number(f: np.ndarray) -> np.ndarray:
    """ERB-rate (Cam) scale [GM]."""
    return 21.4 * np.log10(4.37 * f / 1000.0 + 1.0)


def erb_band_edges(f_lo: float = 20.0, f_hi: float = 16000.0) -> np.ndarray:
    """Edges of consecutive 1-ERB-wide bands covering [f_lo, f_hi]."""
    e_lo, e_hi = erb_number(np.array([f_lo, f_hi]))
    n_bands = int(np.floor(e_hi - e_lo))
    cams = e_lo + np.arange(n_bands + 1)
    return (10.0 ** (cams / 21.4) - 1.0) * 1000.0 / 4.37


# ----------------------------------------------------------------------
# 3. Flexural plate model (cymbal / gong / tam-tam)
# ----------------------------------------------------------------------

@dataclass
class PlateInstrument:
    """
    Circular flexural plate. Valid for cymbals, gongs, tam-tams as a
    first-order model; curvature/dome effects are absorbed empirically
    at low order by the Chladni-law anchor (see below) and flagged as a
    model limitation at high order.

    Parameters
    ----------
    name        : label
    diameter    : plate diameter [m]
    thickness   : effective (mean) thickness [m]
    material    : key into MATERIALS
    chladni     : optional (c, b, p) anchoring low modes to
                  f = c * (m + b*n)**p  [R, Fig. 9.3 / Table 9.1].
                  When absent, low modes are taken from classical
                  free-plate theory scaled by h/d^2 [FR Ch. 3].
    """
    name: str
    diameter: float
    thickness: float
    material: str = "bronze_B20"
    chladni: Optional[Tuple[float, float, float]] = None
    # Instance decay overrides (defaults = literature_derived point fit).
    # Present so Monte Carlo can perturb without changing the public
    # generate_profile call signature; unset behaviour is bit-identical.
    decay_tau_1k: float = 10.0
    decay_alpha: float = 0.84

    # -- asymptotic modal density [C]  (provenance: derived) -----------
    @property
    def area(self) -> float:
        return np.pi * (self.diameter / 2.0) ** 2

    def modal_density(self) -> float:
        """
        Asymptotic modal density n(f) [modes/Hz] of a flexural plate:

            n(f) = sqrt(3) * A / (c_L * h)          [C, sec. 4.3]

        Independent of frequency -- the key structural fact that makes
        the ERB-projected density of a plate GROW with frequency.
        """
        mat = MATERIALS[self.material]
        return float(np.sqrt(3.0) * self.area / (mat.c_L * self.thickness))

    # -- low-order modes ----------------------------------------------
    def low_modes(self, m_max: int = 14, n_max: int = 5) -> np.ndarray:
        """
        Explicit modal frequencies for low (m, n).

        If a Chladni anchor (c, b, p) is provided, use
            f_mn = c * (m + b*n)^p                   [R]
        (provenance: literature_derived, fitted in [R] Table 9.1 /
        Fig. 9.3). Otherwise fall back to the same functional form with
        internal_default parameters scaled from the 46-cm reference
        cymbal by classical plate scaling f ~ h / d^2 [FR]:

            f(d, h) = f_ref * (h / h_ref) * (d_ref / d)^2
        """
        if self.chladni is not None:
            c, b, p = self.chladni
        else:
            # Reference: 46-cm medium crash, n=0 family [R, Fig. 9.3]
            c_ref, b, p = 14.93, 3.0, 1.557          # literature_derived
            d_ref, h_ref = 0.46, 0.0012              # internal_default h_ref
            c = c_ref * (self.thickness / h_ref) * (d_ref / self.diameter) ** 2
        m = np.arange(2, m_max + 1, dtype=float)     # m>=2: free edge
        n = np.arange(0, n_max + 1, dtype=float)
        M, N = np.meshgrid(m, n)
        return np.sort((c * (M + b * N) ** p).ravel())

    # -- decay model [R, Fig. 9.5]  (provenance: literature_derived) ---
    def decay_time(self, f: np.ndarray,
                   tau_1k: Optional[float] = None,
                   alpha: Optional[float] = None) -> np.ndarray:
        """
        60-dB decay time tau(f) ~ tau_1k * (f/1kHz)^-alpha, fitted to the
        log-log trend of [R] Fig. 9.5 (approx. 400 s near 30 Hz down to
        approx. 3 s at 10 kHz). tau_1k, alpha: literature_derived.

        Defaults come from instance fields so Monte Carlo can perturb
        them; explicit kwargs still override (backward compatible).
        """
        t1 = self.decay_tau_1k if tau_1k is None else tau_1k
        al = self.decay_alpha if alpha is None else alpha
        return t1 * (f / 1000.0) ** (-al)


# ----------------------------------------------------------------------
# 4. Membrane model (bass drum)
# ----------------------------------------------------------------------

@dataclass
class MembraneInstrument:
    """
    Circular membrane. First-order model of a bass drum batter head.

    Air loading (which lowers the lowest modes by up to tens of percent)
    and two-head coupling [FR Ch. 18] are NOT modelled; both are flagged
    as validity limits. Either tension or a nominal principal frequency
    f11_nominal may be supplied; if the latter, tension is inferred.
    """
    name: str
    diameter: float                   # [m]
    membrane_thickness: float = 190e-6  # 7.5-mil Mylar (internal_default)
    tension: Optional[float] = None     # [N/m]
    f11_nominal: Optional[float] = None # [Hz] pitch-like (1,1)-mode anchor
    membrane_rho: Optional[float] = None  # None => MYLAR_RHO (MC hook)
    decay_tau_100: float = 0.8
    decay_alpha: float = 1.2

    @property
    def sigma(self) -> float:
        """Areal density [kg/m^2]."""
        rho = MYLAR_RHO if self.membrane_rho is None else self.membrane_rho
        return rho * self.membrane_thickness

    @property
    def radius(self) -> float:
        return self.diameter / 2.0

    def wave_speed(self) -> float:
        if self.tension is not None:
            T = self.tension
        elif self.f11_nominal is not None:
            beta11 = jn_zeros(1, 1)[0]           # 3.8317
            c = 2.0 * np.pi * self.radius * self.f11_nominal / beta11
            return float(c)
        else:
            T = 700.0                            # internal_default [N/m]
        return float(np.sqrt(T / self.sigma))

    def low_modes(self, m_max: int = 10, k_max: int = 8) -> np.ndarray:
        """f_mk = beta_mk * c / (2*pi*a), beta_mk = k-th zero of J_m [FR]."""
        c, a = self.wave_speed(), self.radius
        freqs = []
        for m in range(0, m_max + 1):
            for beta in jn_zeros(m, k_max):
                freqs.append(beta * c / (2.0 * np.pi * a))
        return np.sort(np.asarray(freqs))

    def modal_density(self, f: np.ndarray) -> np.ndarray:
        """
        Asymptotic 2-D membrane modal density (rises linearly with f):

            n(f) = 2*pi*A*f / c^2                  [C]
        """
        c = self.wave_speed()
        return 2.0 * np.pi * (np.pi * self.radius ** 2) * f / c ** 2

    def decay_time(self, f: np.ndarray,
                   tau_100: Optional[float] = None,
                   alpha: Optional[float] = None) -> np.ndarray:
        """Short membrane decays: tau ~ tau_100*(f/100)^-alpha.
        internal_default, order-of-magnitude only [FR Ch. 18]."""
        t0 = self.decay_tau_100 if tau_100 is None else tau_100
        al = self.decay_alpha if alpha is None else alpha
        return t0 * (f / 100.0) ** (-al)


# ----------------------------------------------------------------------
# 5. Spectral-density profile generation
# ----------------------------------------------------------------------

# Canonical temporal phases (cf. [R] sec. 9.3 for plates; membranes use
# strike/decay only). Windows: literature_derived from [R] Figs 9.5-9.6.
PLATE_PHASES = {
    "strike":  (0.000, 0.020),
    "buildup": (0.020, 0.150),
    "shimmer": (0.150, 2.000),
    "residue": (2.000, 6.000),
}
MEMBRANE_PHASES = {
    "strike":  (0.000, 0.050),
    "decay":   (0.050, 1.000),
}


def _band_mode_counts(instr, edges: np.ndarray) -> np.ndarray:
    """Modes per ERB band: explicit low modes below a crossover, asymptotic
    modal density above it (crossover = highest explicitly computed mode)."""
    centres = np.sqrt(edges[:-1] * edges[1:])
    widths = np.diff(edges)
    low = instr.low_modes()
    f_cross = low[-1]
    counts = np.zeros_like(centres)
    # explicit region
    hist, _ = np.histogram(low, bins=edges)
    counts += hist
    # asymptotic region
    if isinstance(instr, PlateInstrument):
        n_of_f = np.full_like(centres, instr.modal_density())
    else:
        n_of_f = instr.modal_density(centres)
    above = centres > f_cross
    counts[above] = n_of_f[above] * widths[above]
    return counts


def _band_energy_weights(
    instr,
    edges: np.ndarray,
    phases: Dict[str, Tuple[float, float]],
    e0: Optional[np.ndarray] = None,
) -> Dict[str, np.ndarray]:
    """
    Relative band energy per phase. Initial excitation defaults to
    equipartition over modes (internal_default). When ``e0`` is supplied
    (AmplitudeLayer measured weights), that replaces equipartition only
    for the initial vector; phase evolution is unchanged.
    """
    centres = np.sqrt(edges[:-1] * edges[1:])
    counts = _band_mode_counts(instr, edges)
    tau = instr.decay_time(centres)
    if e0 is None:
        e0 = counts / max(counts.sum(), 1e-12)  # equipartition per mode
    else:
        e0 = np.asarray(e0, dtype=float)
        s0 = e0.sum()
        e0 = e0 / s0 if s0 > 0 else e0
    out: Dict[str, np.ndarray] = {}
    for phase, (t0, t1) in phases.items():
        tm = 0.5 * (t0 + t1)
        e = e0 * np.exp(-6.91 * tm / tau)  # 60 dB => ln(1e3)=6.91
        if isinstance(instr, PlateInstrument) and phase in ("buildup",
                                                            "shimmer"):
            # 3-5 kHz "shimmer" emphasis [R, Fig. 9.6 obs. 3-4]
            boost = 1.0 + 2.0 * np.exp(
                -0.5 * ((np.log2(centres / 4000.0)) / 0.6) ** 2)
            e = e * boost
        s = e.sum()
        out[phase] = e / s if s > 0 else e
    return out


def power_to_spl_db(power_w: float, distance_m: float,
                    hemisphere: bool = True) -> float:
    """Convert acoustic power [W] at ``distance_m`` to dB SPL re 20 µPa.

    Spreading: hemisphere ``I = P / (2 π r²)`` (internal_default for
    Sivian-style near-source estimates) or full sphere ``4 π r²``.
    """
    if power_w <= 0 or distance_m <= 0:
        return float("nan")
    area = (2.0 if hemisphere else 4.0) * np.pi * distance_m ** 2
    intensity = power_w / area
    p_rms = float(np.sqrt(intensity * _RHO_C))
    return float(20.0 * np.log10(p_rms / _P_REF))


def bars_to_spl_db(bars: float) -> float:
    """Convert Sivian 'bars' (barye) to dB SPL. 1 bar = 0.1 Pa."""
    if bars <= 0:
        return float("nan")
    return float(20.0 * np.log10((bars * 0.1) / _P_REF))


def energy_preserving_remap(
    hist_lo: np.ndarray,
    hist_hi: np.ndarray,
    hist_energy: np.ndarray,
    erb_edges: np.ndarray,
) -> np.ndarray:
    """Map piecewise-constant historical band energies onto an ERB grid.

    Method (documented in README): treat each historical band energy Eᵢ
    as uniform density ρᵢ = Eᵢ / Δfᵢ on [f_lo, f_hi]; each ERB band j
    receives ∫ ρ(f) df over the overlap with band j. Energy is conserved
    up to bands outside the ERB span.
    """
    out = np.zeros(len(erb_edges) - 1, dtype=float)
    for lo, hi, e in zip(hist_lo, hist_hi, hist_energy):
        width = hi - lo
        if width <= 0 or e <= 0:
            continue
        density = e / width
        for j in range(len(out)):
            a, b = erb_edges[j], erb_edges[j + 1]
            ov_lo, ov_hi = max(a, lo), min(b, hi)
            if ov_hi > ov_lo:
                out[j] += density * (ov_hi - ov_lo)
    return out


# Name aliases from prototype instruments -> digitized Sivian keys.
_INSTRUMENT_SOURCE_KEYS = {
    "cymbal_16in_thin": "cymbals_15in",
    "cymbal_18in_medium": "cymbals_15in",
    "cymbal_46cm_medium": "cymbals_15in",
    "bassdrum_32in": "bass_drum_A_36x15",
    "bassdrum_28in": "bass_drum_C_30x12",
    "trumpet": "trumpet",
    "clarinet": "clarinet",
    "flute": "flute",
    "bass_viol": "bass_viol",
    "violin": "violin_soft",
}


@dataclass
class AmplitudeLayer:
    """Absolute-amplitude layer over digitized Sivian/Meyer constants.

    Loads ``data/source_constants.csv``. Where an instrument has peak-band
    power coverage, those energies replace equipartition as the initial
    relative weights (after ERB remap). Otherwise equipartition is kept
    and labelled ``internal_default``.
    """

    csv_path: Path
    table: pd.DataFrame = field(init=False)
    band_edges_hz: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.csv_path = Path(self.csv_path)
        self.table = pd.read_csv(self.csv_path)
        edges = self._load_sivian_edges()
        self.band_edges_hz = edges

    @classmethod
    def default(cls, root: Optional[Path] = None) -> "AmplitudeLayer":
        root = Path(root) if root is not None else Path(__file__).resolve().parent
        return cls(root / "data" / "source_constants.csv")

    def _load_sivian_edges(self) -> np.ndarray:
        sub = self.table[self.table["record_type"] == "sivian_band_edge"]
        pairs = []
        for specimen, g in sub.groupby("specimen", sort=True):
            lo = g.loc[g["parameter"] == "band_edge_lo", "value_si"]
            hi = g.loc[g["parameter"] == "band_edge_hi", "value_si"]
            if len(lo) and len(hi):
                pairs.append((float(lo.iloc[0]), float(hi.iloc[0])))
        pairs.sort(key=lambda x: x[0])
        if not pairs:
            return np.array([20.0, 62.5, 125, 250, 500, 700, 1000, 1400,
                             2000, 2800, 4000, 5600, 8000, 11300])
        edges = [pairs[0][0]]
        for lo, hi in pairs:
            if abs(edges[-1] - lo) > 1e-9:
                edges.append(lo)
            edges.append(hi)
        return np.asarray(edges, dtype=float)

    def source_key(self, instrument_name: str) -> Optional[str]:
        return _INSTRUMENT_SOURCE_KEYS.get(instrument_name)

    def has_coverage(self, instrument_name: str) -> bool:
        key = self.source_key(instrument_name)
        if key is None:
            return False
        sub = self.table[
            (self.table["instrument"] == key)
            & (self.table["parameter"] == "peak_power_band")
            & (self.table["needs_manual_reading"] != 1)
        ]
        return len(sub) > 0

    def historical_band_powers(self, source_key: str
                               ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                                          float, float]:
        """Return (lo, hi, energy_W, whole_W, distance_m) for ``source_key``."""
        meta = self.table[self.table["instrument"] == source_key]
        whole_rows = meta[meta["parameter"] == "peak_power_whole"]
        whole = float(whole_rows["value_si"].iloc[0]) if len(whole_rows) else 0.0
        dist_rows = meta["ref_distance_m"].dropna()
        distance = float(dist_rows.iloc[0]) if len(dist_rows) else 0.9144
        bands = meta[
            (meta["parameter"] == "peak_power_band")
            & (meta["needs_manual_reading"] != 1)
            & meta["value_si"].notna()
        ]
        if len(bands) == 0:
            return (np.array([]), np.array([]), np.array([]), whole, distance)
        lo = bands["band_lo_hz"].to_numpy(float)
        hi = bands["band_hi_hz"].to_numpy(float)
        e = bands["value_si"].to_numpy(float)
        # Residual whole-spectrum power distributed flat over uncovered
        # historical bands (internal_default fill).
        covered = float(e.sum())
        residual = max(whole - covered, 0.0)
        all_lo = self.band_edges_hz[:-1]
        all_hi = self.band_edges_hz[1:]
        density = np.zeros(len(all_lo))
        for a, b, ee in zip(lo, hi, e):
            for i, (x, y) in enumerate(zip(all_lo, all_hi)):
                if abs(x - a) < 1e-9 and abs(y - b) < 1e-9:
                    density[i] += ee
        uncovered = density <= 0
        if residual > 0 and uncovered.any():
            density[uncovered] = residual / uncovered.sum()
        return all_lo, all_hi, density, whole, distance

    def _meyer_hf_factor(self, source_key: str) -> float:
        disc = self.table[
            (self.table["record_type"] == "meyer_hf_discrepancy")
            & (self.table["instrument"] == source_key)
        ]
        if len(disc) == 0:
            return 1.0
        ddb = float(disc["discrepancy_db"].iloc[0])
        return float(10.0 ** (-ddb / 10.0))

    def erb_weights_and_spl(
        self, instrument_name: str, erb_edges: np.ndarray
    ) -> Optional[Tuple[np.ndarray, np.ndarray, float, str]]:
        """Return (relative_w, spl_db_per_band, ref_distance_m, provenance).

        ``None`` if no usable coverage (caller keeps equipartition).
        """
        key = self.source_key(instrument_name)
        if key is None or not self.has_coverage(instrument_name):
            return None
        lo, hi, e, whole, dist = self.historical_band_powers(key)
        if e.size == 0 or e.sum() <= 0:
            return None
        # Prefer Meyer above ~5 kHz: scale historical density in those bands.
        hf = self._meyer_hf_factor(key)
        if hf != 1.0:
            e = e.copy()
            e[lo >= 5000.0] *= hf
        erb_e = energy_preserving_remap(lo, hi, e, erb_edges)
        s = erb_e.sum()
        w = erb_e / s if s > 0 else erb_e
        # Absolute SPL from remapped band powers at reference distance.
        spl = np.array(
            [power_to_spl_db(float(p), dist) if p > 0 else float("nan")
             for p in erb_e],
            dtype=float,
        )
        return w, spl, dist, "primary_source"


@dataclass
class DensityProfile:
    instrument: str
    family: str
    band_edges: np.ndarray
    band_centres: np.ndarray
    modes_per_band: np.ndarray
    energy_weights: Dict[str, np.ndarray]
    notes: List[str] = field(default_factory=list)
    absolute_spl_db: Dict[str, np.ndarray] = field(default_factory=dict)
    energy_provenance: str = "internal_default"
    ref_distance_m: Optional[float] = None

    def composite_index(self, phase: str) -> float:
        """
        Energy-weighted spectral component density (scalar summary):
        effective number of ERB bands times mean modal occupation,
        computed as   D = exp(H) * <modes/band>_w   with H the Shannon
        entropy of the band energy distribution (internal_default
        aggregation; report per-band data as primary output).
        """
        w = self.energy_weights[phase]
        nz = w > 0
        H = -np.sum(w[nz] * np.log(w[nz]))
        mean_occ = float(np.sum(w * self.modes_per_band))
        return float(np.exp(H) * mean_occ) ** 0.5  # geometric compromise

    def to_rows(self) -> List[dict]:
        rows = []
        for i, fc in enumerate(self.band_centres):
            row = {
                "instrument": self.instrument,
                "family": self.family,
                "band_index": i,
                "f_lo_hz": self.band_edges[i],
                "f_hi_hz": self.band_edges[i + 1],
                "f_centre_hz": fc,
                "modes_per_band": self.modes_per_band[i],
                "energy_provenance": self.energy_provenance,
                "ref_distance_m": self.ref_distance_m
                if self.ref_distance_m is not None else "",
            }
            for ph, w in self.energy_weights.items():
                row[f"energy_w_{ph}"] = w[i]
            for ph, spl in self.absolute_spl_db.items():
                row[f"spl_db_{ph}"] = spl[i]
            rows.append(row)
        return rows


def generate_profile(
    instr,
    f_lo: float = 20.0,
    f_hi: float = 16000.0,
    amplitude_layer: Optional[AmplitudeLayer] = None,
) -> DensityProfile:
    """Generate an ERB density profile; optionally apply AmplitudeLayer."""
    edges = erb_band_edges(f_lo, f_hi)
    centres = np.sqrt(edges[:-1] * edges[1:])
    counts = _band_mode_counts(instr, edges)
    if isinstance(instr, PlateInstrument):
        family, phases = "plate", PLATE_PHASES
        notes = [
            "flat-plate approximation: dome/curvature not modelled",
            "linear modal regime: chaotic broadband regime at ff not modelled",
        ]
    else:
        family, phases = "membrane", MEMBRANE_PHASES
        notes = [
            "single membrane in vacuo: air loading and 2-head coupling "
            "not modelled (lowest modes overestimated)",
        ]

    e0 = None
    provenance = "internal_default"
    ref_dist: Optional[float] = None
    abs_spl: Dict[str, np.ndarray] = {}

    if amplitude_layer is not None:
        mapped = amplitude_layer.erb_weights_and_spl(instr.name, edges)
        if mapped is not None:
            e0, spl0, ref_dist, provenance = mapped
            notes.append(
                f"initial energy from AmplitudeLayer ({provenance}); "
                f"ref distance {ref_dist:.4f} m"
            )
            # Phase-evolve absolute SPL with the same relative envelopes.
            rel = _band_energy_weights(instr, edges, phases, e0=e0)
            for ph, w in rel.items():
                # Scale absolute vector so strike matches spl0 energy sum.
                # Use relative shape of phase weight with strike absolute
                # peak as anchor (internal_default absolute evolution).
                scale = np.nansum((10 ** (spl0 / 10.0))[np.isfinite(spl0)])
                lin = w * (scale if scale > 0 else 1.0)
                with np.errstate(divide="ignore"):
                    abs_spl[ph] = 10.0 * np.log10(np.maximum(lin, 1e-30))
            weights = rel
            return DensityProfile(
                instr.name, family, edges, centres, counts, weights, notes,
                absolute_spl_db=abs_spl, energy_provenance=provenance,
                ref_distance_m=ref_dist,
            )

    weights = _band_energy_weights(instr, edges, phases, e0=e0)
    return DensityProfile(
        instr.name, family, edges, centres, counts, weights, notes,
        absolute_spl_db=abs_spl, energy_provenance=provenance,
        ref_distance_m=ref_dist,
    )