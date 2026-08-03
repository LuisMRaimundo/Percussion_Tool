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
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.special import jn_zeros

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
    @staticmethod
    def decay_time(f: np.ndarray,
                   tau_1k: float = 10.0, alpha: float = 0.84) -> np.ndarray:
        """
        60-dB decay time tau(f) ~ tau_1k * (f/1kHz)^-alpha, fitted to the
        log-log trend of [R] Fig. 9.5 (approx. 400 s near 30 Hz down to
        approx. 3 s at 10 kHz). tau_1k, alpha: literature_derived.
        """
        return tau_1k * (f / 1000.0) ** (-alpha)


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

    @property
    def sigma(self) -> float:
        """Areal density [kg/m^2]."""
        return MYLAR_RHO * self.membrane_thickness

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

    @staticmethod
    def decay_time(f: np.ndarray,
                   tau_100: float = 0.8, alpha: float = 1.2) -> np.ndarray:
        """Short membrane decays: tau ~ tau_100*(f/100)^-alpha.
        internal_default, order-of-magnitude only [FR Ch. 18]."""
        return tau_100 * (f / 100.0) ** (-alpha)


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


def _band_energy_weights(instr, edges: np.ndarray,
                         phases: Dict[str, Tuple[float, float]]
                         ) -> Dict[str, np.ndarray]:
    """
    Relative band energy per phase. Initial excitation approximated as
    equipartition over modes (internal_default); phase evolution follows
    exponential decay with tau(f) from the instrument's decay model, plus
    -- for plates -- a nonlinear low->high energy transfer during the
    build-up phase represented as a 3-5 kHz emphasis factor
    (literature_derived from [R] Fig. 9.6, observations 2-4).
    """
    centres = np.sqrt(edges[:-1] * edges[1:])
    counts = _band_mode_counts(instr, edges)
    tau = instr.decay_time(centres)
    e0 = counts / max(counts.sum(), 1e-12)      # equipartition per mode
    out: Dict[str, np.ndarray] = {}
    for phase, (t0, t1) in phases.items():
        tm = 0.5 * (t0 + t1)
        e = e0 * np.exp(-6.91 * tm / tau)       # 60 dB => ln(1e3)=6.91
        if isinstance(instr, PlateInstrument) and phase in ("buildup",
                                                            "shimmer"):
            # 3-5 kHz "shimmer" emphasis [R, Fig. 9.6 obs. 3-4]
            boost = 1.0 + 2.0 * np.exp(
                -0.5 * ((np.log2(centres / 4000.0)) / 0.6) ** 2)
            e = e * boost
        s = e.sum()
        out[phase] = e / s if s > 0 else e
    return out


@dataclass
class DensityProfile:
    instrument: str
    family: str
    band_edges: np.ndarray
    band_centres: np.ndarray
    modes_per_band: np.ndarray
    energy_weights: Dict[str, np.ndarray]
    notes: List[str] = field(default_factory=list)

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
            }
            for ph, w in self.energy_weights.items():
                row[f"energy_w_{ph}"] = w[i]
            rows.append(row)
        return rows


def generate_profile(instr, f_lo: float = 20.0,
                     f_hi: float = 16000.0) -> DensityProfile:
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
    weights = _band_energy_weights(instr, edges, phases)
    return DensityProfile(instr.name, family, edges, centres, counts,
                          weights, notes)