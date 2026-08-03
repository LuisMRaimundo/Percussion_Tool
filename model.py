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
from typing import Dict, List, Optional, Sequence, Tuple

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
    # "cymbal" → PLATE_PHASES; "tamtam" → PLATE_PHASES_TAMTAM (wind gongs
    # use the tam-tam template with a note).
    plate_class: str = "cymbal"

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

_MODE_FREQ_PARAM_RE = None  # compiled lazily


def _mode_freq_param_re():
    import re
    global _MODE_FREQ_PARAM_RE
    if _MODE_FREQ_PARAM_RE is None:
        # mode_freq_01, mode_freq_11b, …
        _MODE_FREQ_PARAM_RE = re.compile(
            r"^mode_freq_(\d)(\d)([a-z]?)$", re.I
        )
    return _MODE_FREQ_PARAM_RE


def membrane_beta(m: int, n: int) -> float:
    """Bessel zero β_mn: n-th zero of J_m (n ≥ 1)."""
    if n < 1:
        raise ValueError("radial index n must be >= 1")
    return float(jn_zeros(m, n)[n - 1])


def fit_effective_wave_speed(
    freqs_hz: np.ndarray,
    mode_indices: Sequence[Tuple[int, int]],
    radius_m: float,
) -> float:
    """Least-squares fit of f = β c / (2π a) → effective c [m/s] (derived).

    Absorbs air-loading and two-head coupling into an effective wave speed
    for the anchored low-mode set, analogous to the Chladni anchor for
    plates (low-order only).
    """
    f = np.asarray(freqs_hz, dtype=float)
    betas = np.array(
        [membrane_beta(m, n) for m, n in mode_indices], dtype=float
    )
    if f.size == 0 or f.size != betas.size:
        raise ValueError("freqs and mode_indices length mismatch")
    # f = β * c / (2πa)  ⇒  c = 2πa * (Σ β f) / (Σ β²)
    c = (
        2.0
        * np.pi
        * radius_m
        * float(np.dot(betas, f) / np.dot(betas, betas))
    )
    return float(c)


def measured_modes_from_csv(
    instrument_key: str,
    specimen: str = "both_heads",
    csv_path: Optional[Path] = None,
) -> Tuple[np.ndarray, Tuple[Tuple[int, int], ...]]:
    """Load Fletcher & Rossing Ch. 18 measured bass-drum modes from CSV.

    Skips ``needs_manual_reading=1`` and blank ``value_si``. Returns
    frequencies sorted ascending with aligned ``(m, n)`` labels. Doublet
    suffixes (e.g. ``mode_freq_01b``) are included when present.
    """
    path = (
        Path(csv_path)
        if csv_path is not None
        else Path(__file__).resolve().parent / "data" / "source_constants.csv"
    )
    df = pd.read_csv(path)
    sub = df[
        (df["record_type"] == "fr_ch18_bassdrum_mode")
        & (df["instrument"] == instrument_key)
        & (df["specimen"] == specimen)
        & (df["needs_manual_reading"] != 1)
        & df["value_si"].notna()
    ]
    if len(sub) == 0:
        return np.array([], dtype=float), tuple()

    rx = _mode_freq_param_re()
    pairs: List[Tuple[float, int, int]] = []
    for _, row in sub.iterrows():
        mobj = rx.match(str(row["parameter"]))
        if not mobj:
            continue
        m_i, n_i = int(mobj.group(1)), int(mobj.group(2))
        pairs.append((float(row["value_si"]), m_i, n_i))
    pairs.sort(key=lambda t: t[0])
    freqs = np.array([p[0] for p in pairs], dtype=float)
    labels = tuple((p[1], p[2]) for p in pairs)
    return freqs, labels


@dataclass
class MembraneInstrument:
    """
    Circular membrane. First-order model of a bass drum batter head.

    Without ``measured_modes``, air loading and two-head coupling
    [FR Ch. 18] are NOT modelled (lowest modes overestimated) — flagged
    as a validity limit. Optional Fletcher & Rossing measured modes
    override the lowest theoretical frequencies (``primary_source``
    anchor); theory continues above with a least-squares fitted effective
    wave speed ``c`` (``derived``), absorbing air-loading / two-head
    effects at low order exactly as the Chladni anchor absorbs curvature
    for plates.

    Either tension, ``f11_nominal``, or ``effective_c`` may set the wave
    speed when no measured-mode fit is active.
    """
    name: str
    diameter: float                   # [m]
    membrane_thickness: float = 190e-6  # 7.5-mil Mylar (internal_default)
    tension: Optional[float] = None     # [N/m]
    f11_nominal: Optional[float] = None # [Hz] pitch-like (1,1)-mode anchor
    membrane_rho: Optional[float] = None  # None => MYLAR_RHO (MC hook)
    decay_tau_100: float = 0.8
    decay_alpha: float = 1.2
    measured_modes: Optional[np.ndarray] = None  # Hz, sorted
    measured_mode_indices: Optional[Tuple[Tuple[int, int], ...]] = None
    effective_c: Optional[float] = None  # [m/s], derived when fitted

    @property
    def sigma(self) -> float:
        """Areal density [kg/m^2]."""
        rho = MYLAR_RHO if self.membrane_rho is None else self.membrane_rho
        return rho * self.membrane_thickness

    @property
    def radius(self) -> float:
        return self.diameter / 2.0

    def has_measured_anchor(self) -> bool:
        return (
            self.measured_modes is not None
            and len(self.measured_modes) > 0
            and self.measured_mode_indices is not None
            and len(self.measured_mode_indices) == len(self.measured_modes)
        )

    def fitted_effective_c(self) -> Optional[float]:
        """LS-squares effective c from measured modes, or stored value."""
        if self.effective_c is not None:
            return float(self.effective_c)
        if not self.has_measured_anchor():
            return None
        return fit_effective_wave_speed(
            self.measured_modes,  # type: ignore[arg-type]
            self.measured_mode_indices,  # type: ignore[arg-type]
            self.radius,
        )

    def wave_speed(self) -> float:
        c_fit = self.fitted_effective_c()
        if c_fit is not None:
            return float(c_fit)
        if self.tension is not None:
            T = self.tension
        elif self.f11_nominal is not None:
            beta11 = jn_zeros(1, 1)[0]           # 3.8317
            c = 2.0 * np.pi * self.radius * self.f11_nominal / beta11
            return float(c)
        else:
            T = 700.0                            # internal_default [N/m]
        return float(np.sqrt(T / self.sigma))

    def theoretical_modes_in_vacuo(
        self, m_max: int = 10, k_max: int = 8
    ) -> np.ndarray:
        """In-vacuo Bessel modes using tension / f11_nominal / default T.

        Ignores ``measured_modes`` and ``effective_c`` so VAL3 can compare
        the unanchored theory against the Fletcher & Rossing measurements.
        """
        if self.tension is not None:
            c = float(np.sqrt(self.tension / self.sigma))
        elif self.f11_nominal is not None:
            beta11 = jn_zeros(1, 1)[0]
            c = 2.0 * np.pi * self.radius * self.f11_nominal / beta11
        else:
            c = float(np.sqrt(700.0 / self.sigma))
        a = self.radius
        freqs = []
        for m in range(0, m_max + 1):
            for beta in jn_zeros(m, k_max):
                freqs.append(beta * c / (2.0 * np.pi * a))
        return np.sort(np.asarray(freqs, dtype=float))

    def low_modes(self, m_max: int = 10, k_max: int = 8) -> np.ndarray:
        """Modal frequencies [Hz].

        With measured anchor: return primary_source measured modes, then
        append theoretical Bessel modes **above** the highest measured
        frequency using the fitted effective wave speed (derived). This
        absorbs air-loading / two-head effects into the low modes only.
        Without anchor: classical ``f_mk = β_mk c / (2π a)`` [FR].
        """
        if self.has_measured_anchor():
            measured = np.sort(np.asarray(self.measured_modes, dtype=float))
            f_max = float(measured[-1])
            c, a = self.wave_speed(), self.radius
            theo: List[float] = []
            for m in range(0, m_max + 1):
                for beta in jn_zeros(m, k_max):
                    f = float(beta * c / (2.0 * np.pi * a))
                    if f > f_max:
                        theo.append(f)
            if theo:
                return np.sort(np.concatenate([measured, np.asarray(theo)]))
            return measured

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


def make_bassdrum_catalogue(
    csv_path: Optional[Path] = None,
) -> List["MembraneInstrument"]:
    """Build bass-drum catalogue: 82 cm anchored + scaled siblings.

    The Fletcher & Rossing Table 18.5 drum is 82 cm; that entry carries
    the measured-mode anchor. Catalogue 32-in / 28-in drums use the
    **fitted effective c** from the 82 cm set with their own diameter
    (``f ∝ 1/a`` at fixed effective c — ``internal_default`` scaling;
    frequencies are not copied across sizes).
    """
    path = csv_path
    freqs, labels = measured_modes_from_csv(
        "bassdrum_82cm", specimen="both_heads", csv_path=path
    )
    d82 = 0.82
    if len(freqs) == 0:
        # Graceful fallback: legacy f11_nominal instruments only.
        return [
            MembraneInstrument("bassdrum_32in", 0.813, f11_nominal=60.0),
            MembraneInstrument("bassdrum_28in", 0.711, f11_nominal=72.0),
        ]
    c_eff = fit_effective_wave_speed(freqs, labels, d82 / 2.0)
    anchored = MembraneInstrument(
        name="bassdrum_82cm",
        diameter=d82,
        membrane_thickness=0.25e-3,  # 0.010 in Mylar [FR Ch.18 prose]
        measured_modes=freqs,
        measured_mode_indices=labels,
        effective_c=c_eff,
    )
    # Sibling sizes: same effective c, own diameter (internal_default).
    siblings = [
        MembraneInstrument(
            "bassdrum_32in",
            0.813,
            effective_c=c_eff,
        ),
        MembraneInstrument(
            "bassdrum_28in",
            0.711,
            effective_c=c_eff,
        ),
    ]
    return [anchored] + siblings


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
# Tam-tam / wind-gong template (literature_derived from [R] tam-tam
# chapter / [FR] Ch. 20): slow nonlinear HF bloom over ~1–2 s, then
# long shimmer/decay out to tens of seconds — NOT the cymbal 0.15 s
# shimmer onset.
PLATE_PHASES_TAMTAM = {
    # Immediate impact transient before the bloom develops.
    "strike":  (0.000, 0.050),
    # Low→high cascade develops over ~1–2 s after the stroke ([R]/[FR]).
    "bloom":   (0.050, 1.500),
    # Established broadband shimmer once HF energy has developed.
    "shimmer": (1.500, 6.000),
    # Long residue of large tam-tams (tens of seconds truncated at 20 s).
    "residue": (6.000, 20.000),
}
MEMBRANE_PHASES = {
    "strike":  (0.000, 0.050),
    "decay":   (0.050, 1.000),
}

# Hertzian-impact contact-time placeholders (s) when CSV has no source row.
# internal_default — NEVER overwrite a source-read CSV value with these.
_CONTACT_TIME_PLACEHOLDERS_S: Dict[str, float] = {
    "stick_tip": 0.15e-3,      # v0.3.5; f_c ≈ 3.33 kHz
    "stick_shoulder": 0.3e-3,  # v0.3.5
    "stick_bell": 0.10e-3,     # v0.3.5; near-impulsive dome strike
    "yarn_mallet": 3.0e-3,     # unchanged v0.3.4 (mallet cohort validated)
    "bass_drum_beater": 6.0e-3,  # unchanged v0.3.4
}
# Dynamic scaling of t_contact (internal_default). Hertz theory predicts
# t ∝ v^(-1/5) (weak force dependence); these factors are deliberately
# conservative / monotone: pp longer, ff shorter.
_DYNAMIC_T_CONTACT_SCALE: Dict[str, float] = {
    "pp": 1.6,
    "mf": 1.0,
    "ff": 0.6,
}
# Shimmer-boost amplitude gate (internal_default). Boost encodes the
# nonlinear low→high cascade of [R] Fig. 9.6 / §9.4 (amplitude-driven).
_DYNAMIC_SHIMMER_GATE: Dict[str, float] = {
    "pp": 0.0,
    "mf": 0.5,
    "ff": 1.0,
}

# Parsed validation stroke → implement label used for contact-time lookup.
_STROKE_TO_IMPLEMENT: Dict[str, str] = {
    "mallet": "yarn_mallet",
    "stick.normal": "stick_tip",
    "stick.shoulder": "stick_shoulder",
    "stick.bell": "stick_bell",
    "stick_tip": "stick_tip",
    "stick_shoulder": "stick_shoulder",
    "stick_bell": "stick_bell",
    "yarn_mallet": "yarn_mallet",
    "bass_drum_beater": "bass_drum_beater",
}


def map_stroke_to_implement(
    stroke: Optional[str],
    *,
    plate_class: Optional[str] = None,
    family: Optional[str] = None,
) -> str:
    """Map a parsed stroke label to an excitation implement key.

    ``unmarked`` → ``yarn_mallet`` for tam-tams / wind gongs (orchestral
    default) and ``stick_tip`` otherwise (``internal_default`` mapping).
    """
    if stroke is None:
        if family == "membrane":
            return "bass_drum_beater"
        if plate_class in {"tamtam", "windgong"}:
            return "yarn_mallet"
        return "stick_tip"
    key = stroke.strip().lower()
    if key in _STROKE_TO_IMPLEMENT:
        return _STROKE_TO_IMPLEMENT[key]
    if key == "unmarked":
        if plate_class in {"tamtam", "windgong"}:
            return "yarn_mallet"
        if family == "membrane":
            return "bass_drum_beater"
        return "stick_tip"
    return "stick_tip"


def load_contact_time_s(
    implement: str,
    csv_path: Optional[Path] = None,
) -> Tuple[float, str, bool]:
    """Return (t_contact_s, provenance, used_placeholder).

    Source-read CSV values (``primary_source`` / ``literature_derived``,
    ``needs_manual_reading != 1``) win over placeholders. Placeholders are
    ``internal_default`` and flagged via the boolean.
    """
    path = (
        Path(csv_path)
        if csv_path is not None
        else Path(__file__).resolve().parent / "data" / "source_constants.csv"
    )
    if path.is_file():
        df = pd.read_csv(path)
        sub = df[
            (df["record_type"] == "excitation_contact_time")
            & (df["instrument"] == implement)
            & (df["needs_manual_reading"] != 1)
            & df["value_si"].notna()
        ]
        # Prefer non-placeholder provenances when both exist.
        if len(sub) > 0:
            ranked = sub.copy()
            ranked["_rank"] = ranked["provenance"].map(
                lambda p: 0
                if p in ("primary_source", "literature_derived")
                else 1
            )
            ranked = ranked.sort_values("_rank")
            row = ranked.iloc[0]
            prov = str(row["provenance"])
            t = float(row["value_si"])
            used_ph = prov == "internal_default" or "placeholder" in str(
                row.get("location", "")
            ).lower()
            return t, prov, used_ph
    if implement not in _CONTACT_TIME_PLACEHOLDERS_S:
        raise KeyError(f"unknown excitation implement: {implement}")
    return (
        float(_CONTACT_TIME_PLACEHOLDERS_S[implement]),
        "internal_default",
        True,
    )


def excitation_cutoff_hz(t_contact_s: float) -> float:
    """Hertzian-impact low-pass cutoff: f_c ≈ 1 / (2 t_contact)."""
    t = max(float(t_contact_s), 1e-9)
    return 1.0 / (2.0 * t)


def excitation_filter_gain(
    f_hz: np.ndarray, f_c_hz: float
) -> np.ndarray:
    """Magnitude-squared low-pass: 1 / (1 + (f/f_c)^4) (internal_default)."""
    f = np.asarray(f_hz, dtype=float)
    fc = max(float(f_c_hz), 1e-9)
    return 1.0 / (1.0 + (f / fc) ** 4)


def apply_excitation_filter(
    e0: np.ndarray,
    centres_hz: np.ndarray,
    t_contact_s: float,
) -> np.ndarray:
    """Filter then renormalize an initial energy vector."""
    e = np.asarray(e0, dtype=float) * excitation_filter_gain(
        centres_hz, excitation_cutoff_hz(t_contact_s)
    )
    s = float(e.sum())
    return e / s if s > 0 else e


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
    shimmer_gate: float = 1.0,
) -> Dict[str, np.ndarray]:
    """
    Relative band energy per phase. Initial excitation defaults to
    equipartition over modes (internal_default). When ``e0`` is supplied
    (AmplitudeLayer measured weights and/or contact-time filter), that
    replaces equipartition only for the initial vector; phase evolution
    is unchanged aside from dynamic-gated shimmer boost.

    ``shimmer_gate`` scales the 3–5 kHz boost amplitude (1.0 = v0.3.3 full
    boost; 0 = off). For tam-tams the boost applies in ``shimmer`` only,
    not ``bloom``; for cymbals it applies in ``buildup`` and ``shimmer``.
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
    is_tamtam = (
        isinstance(instr, PlateInstrument)
        and getattr(instr, "plate_class", "cymbal") == "tamtam"
    )
    boost_phases = ("shimmer",) if is_tamtam else ("buildup", "shimmer")
    out: Dict[str, np.ndarray] = {}
    for phase, (t0, t1) in phases.items():
        tm = 0.5 * (t0 + t1)
        e = e0 * np.exp(-6.91 * tm / tau)  # 60 dB => ln(1e3)=6.91
        if (
            isinstance(instr, PlateInstrument)
            and phase in boost_phases
            and shimmer_gate > 0.0
        ):
            # 3-5 kHz "shimmer" emphasis [R, Fig. 9.6 obs. 3-4]
            boost = 1.0 + (2.0 * shimmer_gate) * np.exp(
                -0.5 * ((np.log2(centres / 4000.0)) / 0.6) ** 2
            )
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
# Cymbal model names map onto Sivian's 15-in. clash PAIR (different
# mechanical system and stroke) — approximation labelled internal_default.
_INSTRUMENT_SOURCE_KEYS = {
    "cymbal_16in_thin": "cymbals_15in",
    "cymbal_18in_medium": "cymbals_15in",
    "cymbal_46cm_medium": "cymbals_15in",
    "bassdrum_32in": "bass_drum_A_36x15",
    # 82 cm FR modal anchor → Sivian 36 in (≈91 cm) bass drum A: size
    # mismatch approximation (internal_default); see data/README.md.
    "bassdrum_82cm": "bass_drum_A_36x15",
    "bassdrum_28in": "bass_drum_C_30x12",
    "trumpet": "trumpet",
    "clarinet": "clarinet",
    "flute": "flute",
    "bass_viol": "bass_viol",
    "violin": "violin_soft",
}

# Fill-fraction thresholds for AmplitudeLayer provenance (internal_default).
# fill_fraction = residual_power / whole_spectrum_peak_power after placing
# digitized peak_power_band rows.
FILL_FRAC_PRIMARY_MAX = 0.10   # <=10% fill → label primary_source
FILL_FRAC_MIXED_MAX = 0.60     # <=60% fill → mixed_primary_and_fill; above → refuse
# ERB band counts as "measured" for calibration if >50% of its remapped
# energy came from measured historical bands (internal_default).
ERB_MEASURED_ENERGY_FRAC = 0.50


@dataclass
class AmplitudeLayer:
    """Absolute-amplitude layer over digitized Sivian/Meyer constants.

    Loads ``data/source_constants.csv``. Where an instrument has peak-band
    power coverage **and** residual fill_fraction ≤ ``FILL_FRAC_MIXED_MAX``,
    those energies replace equipartition as the initial relative weights
    (after ERB remap). Mostly-filled vectors are refused; equipartition
    stays labelled ``internal_default``.
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

    def fill_fraction_for(self, instrument_name: str) -> Optional[float]:
        """Return residual fill fraction for ``instrument_name``, or None."""
        key = self.source_key(instrument_name)
        if key is None:
            return None
        _lo, _hi, e, _meas, fill_frac, _whole, _dist = self.historical_band_powers(
            key
        )
        if e.size == 0:
            return None
        return float(fill_frac)

    def n_measured_bands(self, instrument_name: str) -> int:
        key = self.source_key(instrument_name)
        if key is None:
            return 0
        _lo, _hi, _e, is_measured, _ff, _w, _d = self.historical_band_powers(key)
        if is_measured.size == 0:
            return 0
        return int(np.sum(is_measured))

    def has_coverage(self, instrument_name: str) -> bool:
        """True only when digitized bands exist and fill is not mostly residual.

        Uses the same thresholds as ``erb_weights_and_spl`` so GUI/CLI listings
        agree with runtime behaviour (internal_default thresholds).
        """
        key = self.source_key(instrument_name)
        if key is None:
            return False
        lo, hi, e, is_measured, fill_frac, _whole, _dist = (
            self.historical_band_powers(key)
        )
        if e.size == 0 or not bool(np.any(is_measured)):
            return False
        return float(fill_frac) <= FILL_FRAC_MIXED_MAX

    def historical_band_powers(
        self, source_key: str
    ) -> Tuple[
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float, float
    ]:
        """Return historical band powers for ``source_key``.

        Returns
        -------
        lo, hi, energy_W, is_measured, fill_fraction, whole_W, distance_m

        ``is_measured`` is True only for bands that received a digitized
        ``peak_power_band`` row (not residual fill).

        Residual power (``whole − Σ measured``) is distributed over uncovered
        bands **proportionally to bandwidth** — uniform spectral density in
        W/Hz across the uncovered region (`internal_default` fill,
        uniform-density variant). ``fill_fraction = residual / whole``
        (0 when no fill or whole ≤ 0).
        """
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
            empty = np.array([])
            return empty, empty, empty, empty.astype(bool), 0.0, whole, distance
        lo_m = bands["band_lo_hz"].to_numpy(float)
        hi_m = bands["band_hi_hz"].to_numpy(float)
        e_m = bands["value_si"].to_numpy(float)
        covered = float(e_m.sum())
        residual = max(whole - covered, 0.0)
        fill_fraction = (residual / whole) if whole > 0 else 0.0

        all_lo = self.band_edges_hz[:-1]
        all_hi = self.band_edges_hz[1:]
        energy = np.zeros(len(all_lo), dtype=float)
        is_measured = np.zeros(len(all_lo), dtype=bool)
        for a, b, ee in zip(lo_m, hi_m, e_m):
            for i, (x, y) in enumerate(zip(all_lo, all_hi)):
                if abs(x - a) < 1e-9 and abs(y - b) < 1e-9:
                    energy[i] += ee
                    is_measured[i] = True
        uncovered = ~is_measured
        if residual > 0 and uncovered.any():
            widths = all_hi - all_lo
            uw = widths[uncovered]
            total_w = float(uw.sum())
            if total_w > 0:
                # Uniform spectral density fill: E_i ∝ Δf_i (internal_default).
                energy[uncovered] = residual * (uw / total_w)
        return all_lo, all_hi, energy, is_measured, fill_fraction, whole, distance

    def _meyer_hf_factor(self, source_key: str) -> float:
        disc = self.table[
            (self.table["record_type"] == "meyer_hf_discrepancy")
            & (self.table["instrument"] == source_key)
        ]
        if len(disc) == 0:
            return 1.0
        ddb = float(disc["discrepancy_db"].iloc[0])
        return float(10.0 ** (-ddb / 10.0))

    def provenance_for_fill(self, fill_fraction: float) -> Optional[str]:
        """Map fill fraction → provenance label, or None to refuse coverage."""
        if fill_fraction <= FILL_FRAC_PRIMARY_MAX:
            return "primary_source"
        if fill_fraction <= FILL_FRAC_MIXED_MAX:
            return "mixed_primary_and_fill"
        return None

    def erb_weights_and_spl(
        self, instrument_name: str, erb_edges: np.ndarray
    ) -> Optional[Tuple[np.ndarray, np.ndarray, float, str, float]]:
        """Return (relative_w, spl_db_per_band, ref_distance_m, provenance, fill_fraction).

        ``None`` if no usable coverage (caller keeps equipartition). A vector
        that is mostly residual fill is refused (not labelled as a measurement).
        """
        key = self.source_key(instrument_name)
        if key is None:
            return None
        lo, hi, e, _is_meas, fill_frac, _whole, dist = self.historical_band_powers(
            key
        )
        if e.size == 0 or e.sum() <= 0:
            return None
        provenance = self.provenance_for_fill(float(fill_frac))
        if provenance is None:
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
        return w, spl, dist, provenance, float(fill_frac)


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
    fill_fraction: Optional[float] = None
    stroke: Optional[str] = None
    dynamic: Optional[str] = None
    t_contact_s: Optional[float] = None
    t_contact_provenance: Optional[str] = None
    f_c_hz: Optional[float] = None

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
        note_str = " | ".join(self.notes) if self.notes else ""
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
                "fill_fraction": self.fill_fraction
                if self.fill_fraction is not None else "",
                "stroke": self.stroke if self.stroke is not None else "",
                "dynamic": self.dynamic if self.dynamic is not None else "",
                "t_contact_s": self.t_contact_s
                if self.t_contact_s is not None else "",
                "t_contact_provenance": self.t_contact_provenance
                if self.t_contact_provenance is not None else "",
                "f_c_hz": self.f_c_hz if self.f_c_hz is not None else "",
                "notes": note_str,
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
    stroke: Optional[str] = None,
    dynamic: Optional[str] = None,
    csv_path: Optional[Path] = None,
    t_contact_base_s: Optional[float] = None,
) -> DensityProfile:
    """Generate an ERB density profile; optionally apply AmplitudeLayer.

    When ``stroke`` and ``dynamic`` are **both** ``None``, behaviour is
    bit-identical to v0.3.3 (equipartition / AmplitudeLayer path with full
    shimmer boost). When both are set, a Hertzian contact-time low-pass
    filters the initial energy and the 3–5 kHz shimmer boost is
    dynamic-gated — except at ``dynamic="ff"`` for **plates**, where the
    low-pass is bypassed ([R] §9.4 nonlinear cascade regenerates HF) while
    the boost stays at full amplitude.

    ``t_contact_base_s`` optionally overrides the loaded base contact time
    before dynamic scaling (used by Monte Carlo perturbation).
    """
    edges = erb_band_edges(f_lo, f_hi)
    centres = np.sqrt(edges[:-1] * edges[1:])
    counts = _band_mode_counts(instr, edges)
    if isinstance(instr, PlateInstrument):
        family = "plate"
        pclass = getattr(instr, "plate_class", "cymbal")
        if pclass == "tamtam":
            phases = PLATE_PHASES_TAMTAM
            notes = [
                "flat-plate approximation: dome/curvature not modelled",
                "linear modal regime: chaotic broadband regime at ff not modelled",
                "tam-tam temporal template (literature_derived): slow HF bloom "
                "over ~1-2 s; HF emphasis in shimmer not bloom "
                "([R] tam-tam / [FR] Ch. 20)",
            ]
            if "wind" in instr.name.lower() or "windgong" in instr.name.lower():
                notes.append(
                    "wind gong uses tam-tam phase template with note "
                    "(distinct subtype; internal_default mapping)"
                )
        else:
            phases = PLATE_PHASES
            notes = [
                "flat-plate approximation: dome/curvature not modelled",
                "linear modal regime: chaotic broadband regime at ff not modelled",
            ]
    else:
        family, phases = "membrane", MEMBRANE_PHASES
        pclass = None
        if isinstance(instr, MembraneInstrument) and instr.has_measured_anchor():
            c_eff = instr.fitted_effective_c()
            n_anch = len(instr.measured_modes)  # type: ignore[arg-type]
            notes = [
                f"measured-mode anchor active (Fletcher & Rossing 1998 "
                f"Table 18.5): {n_anch} modes, primary_source; "
                f"fitted effective c = {c_eff:.3f} m/s (derived); "
                f"air-loading / two-head effects absorbed at low order only",
            ]
        elif isinstance(instr, MembraneInstrument) and instr.effective_c is not None:
            notes = [
                f"effective wave speed c = {instr.effective_c:.3f} m/s "
                f"scaled from 82 cm Fletcher & Rossing anchor "
                f"(f ∝ 1/a at fixed c; internal_default); "
                f"in-vacuo bias applies above any local measured range",
            ]
        else:
            notes = [
                "single membrane in vacuo: air loading and 2-head coupling "
                "not modelled (lowest modes overestimated)",
            ]

    e0 = None
    provenance = "internal_default"
    ref_dist: Optional[float] = None
    fill_fraction: Optional[float] = None
    abs_spl: Dict[str, np.ndarray] = {}
    t_contact: Optional[float] = None
    t_prov: Optional[str] = None
    f_c: Optional[float] = None
    # Bit-identity path: both unset → full shimmer boost, no contact filter.
    excitation_active = stroke is not None and dynamic is not None
    shimmer_gate = 1.0
    apply_lp_filter = False
    if excitation_active:
        dyn = str(dynamic).lower()
        shimmer_gate = float(_DYNAMIC_SHIMMER_GATE.get(dyn, 1.0))
        implement = map_stroke_to_implement(
            stroke, plate_class=pclass, family=family
        )
        t0, t_prov, used_ph = load_contact_time_s(implement, csv_path=csv_path)
        if t_contact_base_s is not None:
            t0 = float(t_contact_base_s)
            t_prov = f"{t_prov}+mc_perturbed"
        scale = float(_DYNAMIC_T_CONTACT_SCALE.get(dyn, 1.0))
        t_contact = t0 * scale
        f_c = excitation_cutoff_hz(t_contact)
        # ff plate bypass: nonlinear cascade regenerates HF ([R] §9.4).
        apply_lp_filter = not (dyn == "ff" and family == "plate")
        if apply_lp_filter:
            notes.append(
                f"excitation filter active: stroke={stroke}, dynamic={dyn}, "
                f"implement={implement}, t_contact={t_contact:.4e} s "
                f"(base {t0:.4e} s x dynamic_scale={scale}), "
                f"f_c={f_c:.1f} Hz, t_contact_provenance={t_prov}, "
                f"shimmer_gate={shimmer_gate}"
            )
        else:
            notes.append(
                f"ff excitation low-pass bypassed for plates "
                f"([R] sec. 9.4 nonlinear cascade); shimmer boost at full "
                f"amplitude; stroke={stroke}, implement={implement}, "
                f"t_contact={t_contact:.4e} s (recorded, not applied to E0), "
                f"t_contact_provenance={t_prov}, shimmer_gate={shimmer_gate}"
            )
        if used_ph and t_contact_base_s is None:
            notes.append(
                f"t_contact placeholder in use for {implement} "
                f"(internal_default; awaiting source-read)"
            )

    if amplitude_layer is not None:
        mapped = amplitude_layer.erb_weights_and_spl(instr.name, edges)
        if mapped is not None:
            e0, spl0, ref_dist, provenance, fill_fraction = mapped
            notes.append(
                f"initial energy from AmplitudeLayer ({provenance}); "
                f"fill_fraction={fill_fraction:.3f}; "
                f"ref distance {ref_dist:.4f} m"
            )
            if apply_lp_filter and t_contact is not None:
                e0 = apply_excitation_filter(e0, centres, t_contact)
            rel = _band_energy_weights(
                instr, edges, phases, e0=e0, shimmer_gate=shimmer_gate
            )
            for ph, w in rel.items():
                scale_abs = np.nansum(
                    (10 ** (spl0 / 10.0))[np.isfinite(spl0)]
                )
                lin = w * (scale_abs if scale_abs > 0 else 1.0)
                with np.errstate(divide="ignore"):
                    abs_spl[ph] = 10.0 * np.log10(np.maximum(lin, 1e-30))
            return DensityProfile(
                instr.name, family, edges, centres, counts, rel, notes,
                absolute_spl_db=abs_spl, energy_provenance=provenance,
                ref_distance_m=ref_dist, fill_fraction=fill_fraction,
                stroke=stroke, dynamic=dynamic, t_contact_s=t_contact,
                t_contact_provenance=t_prov, f_c_hz=f_c,
            )
        ff = amplitude_layer.fill_fraction_for(instr.name)
        if ff is not None and ff > FILL_FRAC_MIXED_MAX:
            fill_fraction = ff
            notes.append(
                f"AmplitudeLayer refused: fill_fraction={ff:.3f} > "
                f"{FILL_FRAC_MIXED_MAX} (vector mostly residual fill, not a "
                f"measurement); using equipartition (internal_default)"
            )

    if e0 is None and apply_lp_filter and t_contact is not None:
        e0_eq = counts / max(counts.sum(), 1e-12)
        e0 = apply_excitation_filter(e0_eq, centres, t_contact)
    elif apply_lp_filter and t_contact is not None and e0 is not None:
        e0 = apply_excitation_filter(e0, centres, t_contact)

    weights = _band_energy_weights(
        instr, edges, phases, e0=e0, shimmer_gate=shimmer_gate
    )
    return DensityProfile(
        instr.name, family, edges, centres, counts, weights, notes,
        absolute_spl_db=abs_spl, energy_provenance=provenance,
        ref_distance_m=ref_dist, fill_fraction=fill_fraction,
        stroke=stroke, dynamic=dynamic, t_contact_s=t_contact,
        t_contact_provenance=t_prov, f_c_hz=f_c,
    )