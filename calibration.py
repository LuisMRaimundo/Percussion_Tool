"""
calibration
===========

Scale-commensurability bridge between the model's internal composite
density indices and absolute Sivian/Meyer band data.

Bridge instruments are treated as simple quasi-harmonic fixtures
(partials at n·f0 with band energies from the digitized spectra). This
is not a pitched-instrument model; it exists only to estimate a single
conversion factor and its spread across instruments.

Empirical indices use **measured historical bands only** (residual fill
excluded) so the conversion factor does not trend toward the uniform
maximum merely because of sparse digitization.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from model import (
    ERB_MEASURED_ENERGY_FRAC,
    AmplitudeLayer,
    DensityProfile,
    erb_band_edges,
    energy_preserving_remap,
)

ROOT = Path(__file__).resolve().parent

# Bridge fixtures: (name, f0_Hz, n_partials). internal_default f0 choices
# for typical orchestral tessitura of each instrument.
BRIDGE_INSTRUMENTS: List[Tuple[str, float, int]] = [
    ("trumpet", 233.0, 12),    # Bb3 vicinity
    ("clarinet", 147.0, 14),   # D3 vicinity
    ("flute", 349.0, 10),      # F4 vicinity
    ("bass_viol", 49.0, 16),   # string-family stand-in (violin spectrum absent)
]


@dataclass
class BridgeResult:
    name: str
    model_index: float
    empirical_index: float
    factor: float  # empirical / model
    fill_fraction: float
    n_measured_bands: int


@dataclass
class BridgeExclusion:
    name: str
    reason: str
    fill_fraction: Optional[float] = None
    n_measured_bands: int = 0


def measured_erb_mask(
    lo: np.ndarray,
    hi: np.ndarray,
    energy: np.ndarray,
    is_measured: np.ndarray,
    erb_edges: np.ndarray,
    frac_threshold: float = ERB_MEASURED_ENERGY_FRAC,
) -> np.ndarray:
    """ERB band is measured if > ``frac_threshold`` of its energy is measured.

    Threshold is ``internal_default`` (default 50%).
    """
    e_meas = np.where(is_measured, energy, 0.0)
    erb_total = energy_preserving_remap(lo, hi, energy, erb_edges)
    erb_meas = energy_preserving_remap(lo, hi, e_meas, erb_edges)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(erb_total > 0, erb_meas / erb_total, 0.0)
    return frac > frac_threshold


def _quasi_harmonic_profile(
    name: str,
    f0: float,
    n_partials: int,
    layer: AmplitudeLayer,
) -> DensityProfile:
    """Build a bridge DensityProfile from partials + measured band powers."""
    edges = erb_band_edges(20.0, 16000.0)
    centres = np.sqrt(edges[:-1] * edges[1:])
    partials = f0 * np.arange(1, n_partials + 1, dtype=float)
    modes = np.histogram(partials, bins=edges)[0].astype(float)

    mapped = layer.erb_weights_and_spl(name, edges)
    if mapped is None:
        # Fall back: put equal energy on partial-occupied bands.
        e0 = modes / max(modes.sum(), 1e-12)
        provenance = "internal_default"
        ref_dist = 0.9144
        fill_fraction = layer.fill_fraction_for(name)
        abs_spl: Dict[str, np.ndarray] = {}
    else:
        e0, spl0, ref_dist, provenance, fill_fraction = mapped
        abs_spl = {"strike": spl0}

    # Single "strike" phase for the bridge fixture.
    w = e0 / e0.sum() if e0.sum() > 0 else e0
    return DensityProfile(
        instrument=name,
        family="bridge_quasi_harmonic",
        band_edges=edges,
        band_centres=centres,
        modes_per_band=modes,
        energy_weights={"strike": w},
        notes=[
            f"bridge fixture: partials n*f0, f0={f0} Hz, N={n_partials}",
            f"energy provenance: {provenance}",
            f"fill_fraction: {fill_fraction}",
        ],
        absolute_spl_db=abs_spl,
        energy_provenance=provenance,
        ref_distance_m=ref_dist,
        fill_fraction=fill_fraction,
    )


def _empirical_index_measured_only(
    name: str, layer: AmplitudeLayer
) -> Tuple[float, float, int]:
    """Shannon-entropy index on measured ERB bands only.

    Returns (index, fill_fraction, n_measured_historical_bands).
    Index is NaN when fewer than 2 measured historical bands or fewer
    than 2 measured ERB bands after remap.
    """
    key = layer.source_key(name) or name
    edges = erb_band_edges(20.0, 16000.0)
    lo, hi, e, is_meas, fill_frac, _whole, _dist = layer.historical_band_powers(
        key
    )
    n_meas = int(np.sum(is_meas)) if is_meas.size else 0
    if e.size == 0 or n_meas < 2:
        return float("nan"), float(fill_frac), n_meas

    mask = measured_erb_mask(lo, hi, e, is_meas, edges)
    if int(np.sum(mask)) < 2:
        return float("nan"), float(fill_frac), n_meas

    erb_e = energy_preserving_remap(lo, hi, e, edges)
    # Zero out fill-dominated ERB bands before renormalizing.
    erb_e = np.where(mask, erb_e, 0.0)
    s = float(erb_e.sum())
    if s <= 0:
        return float("nan"), float(fill_frac), n_meas
    w = erb_e / s
    nz = w > 0
    H = -np.sum(w[nz] * np.log(w[nz]))
    mean_occ = float(np.sum(w * (erb_e > 0)))
    idx = float(np.exp(H) * max(mean_occ, 1e-12)) ** 0.5
    return idx, float(fill_frac), n_meas


def run_calibration(
    layer: Optional[AmplitudeLayer] = None,
) -> Tuple[float, float, List[BridgeResult], List[BridgeExclusion]]:
    """Return (factor_mean, factor_spread_stdev, results, exclusions)."""
    layer = layer or AmplitudeLayer.default(ROOT)
    results: List[BridgeResult] = []
    exclusions: List[BridgeExclusion] = []
    for name, f0, npart in BRIDGE_INSTRUMENTS:
        n_meas = layer.n_measured_bands(name)
        ff = layer.fill_fraction_for(name)
        if n_meas < 2:
            exclusions.append(
                BridgeExclusion(
                    name,
                    f"fewer than 2 measured bands (n={n_meas}); "
                    "degenerate empirical index excluded",
                    fill_fraction=ff,
                    n_measured_bands=n_meas,
                )
            )
            continue
        if not layer.has_coverage(name):
            exclusions.append(
                BridgeExclusion(
                    name,
                    "AmplitudeLayer coverage refused "
                    f"(fill_fraction={ff})",
                    fill_fraction=ff,
                    n_measured_bands=n_meas,
                )
            )
            continue
        prof = _quasi_harmonic_profile(name, f0, npart, layer)
        model_idx = prof.composite_index("strike")
        emp_idx, fill_frac, n_meas = _empirical_index_measured_only(name, layer)
        if not np.isfinite(model_idx) or not np.isfinite(emp_idx) or model_idx <= 0:
            exclusions.append(
                BridgeExclusion(
                    name,
                    "non-finite model or measured-only empirical index",
                    fill_fraction=fill_frac,
                    n_measured_bands=n_meas,
                )
            )
            continue
        results.append(
            BridgeResult(
                name,
                model_idx,
                emp_idx,
                emp_idx / model_idx,
                fill_fraction=fill_frac,
                n_measured_bands=n_meas,
            )
        )
    if not results:
        return float("nan"), float("nan"), results, exclusions
    factors = np.array([r.factor for r in results], dtype=float)
    spread = (
        float(np.std(factors, ddof=1)) if len(factors) > 1 else 0.0
    )
    return float(np.mean(factors)), spread, results, exclusions


def write_calibration_report(
    path: Optional[Path] = None,
    layer: Optional[AmplitudeLayer] = None,
) -> Tuple[float, float]:
    """Write ``calibration_report.md`` and return (mean, spread)."""
    path = path or ROOT / "calibration_report.md"
    mean, spread, results, exclusions = run_calibration(layer)
    lines = [
        "# Calibration report — scale commensurability bridge",
        "",
        "Bridge instruments are quasi-harmonic fixtures only (partials at",
        "`n·f0` with band energies from Sivian/Meyer digitized spectra).",
        "They are **not** a pitched-instrument model.",
        "",
        "Empirical indices are computed on **measured bands only** (ERB",
        "bands whose remapped energy is >50% from digitized",
        "`peak_power_band` rows; threshold `internal_default`). Residual",
        "equal-density fill is excluded from the bridge index.",
        "",
        f"**Conversion factor** (empirical index / model index): "
        f"**{mean:.6g}**",
        "",
        f"**Spread across bridge instruments** (sample stdev): "
        f"**{spread:.6g}**",
        "",
        "> The spread IS the uncertainty to be attached to any",
        "> cross-domain ratio claim that converts the model's internal",
        "> composite index into the empirical absolute scale (or the",
        "> reverse).",
        "",
        "> **Sparse-coverage caveat:** this factor is anchored on sparse",
        "> band coverage until the `needs_manual_reading` histogram cells",
        "> listed in `data/README.md` are completed. Until then, treat the",
        "> conversion factor as provisional.",
        "",
        "## Per-instrument bridge results",
        "",
        "| instrument | fill_fraction | n_measured | model index | "
        "empirical index (measured bands) | factor |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for r in results:
        lines.append(
            f"| {r.name} | {r.fill_fraction:.3f} | {r.n_measured_bands} | "
            f"{r.model_index:.6g} | {r.empirical_index:.6g} "
            f"| {r.factor:.6g} |"
        )
    if not results:
        lines.append("| _(none)_ | — | — | — | — | — |")
    lines.extend(
        [
            "",
            "## Excluded from bridge",
            "",
        ]
    )
    if not exclusions:
        lines.append("_None._")
    else:
        for ex in exclusions:
            ff = (
                f"{ex.fill_fraction:.3f}"
                if ex.fill_fraction is not None
                else "—"
            )
            lines.append(
                f"- **{ex.name}**: {ex.reason} "
                f"(fill_fraction={ff}, n_measured={ex.n_measured_bands})"
            )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Violin full-band peak spectrum is not in Sivian et al. (1931);",
            "  `bass_viol` stands in as the string-family bridge member.",
            "- Soft violin average pressure (0.52 bars at 3 ft) is recorded",
            "  in `data/source_constants.csv` but is not used as a spectral",
            "  bridge weight.",
            "- Reference distance for Sivian absolute levels: 3 ft (0.9144 m)",
            "  unless a row states otherwise.",
            "- Instruments with fill_fraction > 0.60 are refused by",
            "  AmplitudeLayer (mostly residual fill) and cannot enter the",
            "  bridge via measured absolute weights.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    return mean, spread


if __name__ == "__main__":
    m, s = write_calibration_report()
    print(f"[CAL] factor={m:.6g}  spread={s:.6g}")
