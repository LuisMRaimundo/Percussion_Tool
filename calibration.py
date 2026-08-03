"""
calibration
===========

Scale-commensurability bridge between the model's internal composite
density indices and absolute Sivian/Meyer band data.

Bridge instruments are treated as simple quasi-harmonic fixtures
(partials at n·f0 with band energies from the digitized spectra). This
is not a pitched-instrument model; it exists only to estimate a single
conversion factor and its spread across instruments.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from model import (
    AmplitudeLayer,
    DensityProfile,
    erb_band_edges,
    energy_preserving_remap,
    power_to_spl_db,
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
        abs_spl: Dict[str, np.ndarray] = {}
    else:
        e0, spl0, ref_dist, provenance = mapped
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
        ],
        absolute_spl_db=abs_spl,
        energy_provenance=provenance,
        ref_distance_m=ref_dist,
    )


def _empirical_index_from_powers(name: str, layer: AmplitudeLayer) -> float:
    """Shannon-entropy style index on absolute band powers (same aggregator)."""
    edges = erb_band_edges(20.0, 16000.0)
    lo, hi, e, whole, dist = layer.historical_band_powers(
        layer.source_key(name) or name
    )
    if e.size == 0:
        return float("nan")
    erb_e = energy_preserving_remap(lo, hi, e, edges)
    w = erb_e / erb_e.sum() if erb_e.sum() > 0 else erb_e
    nz = w > 0
    H = -np.sum(w[nz] * np.log(w[nz]))
    # Occupation proxy: number of historical bands with power, remapped.
    mean_occ = float(np.sum(w * (erb_e > 0)))
    return float(np.exp(H) * max(mean_occ, 1e-12)) ** 0.5


def run_calibration(
    layer: Optional[AmplitudeLayer] = None,
) -> Tuple[float, float, List[BridgeResult]]:
    """Return (factor_mean, factor_spread_stdev, per-instrument results)."""
    layer = layer or AmplitudeLayer.default(ROOT)
    results: List[BridgeResult] = []
    for name, f0, npart in BRIDGE_INSTRUMENTS:
        if not layer.has_coverage(name):
            continue
        prof = _quasi_harmonic_profile(name, f0, npart, layer)
        model_idx = prof.composite_index("strike")
        emp_idx = _empirical_index_from_powers(name, layer)
        if not np.isfinite(model_idx) or not np.isfinite(emp_idx) or model_idx <= 0:
            continue
        results.append(
            BridgeResult(name, model_idx, emp_idx, emp_idx / model_idx)
        )
    if not results:
        return float("nan"), float("nan"), results
    factors = np.array([r.factor for r in results], dtype=float)
    return float(np.mean(factors)), float(np.std(factors, ddof=1) if len(factors) > 1 else 0.0), results


def write_calibration_report(
    path: Optional[Path] = None,
    layer: Optional[AmplitudeLayer] = None,
) -> Tuple[float, float]:
    """Write ``calibration_report.md`` and return (mean, spread)."""
    path = path or ROOT / "calibration_report.md"
    mean, spread, results = run_calibration(layer)
    lines = [
        "# Calibration report — scale commensurability bridge",
        "",
        "Bridge instruments are quasi-harmonic fixtures only (partials at",
        "`n·f0` with band energies from Sivian/Meyer digitized spectra).",
        "They are **not** a pitched-instrument model.",
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
        "## Per-instrument bridge results",
        "",
        "| instrument | model index | empirical index | factor |",
        "|---|---:|---:|---:|",
    ]
    for r in results:
        lines.append(
            f"| {r.name} | {r.model_index:.6g} | {r.empirical_index:.6g} "
            f"| {r.factor:.6g} |"
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
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")
    return mean, spread


if __name__ == "__main__":
    m, s = write_calibration_report()
    print(f"[CAL] factor={m:.6g}  spread={s:.6g}")
