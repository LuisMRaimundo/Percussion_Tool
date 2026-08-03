"""
uncertainty
===========

Monte Carlo uncertainty layer for the idiophone density model.

Perturbs specimen / material / fit parameters and aggregates per-band
and per-phase distributions (median and nested percentiles), restoring
symmetry with the bootstrap-CI convention of the companion empirical
pipeline.

Provenance of perturbation widths is recorded in README / CHANGES;
spans for the Chladni exponent come from Rossing Table 9.1
(``primary_source``).
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model import (
    MATERIALS,
    MYLAR_RHO,
    Material,
    MembraneInstrument,
    PlateInstrument,
    generate_profile,
)

ROOT = Path(__file__).resolve().parent
Instrument = Union[PlateInstrument, MembraneInstrument]

# Lognormal sigma so the 95% multiplicative interval is
# [1/(1+frac), (1+frac)] about the median (asymmetric). internal_default.
def _lognormal_sigma_for_95pct_span(frac: float) -> float:
    return float(np.log(1.0 + frac) / 1.96)


SIGMA_THICKNESS_PLATE = _lognormal_sigma_for_95pct_span(0.25)   # ±25%
SIGMA_THICKNESS_MEMBRANE = _lognormal_sigma_for_95pct_span(0.10)  # ±10%
# Diameter / material / decay: "normal, ±X%" read as 1-sigma = X% of
# nominal (internal_default; tighter than a 95%-span reading).
SIGMA_DIAMETER_FRAC = 0.01
SIGMA_MATERIAL_FRAC = 0.05
SIGMA_DECAY_FRAC = 0.20

DEFAULT_SEED = 20260803
PERCENTILES = (5, 25, 50, 75, 95)
P_LABELS = ("p05", "p25", "p50", "p75", "p95")


# Nearest Table 9.1 class for instruments lacking an exact row.
# 46 cm ≈ 18 in medium (internal_default nearest-class mapping).
_CHLADNI_CLASS_ALIAS = {
    "cymbal_16in_thin": "cymbal_16in_thin",
    "cymbal_18in_medium": "cymbal_18in_medium",
    "cymbal_46cm_medium": "cymbal_18in_medium",
}


def chladni_p_span(
    instrument_name: str,
    csv_path: Optional[Path] = None,
) -> Optional[Tuple[float, float]]:
    """Return (p_lo, p_hi) from Rossing Table 9.1 for the cymbal class.

    Span = min/max of all p1 and p2 rows of that class (primary_source).
    """
    csv_path = csv_path or ROOT / "data" / "source_constants.csv"
    key = _CHLADNI_CLASS_ALIAS.get(instrument_name)
    if key is None or not Path(csv_path).is_file():
        return None
    df = pd.read_csv(csv_path)
    sub = df[
        (df["record_type"] == "rossing_table91")
        & (df["instrument"] == key)
        & (df["parameter"].isin(["p1", "p2"]))
        & df["value_si"].notna()
    ]
    if sub.empty:
        return None
    vals = sub["value_si"].astype(float).to_numpy()
    return float(vals.min()), float(vals.max())


@dataclass
class MonteCarloResult:
    """Aggregated Monte Carlo profiles for one instrument."""

    instrument: str
    family: str
    seed: int
    n_draws: int
    band_edges: np.ndarray
    band_centres: np.ndarray
    # quantiles[metric][label] -> array (n_bands,) or scalar for composites
    modes_quantiles: Dict[str, np.ndarray]
    energy_quantiles: Dict[str, Dict[str, np.ndarray]]  # phase -> label -> arr
    composite_quantiles: Dict[str, Dict[str, float]]    # phase -> label -> float
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_rows(self) -> List[dict]:
        rows = []
        for i, fc in enumerate(self.band_centres):
            row = {
                "instrument": self.instrument,
                "family": self.family,
                "band_index": i,
                "f_lo_hz": float(self.band_edges[i]),
                "f_hi_hz": float(self.band_edges[i + 1]),
                "f_centre_hz": float(fc),
                "mc_seed": self.seed,
                "mc_n_draws": self.n_draws,
            }
            for lab in P_LABELS:
                row[f"modes_per_band_{lab}"] = float(
                    self.modes_quantiles[lab][i]
                )
            # median also as modes_per_band for schema familiarity
            row["modes_per_band"] = float(self.modes_quantiles["p50"][i])
            for ph, qmap in self.energy_quantiles.items():
                for lab in P_LABELS:
                    row[f"energy_w_{ph}_{lab}"] = float(qmap[lab][i])
                row[f"energy_w_{ph}"] = float(qmap["p50"][i])
            rows.append(row)
        return rows


def _perturb_plate(
    instr: PlateInstrument,
    rng: np.random.Generator,
    p_span: Optional[Tuple[float, float]],
    draw_id: int,
) -> Tuple[PlateInstrument, str]:
    """Return perturbed plate + temporary MATERIALS key (may be '')."""
    h = float(
        rng.lognormal(np.log(instr.thickness), SIGMA_THICKNESS_PLATE)
    )
    d = float(rng.normal(instr.diameter, SIGMA_DIAMETER_FRAC * instr.diameter))
    d = max(d, 1e-4)
    h = max(h, 1e-6)

    base = MATERIALS[instr.material]
    E = float(rng.normal(base.E, SIGMA_MATERIAL_FRAC * base.E))
    rho = float(rng.normal(base.rho, SIGMA_MATERIAL_FRAC * base.rho))
    E = max(E, 1e6)
    rho = max(rho, 1.0)
    key = f"_mc_{instr.material}_{draw_id}"
    # Not thread-safe: temporary MATERIALS mutation; relies on finally cleanup.
    MATERIALS[key] = Material(key, E=E, rho=rho, nu=base.nu)

    chladni = instr.chladni
    if chladni is not None and p_span is not None:
        c, b, _p = chladni
        p = float(rng.uniform(p_span[0], p_span[1]))
        chladni = (c, b, p)

    tau = float(
        rng.normal(instr.decay_tau_1k, SIGMA_DECAY_FRAC * instr.decay_tau_1k)
    )
    alpha = float(
        rng.normal(instr.decay_alpha, SIGMA_DECAY_FRAC * instr.decay_alpha)
    )
    tau = max(tau, 1e-3)
    alpha = max(alpha, 1e-3)

    out = replace(
        instr,
        diameter=d,
        thickness=h,
        material=key,
        chladni=chladni,
        decay_tau_1k=tau,
        decay_alpha=alpha,
    )
    return out, key


def _perturb_membrane(
    instr: MembraneInstrument,
    rng: np.random.Generator,
) -> MembraneInstrument:
    h = float(
        rng.lognormal(
            np.log(instr.membrane_thickness), SIGMA_THICKNESS_MEMBRANE
        )
    )
    d = float(rng.normal(instr.diameter, SIGMA_DIAMETER_FRAC * instr.diameter))
    d = max(d, 1e-4)
    h = max(h, 1e-8)
    rho0 = MYLAR_RHO if instr.membrane_rho is None else instr.membrane_rho
    rho = float(rng.normal(rho0, SIGMA_MATERIAL_FRAC * rho0))
    rho = max(rho, 1.0)
    tau = float(
        rng.normal(instr.decay_tau_100, SIGMA_DECAY_FRAC * instr.decay_tau_100)
    )
    alpha = float(
        rng.normal(instr.decay_alpha, SIGMA_DECAY_FRAC * instr.decay_alpha)
    )
    return replace(
        instr,
        diameter=d,
        membrane_thickness=h,
        membrane_rho=rho,
        decay_tau_100=max(tau, 1e-4),
        decay_alpha=max(alpha, 1e-3),
    )


def _quantile_stack(samples: np.ndarray) -> Dict[str, np.ndarray]:
    """samples shape (n_draws, n_bands) -> label -> (n_bands,)."""
    qs = np.percentile(samples, PERCENTILES, axis=0)
    return {lab: qs[i] for i, lab in enumerate(P_LABELS)}


def run_monte_carlo(
    instrument: Instrument,
    n_draws: int = 2000,
    seed: int = DEFAULT_SEED,
    amplitude_layer=None,
    stroke: Optional[str] = None,
    dynamic: Optional[str] = None,
) -> MonteCarloResult:
    """Perturb instrument parameters and aggregate profile distributions.

    Parameters
    ----------
    instrument
        Baseline ``PlateInstrument`` or ``MembraneInstrument``.
    n_draws
        Number of Monte Carlo draws (default 2000).
    seed
        Fixed RNG seed recorded in output metadata.
    amplitude_layer
        Optional ``AmplitudeLayer``; default ``None`` so the MC isolates
        physical-parameter uncertainty (internal_default).
    stroke, dynamic
        Forwarded to ``generate_profile`` (contact-time excitation filter).
        Both ``None`` keeps the v0.3.3 equipartition path.
    """
    rng = np.random.default_rng(seed)
    p_span = None
    if isinstance(instrument, PlateInstrument) and instrument.chladni is not None:
        p_span = chladni_p_span(instrument.name)

    temp_keys: List[str] = []
    modes_list: List[np.ndarray] = []
    energy_lists: Dict[str, List[np.ndarray]] = {}
    composite_lists: Dict[str, List[float]] = {}
    band_edges = None
    band_centres = None
    family = None

    try:
        for i in range(n_draws):
            if isinstance(instrument, PlateInstrument):
                draw, key = _perturb_plate(instrument, rng, p_span, i)
                if key:
                    temp_keys.append(key)
            else:
                draw = _perturb_membrane(instrument, rng)

            prof = generate_profile(
                draw,
                amplitude_layer=amplitude_layer,
                stroke=stroke,
                dynamic=dynamic,
            )
            if band_edges is None:
                band_edges = prof.band_edges.copy()
                band_centres = prof.band_centres.copy()
                family = prof.family
                for ph in prof.energy_weights:
                    energy_lists[ph] = []
                    composite_lists[ph] = []

            modes_list.append(prof.modes_per_band.copy())
            for ph, w in prof.energy_weights.items():
                energy_lists[ph].append(w.copy())
                composite_lists[ph].append(prof.composite_index(ph))
    finally:
        for key in temp_keys:
            MATERIALS.pop(key, None)

    modes_arr = np.stack(modes_list, axis=0)
    modes_q = _quantile_stack(modes_arr)
    energy_q: Dict[str, Dict[str, np.ndarray]] = {}
    for ph, lst in energy_lists.items():
        energy_q[ph] = _quantile_stack(np.stack(lst, axis=0))
    comp_q: Dict[str, Dict[str, float]] = {}
    for ph, lst in composite_lists.items():
        arr = np.asarray(lst, dtype=float)
        qs = np.percentile(arr, PERCENTILES)
        comp_q[ph] = {lab: float(qs[i]) for i, lab in enumerate(P_LABELS)}

    meta = {
        "seed": seed,
        "n_draws": n_draws,
        "thickness_sigma_log": (
            SIGMA_THICKNESS_PLATE
            if isinstance(instrument, PlateInstrument)
            else SIGMA_THICKNESS_MEMBRANE
        ),
        "thickness_95pct_span": (
            "±25%" if isinstance(instrument, PlateInstrument) else "±10%"
        ),
        "diameter_sigma_frac": SIGMA_DIAMETER_FRAC,
        "material_sigma_frac": SIGMA_MATERIAL_FRAC,
        "decay_sigma_frac": SIGMA_DECAY_FRAC,
        "chladni_p_span": list(p_span) if p_span else None,
        "chladni_p_span_provenance": (
            "primary_source" if p_span else None
        ),
        "amplitude_layer": amplitude_layer is not None,
    }
    return MonteCarloResult(
        instrument=instrument.name,
        family=family or "",
        seed=seed,
        n_draws=n_draws,
        band_edges=band_edges,
        band_centres=band_centres,
        modes_quantiles=modes_q,
        energy_quantiles=energy_q,
        composite_quantiles=comp_q,
        metadata=meta,
    )


def export_mc_csv(
    results: Sequence[MonteCarloResult],
    path: Optional[Path] = None,
) -> Path:
    """Write ``density_profiles_mc.csv`` with percentile columns."""
    path = path or ROOT / "density_profiles_mc.csv"
    rows: List[dict] = []
    for r in results:
        rows.extend(r.to_rows())
    if not rows:
        raise ValueError("no MC rows to export")
    fieldnames = sorted(
        set().union(*[set(x.keys()) for x in rows]),
        key=lambda k: (k not in ("instrument", "family"), k),
    )
    for r in rows:
        for fn in fieldnames:
            r.setdefault(fn, "")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    # Sidecar metadata (seed) for the first result / all instruments.
    meta_path = path.with_suffix(".meta.json")
    meta_path.write_text(
        json.dumps(
            {
                "results": [
                    {
                        "instrument": r.instrument,
                        "seed": r.seed,
                        "n_draws": r.n_draws,
                        **r.metadata,
                    }
                    for r in results
                ]
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return path


def plot_fan_chart(
    result: MonteCarloResult,
    phase: str = "shimmer",
    path: Optional[Path] = None,
    quantity: str = "energy",
) -> Path:
    """Fan chart: median + 50% and 90% bands for the 46-cm-style result."""
    path = path or ROOT / f"density_profiles_mc_fan_{result.instrument}.png"
    x = result.band_centres
    if quantity == "modes":
        q = result.modes_quantiles
        ylabel = "modes per ERB band"
        title = f"{result.instrument}: modes/band MC fan (seed={result.seed})"
    else:
        if phase not in result.energy_quantiles:
            phase = next(iter(result.energy_quantiles))
        q = result.energy_quantiles[phase]
        ylabel = f"relative band energy ({phase})"
        title = (
            f"{result.instrument}: {phase} energy MC fan "
            f"(seed={result.seed}, N={result.n_draws})"
        )

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.fill_between(
        x, q["p05"], q["p95"], color="C0", alpha=0.20, label="90% band"
    )
    ax.fill_between(
        x, q["p25"], q["p75"], color="C0", alpha=0.35, label="50% band"
    )
    ax.plot(x, q["p50"], color="C0", lw=1.8, label="median")
    ax.set(
        xscale="log",
        xlabel="ERB-band centre (Hz)",
        ylabel=ylabel,
        title=title,
    )
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def assert_nested_percentiles(result: MonteCarloResult) -> None:
    """Raise AssertionError if percentile nesting fails anywhere."""
    for lab_a, lab_b in zip(P_LABELS, P_LABELS[1:]):
        assert np.all(
            result.modes_quantiles[lab_a] <= result.modes_quantiles[lab_b] + 1e-12
        ), (lab_a, lab_b, "modes")
        for ph, qmap in result.energy_quantiles.items():
            assert np.all(qmap[lab_a] <= qmap[lab_b] + 1e-12), (ph, lab_a, lab_b)
        for ph, qmap in result.composite_quantiles.items():
            assert qmap[lab_a] <= qmap[lab_b] + 1e-12, (ph, lab_a, lab_b)


if __name__ == "__main__":
    cy = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    res = run_monte_carlo(cy, n_draws=200, seed=DEFAULT_SEED)
    export_mc_csv([res])
    plot_fan_chart(res, phase="shimmer")
    print(
        f"[MC] {res.instrument}: seed={res.seed} N={res.n_draws} "
        f"shimmer index p50={res.composite_quantiles['shimmer']['p50']:.4f}"
    )
