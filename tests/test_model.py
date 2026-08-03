"""Pytest suite: VAL1, VAL2, and sanity assertions from the build spec."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from model import (
    AmplitudeLayer,
    MembraneInstrument,
    PlateInstrument,
    erb_band_edges,
    erb_bandwidth,
    generate_profile,
    energy_preserving_remap,
)
from calibration import run_calibration, write_calibration_report

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def cymbal_46() -> PlateInstrument:
    return PlateInstrument(
        "cymbal_46cm_medium",
        0.460,
        0.0012,
        chladni=(14.93, 3.0, 1.557),
    )


def test_val1_modal_density(cymbal_46: PlateInstrument) -> None:
    """VAL1: ~64 modes/kHz; >100 modes below 2 kHz (Wilbur holography)."""
    nd = cymbal_46.modal_density()
    assert 50.0 < nd * 1000.0 < 80.0
    assert abs(nd * 1000.0 - 64.0) < 5.0
    assert nd * 2000.0 > 100.0


def test_val2_chladni_reproduction(cymbal_46: PlateInstrument) -> None:
    """VAL2: n=0 family reproduces Fig. 9.3 fit to ~0% when anchored."""
    low = cymbal_46.low_modes(m_max=7, n_max=0)
    ref = 14.93 * (np.arange(2, 8)) ** 1.557
    err = np.max(np.abs(low[:6] - ref) / ref) * 100.0
    assert err < 0.01


def test_erb_bandwidth_positive() -> None:
    f = np.array([100.0, 1000.0, 5000.0])
    bw = erb_bandwidth(f)
    assert np.all(bw > 0)
    assert bw[0] < bw[1] < bw[2]


def test_erb_band_edges_span() -> None:
    edges = erb_band_edges(20.0, 16000.0)
    assert edges[0] == pytest.approx(20.0, rel=1e-3)
    assert edges[-1] <= 16000.0 + 1.0
    assert np.all(np.diff(edges) > 0)


def test_plate_modes_per_band_grow_with_frequency(cymbal_46: PlateInstrument) -> None:
    """Asymptotic plate n(f) constant ⇒ modes/ERB grow with ERB width."""
    prof = generate_profile(cymbal_46)
    # Compare mid vs high geometric-mean occupation in asymptotic region.
    mid = np.mean(prof.modes_per_band[10:20])
    high = np.mean(prof.modes_per_band[-10:])
    assert high > mid


def test_membrane_modal_density_rises() -> None:
    drum = MembraneInstrument("bassdrum_32in", 0.813, f11_nominal=60.0)
    f = np.array([50.0, 200.0, 800.0])
    n = drum.modal_density(f)
    assert n[0] < n[1] < n[2]


def test_phase_weights_normalize(cymbal_46: PlateInstrument) -> None:
    prof = generate_profile(cymbal_46)
    for ph, w in prof.energy_weights.items():
        assert w.sum() == pytest.approx(1.0, abs=1e-9), ph


def test_csv_schema_columns(cymbal_46: PlateInstrument) -> None:
    rows = generate_profile(cymbal_46).to_rows()
    required = {
        "instrument", "family", "band_index", "f_lo_hz", "f_hi_hz",
        "f_centre_hz", "modes_per_band", "energy_w_strike",
    }
    assert required.issubset(rows[0].keys())


def test_no_coverage_bit_identical_relative_weights() -> None:
    """Instruments without AmplitudeLayer coverage keep equipartition path."""
    gong = PlateInstrument("gong_50cm_bronze", 0.500, 0.0020)
    layer = AmplitudeLayer.default(ROOT)
    assert not layer.has_coverage(gong.name)
    a = generate_profile(gong)
    b = generate_profile(gong, amplitude_layer=layer)
    for ph in a.energy_weights:
        np.testing.assert_allclose(
            a.energy_weights[ph], b.energy_weights[ph], rtol=0, atol=0
        )


def test_coverage_changes_cymbal_weights(cymbal_46: PlateInstrument) -> None:
    layer = AmplitudeLayer.default(ROOT)
    assert layer.has_coverage(cymbal_46.name)
    base = generate_profile(cymbal_46)
    abs_prof = generate_profile(cymbal_46, amplitude_layer=layer)
    assert abs_prof.energy_provenance == "primary_source"
    assert abs_prof.ref_distance_m == pytest.approx(0.9144)
    # Relative shape should differ from equipartition once coverage applies.
    assert not np.allclose(
        base.energy_weights["strike"], abs_prof.energy_weights["strike"]
    )
    assert "strike" in abs_prof.absolute_spl_db


def test_energy_preserving_remap_conserves() -> None:
    hist_lo = np.array([100.0, 200.0])
    hist_hi = np.array([200.0, 400.0])
    hist_e = np.array([2.0, 3.0])
    edges = np.array([100.0, 150.0, 200.0, 300.0, 400.0])
    out = energy_preserving_remap(hist_lo, hist_hi, hist_e, edges)
    assert out.sum() == pytest.approx(5.0, rel=1e-12)


def test_calibration_runs() -> None:
    mean, spread, results = run_calibration(AmplitudeLayer.default(ROOT))
    assert len(results) >= 3
    assert np.isfinite(mean)
    assert np.isfinite(spread)
    assert spread >= 0.0


def test_calibration_report_written(tmp_path: Path) -> None:
    out = tmp_path / "calibration_report.md"
    mean, spread = write_calibration_report(out, AmplitudeLayer.default(ROOT))
    text = out.read_text(encoding="utf-8")
    assert "Conversion factor" in text
    assert "spread IS the uncertainty" in text
    assert np.isfinite(mean)
