"""Pytest suite: VAL1, VAL2, AmplitudeLayer honesty, calibration, sanity."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from model import (
    FILL_FRAC_MIXED_MAX,
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
        "f_centre_hz", "modes_per_band", "energy_w_strike", "fill_fraction",
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


def test_cymbal_high_fill_refuses_primary_label(cymbal_46: PlateInstrument) -> None:
    """cymbals_15in is ~90% fill → no coverage; equipartition + note."""
    layer = AmplitudeLayer.default(ROOT)
    assert not layer.has_coverage(cymbal_46.name)
    ff = layer.fill_fraction_for(cymbal_46.name)
    assert ff is not None and ff > FILL_FRAC_MIXED_MAX

    base = generate_profile(cymbal_46)
    with_layer = generate_profile(cymbal_46, amplitude_layer=layer)
    assert with_layer.energy_provenance == "internal_default"
    assert with_layer.fill_fraction == pytest.approx(ff)
    assert any("fill_fraction=" in n for n in with_layer.notes)
    assert any("refused" in n.lower() for n in with_layer.notes)
    for ph in base.energy_weights:
        np.testing.assert_allclose(
            base.energy_weights[ph], with_layer.energy_weights[ph], rtol=0, atol=0
        )


def test_equal_density_fill_narrow_band_gets_less_than_wide() -> None:
    """With one measured HF band + residual, 20–62.5 Hz < 2000–2800 Hz."""
    layer = AmplitudeLayer.default(ROOT)
    lo, hi, e, is_meas, fill_frac, whole, _d = layer.historical_band_powers(
        "cymbals_15in"
    )
    assert fill_frac == pytest.approx(0.9, rel=1e-6)
    assert whole == pytest.approx(9.5)
    assert int(np.sum(is_meas)) == 1
    i_low = int(np.where(np.isclose(lo, 20.0) & np.isclose(hi, 62.5))[0][0])
    i_mid = int(np.where(np.isclose(lo, 2000.0) & np.isclose(hi, 2800.0))[0][0])
    assert not is_meas[i_low] and not is_meas[i_mid]
    assert e[i_low] < e[i_mid]


def test_is_measured_mask_marks_only_digitized_bands() -> None:
    layer = AmplitudeLayer.default(ROOT)
    lo, hi, e, is_meas, _ff, _w, _d = layer.historical_band_powers("bass_viol")
    assert int(np.sum(is_meas)) == 2
    assert np.all(e[is_meas] > 0)
    # No residual for bass_viol → fill bands should be zero / unmeasured empty
    assert float(np.sum(e[~is_meas])) == pytest.approx(0.0)


def test_energy_preserving_remap_conserves() -> None:
    hist_lo = np.array([100.0, 200.0])
    hist_hi = np.array([200.0, 400.0])
    hist_e = np.array([2.0, 3.0])
    edges = np.array([100.0, 150.0, 200.0, 300.0, 400.0])
    out = energy_preserving_remap(hist_lo, hist_hi, hist_e, edges)
    assert out.sum() == pytest.approx(5.0, rel=1e-12)


def test_calibration_runs() -> None:
    mean, spread, results, exclusions = run_calibration(AmplitudeLayer.default(ROOT))
    assert len(results) + len(exclusions) >= 3
    # Sparse digitization may leave few bridge members; mean finite if any.
    if results:
        assert np.isfinite(mean)
        assert np.isfinite(spread)
        assert spread >= 0.0
        for r in results:
            assert r.n_measured_bands >= 2


def test_calibration_excludes_single_measured_band(tmp_path: Path) -> None:
    """Synthetic instrument: 1 measured band + large residual → excluded."""
    src = ROOT / "data" / "source_constants.csv"
    df = pd.read_csv(src)
    # Minimal edges + one measured band + whole with large residual.
    edge_rows = df[df["record_type"] == "sivian_band_edge"].copy()
    synth = pd.DataFrame(
        [
            {
                "record_type": "sivian_meta",
                "instrument": "synth_sparse",
                "specimen": "",
                "source": "test",
                "location": "test",
                "provenance": "internal_default",
                "parameter": "peak_power_whole",
                "value": 10.0,
                "units_printed": "W",
                "units_si": "W",
                "value_si": 10.0,
                "band_lo_hz": "",
                "band_hi_hz": "",
                "measurement_year": 1931.0,
                "ref_distance_m": 0.9144,
                "discrepancy_db": "",
                "needs_manual_reading": 0,
                "notes": "synthetic sparse",
            },
            {
                "record_type": "sivian_meta",
                "instrument": "synth_sparse",
                "specimen": "",
                "source": "test",
                "location": "test",
                "provenance": "internal_default",
                "parameter": "peak_power_band",
                "value": 1.0,
                "units_printed": "W",
                "units_si": "W",
                "value_si": 1.0,
                "band_lo_hz": 8000.0,
                "band_hi_hz": 11300.0,
                "measurement_year": 1931.0,
                "ref_distance_m": 0.9144,
                "discrepancy_db": "",
                "needs_manual_reading": 0,
                "notes": "one measured HF band",
            },
        ]
    )
    out_csv = tmp_path / "source_constants.csv"
    pd.concat([edge_rows, synth], ignore_index=True).to_csv(out_csv, index=False)

    layer = AmplitudeLayer(out_csv)
    # Patch bridge list via monkeypatch-style: call internals directly.
    from calibration import BridgeExclusion, _empirical_index_measured_only
    from model import _INSTRUMENT_SOURCE_KEYS

    _INSTRUMENT_SOURCE_KEYS["synth_sparse"] = "synth_sparse"
    try:
        assert layer.n_measured_bands("synth_sparse") == 1
        idx, ff, n = _empirical_index_measured_only("synth_sparse", layer)
        assert n == 1
        assert not np.isfinite(idx)
        assert ff > FILL_FRAC_MIXED_MAX
        # Full run_calibration path: inject into BRIDGE_INSTRUMENTS
        import calibration as cal

        old = list(cal.BRIDGE_INSTRUMENTS)
        cal.BRIDGE_INSTRUMENTS = [("synth_sparse", 200.0, 8)]
        try:
            _m, _s, results, exclusions = cal.run_calibration(layer)
            assert results == []
            assert any(ex.name == "synth_sparse" for ex in exclusions)
            assert isinstance(exclusions[0], BridgeExclusion)
        finally:
            cal.BRIDGE_INSTRUMENTS = old
    finally:
        _INSTRUMENT_SOURCE_KEYS.pop("synth_sparse", None)


def test_calibration_report_written(tmp_path: Path) -> None:
    out = tmp_path / "calibration_report.md"
    mean, spread = write_calibration_report(out, AmplitudeLayer.default(ROOT))
    text = out.read_text(encoding="utf-8")
    assert "Conversion factor" in text
    assert "spread IS the uncertainty" in text
    assert "Sparse-coverage caveat" in text
    assert "fill_fraction" in text
    assert "Excluded from bridge" in text


def test_meyer_hf_provenance_is_literature_derived() -> None:
    df = pd.read_csv(ROOT / "data" / "source_constants.csv")
    meyer = df[df["record_type"] == "meyer_hf_discrepancy"]
    assert len(meyer) == 4
    assert set(meyer["provenance"]) == {"literature_derived"}
