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
    fit_effective_wave_speed,
    generate_profile,
    energy_preserving_remap,
    load_contact_time_s,
    make_bassdrum_catalogue,
    measured_modes_from_csv,
    membrane_beta,
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
    # Factor defined only with ≥2 survivors; otherwise nan.
    if len(results) >= 2:
        assert np.isfinite(mean)
        assert np.isfinite(spread)
        assert spread >= 0.0
        for r in results:
            assert r.n_measured_bands >= 2
            # Theory vs measured sides must be independent constructions.
            assert r.model_index != pytest.approx(r.empirical_index, rel=0, abs=0)
    else:
        assert not np.isfinite(mean)
        assert not np.isfinite(spread)


def test_calibration_no_factor_with_single_survivor() -> None:
    """Fewer than 2 survivors → NO CALIBRATION ACHIEVED (nan, nan)."""
    from calibration import NO_CALIBRATION_MSG, format_calibration_cli_line

    mean, spread, results, exclusions = run_calibration(AmplitudeLayer.default(ROOT))
    # Shipping Sivian digitization currently leaves < 2 survivors.
    assert len(results) < 2
    assert not np.isfinite(mean) and not np.isfinite(spread)
    cli = format_calibration_cli_line(mean, spread)
    assert "NO CALIBRATION ACHIEVED" in cli
    assert NO_CALIBRATION_MSG in cli


def test_calibration_theory_index_independent_partial_overlap(tmp_path: Path) -> None:
    """Partials and measured bands only partially overlap → indices differ."""
    import calibration as cal
    from model import _INSTRUMENT_SOURCE_KEYS

    src = ROOT / "data" / "source_constants.csv"
    df = pd.read_csv(src)
    edge_rows = df[df["record_type"] == "sivian_band_edge"].copy()
    # Two measured mid bands; whole = covered ⇒ fill_fraction 0 → coverage ok.
    # f0=80 Hz puts many partials outside 250–700 Hz.
    synth = pd.DataFrame(
        [
            {
                "record_type": "sivian_meta",
                "instrument": "synth_overlap",
                "specimen": "",
                "source": "test",
                "location": "test",
                "provenance": "internal_default",
                "parameter": "peak_power_whole",
                "value": 2.0,
                "units_printed": "W",
                "units_si": "W",
                "value_si": 2.0,
                "band_lo_hz": "",
                "band_hi_hz": "",
                "measurement_year": 1931.0,
                "ref_distance_m": 0.9144,
                "discrepancy_db": "",
                "needs_manual_reading": 0,
                "notes": "synthetic partial-overlap",
            },
            {
                "record_type": "sivian_meta",
                "instrument": "synth_overlap",
                "specimen": "",
                "source": "test",
                "location": "test",
                "provenance": "internal_default",
                "parameter": "peak_power_band",
                "value": 1.2,
                "units_printed": "W",
                "units_si": "W",
                "value_si": 1.2,
                "band_lo_hz": 250.0,
                "band_hi_hz": 500.0,
                "measurement_year": 1931.0,
                "ref_distance_m": 0.9144,
                "discrepancy_db": "",
                "needs_manual_reading": 0,
                "notes": "",
            },
            {
                "record_type": "sivian_meta",
                "instrument": "synth_overlap",
                "specimen": "",
                "source": "test",
                "location": "test",
                "provenance": "internal_default",
                "parameter": "peak_power_band",
                "value": 0.8,
                "units_printed": "W",
                "units_si": "W",
                "value_si": 0.8,
                "band_lo_hz": 500.0,
                "band_hi_hz": 700.0,
                "measurement_year": 1931.0,
                "ref_distance_m": 0.9144,
                "discrepancy_db": "",
                "needs_manual_reading": 0,
                "notes": "",
            },
        ]
    )
    out_csv = tmp_path / "source_constants.csv"
    pd.concat([edge_rows, synth], ignore_index=True).to_csv(out_csv, index=False)
    layer = AmplitudeLayer(out_csv)
    _INSTRUMENT_SOURCE_KEYS["synth_overlap"] = "synth_overlap"
    try:
        f0, npart = 80.0, 16
        model_idx = cal.theory_bridge_model_index(f0, npart)
        emp_idx, ff, n_meas = cal._empirical_index_measured_only(
            "synth_overlap", layer
        )
        assert n_meas == 2
        assert ff == pytest.approx(0.0)
        assert np.isfinite(model_idx) and np.isfinite(emp_idx)
        assert model_idx != pytest.approx(emp_idx, rel=1e-12, abs=1e-12)
    finally:
        _INSTRUMENT_SOURCE_KEYS.pop("synth_overlap", None)


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
    assert "fill_fraction" in text
    assert "Excluded from bridge" in text
    assert "equal-energy-per-partial" in text or "theory only" in text.lower()
    if np.isfinite(mean):
        assert "Conversion factor" in text
        assert "spread IS the uncertainty" in text
    else:
        assert "NO CALIBRATION ACHIEVED" in text
        # Exclusion lines must not duplicate fill_fraction in the reason.
        for line in text.splitlines():
            if line.startswith("- **") and "fill_fraction=" in line:
                # reason body before the trailing "(fill_fraction=...)"
                body = line.split("(fill_fraction=")[0]
                assert "fill_fraction=" not in body


def test_meyer_hf_provenance_is_literature_derived() -> None:
    df = pd.read_csv(ROOT / "data" / "source_constants.csv")
    meyer = df[df["record_type"] == "meyer_hf_discrepancy"]
    assert len(meyer) == 4
    assert set(meyer["provenance"]) == {"literature_derived"}


def test_val3_measured_modes_reproduce_exactly() -> None:
    """VAL3: anchored low modes == FR Table 18.5; vacuo > measured at low order."""
    freqs, labels = measured_modes_from_csv(
        "bassdrum_82cm", specimen="both_heads"
    )
    assert len(freqs) >= 6
    cats = make_bassdrum_catalogue(ROOT / "data" / "source_constants.csv")
    drum = next(i for i in cats if i.name == "bassdrum_82cm")
    assert drum.has_measured_anchor()
    anchored = drum.low_modes()[: len(freqs)]
    np.testing.assert_allclose(anchored, freqs, rtol=0, atol=0)

    # In-vacuo from m ≥ 2: air-loading-dominated (0,1) and (1,1) overestimated
    hi_f = np.array(
        [float(f) for (m, n), f in zip(labels, freqs) if m >= 2]
    )
    hi_lab = tuple((m, n) for (m, n) in labels if m >= 2)
    c_vac = fit_effective_wave_speed(hi_f, hi_lab, drum.radius)
    for (m, n), f_m in zip(labels, freqs):
        if m > 1:
            continue
        f_v = membrane_beta(m, n) * c_vac / (2.0 * np.pi * drum.radius)
        assert f_v > f_m, f"({m}{n}): vacuo {f_v} should exceed measured {f_m}"

    prof = generate_profile(drum)
    assert any("measured-mode anchor active" in n for n in prof.notes)
    assert any("fitted effective c" in n for n in prof.notes)


def test_membrane_no_anchor_fallback_bit_identical_v031(tmp_path: Path) -> None:
    """No usable measured rows → catalogue matches pre-anchor f11_nominal path."""
    src = pd.read_csv(ROOT / "data" / "source_constants.csv")
    # Flag every bass-drum mode row so the loader returns empty.
    mask = src["record_type"] == "fr_ch18_bassdrum_mode"
    src.loc[mask, "needs_manual_reading"] = 1
    src.loc[mask, "value_si"] = np.nan
    csv_path = tmp_path / "source_constants.csv"
    src.to_csv(csv_path, index=False)

    cats = make_bassdrum_catalogue(csv_path)
    assert [c.name for c in cats] == ["bassdrum_32in", "bassdrum_28in"]
    legacy = [
        MembraneInstrument("bassdrum_32in", 0.813, f11_nominal=60.0),
        MembraneInstrument("bassdrum_28in", 0.711, f11_nominal=72.0),
    ]
    for a, b in zip(cats, legacy):
        np.testing.assert_array_equal(a.low_modes(), b.low_modes())
        pa = generate_profile(a)
        pb = generate_profile(b)
        for ph in pa.energy_weights:
            np.testing.assert_allclose(
                pa.energy_weights[ph], pb.energy_weights[ph], rtol=0, atol=0
            )


def test_bassdrum_size_scaling_uses_fitted_c_not_copied_freqs() -> None:
    """32-in / 28-in inherit effective c; do not copy 82 cm measured Hz."""
    cats = make_bassdrum_catalogue(ROOT / "data" / "source_constants.csv")
    by_name = {c.name: c for c in cats}
    d82 = by_name["bassdrum_82cm"]
    d32 = by_name["bassdrum_32in"]
    assert d82.has_measured_anchor()
    assert not d32.has_measured_anchor()
    assert d32.effective_c == pytest.approx(d82.fitted_effective_c())
    # At fixed c, f_mk ∝ 1/a ⇒ ratio of (1,1)-like scales = diameter ratio
    scale_82 = d82.wave_speed() / d82.diameter
    scale_32 = d32.wave_speed() / d32.diameter
    assert scale_32 / scale_82 == pytest.approx(
        d82.diameter / d32.diameter, rel=1e-9
    )
    # Must not equal the 82 cm measured (01) frequency
    assert float(d32.low_modes()[0]) != pytest.approx(
        float(d82.measured_modes[0]), rel=0, abs=0  # type: ignore[index]
    )
    prof = generate_profile(d32)
    assert any("internal_default" in n for n in prof.notes)


def _hf_share(prof, phase: str = "strike", fmin: float = 3000.0) -> float:
    w = prof.energy_weights[phase]
    m = prof.band_centres >= fmin
    return float(w[m].sum())


def test_generate_profile_stroke_none_bit_identical_v033() -> None:
    """stroke=dynamic=None reproduces v0.3.3 path for cymbal and anchored drum."""
    cy = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    a = generate_profile(cy)
    b = generate_profile(cy, stroke=None, dynamic=None)
    for ph in a.energy_weights:
        np.testing.assert_allclose(
            a.energy_weights[ph], b.energy_weights[ph], rtol=0, atol=0
        )
    drum = next(
        i
        for i in make_bassdrum_catalogue(ROOT / "data" / "source_constants.csv")
        if i.name == "bassdrum_82cm"
    )
    da = generate_profile(drum)
    db = generate_profile(drum, stroke=None, dynamic=None)
    for ph in da.energy_weights:
        np.testing.assert_allclose(
            da.energy_weights[ph], db.energy_weights[ph], rtol=0, atol=0
        )


def test_excitation_hf_monotonic_dynamic_and_stroke() -> None:
    cy = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    shares = [
        _hf_share(generate_profile(cy, stroke="stick.normal", dynamic=d))
        for d in ("pp", "mf", "ff")
    ]
    assert shares[0] <= shares[1] <= shares[2]
    mal = _hf_share(generate_profile(cy, stroke="mallet", dynamic="mf"))
    tip = _hf_share(generate_profile(cy, stroke="stick.normal", dynamic="mf"))
    assert mal < tip


def test_tamtam_template_early_hf_and_shimmer_gate() -> None:
    tt = PlateInstrument(
        "tamtam_80cm_bronze", 0.800, 0.0015, plate_class="tamtam"
    )
    cy = PlateInstrument(
        "cymbal_same_size", 0.800, 0.0015, plate_class="cymbal"
    )
    pt = generate_profile(tt, stroke="mallet", dynamic="mf")
    pc = generate_profile(cy, stroke="mallet", dynamic="mf")
    assert "bloom" in pt.energy_weights
    assert "buildup" in pc.energy_weights
    # ~0.2 s: tam-tam still in bloom (no HF boost); cymbal already in shimmer
    assert _hf_share(pt, "bloom") < _hf_share(pc, "shimmer")
    # Shimmer emphasis active and dynamic-gated on tam-tams
    pp = generate_profile(tt, stroke="mallet", dynamic="pp")
    ff = generate_profile(tt, stroke="mallet", dynamic="ff")
    assert _hf_share(ff, "shimmer") > _hf_share(pp, "shimmer")


def test_contact_time_source_wins_over_placeholder(tmp_path: Path) -> None:
    src = pd.read_csv(ROOT / "data" / "source_constants.csv")
    # Inject a primary_source stick_tip that differs from the placeholder.
    src = src[src["instrument"] != "stick_tip"]
    row = {
        "record_type": "excitation_contact_time",
        "instrument": "stick_tip",
        "specimen": "",
        "source": "Fletcher_Rossing_1998",
        "location": "Ch19_test_inject",
        "provenance": "primary_source",
        "parameter": "t_contact_s",
        "value": 0.00055,
        "units_printed": "s",
        "units_si": "s",
        "value_si": 0.00055,
        "band_lo_hz": "",
        "band_hi_hz": "",
        "measurement_year": 1998.0,
        "ref_distance_m": "",
        "discrepancy_db": "",
        "needs_manual_reading": 0,
        "notes": "test inject",
    }
    src = pd.concat([src, pd.DataFrame([row])], ignore_index=True)
    csv_path = tmp_path / "source_constants.csv"
    src.to_csv(csv_path, index=False)
    t, prov, used_ph = load_contact_time_s("stick_tip", csv_path=csv_path)
    assert t == pytest.approx(0.00055)
    assert prov == "primary_source"
    assert used_ph is False
    cy = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    # Placeholder path flags notes
    ph = generate_profile(
        cy, stroke="stick.normal", dynamic="mf",
        csv_path=ROOT / "data" / "source_constants.csv",
    )
    assert any("placeholder" in n for n in ph.notes)


def test_bassdrum_82cm_amplitude_alias() -> None:
    layer = AmplitudeLayer.default(ROOT)
    assert layer.has_coverage("bassdrum_82cm")
    ff = layer.fill_fraction_for("bassdrum_82cm")
    assert ff is not None and ff == pytest.approx(0.5934959349593496, rel=1e-6)


def test_size_sweep_mc_medians_monotone_intervals_nested() -> None:
    from uncertainty import run_monte_carlo

    sizes = np.linspace(0.30, 0.60, 7)
    rows = []
    for d in sizes:
        mc = run_monte_carlo(
            PlateInstrument(f"cymbal_{d:.2f}", d, 0.0012),
            n_draws=80,
            seed=20260803,
        )
        c = mc.composite_quantiles["shimmer"]
        rows.append(c)
    meds = [r["p50"] for r in rows]
    # Modal density grows with area ⇒ composite index non-decreasing in diameter
    assert all(meds[i] <= meds[i + 1] + 1e-9 for i in range(len(meds) - 1))
    for r in rows:
        assert r["p05"] <= r["p25"] <= r["p50"] <= r["p75"] <= r["p95"]


def test_v035_stick_hf_above_mallet_and_ff_bypass() -> None:
    cy = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    tip = _hf_share(generate_profile(cy, stroke="stick.normal", dynamic="mf"))
    mal = _hf_share(generate_profile(cy, stroke="mallet", dynamic="mf"))
    assert tip > mal
    ff = _hf_share(generate_profile(cy, stroke="stick.normal", dynamic="ff"))
    mf = _hf_share(generate_profile(cy, stroke="stick.normal", dynamic="mf"))
    assert ff >= mf


def test_mc_meta_contains_t_contact_sigma() -> None:
    from uncertainty import SIGMA_T_CONTACT, run_monte_carlo

    cy = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    mc = run_monte_carlo(
        cy, n_draws=30, seed=1, stroke="stick.normal", dynamic="mf"
    )
    assert mc.metadata["t_contact_perturbed"] is True
    assert mc.metadata["t_contact_sigma_log"] == pytest.approx(SIGMA_T_CONTACT)
    assert mc.metadata["t_contact_95pct_span"] == "±50%"
    assert mc.metadata["t_contact_nominal_s"] == pytest.approx(0.15e-3)


def test_mallet_path_bit_identical_v034() -> None:
    """yarn_mallet @ pp/mf unchanged from v0.3.4 (contact time + filter path)."""
    from model import (
        PLATE_PHASES,
        apply_excitation_filter,
        erb_band_edges,
        _band_energy_weights,
        _band_mode_counts,
    )

    t_yarn, prov, used_ph = load_contact_time_s("yarn_mallet")
    assert t_yarn == pytest.approx(3.0e-3)
    assert used_ph is True
    t_beater, _, _ = load_contact_time_s("bass_drum_beater")
    assert t_beater == pytest.approx(6.0e-3)

    cy = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    for dyn, gate, scale in (("mf", 0.5, 1.0), ("pp", 0.0, 1.6)):
        got = generate_profile(cy, stroke="mallet", dynamic=dyn)
        edges = erb_band_edges(20.0, 16000.0)
        centres = np.sqrt(edges[:-1] * edges[1:])
        counts = _band_mode_counts(cy, edges)
        e0 = apply_excitation_filter(
            counts / counts.sum(), centres, t_yarn * scale
        )
        ref = _band_energy_weights(
            cy, edges, PLATE_PHASES, e0=e0, shimmer_gate=gate
        )
        for ph in got.energy_weights:
            np.testing.assert_allclose(
                got.energy_weights[ph], ref[ph], rtol=0, atol=0
            )
        assert got.t_contact_s == pytest.approx(t_yarn * scale)
