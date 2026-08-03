"""Metadata-only sample grouping tests (no audio parameter fitting)."""

from __future__ import annotations

from pathlib import Path

from sample_metadata import parse_sample_path, resolve_model


def test_parse_crash_stick_dynamic() -> None:
    p = Path("Samples/Cymbals/17crash/17crash.stick.bell.ff.aif")
    m = parse_sample_path(p)
    assert m.status == "ok"
    assert m.diameter_in == 17
    assert m.subtype == "crash"
    assert m.stroke == "stick.bell"
    assert m.dynamic == "ff"
    assert m.group_key == ("cymbal_17in_crash", "stick.bell", "ff")
    assert not m.contributes_to_aggregate  # ff


def test_parse_mallet_pp_primary() -> None:
    p = Path("Samples/Cymbals/mallet.pp/18crash.mallet.pp.aif")
    m = parse_sample_path(p)
    assert m.status == "ok"
    assert m.stroke == "mallet"
    assert m.dynamic == "pp"
    assert m.contributes_to_aggregate


def test_parse_chinese_transfer_caution() -> None:
    p = Path("Samples/Cymbals/mallet.mf/16chinese.mallet.mf.aif")
    m = parse_sample_path(p)
    assert m.status == "ok"
    assert m.is_transfer_caution
    assert not m.contributes_to_aggregate


def test_parse_splash_without_size_unparseable() -> None:
    p = Path("Samples/Cymbals/mallet.mf/splash.mallet.mf.aif")
    m = parse_sample_path(p)
    assert m.status == "unparseable"
    assert "size" in (m.skip_reason or "")


def test_skip_thai_gong_note() -> None:
    p = Path("Samples/Gong/thaigong.mf/thaigong.A4.mf.aif")
    m = parse_sample_path(p)
    assert m.status == "skip_pitched"
    assert m.group_key is None


def test_parse_tamtam_and_windgong() -> None:
    tt = parse_sample_path(Path("Samples/Tam-Tam/tamtam.mf/22tamtam.mf.aif"))
    assert tt.status == "ok"
    assert tt.plate_class == "tamtam"
    assert tt.stroke == "unmarked"
    assert tt.dynamic == "mf"

    wg = parse_sample_path(Path("Samples/Tam-Tam/tamtam.pp/20windgong.pp.aif"))
    assert wg.status == "ok"
    assert wg.plate_class == "windgong"
    assert wg.subtype == "windgong"


def test_exact_catalogue_inch_and_parametric() -> None:
    m18 = parse_sample_path(
        Path("Samples/Cymbals/18crash/18crash.stick.normal.mf.aif")
    )
    map18 = resolve_model(m18)
    assert map18.catalogue_match == "cymbal_18in_medium"
    assert "exact_catalogue" in map18.provenance

    m17 = parse_sample_path(
        Path("Samples/Cymbals/17crash/17crash.stick.normal.mf.aif")
    )
    map17 = resolve_model(m17)
    assert map17.catalogue_match is None
    assert map17.provenance == "filename_metadata + internal_default_thickness"
    assert abs(map17.diameter_m - 17 * 0.0254) < 1e-9


def test_tamtam_parametric_bronze() -> None:
    m = parse_sample_path(Path("Samples/Tam-Tam/tamtam.mf/40tamtam.mf.aif"))
    mapping = resolve_model(m)
    assert mapping.instrument.material == "bronze_B20"
    assert mapping.thickness_m == 0.0015
    assert mapping.provenance == "filename_metadata + internal_default_thickness"
