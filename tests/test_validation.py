"""Synthetic-WAV tests for validate_against_recordings.py."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import soundfile as sf
from scipy.stats import spearmanr

from model import PlateInstrument, erb_band_edges
from validate_against_recordings import (
    analyse_file,
    aggregate_weights,
    run_validation,
)


def _planted_target(instr: PlateInstrument, edges: np.ndarray) -> np.ndarray:
    """Target relative weights = normalized ERB mode histogram."""
    modes = instr.low_modes(m_max=12, n_max=3)
    modes = modes[(modes >= edges[0]) & (modes < min(edges[-1], 8000.0))]
    hist, _ = np.histogram(modes, bins=edges)
    w = hist.astype(float)
    # Ensure a few empty high bands stay near-zero for rank contrast.
    return w / max(w.sum(), 1e-12)


def _synthetic_cymbal_wav(
    path: Path,
    instr: PlateInstrument,
    edges: np.ndarray,
    target: np.ndarray,
    sr: int = 44100,
    duration: float = 4.0,
    rng: np.random.Generator | None = None,
) -> None:
    """FFT-domain ERB-band planting with slow decay (test fixture).

    Assigns random-phase FFT bins inside each ERB band with total power
    equal to ``target[i]``, then applies a slow exponential envelope so
    the shimmer window still carries the planted spectrum. Modal
    frequencies from ``instr`` define ``target`` only; they are not
    re-injected as sinusoids (keeps the planted rank order clean).
    """
    del instr  # target already encodes the instrument's mode histogram
    rng = rng or np.random.default_rng(0)
    n = int(sr * duration)
    t = np.arange(n) / sr
    env = np.exp(-t / 2.0)

    spec = np.zeros(n // 2 + 1, dtype=np.complex128)
    freqs = np.fft.rfftfreq(n, d=1.0 / sr)
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        mask = (freqs >= lo) & (freqs < hi) & (freqs < 0.45 * sr)
        k = int(mask.sum())
        if k == 0 or target[i] <= 0:
            continue
        amp = np.sqrt(target[i] / k)
        phases = rng.uniform(0.0, 2.0 * np.pi, size=k)
        spec[mask] = amp * np.exp(1j * phases)
    y = np.fft.irfft(spec, n=n).astype(float) * env

    pad = np.zeros(int(0.05 * sr))
    y = np.concatenate([pad, y])
    y = y / (np.max(np.abs(y)) + 1e-12) * 0.8
    sf.write(str(path), y.astype(np.float32), sr)


@pytest.fixture
def synthetic_bundle(tmp_path: Path):
    instr = PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    )
    edges = erb_band_edges(20.0, 16000.0)
    target = _planted_target(instr, edges)
    d = tmp_path / "wavs"
    d.mkdir()
    for i in range(3):
        _synthetic_cymbal_wav(
            d / f"synth_{i}.wav",
            instr,
            edges,
            target,
            rng=np.random.default_rng(100 + i),
        )
    return d, target, edges, instr


def test_synthetic_recovers_band_structure(synthetic_bundle) -> None:
    """Pipeline recovers planted ERB band-power rank order (ρ > 0.9)."""
    wav_dir, target, edges, _instr = synthetic_bundle
    file_results = [
        analyse_file(p, edges) for p in sorted(wav_dir.glob("*.wav"))
    ]
    meas_med, _ = aggregate_weights(file_results, "shimmer")
    # Rank correlation over bands with meaningful planted mass.
    mask = target > (target.max() * 0.05)
    rho, _ = spearmanr(meas_med[mask], target[mask])
    assert rho > 0.9, f"expected Spearman > 0.9 on planted bands, got {rho}"

    summary = run_validation(
        wav_dir,
        instrument_name="cymbal_46cm_medium",
        out_dir=wav_dir.parent / "report",
        mc_draws=40,
        mc_seed=20260803,
    )
    assert (wav_dir.parent / "report" / "validation_report.md").is_file()
    assert summary["n_files"] == 3


def test_analyse_file_phase_windows(synthetic_bundle) -> None:
    wav_dir, _target, edges, _instr = synthetic_bundle
    wav = next(wav_dir.glob("*.wav"))
    res = analyse_file(wav, edges)
    assert "shimmer" in res["phases"]
    w = res["phases"]["shimmer"]["weights"]
    assert w.sum() == pytest.approx(1.0, abs=1e-6)
    assert "9.0" in res["peak_counts"] and "15.0" in res["peak_counts"]


def test_no_writeback_to_source_constants(synthetic_bundle, tmp_path: Path) -> None:
    """Validation must not alter data/source_constants.csv."""
    wav_dir, *_ = synthetic_bundle
    src = Path(__file__).resolve().parents[1] / "data" / "source_constants.csv"
    before = src.read_bytes()
    run_validation(wav_dir, out_dir=tmp_path / "rep", mc_draws=30)
    after = src.read_bytes()
    assert before == after
