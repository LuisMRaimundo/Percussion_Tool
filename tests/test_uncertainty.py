"""Tests for the Monte Carlo uncertainty layer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from model import MembraneInstrument, PlateInstrument
from uncertainty import (
    DEFAULT_SEED,
    assert_nested_percentiles,
    export_mc_csv,
    plot_fan_chart,
    run_monte_carlo,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def cymbal_46() -> PlateInstrument:
    return PlateInstrument(
        "cymbal_46cm_medium",
        0.460,
        0.0012,
        chladni=(14.93, 3.0, 1.557),
    )


def test_mc_reproducibility_medians(cymbal_46: PlateInstrument) -> None:
    """Same seed ⇒ medians agree to 6 decimals."""
    a = run_monte_carlo(cymbal_46, n_draws=40, seed=DEFAULT_SEED)
    b = run_monte_carlo(cymbal_46, n_draws=40, seed=DEFAULT_SEED)
    np.testing.assert_allclose(
        a.modes_quantiles["p50"], b.modes_quantiles["p50"], rtol=0, atol=1e-6
    )
    for ph in a.energy_quantiles:
        np.testing.assert_allclose(
            a.energy_quantiles[ph]["p50"],
            b.energy_quantiles[ph]["p50"],
            rtol=0,
            atol=1e-6,
        )
        assert abs(
            a.composite_quantiles[ph]["p50"]
            - b.composite_quantiles[ph]["p50"]
        ) < 1e-6


def test_mc_nested_nondegenerate(cymbal_46: PlateInstrument) -> None:
    res = run_monte_carlo(cymbal_46, n_draws=60, seed=DEFAULT_SEED)
    assert_nested_percentiles(res)
    # Non-degenerate 90% band somewhere in the spectrum.
    width = res.modes_quantiles["p95"] - res.modes_quantiles["p05"]
    assert np.max(width) > 0.0
    ewidth = (
        res.energy_quantiles["shimmer"]["p95"]
        - res.energy_quantiles["shimmer"]["p05"]
    )
    assert np.max(ewidth) > 0.0


def test_mc_seed_in_metadata(cymbal_46: PlateInstrument) -> None:
    res = run_monte_carlo(cymbal_46, n_draws=20, seed=12345)
    assert res.seed == 12345
    assert res.metadata["seed"] == 12345
    assert res.metadata["chladni_p_span"] is not None


def test_mc_export_and_fan(cymbal_46: PlateInstrument, tmp_path: Path) -> None:
    res = run_monte_carlo(cymbal_46, n_draws=30, seed=DEFAULT_SEED)
    csv_path = export_mc_csv([res], tmp_path / "density_profiles_mc.csv")
    assert csv_path.is_file()
    text = csv_path.read_text(encoding="utf-8")
    assert "modes_per_band_p50" in text
    assert "mc_seed" in text
    fig = plot_fan_chart(
        res, phase="shimmer", path=tmp_path / "fan.png"
    )
    assert fig.is_file()


def test_mc_membrane_runs() -> None:
    drum = MembraneInstrument("bassdrum_32in", 0.813, f11_nominal=60.0)
    res = run_monte_carlo(drum, n_draws=25, seed=7)
    assert "decay" in res.energy_quantiles
    assert_nested_percentiles(res)
