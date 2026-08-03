# CHANGES

## v0.3 — Monte Carlo uncertainty + recording validation

- Added `QUICK_REFERENCE.md` (GUI/options/outputs) and `TECHNICAL_MANUAL.md` (architecture, equations, MC, validation, validity).
- Cleared regenerable outputs (`density_profiles*.csv/png`, MC fan, `size_sweep.png`, `calibration_report.md`) and `__pycache__` / `.pytest_cache`; re-create via NonTunPerc. Kept PDFs, `data/`, and source code.
- Renamed user-facing app to **NonTunPerc** (`nontunperc.py`): tkinter GUI (stages, instruments, MC draws/seed, output folder) plus `--cli` headless pipeline; `demo.py` kept as thin compatibility wrapper.
- Added `run_nontunperc.bat` (GUI by default; `--cli` for headless); `run_demo.bat` now forwards to it.
- Added validation GUI (`validate_against_recordings.py --gui`): browse for WAV sample folder, instrument, report dir; `run_validate.bat` opens GUI by default (`--cli` keeps headless path).
- Added `run_validate.bat` Windows launcher (WAV validation with optional folder/instrument/out args, default `wavs\`).
- Added `uncertainty.py` with `run_monte_carlo` (default N=2000, fixed seed) aggregating p5/p25/p50/p75/p95 for modes, energy weights, and composite indices.
- Thickness lognormal sigmas set so the 95% span is ±25% (plates) / ±10% (membranes) (`internal_default`: hammering/taper vs membrane film tolerance).
- Diameter / E,ρ / decay widths taken as 1-σ = 1% / 5% / 20% of nominal (`internal_default` reading of “normal, ±X%”).
- Chladni *p* drawn uniform on Rossing Table 9.1 p1–p2 span of the instrument class (`primary_source` span); 46 cm mapped to 18-in medium as nearest class (`internal_default`).
- Added optional instance decay / `membrane_rho` fields on plate/membrane dataclasses so MC can perturb without changing `generate_profile` call sites; default values keep deterministic numerics.
- Exported `density_profiles_mc.csv` (+ seed sidecar) and 46-cm shimmer fan chart from `demo.py`.
- README: deterministic point estimates deprecated for reporting; MC median + intervals are the citable output.
- Added `validate_against_recordings.py` (numpy/scipy/soundfile + local model only): onset, phase windows, Welch ERB weights, peak-count prominence sweep, decay proxy, MC comparison report.
- Validation writes only to `validation_report.md` / figures — never to `data/source_constants.csv` or any `primary_source` field.
- Tests: MC reproducibility (6 d.p.), nested non-degenerate percentiles; synthetic-WAV Spearman > 0.9 and no source-constants write-back.
- Welch windows: ~100 ms (long shimmer/residue), ~40 ms (buildup), short `n//2` (strike); 50% overlap; peak prominence default 12 dB (`internal_default`).
- WAV discovery de-duplicates case-insensitive `*.wav`/`*.WAV` matches (Windows).

## v0.2 — Audit / absolute amplitude / calibration

Audit / extension log (prototype v0.1 → v0.2). One-line rationale each.

- Restored canonical filenames `model.py` / `demo.py` / `README.md` from `Model.py` / `Demo.py` / `README_1.md` so imports and the build-spec names match.
- Noted absence of shipping `cursor_prompt_idiophone_density.md`; reconstructed a normative checklist from the shipping README so Task-1 audits have a local reference (`cursor_prompt_idiophone_density.md`).
- Added `tests/test_model.py` implementing VAL1, VAL2, ERB sanity, phase normalization, coverage bit-identity, remap conservation, and calibration checks (spec required pytest suite was missing).
- Digitized Rossing Table 9.1 and Fig. 9.3 fits into `data/source_constants.csv` (primary_source / literature_derived) without altering in-code Chladni anchors already matching those rows.
- Digitized Sivian–Dunn–White (1931) band edges, peak-power anchors, and average total pressures; left full figure-curve ordinates blank with `needs_manual_reading=1` rather than guessing.
- Digitized Meyer (2009) SPL / sound-power corroboration and recorded HF `discrepancy_db` where Meyer is preferred above ~5 kHz.
- Added `AmplitudeLayer` to `model.py`: loads CSV, energy-preserving ERB remap, absolute SPL at source distances; equipartition path unchanged when coverage is absent (bit-identical relative weights).
- Extended `DensityProfile` / CSV export with `energy_provenance`, `ref_distance_m`, and optional `spl_db_*` columns (additive schema; prior columns preserved).
- Added `calibration.py` + `calibration_report.md` writer: quasi-harmonic bridge on trumpet / clarinet / flute / bass_viol; spread reported as the cross-domain uncertainty.
- Updated `demo.py` to apply AmplitudeLayer and print the calibration factor at the end.
- Updated `README.md` with AmplitudeLayer, calibration protocol, and validity limits (6) and (7); kept limits (1)–(5).
- Chose hemispherical spreading `I = P/(2πr²)` for power→SPL conversion as `internal_default` matching Sivian near-source estimates.
- Chose bridge `f0` values (trumpet 233 Hz, clarinet 147 Hz, flute 349 Hz, bass_viol 49 Hz) as `internal_default` orchestral-tessitura anchors.
- Aliased prototype cymbal/bass-drum names to Sivian specimen keys; deliberately did **not** alias gong/tam-tam (no coverage → bit-identical equipartition).
- Residual whole-spectrum peak power after named bands is spread flat across uncovered historical bands (`internal_default` fill) so total power is conserved when only dominant bands are textual.
