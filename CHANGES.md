# CHANGES

## v0.3.5 — Stick contact times, ff filter bypass, t_contact in MC

- **Stick contact-time revision** (`internal_default` placeholders; FR /
  Rossing 2000 discuss stick–cymbal interaction qualitatively but do not
  print unambiguous absolute `t_contact` for these implement classes —
  no `primary_source` overwrite). New placeholders:

  | implement | t_contact (mf) | notes |
  |---|---:|---|
  | `stick_tip` | 0.15 ms | was 0.4 ms |
  | `stick_shoulder` | 0.3 ms | was 0.7 ms |
  | `stick_bell` | 0.10 ms | new; near-impulsive dome/bell |
  | `yarn_mallet` | 3.0 ms | **UNCHANGED** — validated by mallet-cohort improvement in the Iowa re-gate |
  | `bass_drum_beater` | 6.0 ms | **UNCHANGED** — same rationale |

  `stick.bell` now maps to `stick_bell` (was aliased to `stick_tip`).
- **ff plate bypass:** at `dynamic="ff"` the excitation low-pass is
  disabled for plate instruments (nonlinear cascade regenerates HF
  independently of excitation content, Rossing 2000 §9.4); the 3–5 kHz
  boost remains at full amplitude. pp/mf behaviour unchanged. Membranes
  keep the low-pass at ff.
- **MC:** when `stroke` is specified, `t_contact` is perturbed lognormal
  with 95% span ±50% (`internal_default`); recorded in MC meta as
  `t_contact_sigma_log` / `t_contact_95pct_span`. Expected: wider
  intervals, better inside90 coverage.
- **Validation re-run:** pass `--baseline <prior validation_summary.json>`
  so the report fills previous-ρ vs new-ρ (see README).
- Package version → `0.3.5`. Mallet-path (yarn_mallet @ pp/mf) remains
  bit-identical to v0.3.4. Docs synced: README (stick times, ff bypass,
  `--baseline` re-run); TECHNICAL_MANUAL (§4.4, MC `t_contact`);
  QUICK_REFERENCE; `data/README.md`; `cursor_prompt_idiophone_density.md`.

## v0.3.4 — (stroke, dynamic) excitation filter + tam-tam template

Repairs the falsified universal equipartition / ff-shimmer convention
with bibliography-anchored Hertzian contact-time filtering. Validation
recordings exposed the defect but **did not fit** the repair.

- **Excitation filter:** `E0(f) ∝ 1/(1+(f/f_c)^4)`, `f_c = 1/(2 t_contact)`.
  API: `generate_profile(..., stroke=, dynamic=)`. Both unset →
  bit-identical to v0.3.3.
- **Contact times** (`excitation_contact_time` rows) — all currently
  `internal_default` placeholders awaiting source-read values (FR Ch. 19.7
  / Fig. 19.12 give qualitative / f_max curves, not unambiguous printed
  `t_contact` in ms for these implement classes):

  | implement | t_contact (mf) | notes |
  |---|---:|---|
  | `stick_tip` | 0.4 ms | f_c ≈ 1.25 kHz |
  | `stick_shoulder` | 0.7 ms | |
  | `yarn_mallet` | 3.0 ms | f_c ≈ 167 Hz |
  | `bass_drum_beater` | 6.0 ms | |

  Dynamic scale of `t_contact`: pp×1.6, mf×1.0, ff×0.6 (`internal_default`;
  Hertz predicts only `t ∝ v^(-1/5)` — factors deliberately conservative).
  Shimmer-boost gate: pp 0, mf ×0.5, ff full (`internal_default`).
- **Tam-tam template** `PLATE_PHASES_TAMTAM` + `plate_class="tamtam"`;
  HF boost in shimmer not bloom; wind gongs share the template.
- **Validation:** passes (stroke, dynamic) into MC; stick/mallet cohorts
  separate; optional `--baseline validation_summary.json`.
- **MC size sweep** replaces deterministic `size_sweep.png`
  (`size_sweep_mc.csv` + fan). Deterministic sweep deprecated for citation.
- **Alias:** `bassdrum_82cm` → Sivian `bass_drum_A_36x15` (36 in ≈91 cm vs
  82 cm; size-mismatch `internal_default`).
- Package version → `0.3.4`. Docs synced: README (excitation, tam-tam,
  validity limit 4, MC size sweep); TECHNICAL_MANUAL (§4.4–4.5, API);
  QUICK_REFERENCE (instruments, phases, outputs); `data/README.md`
  (contact times + 82 cm Sivian alias); `cursor_prompt_idiophone_density.md`.

## v0.3.3 — Membrane measured-mode anchor (Fletcher & Rossing Ch. 18)

- Digitized Fletcher & Rossing (1998) Table 18.5 concert bass-drum modes
  into `data/source_constants.csv` as `fr_ch18_bassdrum_mode` /
  `bassdrum_82cm` (`primary_source`). Default specimen `both_heads`
  (carry-head lower): 39, 80, 121, 162, 204, 248 Hz for (01)…(51).
  Equal-tension doublets recorded under `both_heads_equal`;
  `batter_only` left blank with `needs_manual_reading=1`. **No existing
  `primary_source` numeric value was altered.**
- `MembraneInstrument.measured_modes` + `measured_modes_from_csv` /
  `make_bassdrum_catalogue`: 82 cm entry carries the anchor; 32-in /
  28-in reuse fitted effective `c` with own diameter (`f ∝ 1/a`,
  `internal_default`). Fitted `c_eff = 67.036 m/s` (derived) from the
  both-heads set.
- `low_modes()` returns measured frequencies then Bessel theory above
  the highest measured mode (air-loading / two-head effects absorbed at
  low order, Chladni-analogue).
- **VAL3** in CLI / tests: anchored modes reproduce Table 18.5 exactly;
  in-vacuo (`c` from m>=2 = 68.725 m/s) vs measured side-by-side:

  | mode | measured [Hz] | in-vacuo [Hz] | pct_dev |
  |---|---:|---:|---:|
  | (01) | 39 | 64.16 | +64.50% |
  | (11) | 80 | 102.22 | +27.78% |
  | (21) | 121 | 137.01 | +13.23% |
  | (31) | 162 | 170.21 | +5.07% |
  | (41) | 204 | 202.44 | -0.76% |
  | (51) | 248 | 234.00 | -5.64% |

- Graceful fallback: with no usable measured rows, catalogue is
  bit-identical to the v0.3.1 `f11_nominal` path.
- Docs synced to **0.3.3**: README (§2.2, VAL3, validity limit 3, `notes`
  export); TECHNICAL_MANUAL (§5.1, API, packaging); QUICK_REFERENCE
  (instruments, VAL3, validity); `data/README.md` Ch. 18 extraction +
  `needs_manual_reading` batter-only; `cursor_prompt_idiophone_density.md`.
  Package version → `0.3.3`.

## v0.3.2 — Calibration bridge: independent theory side + no-factor verdict

- If fewer than 2 bridge instruments survive exclusion, the conversion
  factor is **undefined** (`nan`, `nan`); report and CLI print
  `NO CALIBRATION ACHIEVED - factor undefined until the
  needs_manual_reading Sivian histograms are completed`.
- Model-side bridge index no longer reuses AmplitudeLayer measured
  weights (structural self-comparison removed). It is built from theory
  only: partial histogram + equal-energy-per-partial weighting
  (`internal_default`).
- Empirical side remains measured-bands-only (unchanged).
- Exclusion lines no longer embed `fill_fraction` twice (reason + trailer).
- Package version → `0.3.2`.
- Docs updated (README, TECHNICAL_MANUAL, QUICK_REFERENCE, data/README,
  cursor_prompt) for the independent theory-side index and no-calibration
  verdict.

## v0.3.1 — Honest AmplitudeLayer labelling, density fill, calibration, packaging

- **Fix 1 — provenance / refusal:** `historical_band_powers` now returns an
  `is_measured` mask; `fill_fraction = residual/whole` drives provenance
  (`primary_source` ≤0.10, `mixed_primary_and_fill` ≤0.60, else refuse and
  keep equipartition). Thresholds are module constants (`internal_default`).
  `DensityProfile.fill_fraction` exported in `to_rows()`; notes state the
  fraction. `has_coverage` uses the same rule. Rationale: a vector that is
  ~90% residual fill must not be stamped `primary_source`.
- **Fix 2 — equal-density fill:** residual power is spread over uncovered
  historical bands ∝ bandwidth (uniform W/Hz), not equal energy per band.
  For `cymbals_15in`, mass below 250 Hz drops from ~2.14 W (22.5% of 9.5 W)
  under equal-energy fill to ~0.25 W (2.6%) under equal-density fill; the
  20–62.5 Hz band receives less energy than 2000–2800 Hz. Remap mechanics
  unchanged.
- **Fix 3 — calibration bridge:** empirical index uses measured ERB bands
  only (>50% energy from measured historical bands); instruments with
  fewer than 2 measured bands are excluded and listed in
  `calibration_report.md`, which also reports fill_fraction / n_measured
  and carries a sparse-coverage caveat.
- **Fix 4 — label corrections (no value changes):** four
  `meyer_hf_discrepancy` rows re-labelled `literature_derived` (were
  mis-tagged `primary_source`); clash-pair alias documented in CSV notes,
  `data/README.md`, and README validity limit (8). **No `primary_source`
  numeric value was altered.**
- **Fix 5 — packaging:** `pyproject.toml` (`nontunperc` 0.3.1,
  `requires-python >=3.10`) + `requirements.txt`; uncertainty docstring
  corrected to `[1/(1+frac), (1+frac)]`; `MATERIALS` mutation noted as
  not thread-safe.
- **Docs sync:** README / TECHNICAL_MANUAL / QUICK_REFERENCE /
  `data/README.md` / `cursor_prompt_idiophone_density.md` updated for
  fill refusal, equal-density fill, measured-bands calibration, packaging,
  and validation auto-group defaults.
- **Local folder rename:** Desktop project directory corrected from
  `Percussion Tôol` → `Percussion Tool` (misspelling; no code-path change).

## v0.3 — Monte Carlo uncertainty + recording validation

- Added `QUICK_REFERENCE.md` (GUI/options/outputs) and `TECHNICAL_MANUAL.md` (architecture, equations, MC, validation, validity).
- Cleared regenerable outputs (`density_profiles*.csv/png`, MC fan, `size_sweep.png`, `calibration_report.md`) and `__pycache__` / `.pytest_cache`; re-create via NonTunPerc. Kept PDFs, `data/`, and source code.
- Renamed user-facing app to **NonTunPerc** (`nontunperc.py`): tkinter GUI (stages, instruments, MC draws/seed, output folder) plus `--cli` headless pipeline; `demo.py` kept as thin compatibility wrapper.
- Added `run_nontunperc.bat` (GUI by default; `--cli` for headless); `run_demo.bat` now forwards to it.
- Added validation GUI (`validate_against_recordings.py --gui`): browse for sample folder, instrument, report dir; `run_validate.bat` opens GUI by default (`--cli` keeps headless path).
- Validation sample discovery walks **subfolders by default** (`.wav` / `.aif` / `.aiff` / `.flac`); skips `__MACOSX` and `._*` junk; GUI checkbox + CLI `--no-recursive`.
- Metadata-only **auto-group** (`sample_metadata.py`, `--auto-group`): group key `(instrument, stroke, dynamic)`; exact inch→catalogue else parametric `filename_metadata + internal_default_thickness`; skip pitched Thai gongs; skip unparseable (no guessing); chinese/splash/ride band-profile only (out of aggregate); `ff` labelled nonlinear-regime probe; circularity refusal stated in report; full file→model mapping table.
- Added `run_validate.bat` Windows launcher (folder-only ⇒ auto-group; optional catalogue instrument ⇒ manual).
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
- Residual whole-spectrum peak power after named bands is spread across uncovered historical bands (`internal_default` fill) so total power is conserved when only dominant bands are textual. *(Superseded in v0.3.1 by equal-density / bandwidth-proportional fill.)*
