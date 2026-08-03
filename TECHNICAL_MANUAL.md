# NonTunPerc — Technical manual

Detailed description of architecture, physics, numerical methods, data
provenance, and interfaces. For day-to-day GUI use see
`QUICK_REFERENCE.md`.

---

## 1. Purpose and scope

**NonTunPerc** (non-tuned percussion) generates theoretical
**ERB-band spectral density profiles** for unpitched percussion
(cymbals, gongs, tam-tams, bass drums) from physical parameters, without
requiring audio analysis.

It is designed as an autonomous companion to empirical spectral-density
pipelines (e.g. *Spectral_Analyser*): schema-aligned metadata out, no
shared code, no write-back into empirical corpora.

**In scope**
- Modal occupation of the auditory scale (modes per 1-ERB band)
- Relative band-energy envelopes by temporal phase
- Optional absolute SPL anchors from Sivian–Dunn–White (1931) / Meyer (2009)
- Monte Carlo specimen/parameter uncertainty
- Optional offline validation against WAV strokes

**Out of scope**
- Pitched-instrument synthesis (bridge fixtures only)
- Chaotic / nonlinear ff broadband regime
- Dome curvature; air loading / two-head coupling except via the optional
  FR Table 18.5 measured-mode anchor at low order (flagged as limits)

---

## 2. Software architecture

```text
nontunperc.py          GUI + pipeline orchestration (entry point)
demo.py                Thin headless wrapper → run_pipeline()
model.py               Physics: plates, membranes, ERB, AmplitudeLayer
uncertainty.py         Monte Carlo layer
calibration.py         Scale-commensurability bridge
sample_metadata.py     Filename/folder parse for validation auto-group
validate_against_recordings.py   Standalone recording check
data/source_constants.csv        Digitized literature constants
pyproject.toml / requirements.txt   Packaging (nontunperc 0.3.4)
tests/                           pytest (VAL1/VAL2/VAL3, excitation, MC, WAV)
```

### 2.1 Public API (stable)

| Symbol | Module | Role |
|---|---|---|
| `PlateInstrument` / `MembraneInstrument` | `model` | Instrument dataclasses (membranes: optional `measured_modes`) |
| `measured_modes_from_csv` / `make_bassdrum_catalogue` | `model` | FR Ch. 18 loader + 82 cm / scaled siblings |
| `generate_profile(instr, …, stroke=, dynamic=)` | `model` | Deterministic `DensityProfile` (+ optional excitation filter) |
| `map_stroke_to_implement` / `load_contact_time_s` | `model` | Stroke→implement + CSV contact times |
| `AmplitudeLayer` | `model` | Measured/mixed band weights + SPL (with fill refusal) |
| `FILL_FRAC_*` / `ERB_MEASURED_ENERGY_FRAC` | `model` | Provenance & calibration thresholds (`internal_default`) |
| `run_monte_carlo(instr, …, stroke=, dynamic=)` | `uncertainty` | Distributional profiles |
| `run_pipeline(options, log)` | `nontunperc` | Full staged run |
| `launch_gui()` | `nontunperc` | Desktop UI |
| `parse_sample_path` / `resolve_model` | `sample_metadata` | Metadata-only validation grouping |

Do not silently alter constants labelled `primary_source` in
`data/source_constants.csv`.

### 2.2 Provenance classes

Every numeric claim carries one of:

| Class | Definition |
|---|---|
| `primary_source` | Copied from a cited table/prose value |
| `derived` | Computed from primary values via stated theory |
| `literature_derived` | Fitted, estimated, or read from a published figure |
| `mixed_primary_and_fill` | AmplitudeLayer initial weights: measured bands + residual fill (10–60%) |
| `internal_default` | Tunable engineering choice, documented |

---

## 3. Auditory grid (ERB)

Glasberg & Moore (1990):

```text
ERB(f) = 24.7 · (4.37·f/1000 + 1)          [Hz]
E(f)   = 21.4 · log10(4.37·f/1000 + 1)     [Cam]
```

Band edges: consecutive 1-ERB steps from 20 Hz to 16 kHz
(`erb_band_edges`). Band centres use geometric means
`√(f_lo·f_hi)`.

---

## 4. Plate model (cymbal / gong / tam-tam)

### 4.1 Low-order modes — Chladni law

```text
f_mn = c · (m + b·n)^p
```

- Anchors `(c, b, p)` from Rossing (2000) Table 9.1 / Fig. 9.3 when
  supplied (`literature_derived` / table `primary_source`).
- Else classical scaling from the 46 cm reference:
  `f ∝ h / d²` (Fletcher & Rossing).

Fig. 9.3 (46 cm) uses `b = 3` and six family fits `f0…f5`. Table 9.1
uses `b = 2` form for the `n = 0` piecewise `(p1,c1)/(p2,c2)` fits.

### 4.2 Asymptotic modal density

Cremer, Heckl & Petersson:

```text
n(f) = √3 · A / (c_L · h)     [modes/Hz, frequency-independent]
c_L  = √( E / (ρ·(1−ν²)) )
```

Above the highest explicit mode, ERB occupation is
`n(f)·Δf_ERB`. Because ERB width grows with `f` while `n(f)` is flat,
**modes per ERB grow with frequency**.

### 4.3 Decay

```text
τ(f) = τ_1k · (f / 1 kHz)^(−α)
```

Default `τ_1k ≈ 10 s`, `α ≈ 0.84` (`literature_derived` from Rossing
Fig. 9.5).

### 4.4 Phase energy

1. Initial weights `e0`: equipartition over modes
   (`internal_default`), **or** AmplitudeLayer mapped weights when
   coverage is **accepted** (fill_fraction ≤ 0.60).
2. Optional Hertzian contact-time filter (when `stroke` and `dynamic`
   both set): `E0(f) ← E0(f) / (1+(f/f_c)^4)`, `f_c = 1/(2 t_contact)`.
   `t_contact` from `excitation_contact_time` CSV (source-read wins) else
   placeholders (`internal_default`). Dynamic scale pp/mf/ff =
   1.6/1.0/0.6 (`internal_default`; Hertz `t ∝ v^(-1/5)` is weaker —
   factors are conservative).
3. Time evolution: `e(t) ∝ e0 · exp(−6.91·t / τ(f))` (60 dB convention).
4. Plate buildup/shimmer (cymbal) or shimmer only (tam-tam): Gaussian
   boost around 4 kHz (Rossing Fig. 9.6), amplitude-gated by dynamic
   (pp 0, mf 0.5, ff 1.0; `internal_default`). With stroke/dynamic
   unset, gate = 1 (bit-identical to v0.3.3).

Weights are renormalized per phase so they sum to 1.

### 4.5 Tam-tam temporal template

`PLATE_PHASES_TAMTAM` (`literature_derived`): strike 0–50 ms, bloom
50 ms–1.5 s (slow HF cascade), shimmer 1.5–6 s, residue 6–20 s.
Selected when `PlateInstrument.plate_class == "tamtam"` (wind gongs
included with a note).

### 4.6 Phase windows (cymbal class)

| Phase | t₀–t₁ |
|---|---|
| strike | 0–0.020 s |
| buildup | 0.020–0.150 s |
| shimmer | 0.150–2.000 s |
| residue | 2.000–6.000 s |

### 4.7 Composite index

Scalar summary (not the primary datum):

```text
H = −Σ w_i log w_i          (Shannon entropy of band weights)
D = [ exp(H) · ⟨modes/band⟩_w ]^(1/2)
```

Report **per-band vectors** in publications; quote `D` only as a
summary.

---

## 5. Membrane model (bass drum)

```text
f_mk = β_mk · c / (2π a),   β_mk = k-th zero of J_m
c    = √(T/σ)   or inferred from nominal f_11
n(f) = 2π A f / c²          (rises linearly with f)
τ(f) = τ_100 · (f/100)^(−α)  (order-of-magnitude, internal_default)
```

Phases: strike (0–50 ms), decay (50–1000 ms).

### 5.1 Measured-mode anchor and effective wave speed

Optional `MembraneInstrument.measured_modes` (Hz) from Fletcher & Rossing
(1998) Table 18.5 (`fr_ch18_bassdrum_mode` rows; default specimen
`both_heads`). Loader: `measured_modes_from_csv` (skips
`needs_manual_reading=1`).

When present, `low_modes()` returns the measured frequencies for the
anchored range, then appends theoretical Bessel modes **only above**
the highest measured frequency. The wave speed used for that
continuation is a least-squares fit of `f = β c / (2π a)` to the
measured set:

```text
c_eff = 2π a · (Σ β_i f_i) / (Σ β_i²)     [derived]
```

This absorbs air-loading and two-head coupling into the low modes
(analogous to the Chladni curvature absorption for plates) and only at
low order. Provenance of `c_eff`: `derived`.

**Size scaling.** Table 18.5 is an 82 cm drum. That size carries the
frequency anchor. Other catalogue diameters reuse the same fitted
`c_eff` with their own radius (`f ∝ 1/a` at fixed effective `c`) —
`internal_default`. Measured Hz are never copied across sizes.

Without usable measured rows, behaviour falls back to the pre-anchor
`f11_nominal` catalogue (bit-identical to v0.3.1).

---

## 6. Absolute amplitude layer

`AmplitudeLayer` loads `data/source_constants.csv`.

### 6.1 Coverage and fill refusal

Prototype names map to Sivian specimen keys via `_INSTRUMENT_SOURCE_KEYS`
(e.g. suspended cymbals → `cymbals_15in` clash PAIR — different mechanical
system/stroke, `internal_default` alias). Gong/tam-tam have **no** key →
bit-identical equipartition.

`historical_band_powers` returns `(lo, hi, energy, is_measured,
fill_fraction, whole, distance)`. Digitized `peak_power_band` rows set
`is_measured=True`. Residual `whole − Σ measured` is filled with **uniform
spectral density** (energy ∝ bandwidth) over uncovered bands
(`internal_default`).

| fill_fraction | Provenance / action |
|---|---|
| ≤ 0.10 (`FILL_FRAC_PRIMARY_MAX`) | `primary_source` |
| ≤ 0.60 (`FILL_FRAC_MIXED_MAX`) | `mixed_primary_and_fill` |
| > 0.60 | refuse → equipartition; profile notes state fill_fraction |

`has_coverage` and `erb_weights_and_spl` share these thresholds.
`DensityProfile.fill_fraction` is exported in `to_rows()`.

### 6.2 ERB mapping (energy-preserving)

Historical band energy `E_i` on `[f_lo, f_hi]` → uniform density
`ρ_i = E_i / Δf_i`. Each ERB band receives `∫ ρ(f) df` over the
overlap. Band edges are never resampled before this step. Remap itself
is unchanged by the fill-policy fixes; only residual placement and
labelling changed.

### 6.3 High-frequency policy

Above ~5 kHz, Meyer is preferred where Sivian is weak; corrections live
in `meyer_hf_discrepancy` rows labelled `literature_derived` (conservative
offsets, not figure readings). Reference distance for Sivian levels is
typically **3 ft (0.9144 m)**.

### 6.4 Power → SPL

Hemispherical spreading (`internal_default`):

```text
I = P / (2π r²),   p_rms = √(I · ρc),   L_p = 20 log10(p_rms / 20 µPa)
```

with `ρc ≈ 413 Pa·s/m`. Sivian “bars” (barye): `1 bar = 0.1 Pa`.

---

## 7. Monte Carlo uncertainty

`uncertainty.run_monte_carlo(instrument, n_draws=2000, seed=…)`

### 7.1 Perturbations (per draw)

| Parameter | Law | Width | Provenance |
|---|---|---|---|
| Thickness | lognormal | 95% multiplicative interval `[1/(1+f), 1+f]` (f=0.25 plates / 0.10 membranes) | `internal_default` (taper/hammering vs film) |
| Diameter | normal | σ = 1% of nominal | `internal_default` |
| E, ρ | normal | σ = 5% | `internal_default` (alloy/film) |
| Chladni *p* | uniform on Table 9.1 p1–p2 class span | class span | span: `primary_source`; 46 cm → 18″ medium nearest class (`internal_default`) |
| Decay τ, α | normal | σ = 20% | `literature_derived` fit uncertainty |

Temporary material keys are registered in `MATERIALS` during a draw and
removed in a `finally` block. **Not thread-safe** (documented limitation).

### 7.2 Aggregation

For each band and phase: percentiles **p5, p25, p50, p75, p95** of
`modes_per_band`, `energy_w_<phase>`, and composite indices. Nesting
`p5 ≤ p25 ≤ p50 ≤ p75 ≤ p95` is tested.

### 7.3 Reporting rule

Deterministic `generate_profile` point estimates are **deprecated for
citation**. Cite MC **p50** with intervals. Seed is stored in
`density_profiles_mc.meta.json`.

Default seed: `20260803`.

---

## 8. Calibration bridge

`calibration.py` builds quasi-harmonic fixtures (partials at `n·f0`) for
trumpet, clarinet, flute, and bass viol (string stand-in; violin
full-band spectrum absent from Sivian 1931).

**Model index (theory only):** partial histogram with
equal-energy-per-partial weighting (`internal_default`). Must **not**
reuse AmplitudeLayer measured band weights (avoids structural
self-comparison).

**Empirical index:** measured ERB bands only — an ERB band counts as
measured if >50% of its remapped energy came from `is_measured`
historical bands (`ERB_MEASURED_ENERGY_FRAC`, `internal_default`).
Instruments with fewer than 2 digitized measured bands, or refused
AmplitudeLayer coverage, are **excluded** and listed in
`calibration_report.md` with fill_fraction and n_measured (once each).

If fewer than **two** survivors remain, return `(nan, nan)` and write /
print `NO CALIBRATION ACHIEVED - factor undefined until the
needs_manual_reading Sivian histograms are completed`.

When ≥2 survivors exist, **the sample standard deviation of the
conversion factor IS the uncertainty** for any cross-domain ratio claim.
Treat the factor as provisional until `needs_manual_reading` cells are
completed (`data/README.md`).

Bridge `f0` choices are `internal_default` tessitura anchors.

---

## 9. Empirical validation (recordings)

`validate_against_recordings.py` is self-contained (numpy / scipy /
soundfile + local `model` / `uncertainty` / `sample_metadata` only).

### 9.1 CLI / GUI

```text
python validate_against_recordings.py --gui
python validate_against_recordings.py --cli --auto-group --wav-dir <folder>
python validate_against_recordings.py --cli --no-auto-group --wav-dir <folder> \
    --instrument cymbal_46cm_medium --out <report_dir>
```

**Auto-group** (metadata only): group key = `(instrument, stroke, dynamic)`
from filenames/folders. Never estimate physical parameters from audio
(circularity refusal). Exact inch → catalogue; else parametric diameter +
`internal_default` thickness. Skip pitched Thai gongs; skip unparseable
names; chinese/splash/ride band-profile only (out of aggregate); `ff`
labelled nonlinear-regime probe.

### 9.2 Per-file processing

1. Mono-mix; resample to 44.1 kHz if needed  
2. Onset: first frame with energy > 10× noise-floor median  
3. Phase cut using model plate windows (truncated to file length)  
4. Welch PSD → integrate on model ERB edges → relative weights  
   - Long segments (≳0.5 s): ~100 ms / 50% overlap  
   - Buildup: ~40 ms  
   - Strike: `n//2`  
5. Shimmer peak counts (Hann, zero-padded FFT) at prominence 9 / 12 / 15 dB  
6. Octave-band 60 dB decay proxy via linear fit to dB envelope  

### 9.3 Comparison metrics

- Spearman ρ (measured vs model MC median, shimmer)  
- Fraction of bands inside model 90% interval  
- Log-ratio bias; separate low-f (<500 Hz) and >10 kHz summaries  

### 9.4 Hard rule

Measured values **never** write back into `data/source_constants.csv` or
any `primary_source` field. Findings go to `validation_report.md` (and
figures) only.

---

## 10. Analytic validation (built-in)

| Check | Criterion |
|---|---|
| VAL1 | 46 cm bronze plate: `n(f)·1000 ≈ 64` modes/kHz; `n·2000 > 100` modes below 2 kHz (Wilbur holography in Rossing §9.2) |
| VAL2 | Anchored `n=0` family reproduces Fig. 9.3 to ~0% for `m = 2…7` |
| VAL3 | `bassdrum_82cm` anchored modes == FR Table 18.5; in-vacuo (c from m≥2) > measured for (0,1) and (1,1); % deviations printed |

Implemented in `nontunperc` pipeline, `demo` wrapper, and
`tests/test_model.py`.

---

## 11. Materials and defaults

| Key | E [Pa] | ρ [kg/m³] | ν |
|---|---:|---:|---:|
| `bronze_B20` | 110e9 | 8700 | 0.33 |
| `brass` | 100e9 | 8500 | 0.35 |
| `steel` | 200e9 | 7850 | 0.30 |

Mylar membrane: `ρ = 1390 kg/m³`, default thickness 190 µm
(`internal_default`).

---

## 12. File formats

### 12.1 `density_profiles.csv`

Minimum columns: `instrument`, `family`, `band_index`, `f_lo_hz`,
`f_hi_hz`, `f_centre_hz`, `modes_per_band`, `energy_w_<phase>`.  
Also: `energy_provenance`, `ref_distance_m`, `fill_fraction` (empty when
N/A), `notes` (profile notes joined with ` | `, including measured-mode
anchor / fitted `c` when active), optional `spl_db_<phase>`.

### 12.2 `density_profiles_mc.csv`

Same base columns plus `modes_per_band_p05…p95`,
`energy_w_<phase>_p05…p95`, `mc_seed`, `mc_n_draws`.  
Plain `modes_per_band` / `energy_w_*` duplicate the **p50** values for
schema familiarity.

### 12.3 `data/source_constants.csv`

Long table: Rossing Table 9.1 & Fig. 9.3, Sivian band edges / powers /
pressures, Meyer levels, HF discrepancies. Extraction notes and
`needs_manual_reading` cells: `data/README.md`.

---

## 13. CLI reference

```text
pip install -e ".[dev]"           # or pip install -r requirements.txt
python nontunperc.py              # GUI
python nontunperc.py --cli        # full headless pipeline
python nontunperc.py --cli --no-mc
python nontunperc.py --cli --no-plots --no-calibration
python nontunperc.py --cli --mc-draws 500 --seed 20260803 --out <dir>

python validate_against_recordings.py --gui
# CLI defaults to metadata auto-group; --auto-group is optional
python validate_against_recordings.py --cli --wav-dir <dir>
python validate_against_recordings.py --cli --no-auto-group --wav-dir <dir> \
    --instrument cymbal_46cm_medium --out <report_dir>
# or: run_validate.bat --cli "D:\Samples"

python -m pytest tests/ -q
```

---

## 14. Dependencies

Declared in `pyproject.toml` (`nontunperc` 0.3.4, `requires-python >=3.10`)
and mirrored in `requirements.txt`:

Python ≥ 3.10 · numpy · scipy · matplotlib · pandas · soundfile  
Optional `[dev]`: pytest · GUI: tkinter (stdlib)

---

## 15. Validity limits (must travel with exports)

1. **Flat-plate approximation** — dome/bow curvature not modelled.  
2. **Linear regime** — chaotic ff broadband outside scope.  
3. **Membranes in vacuo** — with measured anchor, low modes are
   `primary_source` (FR 1998); in-vacuo bias applies only above the
   anchored range. Without anchor, air loading / two-head coupling omitted.  
4. **Excitation** — stroke/dynamic ARE parameters via the contact-time
   filter when both are set; striking position remains unmodelled.
   Unset stroke/dynamic → v0.3.3 equipartition path. AmplitudeLayer
   still requires fill ≤ 0.60 to accept coverage.  
5. **Scale commensurability** — use calibration factor when >=2 bridge survivors exist; otherwise NO CALIBRATION ACHIEVED; spread = uncertainty; factor provisional under sparse digitization.
6. **1931 HF chain** — above ~5 kHz prefer Meyer; `discrepancy_db` (`literature_derived`) records corrections.  
7. **Specimen anchoring** — absolute levels = single specimens; variance via MC.  
8. **Clash-pair alias** — suspended-cymbal model names → Sivian 15-in. clash PAIR (`internal_default`).

---

## 16. References

- Cremer, Heckl & Petersson (2005). *Structure-borne sound* (3rd ed.).  
- Fletcher & Rossing (1998). *The physics of musical instruments* (2nd ed.).  
- Glasberg & Moore (1990). *Hearing Research*, 47, 103–138.  
- Meyer (2009). *Acoustics and the performance of music* (5th ed.).  
- Rossing (2000). *Science of percussion instruments*.  
- Sivian, Dunn & White (1931). *JASA*, 2, 330–371.

---

## 17. Change control

All behavioural and constant changes are logged in `CHANGES.md`.
Parameter edits that touch `primary_source` rows require explicit human
approval; validation findings alone do not alter those rows.
