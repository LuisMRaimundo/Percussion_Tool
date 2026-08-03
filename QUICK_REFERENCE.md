# NonTunPerc — Quick reference

One-page guide to the GUI options, outputs, and everyday use.

---

## Start

| Action | How |
|---|---|
| Open the app | Double-click `run_nontunperc.bat` |
| Headless full run | `run_nontunperc.bat --cli` or `python nontunperc.py --cli` |
| Validate samples | `run_validate.bat` → GUI (**subfolders** searched by default; `.wav`/`.aif`/`.flac`); or `run_validate.bat --cli "D:\samples"` |

**Citable result:** Monte Carlo **median** (`p50`) with percentile bands — not the deterministic single-run numbers alone.

---

## What the model does (in one paragraph)

From physical parameters (size, thickness/tension, material), NonTunPerc predicts how densely each auditory ERB band is occupied by vibrational modes, and how relative energy is distributed across time phases after a stroke. It does **not** analyse audio unless you run the separate validation tool.

---

## GUI — Pipeline stages

| Option | Meaning |
|---|---|
| **VAL1 / VAL2 checks** | Sanity checks on the 46 cm cymbal: modal density ≈ 64 modes/kHz (VAL1); Chladni fit reproduces Rossing Fig. 9.3 (VAL2). |
| **Density profiles (CSV)** | Writes `density_profiles.csv`: one row per instrument × ERB band with mode counts and relative energy weights. |
| **Plots (PNG)** | Writes `density_profiles.png` (modes + phase energies) and `size_sweep.png` (index vs cymbal diameter). |
| **Calibration bridge** | Compares model scale to Sivian/Meyer pitched-instrument anchors → `calibration_report.md` (conversion factor + spread). |
| **Monte Carlo uncertainty** | Re-runs the model with specimen/parameter noise; exports medians and intervals (`density_profiles_mc.csv` + fan chart). |
| **Use AmplitudeLayer** | When on, cymbals/bass drums with Sivian–Meyer coverage use measured band weights instead of equipartition. Gong/tam-tam stay on equipartition either way. |

---

## GUI — Monte Carlo controls

| Control | Meaning | Typical |
|---|---|---|
| **Draws (46 cm)** | Number of random specimens for the focus cymbal | 2000 (publication); 200–400 for a quick look |
| **Draws (others)** | Draws for every other selected instrument | 400 |
| **RNG seed** | Fixes the random stream so results are reproducible | `20260803` (default) |

Higher draws → smoother intervals, longer runtime.

---

## GUI — Instruments

| Name | Family | Notes |
|---|---|---|
| `cymbal_16in_thin` | plate | Table 9.1 Chladni anchor |
| `cymbal_18in_medium` | plate | Table 9.1 Chladni anchor |
| `cymbal_46cm_medium` | plate | Fig. 9.3 anchor; main validation / fan-chart target |
| `gong_50cm_bronze` | plate | Scaled (no Table 9.1 row) |
| `tamtam_80cm_bronze` | plate | Scaled |
| `bassdrum_32in` / `28in` | membrane | Tension from nominal (1,1) frequency |

Uncheck instruments you do not need to speed up the run.

---

## GUI — Buttons

| Button | Action |
|---|---|
| **Run pipeline** | Runs the checked stages (log on the right; UI stays responsive). |
| **Open output folder** | Opens the folder where CSVs/PNGs/reports are written. |
| **Validate against WAVs…** | Launches `run_validate.bat` for empirical checks on recordings. |

**Output folder:** where files are saved (defaults to the NonTunPerc program folder).

---

## Key outputs

| File | Contents |
|---|---|
| `density_profiles.csv` | Deterministic per-band profile (modes + `energy_w_*`) |
| `density_profiles_mc.csv` | Same schema + `_p05`…`_p95` columns; use **`_p50`** to cite |
| `density_profiles_mc.meta.json` | Seed, draw counts, perturbation metadata |
| `density_profiles_mc_fan_*.png` | Median line + 50% / 90% bands (fan chart) |
| `calibration_report.md` | Scale factor and **spread = uncertainty** for cross-domain ratios |
| `validation_report.md` | Only after WAV validation — measured vs model comparison |

---

## Time phases (plates)

| Phase | Window | What it represents |
|---|---|---|
| strike | 0–20 ms | Immediate strike sound |
| buildup | 20–150 ms | Early spectral build |
| shimmer | 0.15–2 s | 3–5 kHz aftersound (“shimmer”) |
| residue | 2–6 s | Late decay |

Membranes use only **strike** / **decay**.

---

## How to read a result quickly

1. Open `density_profiles_mc.csv`.
2. Filter to your instrument.
3. Cite `modes_per_band_p50` and `energy_w_shimmer_p50` (or the phase you need).
4. Report the interval: `p05`–`p95` (90% band) or `p25`–`p75` (50% band).
5. If comparing to pitched-instrument metadata, multiply/divide by the calibration factor and attach the **spread** as uncertainty.

---

## Provenance labels (short)

| Label | Meaning |
|---|---|
| `primary_source` | Taken from a cited table/text |
| `derived` | Computed from primary values by stated theory |
| `literature_derived` | Fit / read from a published figure |
| `internal_default` | Engineering choice; documented, not a measurement |

---

## Validity (do not ignore)

1. Flat-plate / single-membrane idealisations  
2. Linear regime only (not chaotic ff crashes)  
3. Equipartition at strike is a convention (unless AmplitudeLayer applies)  
4. Absolute 1931 levels weak above ~5 kHz (Meyer corrections recorded)  
5. Absolute levels = single specimens, not variance (variance → MC layer)

Full equations and file map: see **TECHNICAL_MANUAL.md**.
