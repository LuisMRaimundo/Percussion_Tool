# NonTunPerc — Quick reference

One-page guide to the GUI options, outputs, and everyday use.
Version **0.3.4**.

---

## Start

| Action | How |
|---|---|
| Project folder | `…\Desktop\Percussion Tool` (not `Tôol`) |
| Install deps | `pip install -e ".[dev]"` or `pip install -r requirements.txt` |
| Open the app | Double-click `run_nontunperc.bat` |
| Headless full run | `run_nontunperc.bat --cli` or `python nontunperc.py --cli` |
| Validate samples | `run_validate.bat` → GUI (**auto-group** from filenames; subfolders on); or `run_validate.bat --cli "D:\Samples"`; manual: `run_validate.bat --cli "D:\folder" cymbal_18in_medium` |

**Citable result:** Monte Carlo **median** (`p50`) with percentile bands — not the deterministic single-run numbers alone.

---

## What the model does (in one paragraph)

From physical parameters (size, thickness/tension, material), NonTunPerc predicts how densely each auditory ERB band is occupied by vibrational modes, and how relative energy is distributed across time phases after a stroke. It does **not** analyse audio unless you run the separate validation tool.

---

## GUI — Pipeline stages

| Option | Meaning |
|---|---|
| **VAL1 / VAL2 / VAL3 checks** | VAL1/VAL2: 46 cm cymbal. VAL3: `bassdrum_82cm` vs FR Table 18.5. CLI also prints `[EXC ]` excitation state per instrument. |
| **Density profiles (CSV)** | Writes `density_profiles.csv`: one row per instrument × ERB band with mode counts, relative energy weights, `energy_provenance`, `fill_fraction`. |
| **Plots (PNG)** | Writes `density_profiles.png` and MC `size_sweep.png` / `size_sweep_mc.csv` (fan along diameter; deterministic sweep deprecated). |
| **Calibration bridge** | Theory-side partials vs measured-bands-only empirical index → `calibration_report.md`. Needs >=2 survivors; otherwise **NO CALIBRATION ACHIEVED**. |
| **Monte Carlo uncertainty** | Re-runs the model with specimen/parameter noise; exports medians and intervals (`density_profiles_mc.csv` + fan chart). |
| **Use AmplitudeLayer** | When on, instruments with Sivian–Meyer coverage **and** fill_fraction ≤ 0.60 use measured/mixed weights; mostly-filled vectors (e.g. cymbals ≈ 0.90) refuse and keep equipartition. Gong/tam-tam stay equipartition. |

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
| `cymbal_16in_thin` | plate | Table 9.1 Chladni anchor; AmplitudeLayer currently **refuses** (high fill) |
| `cymbal_18in_medium` | plate | Table 9.1 Chladni anchor; AmplitudeLayer currently **refuses** (high fill) |
| `cymbal_46cm_medium` | plate | Fig. 9.3 anchor; main validation / fan-chart target; AmplitudeLayer **refuses** |
| `gong_50cm_bronze` | plate | Scaled (no Table 9.1 row); no Sivian alias |
| `tamtam_80cm_bronze` | plate | Scaled; no Sivian alias |
| `bassdrum_82cm` | membrane | FR Table 18.5 mode anchor; Sivian alias → `bass_drum_A_36x15` (size-mismatch `internal_default`) |
| `bassdrum_32in` | membrane | Fitted `c` from 82 cm + own diameter; may accept mixed fill ≈ 0.59 |
| `bassdrum_28in` | membrane | Fitted `c` from 82 cm + own diameter; currently refuses (fill ≈ 0.87) |

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
| `density_profiles.csv` | Deterministic per-band profile (modes + `energy_w_*` + `fill_fraction` + excitation notes) |
| `density_profiles_mc.csv` | Same schema + `_p05`…`_p95` columns; use **`_p50`** to cite |
| `density_profiles_mc.meta.json` | Seed, draw counts, perturbation metadata |
| `density_profiles_mc_fan_*.png` | Median line + 50% / 90% bands (fan chart) |
| `size_sweep_mc.csv` / `size_sweep.png` | MC composite index vs diameter (cite median + bands) |
| `calibration_report.md` | Factor + spread if >=2 survivors; else **NO CALIBRATION ACHIEVED**; fill_fraction / exclusions |
| `validation_report.md` | Auto-group: stick/mallet cohorts, optional baseline ρ deltas |

---

## Time phases (plates)

**Cymbal class**

| Phase | Window | What it represents |
|---|---|---|
| strike | 0–20 ms | Immediate strike sound |
| buildup | 20–150 ms | Early spectral build |
| shimmer | 0.15–2 s | 3–5 kHz aftersound (“shimmer”) |
| residue | 2–6 s | Late decay |

**Tam-tam class** (`plate_class="tamtam"`; wind gongs share this)

| Phase | Window | What it represents |
|---|---|---|
| strike | 0–50 ms | Impact transient |
| bloom | 50 ms–1.5 s | Slow low→high cascade |
| shimmer | 1.5–6 s | Established HF shimmer |
| residue | 6–20 s | Long decay |

Membranes use only **strike** / **decay**. When validation supplies
`stroke` + `dynamic`, a Hertzian contact-time filter shapes `E0`
(placeholders listed in `CHANGES.md` until source-read).

---

## How to read a result quickly

1. Open `density_profiles_mc.csv`.
2. Filter to your instrument.
3. Cite `modes_per_band_p50` and `energy_w_shimmer_p50` (or the phase you need).
4. Report the interval: `p05`–`p95` (90% band) or `p25`–`p75` (50% band).
5. Check `energy_provenance` / `fill_fraction` on the deterministic CSV before treating absolute SPL columns as measurements.
6. If comparing to pitched-instrument metadata, use the calibration factor only when the report defines one (>=2 survivors); otherwise there is **no** conversion factor yet. Attach the **spread** as uncertainty when a factor exists (provisional under sparse digitization).

---

## Provenance labels (short)

| Label | Meaning |
|---|---|
| `primary_source` | Taken from a cited table/text (or AmplitudeLayer with ≤10% fill) |
| `mixed_primary_and_fill` | AmplitudeLayer: measured bands + residual fill (10–60%) |
| `derived` | Computed from primary values by stated theory |
| `literature_derived` | Fit / estimate / read from a published figure |
| `internal_default` | Engineering choice; documented, not a measurement |

---

## Validity (do not ignore)

1. Flat-plate / single-membrane idealisations  
2. Linear regime only (not chaotic ff crashes)  
3. Membrane in-vacuo bias: with FR Table 18.5 anchor, only **above** the
   measured range; without anchor, lowest modes overestimated  
4. Stroke/dynamic ARE excitation parameters (contact-time filter); striking
   position unmodelled; AmplitudeLayer still needs fill ≤ 0.60 to accept  
5. Absolute 1931 levels weak above ~5 kHz (Meyer corrections = `literature_derived`)  
6. Absolute levels = single specimens, not variance (variance → MC layer)  
7. Cymbal model names → Sivian clash PAIR alias (`internal_default`)  
8. Calibration factor provisional until `needs_manual_reading` cells are filled  

Full equations and file map: see **TECHNICAL_MANUAL.md**.
