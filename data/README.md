# Source-constant extraction notes

Digitized from the three PDFs in the project root. Rows live in
`source_constants.csv`. Every numeric cell carries a provenance class and
unit pair (as printed → SI / dB re 20 µPa or dB re 1 pW).

## Needs manual reading (`needs_manual_reading=1`)

Complete these before treating the corresponding absolute band curves as
primary:

1. **Sivian per-band peak histograms** for bass drum A (Fig. 33), 15″
   cymbals (Fig. 38), trumpet (Fig. 44), clarinet (Fig. 46), flute
   (Fig. 47): ordinates (% of observations / 1 dB zone) per printed band.
2. **Sivian long-average pressure-per-cycle curves** (Figs. 7–32 family)
   for the same instruments: dB re total average pressure, per band.
3. **Violin band spectrum** — not present as a named peak-power instrument
   in Sivian et al. (1931); only the soft average pressure (0.52 bars at
   3 ft) is textual. Full band levels require another source or human
   figure work if a different edition contains them.
4. **Fletcher & Rossing Table 18.5 — batter-only / single-head column:**
   prose (book p. 513) states that removing the carry head changes modal
   frequencies but little from the carry-lower column; **no separate
   tabulated single-head frequencies** are printed. Six
   `fr_ch18_bassdrum_mode` / `specimen=batter_only` rows are blank with
   `needs_manual_reading=1` (do not invent values).

Blank cells for **15 in. thick** cymbal `p2` / `c2` / `m_c` are blank in
the printed Rossing Table 9.1 (not illegible) and are *not* flagged.

## Extraction method

- **Rossing**: scan OCR + direct reading of Table 9.1 and Fig. 9.3 equation
  labels on PDF pages 101–102 (book pp. 91–92).
- **Sivian, Dunn & White (1931)**: text-layer extraction for narrative
  numbers; Fig. 59 OCR cross-checked against the prose; band edges taken
  from figure axis labels (not resampled). Pressure unit: 1 bar (barye) =
  0.1 Pa; SPL ≈ 20·log₁₀(p/20µPa). Reference distance: **3 ft (0.9144 m)**
  unless a row states otherwise (snare drum 4 ft in prose; piano 10 ft).
- **Meyer (2009)**: text-layer SPL / sound-power-level statements for the
  same instrument set. Where Meyer and Sivian disagree above ~5 kHz,
  Meyer is preferred and `discrepancy_db` records the correction.
- **Fletcher & Rossing (1998) Ch. 18 — concert bass drum modes:**
  Table 18.5 (location tag `Ch18_Table18.5_pp512-513`; measurement year
  1987 via Rossing, reprinted in FR 1998). Instrument key
  `bassdrum_82cm` records the **printed 82 cm** diameter — not equated
  with the catalogue 32-in entry. Columns digitized:
  - `specimen=both_heads` — carry-head lower tension (orchestral /
    preferred anchor default): (01)=39, (11)=80, (21)=121, (31)=162,
    (41)=204, (51)=248 Hz.
  - `specimen=both_heads_equal` — heads at same tension, including
    doublet members tagged `mode_freq_*b` where printed.
  - `specimen=batter_only` — blank / `needs_manual_reading=1` (see list
    item 4 above).
  Head thickness in FR prose: 0.010 in Mylar (used as catalogue thickness
  for the 82 cm entry). Model loader default specimen: `both_heads`.

## Mapping note (used by `AmplitudeLayer`)

Historical band power densities are treated as piecewise-constant in
frequency. Mapping onto the model ERB grid is **energy-preserving
overlap integration**: density ρᵢ = Eᵢ/Δfᵢ inside each historical band;
each ERB band receives ∫ ρ(f) df over its overlap. See README § Absolute
amplitude layer.

Residual whole-spectrum power after named bands is filled with **uniform
spectral density** (W/Hz) over uncovered historical bands
(`internal_default` fill, bandwidth-proportional). Provenance thresholds
(`internal_default`): fill_fraction ≤ 0.10 → `primary_source`; ≤ 0.60 →
`mixed_primary_and_fill`; > 0.60 → AmplitudeLayer **refuses** coverage
(equipartition kept; not labelled as a measurement).

### Instrument-name aliases (`_INSTRUMENT_SOURCE_KEYS`)

Suspended-cymbal model entries (`cymbal_16in_thin`, `cymbal_18in_medium`,
`cymbal_46cm_medium`) map onto Sivian's **15-in. clash PAIR** specimen
(`cymbals_15in`) — a different mechanical system and stroke. This alias
is an approximation labelled `internal_default` (see the
`peak_power_whole` notes cell for `cymbals_15in`).

`bassdrum_82cm` aliases Sivian `bass_drum_A_36x15` for AmplitudeLayer
(v0.3.4): the Sivian drum is 36 in (≈91 cm) vs the FR modal anchor at
82 cm — size-mismatch approximation (`internal_default`). Expected
fill ≈ 0.59 → `mixed_primary_and_fill`, same machinery as `bassdrum_32in`.
`bassdrum_32in` / `bassdrum_28in` keep their keys; they inherit the
fitted effective wave speed from the 82 cm modal anchor (`f ∝ 1/a` at
fixed `c`, `internal_default`) rather than copying measured frequencies.

### Excitation contact times (v0.3.4)

`record_type=excitation_contact_time` rows store Hertzian-impact
`t_contact_s` per implement (`stick_tip`, `stick_shoulder`, `yarn_mallet`,
`bass_drum_beater`). FR 1998 §19.7 / Fig. 19.12 discuss contact and
`f_max(v)` but do not print unambiguous absolute contact times for these
implement classes — current rows are `internal_default` placeholders
(see CHANGES.md). Source-read values must never be overwritten by
placeholders.

### Provenance label correction (values unchanged)

The four `meyer_hf_discrepancy` rows were re-labelled from
`primary_source` to `literature_derived` because their notes already
state the offsets are conservative estimates, not figure readings. **No
numeric values were changed.**

### Calibration bridge (v0.3.2+)

Until the `needs_manual_reading` histogram cells above are completed,
the scale-commensurability bridge typically has fewer than two surviving
instruments and reports **NO CALIBRATION ACHIEVED** (factor undefined).
The model side of the bridge is theory-only (partials); it does not reuse
AmplitudeLayer measured weights.

### Membrane mode anchor (v0.3.3)

`fr_ch18_bassdrum_mode` rows feed `MembraneInstrument.measured_modes`
(`both_heads` default). They are **modal frequencies**, not
AmplitudeLayer band powers. See README §2.2 and TECHNICAL_MANUAL §5.1.
