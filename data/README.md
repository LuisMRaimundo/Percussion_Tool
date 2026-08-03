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

### Provenance label correction (values unchanged)

The four `meyer_hf_discrepancy` rows were re-labelled from
`primary_source` to `literature_derived` because their notes already
state the offsets are conservative estimates, not figure readings. **No
numeric values were changed.**

### Calibration bridge (v0.3.2)

Until the `needs_manual_reading` histogram cells above are completed,
the scale-commensurability bridge typically has fewer than two surviving
instruments and reports **NO CALIBRATION ACHIEVED** (factor undefined).
The model side of the bridge is theory-only (partials); it does not reuse
AmplitudeLayer measured weights.
