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
