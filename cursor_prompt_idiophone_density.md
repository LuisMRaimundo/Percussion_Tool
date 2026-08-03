# Normative build specification (v0.3.5)

In-repo checklist aligned with current `README.md`, `model.py`, and
`tests/`. Do not treat this file as a primary source for physical
constants — those live in `data/source_constants.csv` with provenance
labels.

## Provenance classes

`primary_source` | `derived` | `literature_derived` |
`mixed_primary_and_fill` | `internal_default`

AmplitudeLayer: fill_fraction ≤ 0.10 → `primary_source`; ≤ 0.60 →
`mixed_primary_and_fill`; > 0.60 → refuse (equipartition). Residual fill
is equal-density (∝ bandwidth), not equal energy per band.

## Equations (must hold)

- ERB bandwidth: `ERB(f) = 24.7 (4.37 f/1000 + 1)` [Glasberg & Moore]
- Plate asymptotic density: `n(f) = √3 · A / (c_L · h)` [Cremer et al.]
- Chladni: `f = c (m + b n)^p` with Table 9.1 / Fig. 9.3 anchors
- Plate scaling fallback: `f ∝ h / d²`
- Membrane: `f_mk = β_mk c / (2π a)`, `n(f) = 2π A f / c²`;
  optional FR Table 18.5 `measured_modes` override low modes;
  `c_eff` LS-fit (`derived`); siblings scale `f ∝ 1/a` at fixed `c`
  (`internal_default`)
- Decay plates: `τ(f) = τ_1k (f/1 kHz)^−α`, `τ_1k≈10 s`, `α≈0.84`
- Phase windows (plates): strike 0–20 ms; buildup 20–150 ms;
  shimmer 0.15–2 s; residue 2–6 s
- Phase windows (membranes): strike 0–50 ms; decay 50–1000 ms
- Excitation: Hertzian `E0 ∝ 1/(1+(f/f_c)^4)` when stroke+dynamic set;
  shimmer boost dynamic-gated; ff plates bypass low-pass ([R] §9.4);
  MC perturbs t_contact ±50% when stroke set; bit-identical if both unset
- Tam-tam: `PLATE_PHASES_TAMTAM`; HF boost in shimmer not bloom

## Validation

- **VAL1:** 46 cm cymbal ≈ 64 modes/kHz; >100 modes below 2 kHz
- **VAL2:** anchored n=0 family reproduces Fig. 9.3 to ~0%
- **VAL3:** bassdrum_82cm anchored modes reproduce FR Table 18.5; in-vacuo bias documented for lowest modes
- Recording validation: metadata-only auto-group; never fit parameters
  from audio; no write-back to `primary_source` rows
- Calibration bridge: model index = theory partials + equal energy per
  partial (never measured weights); empirical = measured bands only;
  <2 survivors → `NO CALIBRATION ACHIEVED` (factor undefined)

## Export schema (minimum columns)

`instrument, family, band_index, f_lo_hz, f_hi_hz, f_centre_hz,
modes_per_band, energy_w_<phase>, energy_provenance, fill_fraction, notes`

## Sanity assertions (tests)

- ERB bandwidth increases with frequency
- ERB edges span ~20 Hz–16 kHz, strictly increasing
- Plate modes/ERB grow toward high frequency
- Membrane modal density rises with f
- Phase weights normalize to 1
- High fill_fraction refuses AmplitudeLayer (`internal_default`)
- Equal-density residual: narrow uncovered band gets less energy than wide
- Calibration excludes instruments with fewer than 2 measured bands
- VAL1 / VAL2 / VAL3 as above
