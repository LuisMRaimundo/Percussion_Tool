# Normative build specification (reconstructed)

> This file was **not present** in the shipping folder at audit time.
> It is reconstructed from `README.md` (v0.1) and the working
> `model.py` / `demo.py` so that Task-1 audits have a local normative
> checklist. Do not treat reconstructed prose as a primary source for
> physical constants.

## Provenance classes

`primary_source` | `derived` | `literature_derived` | `internal_default`

## Equations (must hold)

- ERB bandwidth: `ERB(f) = 24.7 (4.37 f/1000 + 1)` [Glasberg & Moore]
- Plate asymptotic density: `n(f) = √3 · A / (c_L · h)` [Cremer et al.]
- Chladni: `f = c (m + b n)^p` with Table 9.1 / Fig. 9.3 anchors
- Plate scaling fallback: `f ∝ h / d²`
- Membrane: `f_mk = β_mk c / (2π a)`, `n(f) = 2π A f / c²`
- Decay plates: `τ(f) = τ_1k (f/1 kHz)^−α`, `τ_1k≈10 s`, `α≈0.84`
- Phase windows (plates): strike 0–20 ms; buildup 20–150 ms;
  shimmer 0.15–2 s; residue 2–6 s
- Phase windows (membranes): strike 0–50 ms; decay 50–1000 ms

## Validation

- **VAL1:** 46 cm cymbal ≈ 64 modes/kHz; >100 modes below 2 kHz
- **VAL2:** anchored n=0 family reproduces Fig. 9.3 to ~0%

## Export schema (minimum columns)

`instrument, family, band_index, f_lo_hz, f_hi_hz, f_centre_hz,
modes_per_band, energy_w_<phase>`

## Sanity assertions (tests)

- ERB bandwidth increases with frequency
- ERB edges span ~20 Hz–16 kHz, strictly increasing
- Plate modes/ERB grow toward high frequency
- Membrane modal density rises with f
- Phase weights normalize to 1
- VAL1 / VAL2 as above
