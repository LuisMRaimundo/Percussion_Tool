# Idiophone / Membranophone Theoretical Density Model (prototype v0.1)

A bibliography-based, parametric model that generates **ERB-band spectral
density profiles** for unpitched percussion from physical input parameters
(diameter, thickness or tension, material), without audio analysis. It is
an autonomous companion to *Spectral_Analyser*: schema-aligned metadata
out, no shared code, no modification to the existing pipeline.

## 1. Input → output

**Input:** instrument family (`plate` | `membrane`), diameter, effective
thickness (plates) or tension / nominal (1,1)-mode frequency (membranes),
material preset, optional Chladni anchor `(c, b, p)`.

**Output (per instrument):**
- `modes_per_band` — modal count per 1-ERB band, 20 Hz–16 kHz;
- `energy_w_<phase>` — relative band-energy weights per temporal phase
  (plates: strike / buildup / shimmer / residue; membranes: strike / decay);
- a composite scalar index per phase (explicitly labelled as a summary;
  the per-band vector is the primary datum);
- CSV export (`density_profiles.csv`) with one row per (instrument, band).

## 2. Model core

### 2.1 Plates (cymbal, gong, tam-tam)

Low-order modes follow Chladni's law as fitted by Rossing:
`f_mn = c (m + b n)^p`, with `(c, b, p)` taken from Rossing (2000,
Table 9.1 / Fig. 9.3) when available (provenance: `literature_derived`),
or scaled from the 46-cm reference by classical plate scaling
`f ∝ h/d²` (Fletcher & Rossing, 1998) otherwise.

Above the highest explicitly enumerated mode, the **asymptotic modal
density of a flexural plate** is used (Cremer, Heckl & Petersson, 2005):

    n(f) = √3 · A / (c_L · h)      [modes/Hz, frequency-independent]

with `A` plate area, `h` thickness, `c_L = √(E/ρ(1−ν²))`. Because ERB
bandwidth grows with frequency while `n(f)` is constant, **modes per ERB
grow monotonically** — the theoretical form of the cymbal's saturation of
the upper auditory range.

Decay: `τ(f) = τ_1k (f/1 kHz)^−α`, with `τ_1k ≈ 10 s`, `α ≈ 0.84` fitted
to the log–log trend of Rossing (2000, Fig. 9.5)
(provenance: `literature_derived`).

Phase energy: equipartition over modes at excitation
(`internal_default`), evolved by exponential decay per band, with a
3–5 kHz emphasis in buildup/shimmer representing the nonlinear low→high
energy transfer (Rossing, 2000, Fig. 9.6, observations 2–4).

### 2.2 Membranes (bass drum)

Modes `f_mk = β_mk c / (2πa)` with `β_mk` Bessel zeros and
`c = √(T/σ)`; asymptotic modal density `n(f) = 2πA f / c²` (rises
linearly). Tension may be inferred from a nominal (1,1)-mode frequency.

## 3. Validation (built into `demo.py`)

- **VAL1** — 46-cm cymbal: predicted 64 modes/kHz, i.e. ≈127 modes below
  2 kHz; consistent with the >100 modes recorded holographically by
  Wilbur (in Rossing, 2000, §9.2).
- **VAL2** — n = 0 mode family reproduces the Rossing Fig. 9.3 fit to
  0.00% by construction when the anchor is supplied; the h/d² scaling
  path is the extrapolative (and weaker) branch.

## 4. Validity limits (must accompany any exported value)

1. **Flat-plate approximation.** Dome and bow curvature are not modelled;
   the Chladni anchor absorbs them at low order only.
2. **Linear regime.** The chaotic broadband regime at ff dynamics
   (Rossing, 2000, §9.4) is outside the model; at high dynamic levels the
   discrete-mode picture itself dissolves and profiles should be read as
   lower bounds on spectral occupation.
3. **Membranes in vacuo.** Air loading and two-head coupling are not
   modelled; the lowest bass-drum modes are overestimated in frequency.
4. **Equipartition at excitation** is a stated convention, not a
   measurement; stroke type and striking point are not yet parameters.
5. **Scale commensurability.** Indices are on the model's own scale.
   Ratio comparisons against empirically derived (Iowa/OrchideaSOL)
   pitched-instrument metadata require the calibration bridge:
   run the same derivation for 1–2 pitched instruments with known
   empirical values and use the discrepancy as the conversion factor and
   its magnitude as the reported uncertainty.

## 5. References

Cremer, L., Heckl, M., & Petersson, B. A. T. (2005). *Structure-borne
sound* (3rd ed.). Springer.

Fletcher, N. H., & Rossing, T. D. (1998). *The physics of musical
instruments* (2nd ed.). Springer.

Glasberg, B. R., & Moore, B. C. J. (1990). Derivation of auditory filter
shapes from notched-noise data. *Hearing Research, 47*(1–2), 103–138.

Meyer, J. (2009). *Acoustics and the performance of music* (5th ed.).
Springer. [comparative SPL layer; not yet wired into v0.1]

Rossing, T. D. (2000). *Science of percussion instruments*. World
Scientific.
