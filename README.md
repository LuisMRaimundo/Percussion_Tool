# NonTunPerc — non-tuned percussion density model (v0.3)

A bibliography-based, parametric model that generates **ERB-band spectral
density profiles** for unpitched percussion from physical input parameters
(diameter, thickness or tension, material), without audio analysis. It is
an autonomous companion to *Spectral_Analyser*: schema-aligned metadata
out, no shared code, no modification to the existing pipeline.

## Quick start

| Launcher | What it does |
|---|---|
| `run_nontunperc.bat` | Opens the **NonTunPerc GUI** |
| `run_nontunperc.bat --cli` | Headless full pipeline (profiles + MC + calibration) |
| `run_validate.bat [wav-dir]` | Validate against cymbal WAV recordings |

```text
python nontunperc.py           # GUI
python nontunperc.py --cli     # headless (replaces demo.py)
```

`demo.py` remains as a thin compatibility wrapper around the headless pipeline.

**Documentation**
- [`QUICK_REFERENCE.md`](QUICK_REFERENCE.md) — GUI options, outputs, everyday use  
- [`TECHNICAL_MANUAL.md`](TECHNICAL_MANUAL.md) — equations, architecture, methods, validity  

## 1. Input → output

**Input:** instrument family (`plate` | `membrane`), diameter, effective
thickness (plates) or tension / nominal (1,1)-mode frequency (membranes),
material preset, optional Chladni anchor `(c, b, p)`.

**Output (per instrument):**
- `modes_per_band` — modal count per 1-ERB band, 20 Hz–16 kHz;
- `energy_w_<phase>` — relative band-energy weights per temporal phase
  (plates: strike / buildup / shimmer / residue; membranes: strike / decay);
- optional `spl_db_<phase>` — absolute band levels (dB SPL) when the
  AmplitudeLayer has Sivian/Meyer coverage, at the source reference
  distance (typically **3 ft / 0.9144 m**);
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
energy transfer (Rossing, 2000, Fig. 9.6, observations 2–4). Where the
AmplitudeLayer has coverage, equipartition is replaced by measured band
weights (see §3).

### 2.2 Membranes (bass drum)

Modes `f_mk = β_mk c / (2πa)` with `β_mk` Bessel zeros and
`c = √(T/σ)`; asymptotic modal density `n(f) = 2πA f / c²` (rises
linearly). Tension may be inferred from a nominal (1,1)-mode frequency.

## 3. Absolute amplitude layer

`AmplitudeLayer` (`model.py`) loads `data/source_constants.csv`
(Rossing Table 9.1 / Fig. 9.3; Sivian–Dunn–White 1931 band edges and
peak/average anchors; Meyer 2009 corroboration).

**ERB mapping (energy-preserving overlap integration):** each historical
band energy `Eᵢ` is treated as uniform density `ρᵢ = Eᵢ/Δfᵢ` on its
printed `[f_lo, f_hi]`; each ERB band receives `∫ ρ(f) df` over the
overlap. Band edges are never resampled before this step. Above ~5 kHz,
Meyer corrections recorded in `discrepancy_db` are applied.

Relative outputs remain **bit-identical** to the equipartition path for
instruments without coverage (e.g. gong, tam-tam).

## 4. Calibration bridge

`calibration.py` runs quasi-harmonic bridge fixtures for trumpet,
clarinet, flute, and bass viol (string-family stand-in; violin full-band
spectrum is absent from Sivian 1931). It compares model composite
indices to empirical band-power indices and writes
`calibration_report.md` with the conversion factor and its spread.

**The spread IS the uncertainty to be attached to any cross-domain ratio
claim.**

## 5. Monte Carlo uncertainty (`uncertainty.py`)

`run_monte_carlo(instrument, n_draws=2000, seed=...)` perturbs, per draw:

| parameter | distribution | width | provenance |
|---|---|---|---|
| thickness | lognormal | 95% span ±25% (plates) / ±10% (membranes) | `internal_default` — manufacturing taper and hammering on cymbals/gongs; tighter membrane film tolerance |
| diameter | normal | σ = 1% of nominal | `internal_default` — nominal sizes are tight |
| E, ρ | normal | σ = 5% | `internal_default` — alloy / film variation |
| Chladni *p* | uniform on Table 9.1 p1–p2 span of the same class | class span | span itself: `primary_source`; 46 cm uses 18-in medium as nearest class (`internal_default`) |
| τ₁ₖ / α (plates) or τ₁₀₀ / α (membranes) | normal | σ = 20% | `literature_derived` fit uncertainty |

Aggregates per band / phase: p5, p25, **p50 (median)**, p75, p95 of
`modes_per_band`, `energy_w_<phase>`, and composite indices. Exports
`density_profiles_mc.csv` (+ `.meta.json` with the RNG seed) and a fan
chart for the 46-cm cymbal.

**Reporting rule:** point estimates from the deterministic
`generate_profile` path are **deprecated for citation**. The citable
output is the MC median with nested percentile intervals (symmetry with
the companion empirical pipeline's bootstrap-CI convention).

## 6. Empirical validation against recordings

`validate_against_recordings.py` is a standalone one-time check:

```text
python validate_against_recordings.py --wav-dir <folder> \
    --instrument cymbal_46cm_medium [--out report_dir]
```

It mono-mixes / resamples WAVs to 44.1 kHz, detects onsets (10× noise
floor), segments with the model's plate phase windows, computes Welch
band energies on the ERB grid, counts shimmer-phase spectral peaks at
prominence 9/12/15 dB, and compares the measured shimmer profile to the
Task-6 MC fan (Spearman ρ, 90% coverage, log-ratio bias). Results go to
`validation_report.md` only — **never** into `data/source_constants.csv`.

Welch windows: ~20 ms, 50% overlap when the phase segment is long
enough; shorter strike segments use `n//4` (documented in the script).

## 7. Built-in analytic validation (`demo.py` / tests)

- **VAL1** — 46-cm cymbal: predicted 64 modes/kHz, i.e. ≈127 modes below
  2 kHz; consistent with the >100 modes recorded holographically by
  Wilbur (in Rossing, 2000, §9.2).
- **VAL2** — n = 0 mode family reproduces the Rossing Fig. 9.3 fit to
  0.00% by construction when the anchor is supplied; the h/d² scaling
  path is the extrapolative (and weaker) branch.

## 8. Validity limits (must accompany any exported value)

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
   (AmplitudeLayer replaces this only where Sivian/Meyer coverage exists.)
5. **Scale commensurability.** Indices are on the model's own scale.
   Ratio comparisons against empirically derived (Iowa/OrchideaSOL)
   pitched-instrument metadata require the calibration bridge:
   run the same derivation for 1–2 pitched instruments with known
   empirical values and use the discrepancy as the conversion factor and
   its magnitude as the reported uncertainty.
6. **1931 high-frequency chain.** The 1931 measurement chain has limited
   high-frequency accuracy; values above ~5 kHz are corroborated or
   corrected by Meyer, and the correction is recorded (`discrepancy_db`).
7. **Specimen anchoring.** Absolute levels describe single
   specimens/players and anchor the scale, not specimen variance.
   Specimen variance is instead carried by the Monte Carlo layer (§5).

## 9. References

Cremer, L., Heckl, M., & Petersson, B. A. T. (2005). *Structure-borne
sound* (3rd ed.). Springer.

Fletcher, N. H., & Rossing, T. D. (1998). *The physics of musical
instruments* (2nd ed.). Springer.

Glasberg, B. R., & Moore, B. C. J. (1990). Derivation of auditory filter
shapes from notched-noise data. *Hearing Research, 47*(1–2), 103–138.

Meyer, J. (2009). *Acoustics and the performance of music* (5th ed.).
Springer.

Rossing, T. D. (2000). *Science of percussion instruments*. World
Scientific.

Sivian, L. J., Dunn, H. K., & White, S. D. (1931). Absolute amplitudes
and spectra of certain musical instruments and orchestras. *Journal of
the Acoustical Society of America, 2*, 330–371.
