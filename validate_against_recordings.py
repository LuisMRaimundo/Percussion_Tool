"""
validate_against_recordings
===========================

Standalone empirical check of the idiophone density model against real
(or synthetic) cymbal stroke recordings.

Self-contained: imports only numpy, scipy, soundfile (librosa optional
fallback for resample), matplotlib, and local ``model.py`` /
``uncertainty.py``. Does **not** write measured values into
``data/source_constants.csv`` or any ``primary_source`` field.

CLI
---
python validate_against_recordings.py --wav-dir <folder> \\
    --instrument cymbal_46cm_medium [--out report_dir]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import soundfile as sf
from scipy import signal
from scipy.stats import spearmanr

from model import (
    PLATE_PHASES,
    PlateInstrument,
    erb_band_edges,
    generate_profile,
)
from uncertainty import DEFAULT_SEED, run_monte_carlo

ROOT = Path(__file__).resolve().parent
TARGET_SR = 44100

# Peak-picking prominence defaults (internal_default).
DEFAULT_PROMINENCE_DB = 12.0
PROMINENCE_SWEEP_DB = (9.0, 12.0, 15.0)

# Instrument catalogue used when --instrument matches a known name.
_CATALOGUE = {
    "cymbal_16in_thin": PlateInstrument(
        "cymbal_16in_thin", 0.406, 0.0008, chladni=(10.8, 2.0, 1.81)
    ),
    "cymbal_18in_medium": PlateInstrument(
        "cymbal_18in_medium", 0.457, 0.0012, chladni=(13.4, 2.0, 1.65)
    ),
    "cymbal_46cm_medium": PlateInstrument(
        "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
    ),
}


# ----------------------------------------------------------------------
# Audio I/O
# ----------------------------------------------------------------------

_AUDIO_SUFFIXES = {".wav", ".aif", ".aiff", ".flac"}


def find_sample_files(
    folder: Path,
    recursive: bool = True,
) -> List[Path]:
    """List audio samples under ``folder``.

    When ``recursive`` is True (default), walks all subfolders. Skips
    ``__MACOSX`` directories and AppleDouble ``._*`` junk files. Case of
    the suffix is ignored; duplicates are de-duped by resolved path.
    """
    folder = Path(folder)
    if not folder.is_dir():
        return []
    found: set[Path] = set()
    iterator = folder.rglob("*") if recursive else folder.iterdir()
    for p in iterator:
        if not p.is_file():
            continue
        if p.name.startswith("._"):
            continue
        if "__MACOSX" in p.parts:
            continue
        if p.suffix.lower() in _AUDIO_SUFFIXES:
            found.add(p.resolve())
    return sorted(found)


def load_mono(path: Path, target_sr: int = TARGET_SR) -> Tuple[np.ndarray, int]:
    """Load WAV/AIFF/FLAC, mono-mix if stereo, resample to target_sr if needed."""
    y, sr = sf.read(str(path), always_2d=True)
    y = np.mean(y, axis=1).astype(np.float64)
    if sr != target_sr:
        # Polyphase resample (scipy); avoid librosa dependency for core path.
        g = np.gcd(sr, target_sr)
        up, down = target_sr // g, sr // g
        y = signal.resample_poly(y, up, down)
        sr = target_sr
    return y, sr


def detect_onset(y: np.ndarray, sr: int, frame_ms: float = 2.0) -> int:
    """Energy-based onset: first frame exceeding 10× noise-floor energy.

    Noise floor = median frame energy of the first 50 ms (or 5% of file).
    """
    frame = max(1, int(sr * frame_ms / 1000.0))
    n = len(y)
    if n < frame * 4:
        return 0
    energies = np.array(
        [np.mean(y[i : i + frame] ** 2) for i in range(0, n - frame, frame)]
    )
    n_floor = max(1, int(0.05 * len(energies)))
    floor = float(np.median(energies[:n_floor])) + 1e-20
    thr = 10.0 * floor
    hits = np.where(energies > thr)[0]
    if hits.size == 0:
        # Fallback: peak absolute sample.
        return int(np.argmax(np.abs(y)))
    return int(hits[0] * frame)


def phase_segments(
    y: np.ndarray,
    sr: int,
    onset: int,
    phases: Dict[str, Tuple[float, float]] = PLATE_PHASES,
) -> Dict[str, Tuple[np.ndarray, Tuple[float, float]]]:
    """Slice ``y`` into model phase windows truncated to file length."""
    out = {}
    n = len(y)
    for name, (t0, t1) in phases.items():
        i0 = onset + int(round(t0 * sr))
        i1 = onset + int(round(t1 * sr))
        i0 = max(0, min(i0, n))
        i1 = max(i0, min(i1, n))
        seg = y[i0:i1]
        actual = (i0 / sr, i1 / sr)
        out[name] = (seg, actual)
    return out


# ----------------------------------------------------------------------
# Spectral features
# ----------------------------------------------------------------------

def _welch_params(n: int, sr: int) -> Tuple[int, int]:
    """Choose nperseg / noverlap from segment length (documented defaults).

    Frequency resolution must resolve modal spacing for long segments:
    - shimmer / residue (≳0.5 s): ~100 ms windows (or 8192 samples), 50% overlap;
    - buildup: ~40 ms;
    - strike (short): n//2, no less than 8 samples;
    - always ``nperseg <= n``.
    """
    dur = n / float(sr)
    if n < 16:
        nperseg = max(4, n)
    elif dur >= 0.5:
        nperseg = min(max(int(0.100 * sr), 4096), 8192)
    elif dur >= 0.1:
        nperseg = int(0.040 * sr)
    else:
        nperseg = max(8, n // 2)
    nperseg = int(min(nperseg, n))
    noverlap = nperseg // 2
    return nperseg, noverlap


def band_energies_welch(
    seg: np.ndarray,
    sr: int,
    edges: np.ndarray,
) -> np.ndarray:
    """Welch PSD integrated over ERB bands; relative weights (sum=1)."""
    n = len(seg)
    n_bands = len(edges) - 1
    if n < 8:
        return np.zeros(n_bands)
    nperseg, noverlap = _welch_params(n, sr)
    freqs, psd = signal.welch(
        seg, fs=sr, nperseg=nperseg, noverlap=noverlap, scaling="density"
    )
    e = np.zeros(n_bands)
    for i in range(n_bands):
        lo, hi = edges[i], edges[i + 1]
        mask = (freqs >= lo) & (freqs < hi)
        if not np.any(mask):
            continue
        # Trapezoidal integrate PSD over band (V^2).
        integ = getattr(np, "trapezoid", None) or np.trapz
        e[i] = float(integ(psd[mask], freqs[mask]))
    s = e.sum()
    return e / s if s > 0 else e


def count_peaks_per_band(
    seg: np.ndarray,
    sr: int,
    edges: np.ndarray,
    prominence_db: float = DEFAULT_PROMINENCE_DB,
) -> np.ndarray:
    """High-res Hann FFT peak count per ERB band (shimmer phase)."""
    n_bands = len(edges) - 1
    n = len(seg)
    if n < 32:
        return np.zeros(n_bands)
    # Zero-pad to at least 4× length, next pow2 (internal_default).
    nfft = 1 << int(np.ceil(np.log2(max(4 * n, 4096))))
    window = np.hanning(n)
    spec = np.fft.rfft(seg * window, n=nfft)
    mag = np.abs(spec)
    freqs = np.fft.rfftfreq(nfft, d=1.0 / sr)
    mag_db = 20.0 * np.log10(mag + 1e-20)
    med = float(np.median(mag_db))
    prom = prominence_db  # dB above median spectral level
    # find_peaks prominence is linear-magnitude; convert approx via ratio.
    # Use dB-threshold on mag_db directly with a simple local-max picker.
    peaks = []
    for i in range(2, len(mag_db) - 2):
        if mag_db[i] >= mag_db[i - 1] and mag_db[i] >= mag_db[i + 1]:
            if mag_db[i] >= med + prom:
                peaks.append(i)
    counts = np.zeros(n_bands)
    for idx in peaks:
        f = freqs[idx]
        if f < edges[0] or f >= edges[-1]:
            continue
        b = int(np.searchsorted(edges, f, side="right") - 1)
        if 0 <= b < n_bands:
            counts[b] += 1
    return counts


def octave_decay_times(
    y: np.ndarray,
    sr: int,
    onset: int,
    f_lo: float = 62.5,
    f_hi: float = 8000.0,
) -> Dict[float, float]:
    """Per-octave-band 60-dB decay proxy via linear fit to dB envelope.

    Returns {band_centre_hz: tau_60_s}.
    """
    # Build octave edges from f_lo.
    edges = [f_lo]
    while edges[-1] * 2 < f_hi:
        edges.append(edges[-1] * 2)
    edges.append(f_hi)
    edges = np.asarray(edges, dtype=float)

    # Analytic envelopes via bandpass + Hilbert on post-onset audio.
    seg = y[onset:]
    if len(seg) < sr // 10:
        return {}
    # Frame the decay: hop 10 ms.
    hop = max(1, int(0.010 * sr))
    out = {}
    for lo, hi in zip(edges[:-1], edges[1:]):
        # Butterworth bandpass
        wn = np.array([lo, hi]) / (sr / 2.0)
        if wn[1] >= 1.0:
            wn[1] = 0.999
        if wn[0] <= 0:
            continue
        try:
            b, a = signal.butter(2, wn, btype="band")
            filt = signal.filtfilt(b, a, seg)
        except Exception:
            continue
        # Frame RMS in dB
        n_frames = max(1, (len(filt) - hop) // hop)
        env = []
        times = []
        for i in range(n_frames):
            sl = filt[i * hop : (i + 1) * hop]
            rms = np.sqrt(np.mean(sl ** 2) + 1e-20)
            env.append(20.0 * np.log10(rms))
            times.append((i + 0.5) * hop / sr)
        env = np.asarray(env)
        times = np.asarray(times)
        # Fit from peak to −30 dB relative (or end); extrapolate to −60.
        peak_i = int(np.argmax(env))
        peak = env[peak_i]
        mask = np.arange(len(env)) >= peak_i
        usable = mask & (env > peak - 30.0)
        if usable.sum() < 5:
            continue
        # Linear fit env ≈ a + b t ; tau_60 = -60 / b
        coef = np.polyfit(times[usable], env[usable], 1)
        slope = coef[0]
        if slope >= -1e-6:
            continue
        tau60 = -60.0 / slope
        if 0.01 < tau60 < 1000:
            out[float(np.sqrt(lo * hi))] = float(tau60)
    return out


# ----------------------------------------------------------------------
# Per-file analysis
# ----------------------------------------------------------------------

def analyse_file(
    path: Path,
    edges: np.ndarray,
    prominence_db: float = DEFAULT_PROMINENCE_DB,
) -> dict:
    y, sr = load_mono(path)
    onset = detect_onset(y, sr)
    segs = phase_segments(y, sr, onset)
    result = {
        "file": path.name,
        "sr": sr,
        "onset_s": onset / sr,
        "phases": {},
        "peak_counts": {},
        "decay_tau": octave_decay_times(y, sr, onset),
    }
    for ph, (seg, actual) in segs.items():
        w = band_energies_welch(seg, sr, edges)
        result["phases"][ph] = {
            "weights": w,
            "window_s": actual,
            "n_samples": len(seg),
        }
    # Peak counts on shimmer only, with prominence sweep.
    shimmer_seg = segs["shimmer"][0]
    for pdb in PROMINENCE_SWEEP_DB:
        result["peak_counts"][str(pdb)] = count_peaks_per_band(
            shimmer_seg, sr, edges, prominence_db=pdb
        )
    result["peak_counts"]["default"] = result["peak_counts"][
        str(prominence_db)
    ]
    return result


# ----------------------------------------------------------------------
# Aggregation / comparison
# ----------------------------------------------------------------------

def aggregate_weights(
    file_results: Sequence[dict], phase: str
) -> Tuple[np.ndarray, np.ndarray]:
    """Median and MAD-like spread across files for one phase."""
    stack = np.stack(
        [r["phases"][phase]["weights"] for r in file_results], axis=0
    )
    med = np.median(stack, axis=0)
    spread = np.median(np.abs(stack - med), axis=0)
    return med, spread


def compare_to_model(
    measured_med: np.ndarray,
    model_result,
    phase: str = "shimmer",
) -> dict:
    """Spearman, interval coverage, log-ratio bias vs MC fan."""
    mq = model_result.energy_quantiles[phase]
    model_med = mq["p50"]
    # Align lengths
    n = min(len(measured_med), len(model_med))
    m = measured_med[:n]
    mod = model_med[:n]
    p05, p95 = mq["p05"][:n], mq["p95"][:n]

    # Avoid zeros in rank / log
    m_eps = m + 1e-12
    mod_eps = mod + 1e-12
    rho, pval = spearmanr(m_eps, mod_eps)
    inside = (m >= p05) & (m <= p95)
    frac_inside = float(np.mean(inside))
    log_bias = np.log(m_eps) - np.log(mod_eps)

    centres = model_result.band_centres[:n]
    # Weak-point masks
    low_mask = centres < 500.0
    hf_mask = centres > 10000.0
    return {
        "spearman_rho": float(rho) if np.isfinite(rho) else float("nan"),
        "spearman_p": float(pval) if np.isfinite(pval) else float("nan"),
        "frac_inside_90": frac_inside,
        "log_ratio_bias": log_bias,
        "log_ratio_bias_mean": float(np.mean(log_bias)),
        "log_ratio_bias_low_f": float(np.mean(log_bias[low_mask]))
        if low_mask.any()
        else float("nan"),
        "log_ratio_bias_above_10k": float(np.mean(log_bias[hf_mask]))
        if hf_mask.any()
        else float("nan"),
        "n_bands": n,
    }


def plot_comparison(
    measured_med: np.ndarray,
    measured_spread: np.ndarray,
    model_result,
    phase: str,
    path: Path,
) -> None:
    x = model_result.band_centres
    n = min(len(x), len(measured_med))
    x = x[:n]
    mq = model_result.energy_quantiles[phase]
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    ax.fill_between(
        x, mq["p05"][:n], mq["p95"][:n], color="C0", alpha=0.20, label="model 90%"
    )
    ax.fill_between(
        x, mq["p25"][:n], mq["p75"][:n], color="C0", alpha=0.35, label="model 50%"
    )
    ax.plot(x, mq["p50"][:n], color="C0", lw=1.6, label="model median")
    ax.plot(x, measured_med[:n], color="C3", lw=1.6, label="measured median")
    ax.fill_between(
        x,
        np.maximum(measured_med[:n] - measured_spread[:n], 0),
        measured_med[:n] + measured_spread[:n],
        color="C3",
        alpha=0.25,
        label="measured ± MAD",
    )
    ax.set(
        xscale="log",
        xlabel="ERB-band centre (Hz)",
        ylabel=f"relative band energy ({phase})",
        title=f"Validation: measured vs model ({model_result.instrument})",
    )
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_report(
    out_dir: Path,
    instrument_name: str,
    metrics: dict,
    file_results: Sequence[dict],
    model_meta: dict,
) -> Path:
    path = out_dir / "validation_report.md"
    lines = [
        f"# Validation report — `{instrument_name}`",
        "",
        f"Files analysed: **{len(file_results)}**",
        f"Model MC seed: `{model_meta.get('seed')}`, "
        f"N={model_meta.get('n_draws')}",
        "",
        "## Metrics (shimmer phase)",
        "",
        f"- Spearman rank correlation (measured vs model median): "
        f"**{metrics['spearman_rho']:.4f}** (p={metrics['spearman_p']:.3g})",
        f"- Fraction of measured bands inside model 90% interval: "
        f"**{metrics['frac_inside_90']:.3f}**",
        f"- Mean log-ratio bias ln(meas/model): "
        f"**{metrics['log_ratio_bias_mean']:.4f}**",
        f"- Mean log-ratio bias below 500 Hz: "
        f"**{metrics['log_ratio_bias_low_f']:.4f}**",
        f"- Mean log-ratio bias above 10 kHz: "
        f"**{metrics['log_ratio_bias_above_10k']:.4f}**",
        "",
        "## Peak-count sensitivity (shimmer, prominence sweep)",
        "",
    ]
    # Aggregate peak totals
    for pdb in PROMINENCE_SWEEP_DB:
        totals = [
            float(np.sum(r["peak_counts"][str(pdb)])) for r in file_results
        ]
        lines.append(
            f"- prominence {pdb:.0f} dB: median total peaks = "
            f"{np.median(totals):.1f} "
            f"(range {np.min(totals):.0f}–{np.max(totals):.0f})"
        )
    lines += [
        "",
        "The prominence sensitivity band is part of the result: absolute",
        "peak counts are threshold-dependent; rank/shape comparisons are",
        "the durable claim.",
        "",
        "## Verdict template",
        "",
        "### Claims this validation supports",
        "",
        "- **Profile shape / rank order** of shimmer-phase relative band",
        "  energies, if Spearman ρ is high and most bands fall inside the",
        "  model 90% fan.",
        "- **Qualitative high-frequency occupation** of the auditory scale",
        "  for plates (modes-per-ERB growing with frequency), when peak",
        "  counts and energy weights both rise toward the upper midrange.",
        "",
        "### Claims this validation does **not** support",
        "",
        "- **Absolute levels** (dB SPL): the script compares relative",
        "  weights only; recording gain and mic distance are uncontrolled.",
        "- **Loud-dynamics / chaotic broadband behaviour** (Rossing §9.4):",
        "  outside the model's linear discrete-mode regime.",
        "- **Low-frequency strike/buildup phases**: dominated by",
        "  nonlinearity and transient radiation not captured by the",
        "  equipartition + linear decay picture.",
        "- **Levels above ~10 kHz**: recording-chain and microphone",
        "  dependent; treat disagreements here as inconclusive for the",
        "  physical model.",
        "",
        "### Bands / regimes outside stated validity",
        "",
        "- Strike and early buildup (nonlinearity).",
        "- Above ~10 kHz (measurement chain).",
        "- Any ff / crash stroke that leaves the linear modal regime.",
        "",
        "## Per-file windows (actual)",
        "",
    ]
    for r in file_results:
        lines.append(f"- `{r['file']}` onset={r['onset_s']:.4f}s")
        for ph, info in r["phases"].items():
            t0, t1 = info["window_s"]
            lines.append(
                f"  - {ph}: [{t0:.4f}, {t1:.4f}] s "
                f"({info['n_samples']} samples)"
            )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


# ----------------------------------------------------------------------
# Public entry
# ----------------------------------------------------------------------

def run_validation(
    wav_dir: Path,
    instrument_name: str = "cymbal_46cm_medium",
    out_dir: Optional[Path] = None,
    mc_draws: int = 400,
    mc_seed: int = DEFAULT_SEED,
    recursive: bool = True,
) -> dict:
    """Run the full validation pipeline; write report + figure."""
    wav_dir = Path(wav_dir)
    out_dir = Path(out_dir) if out_dir else ROOT / "validation_out"
    out_dir.mkdir(parents=True, exist_ok=True)

    if instrument_name not in _CATALOGUE:
        raise ValueError(
            f"unknown instrument {instrument_name!r}; "
            f"known: {sorted(_CATALOGUE)}"
        )
    instr = _CATALOGUE[instrument_name]
    edges = erb_band_edges(20.0, 16000.0)

    wavs = find_sample_files(wav_dir, recursive=recursive)
    if not wavs:
        scope = "recursively under" if recursive else "in"
        raise FileNotFoundError(
            f"no audio samples (.wav/.aif/.aiff/.flac) {scope} {wav_dir}"
        )

    file_results = [analyse_file(p, edges) for p in wavs]
    meas_med, meas_spread = aggregate_weights(file_results, "shimmer")

    mc = run_monte_carlo(instr, n_draws=mc_draws, seed=mc_seed)
    metrics = compare_to_model(meas_med, mc, phase="shimmer")

    fig_path = out_dir / "validation_comparison.png"
    plot_comparison(meas_med, meas_spread, mc, "shimmer", fig_path)
    report = write_report(
        out_dir, instrument_name, metrics, file_results, mc.metadata
    )

    # JSON sidecar (no write-back to source_constants)
    summary = {
        "instrument": instrument_name,
        "n_files": len(file_results),
        "recursive": recursive,
        "sample_dir": str(wav_dir.resolve()),
        "metrics": {
            k: (v.tolist() if isinstance(v, np.ndarray) else v)
            for k, v in metrics.items()
        },
        "report": str(report),
        "figure": str(fig_path),
    }
    (out_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def launch_validate_gui() -> None:
    """Desktop UI to pick the WAV sample folder and run validation."""
    import os
    import threading
    import traceback
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk

    root = tk.Tk()
    root.title("NonTunPerc — Validate against recordings")
    root.geometry("720x520")
    root.minsize(640, 440)

    bg = "#1e2a24"
    panel = "#2a3a32"
    accent = "#c4a35a"
    text = "#f2efe6"
    muted = "#a8b5ad"
    root.configure(bg=bg)

    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    style.configure("TFrame", background=bg)
    style.configure("Card.TFrame", background=panel)
    style.configure("TLabel", background=bg, foreground=text, font=("Segoe UI", 10))
    style.configure(
        "Title.TLabel", background=bg, foreground=accent,
        font=("Georgia", 16, "bold"),
    )
    style.configure(
        "Sub.TLabel", background=bg, foreground=muted, font=("Segoe UI", 9)
    )
    style.configure(
        "Card.TLabel", background=panel, foreground=text, font=("Segoe UI", 10)
    )
    style.configure("TButton", font=("Segoe UI", 10), padding=6)
    style.configure(
        "Accent.TButton", background=accent, foreground="#1a1510",
        font=("Segoe UI", 10, "bold"), padding=8,
    )
    style.map("Accent.TButton", background=[("active", "#d4b56a")])
    style.configure("TCombobox", padding=4)

    hdr = ttk.Frame(root)
    hdr.pack(fill="x", padx=16, pady=(14, 6))
    ttk.Label(hdr, text="Validate against WAVs", style="Title.TLabel").pack(
        anchor="w"
    )
    ttk.Label(
        hdr,
        text="Choose the folder that holds single-stroke samples "
        "(.wav / .aif / .flac). Subfolders are searched by default. "
        "Report only — no write-back to source constants.",
        style="Sub.TLabel",
    ).pack(anchor="w")

    card = ttk.Frame(root, style="Card.TFrame")
    card.pack(fill="x", padx=16, pady=8)
    card.configure(padding=12)

    wav_var = tk.StringVar(value=str(ROOT / "wavs"))
    out_var = tk.StringVar(value=str(ROOT / "validation_out"))
    instr_var = tk.StringVar(value="cymbal_46cm_medium")
    draws_var = tk.StringVar(value="400")
    recursive_var = tk.BooleanVar(value=True)
    status_var = tk.StringVar(value="Select a sample folder, then Run.")

    def browse_wav() -> None:
        initial = wav_var.get() if Path(wav_var.get()).is_dir() else str(ROOT)
        d = filedialog.askdirectory(
            title="Select folder with audio samples (searches subfolders)",
            initialdir=initial,
        )
        if d:
            wav_var.set(d)
            refresh_count()

    def browse_out() -> None:
        initial = out_var.get() if Path(out_var.get()).is_dir() else str(ROOT)
        d = filedialog.askdirectory(
            title="Select report output folder",
            initialdir=initial,
        )
        if d:
            out_var.set(d)

    def refresh_count(*_args) -> None:
        folder = Path(wav_var.get())
        files = find_sample_files(folder, recursive=recursive_var.get())
        n = len(files)
        scope = "including subfolders" if recursive_var.get() else "top level only"
        status_var.set(f"{n} sample file(s) found ({scope}).")

    row = 0
    ttk.Label(card, text="Sample folder", style="Card.TLabel").grid(
        row=row, column=0, sticky="w", pady=4
    )
    row += 1
    wav_row = ttk.Frame(card, style="Card.TFrame")
    wav_row.grid(row=row, column=0, sticky="ew", pady=2)
    card.columnconfigure(0, weight=1)
    tk.Entry(wav_row, textvariable=wav_var, width=56).pack(
        side="left", fill="x", expand=True
    )
    ttk.Button(wav_row, text="Browse…", command=browse_wav).pack(
        side="left", padx=(6, 0)
    )
    row += 1

    style.configure(
        "TCheckbutton", background=panel, foreground=text, font=("Segoe UI", 9)
    )
    ttk.Checkbutton(
        card,
        text="Search inside subfolders",
        variable=recursive_var,
        command=refresh_count,
    ).grid(row=row, column=0, sticky="w", pady=(6, 2))
    row += 1

    ttk.Label(card, text="Instrument model", style="Card.TLabel").grid(
        row=row, column=0, sticky="w", pady=(10, 2)
    )
    row += 1
    ttk.Combobox(
        card,
        textvariable=instr_var,
        values=sorted(_CATALOGUE),
        state="readonly",
        width=40,
    ).grid(row=row, column=0, sticky="w")
    row += 1

    ttk.Label(card, text="Report output folder", style="Card.TLabel").grid(
        row=row, column=0, sticky="w", pady=(10, 2)
    )
    row += 1
    out_row = ttk.Frame(card, style="Card.TFrame")
    out_row.grid(row=row, column=0, sticky="ew", pady=2)
    tk.Entry(out_row, textvariable=out_var, width=56).pack(
        side="left", fill="x", expand=True
    )
    ttk.Button(out_row, text="Browse…", command=browse_out).pack(
        side="left", padx=(6, 0)
    )
    row += 1

    ttk.Label(card, text="MC draws (model fan)", style="Card.TLabel").grid(
        row=row, column=0, sticky="w", pady=(10, 2)
    )
    row += 1
    tk.Spinbox(
        card, from_=40, to=5000, increment=20, width=10, textvariable=draws_var
    ).grid(row=row, column=0, sticky="w")

    log_box = scrolledtext.ScrolledText(
        root,
        wrap="word",
        height=14,
        bg="#121a16",
        fg=text,
        insertbackground=text,
        font=("Consolas", 9),
        relief="flat",
    )
    log_box.pack(fill="both", expand=True, padx=16, pady=(4, 8))

    btn_row = ttk.Frame(root)
    btn_row.pack(fill="x", padx=16, pady=(0, 12))
    ttk.Label(btn_row, textvariable=status_var, style="Sub.TLabel").pack(
        side="left"
    )

    running = {"flag": False}

    def append_log(msg: str) -> None:
        root.after(0, lambda: (log_box.insert("end", msg + "\n"), log_box.see("end")))

    def open_report_dir() -> None:
        folder = Path(out_var.get())
        folder.mkdir(parents=True, exist_ok=True)
        os.startfile(str(folder))  # noqa: S606

    def on_run() -> None:
        if running["flag"]:
            return
        wav_dir = Path(wav_var.get())
        if not wav_dir.is_dir():
            messagebox.showwarning(
                "Validate", f"Sample folder not found:\n{wav_dir}"
            )
            return
        recursive = bool(recursive_var.get())
        n = len(find_sample_files(wav_dir, recursive=recursive))
        if n == 0:
            messagebox.showwarning(
                "Validate",
                f"No audio samples (.wav/.aif/.flac) found "
                f"{'under' if recursive else 'in'}:\n{wav_dir}",
            )
            return
        running["flag"] = True
        run_btn.state(["disabled"])
        status_var.set("Running validation…")
        append_log("=== Validation start ===")
        append_log(f"Sample dir: {wav_dir}")
        append_log(f"Recursive : {recursive}")
        append_log(f"Instrument: {instr_var.get()}")
        append_log(f"Out dir   : {out_var.get()}")
        append_log(f"Files     : {n}")

        def worker() -> None:
            try:
                summary = run_validation(
                    wav_dir,
                    instrument_name=instr_var.get(),
                    out_dir=Path(out_var.get()),
                    mc_draws=int(draws_var.get()),
                    recursive=recursive,
                )
                rho = summary["metrics"]["spearman_rho"]
                inside = summary["metrics"]["frac_inside_90"]
                append_log(
                    f"[VAL] files={summary['n_files']}  "
                    f"Spearman={rho:.4f}  inside90={inside:.3f}"
                )
                append_log(f"[OUT] {summary['report']}")
                root.after(
                    0,
                    lambda: (
                        status_var.set("Finished."),
                        messagebox.showinfo(
                            "Validate",
                            f"Done.\n\nSpearman ρ = {rho:.4f}\n"
                            f"Inside 90% = {inside:.3f}\n\n"
                            f"{summary['report']}",
                        ),
                    ),
                )
            except Exception as exc:
                append_log(traceback.format_exc())
                root.after(
                    0,
                    lambda: (
                        status_var.set("Failed."),
                        messagebox.showerror("Validate", str(exc)),
                    ),
                )
            finally:
                running["flag"] = False
                root.after(0, lambda: run_btn.state(["!disabled"]))

        threading.Thread(target=worker, daemon=True).start()

    run_btn = ttk.Button(
        btn_row, text="Run validation", style="Accent.TButton", command=on_run
    )
    run_btn.pack(side="right", padx=(6, 0))
    ttk.Button(btn_row, text="Open report folder", command=open_report_dir).pack(
        side="right"
    )

    wav_var.trace_add("write", refresh_count)
    refresh_count()
    append_log(
        "Ready. Browse to a folder of samples "
        "(.wav / .aif / .flac; subfolders included by default)."
    )
    root.mainloop()


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Validate idiophone density model against cymbal WAVs"
    )
    ap.add_argument(
        "--wav-dir", type=Path, default=None,
        help="folder of audio samples (omit to open the GUI)",
    )
    ap.add_argument(
        "--instrument", default="cymbal_46cm_medium",
        choices=sorted(_CATALOGUE),
    )
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--mc-draws", type=int, default=400)
    ap.add_argument("--mc-seed", type=int, default=DEFAULT_SEED)
    ap.add_argument(
        "--no-recursive", action="store_true",
        help="do not search subfolders (top level only)",
    )
    ap.add_argument("--gui", action="store_true", help="force GUI")
    ap.add_argument("--cli", action="store_true", help="require CLI (needs --wav-dir)")
    args = ap.parse_args(argv)

    if args.gui or (args.wav_dir is None and not args.cli):
        launch_validate_gui()
        return 0

    if args.wav_dir is None:
        ap.error("--wav-dir is required with --cli")

    summary = run_validation(
        args.wav_dir,
        instrument_name=args.instrument,
        out_dir=args.out,
        mc_draws=args.mc_draws,
        mc_seed=args.mc_seed,
        recursive=not args.no_recursive,
    )
    print(
        f"[VAL] files={summary['n_files']}  "
        f"Spearman={summary['metrics']['spearman_rho']:.4f}  "
        f"inside90={summary['metrics']['frac_inside_90']:.3f}"
    )
    print(f"[OUT] {summary['report']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
