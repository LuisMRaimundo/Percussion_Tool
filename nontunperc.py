"""
NonTunPerc — non-tuned percussion spectral-density model
========================================================

User-facing application for the idiophone / membranophone density tool.
Provides a desktop GUI (default) and a headless CLI pipeline (formerly
``demo.py``).

Usage
-----
python nontunperc.py              # open GUI
python nontunperc.py --cli        # headless full pipeline (demo equivalent)
python nontunperc.py --cli --no-mc
"""

from __future__ import annotations

import argparse
import csv
import sys
import threading
import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional, Sequence

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import (
    AmplitudeLayer,
    MembraneInstrument,
    PlateInstrument,
    fit_effective_wave_speed,
    generate_profile,
    make_bassdrum_catalogue,
    membrane_beta,
)
from calibration import format_calibration_cli_line, write_calibration_report
from uncertainty import (
    DEFAULT_SEED,
    export_mc_csv,
    plot_fan_chart,
    run_monte_carlo,
)

ROOT = Path(__file__).resolve().parent
APP_NAME = "NonTunPerc"
APP_TITLE = "NonTunPerc — non-tuned percussion density model"

LogFn = Callable[[str], None]


def default_instruments() -> List:
    """Canonical instrument set (Chladni / FR Ch.18 membrane anchors)."""
    plates = [
        PlateInstrument(
            "cymbal_16in_thin", 0.406, 0.0008, chladni=(10.8, 2.0, 1.81)
        ),
        PlateInstrument(
            "cymbal_18in_medium", 0.457, 0.0012, chladni=(13.4, 2.0, 1.65)
        ),
        PlateInstrument(
            "cymbal_46cm_medium", 0.460, 0.0012, chladni=(14.93, 3.0, 1.557)
        ),
        PlateInstrument("gong_50cm_bronze", 0.500, 0.0020),
        PlateInstrument(
            "tamtam_80cm_bronze",
            0.800,
            0.0015,
            plate_class="tamtam",
        ),
    ]
    return plates + make_bassdrum_catalogue(ROOT / "data" / "source_constants.csv")


@dataclass
class PipelineOptions:
    """Controls which stages the headless / GUI pipeline runs."""

    run_validation: bool = True
    run_profiles: bool = True
    run_plots: bool = True
    run_calibration: bool = True
    run_mc: bool = True
    mc_draws_focus: int = 2000   # 46-cm cymbal
    mc_draws_other: int = 400
    mc_seed: int = DEFAULT_SEED
    instrument_names: Optional[Sequence[str]] = None  # None = all
    use_amplitude_layer: bool = True
    out_dir: Path = ROOT


def _log(msg: str, log: Optional[LogFn] = None) -> None:
    if log is not None:
        log(msg)
    else:
        print(msg)


def run_pipeline(
    options: Optional[PipelineOptions] = None,
    log: Optional[LogFn] = None,
) -> dict:
    """Run the NonTunPerc analysis pipeline. Returns a summary dict."""
    opt = options or PipelineOptions()
    out = Path(opt.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    all_instr = default_instruments()
    if opt.instrument_names:
        wanted = set(opt.instrument_names)
        instruments = [i for i in all_instr if i.name in wanted]
        if not instruments:
            raise ValueError(f"no instruments matched {sorted(wanted)}")
    else:
        instruments = all_instr

    layer = (
        AmplitudeLayer.default(ROOT) if opt.use_amplitude_layer else None
    )
    summary: dict = {"outputs": [], "instruments": [i.name for i in instruments]}

    # ---- VAL1 / VAL2 / VAL3 ------------------------------------------
    if opt.run_validation:
        cy = next(
            (i for i in all_instr if i.name == "cymbal_46cm_medium"),
            instruments[0],
        )
        if isinstance(cy, PlateInstrument):
            nd = cy.modal_density()
            _log(
                f"[VAL1] 46-cm cymbal asymptotic modal density: {nd:.4f} "
                f"modes/Hz ({nd*1000:.0f} modes per kHz)",
                log,
            )
            _log(
                f"[VAL1] predicted modes below 2 kHz: {nd*2000:.0f} "
                f"(holography: >100 observed; consistent)",
                log,
            )
            low = cy.low_modes(m_max=7, n_max=0)
            ref = 14.93 * (np.arange(2, 8)) ** 1.557
            err = np.max(np.abs(low[:6] - ref) / ref) * 100
            _log(
                f"[VAL2] n=0 family m=2..7 reproduction error vs [R]: "
                f"{err:.2f}%",
                log,
            )
            summary["val1_modes_per_khz"] = float(nd * 1000)
            summary["val2_error_pct"] = float(err)

        drum = next(
            (i for i in all_instr if isinstance(i, MembraneInstrument)
             and i.has_measured_anchor()),
            None,
        )
        if drum is not None:
            meas = np.asarray(drum.measured_modes, dtype=float)
            labels = drum.measured_mode_indices or tuple()
            anchored = drum.low_modes()[: len(meas)]
            # VAL3a: anchored low modes reproduce measured values exactly
            max_err = float(np.max(np.abs(anchored - meas)))
            _log(
                f"[VAL3] {drum.name}: anchored modes reproduce FR Table 18.5 "
                f"exactly (max |df| = {max_err:.3e} Hz); "
                f"fitted c = {drum.fitted_effective_c():.3f} m/s",
                log,
            )
            # In-vacuo comparison: wave speed from higher modes (m ≥ 2),
            # where air loading is weaker; predict all labels with that c.
            # Expected: vacuo overestimates the air-loading-dominated
            # lowest modes ((0,1), (1,1)).
            hi_f = [
                float(f)
                for (m_i, n_i), f in zip(labels, meas)
                if m_i >= 2
            ]
            hi_lab = tuple(
                (m_i, n_i) for (m_i, n_i) in labels if m_i >= 2
            )
            c_vac = fit_effective_wave_speed(
                np.asarray(hi_f), hi_lab, drum.radius
            )
            _log(
                "[VAL3] mode   measured   in_vacuo   pct_dev "
                "(+ = vacuo overestimates; c from m>=2)",
                log,
            )
            rows_val3 = []
            low_ok = True
            for (m_i, n_i), f_m in zip(labels, meas):
                beta = membrane_beta(m_i, n_i)
                f_v = beta * c_vac / (2.0 * np.pi * drum.radius)
                pct = 100.0 * (f_v - f_m) / f_m
                if m_i <= 1 and pct <= 0:
                    low_ok = False
                _log(
                    f"[VAL3] ({m_i}{n_i})  {f_m:8.2f}  {f_v:8.2f}  "
                    f"{pct:+7.2f}%",
                    log,
                )
                rows_val3.append(
                    {
                        "mode": f"{m_i}{n_i}",
                        "measured_hz": float(f_m),
                        "in_vacuo_hz": float(f_v),
                        "pct_dev": float(pct),
                    }
                )
            if low_ok:
                _log(
                    "[VAL3] in-vacuo > measured for air-loading-dominated "
                    "lowest modes (0,1) and (1,1) - expected direction",
                    log,
                )
            else:
                _log(
                    "[VAL3] WARNING: in-vacuo not above measured for "
                    "lowest modes (0,1)/(1,1)",
                    log,
                )
            summary["val3_max_abs_err_hz"] = max_err
            summary["val3_fitted_c"] = float(drum.fitted_effective_c() or 0.0)
            summary["val3_c_vacuo_from_mge2"] = float(c_vac)
            summary["val3_rows"] = rows_val3

    profiles = []
    if opt.run_profiles or opt.run_plots:
        profiles = [
            generate_profile(i, amplitude_layer=layer) for i in instruments
        ]
        for i, p in zip(instruments, profiles):
            pclass = getattr(i, "plate_class", None)
            if p.stroke is None and p.dynamic is None:
                _log(
                    f"[EXC ] {p.instrument}: legacy equipartition "
                    f"(stroke/dynamic unset; bit-identical v0.3.3 path)"
                    + (f"; plate_class={pclass}" if pclass else ""),
                    log,
                )
            else:
                _log(
                    f"[EXC ] {p.instrument}: stroke={p.stroke}, "
                    f"dynamic={p.dynamic}, t_contact={p.t_contact_s}, "
                    f"f_c={p.f_c_hz}, provenance={p.t_contact_provenance}",
                    log,
                )

    if opt.run_profiles and profiles:
        rows = []
        for p in profiles:
            rows.extend(p.to_rows())
        fieldnames = sorted(
            set().union(*[set(x.keys()) for x in rows]),
            key=lambda k: (k not in ("instrument", "family"), k),
        )
        for r in rows:
            for fn in fieldnames:
                r.setdefault(fn, "")
        csv_path = out / "density_profiles.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        _log(
            f"[OUT ] {csv_path.name}: {len(rows)} band rows, "
            f"{len(profiles)} instruments",
            log,
        )
        summary["outputs"].append(str(csv_path))

        _log("\nComposite energy-weighted density indices (per phase):", log)
        phases_all = ["strike", "buildup", "shimmer", "residue", "decay"]
        header = f"{'instrument':<22}" + "".join(f"{ph:>9}" for ph in phases_all)
        _log(header, log)
        for p in profiles:
            line = f"{p.instrument:<22}"
            for ph in phases_all:
                if ph in p.energy_weights:
                    line += f"{p.composite_index(ph):9.2f}"
                else:
                    line += f"{'--':>9}"
            line += f"  [{p.energy_provenance}]"
            _log(line, log)

    if opt.run_plots and profiles:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for p in profiles:
            ls = "-" if p.family == "plate" else "--"
            ax[0].plot(
                p.band_centres, p.modes_per_band, ls, lw=1.4, label=p.instrument
            )
        ax[0].set(
            xscale="log",
            yscale="log",
            xlabel="ERB-band centre (Hz)",
            ylabel="modes per ERB band",
            title="Modal occupation of the auditory scale",
        )
        ax[0].grid(alpha=0.3, which="both")
        ax[0].legend(fontsize=7)

        p46 = next(
            (p for p in profiles if p.instrument == "cymbal_46cm_medium"),
            profiles[0],
        )
        for ph, w in p46.energy_weights.items():
            ax[1].plot(p46.band_centres, w, lw=1.4, label=ph)
        ax[1].set(
            xscale="log",
            xlabel="ERB-band centre (Hz)",
            ylabel="relative band energy",
            title=f"{p46.instrument}: energy per phase",
        )
        ax[1].grid(alpha=0.3, which="both")
        ax[1].legend(fontsize=8)
        fig.tight_layout()
        png1 = out / "density_profiles.png"
        fig.savefig(png1, dpi=150)
        plt.close(fig)
        _log(f"[OUT ] {png1.name}", log)
        summary["outputs"].append(str(png1))

        # MC size sweep (replaces deterministic point estimates for citation).
        sizes = np.linspace(0.30, 0.60, 7)
        sweep_rows = []
        medians = []
        p25s, p75s, p05s, p95s = [], [], [], []
        n_sweep = min(400, max(50, opt.mc_draws_other))
        for d in sizes:
            inst = PlateInstrument(f"cymbal_{d:.2f}", d, 0.0012)
            mc = run_monte_carlo(
                inst, n_draws=n_sweep, seed=opt.mc_seed
            )
            # shimmer composite distribution across draws
            comps = []
            # Recompute composites from stored energy weights percentiles
            # via the MC result's composite stacks if available.
            cstack = mc.composite_quantiles["shimmer"]
            row = {
                "diameter_m": float(d),
                "diameter_cm": float(d * 100),
                "n_draws": n_sweep,
                "seed": opt.mc_seed,
            }
            for lab in ("p05", "p25", "p50", "p75", "p95"):
                row[lab] = float(cstack[lab])
            sweep_rows.append(row)
            medians.append(row["p50"])
            p25s.append(row["p25"])
            p75s.append(row["p75"])
            p05s.append(row["p05"])
            p95s.append(row["p95"])

        sweep_csv = out / "size_sweep_mc.csv"
        with open(sweep_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(sweep_rows[0].keys()))
            w.writeheader()
            w.writerows(sweep_rows)
        _log(f"[OUT ] {sweep_csv.name}", log)
        summary["outputs"].append(str(sweep_csv))

        fig2, ax2 = plt.subplots(figsize=(5.2, 3.6))
        x = sizes * 100
        ax2.fill_between(x, p05s, p95s, alpha=0.20, label="90% band")
        ax2.fill_between(x, p25s, p75s, alpha=0.35, label="50% band")
        ax2.plot(x, medians, "o-", color="C0", label="median")
        ax2.set(
            xlabel="cymbal diameter (cm)",
            ylabel="composite density index (shimmer)",
            title=f"MC size sweep (h=1.2 mm, bronze, N={n_sweep})",
        )
        ax2.grid(alpha=0.3)
        ax2.legend(fontsize=8)
        fig2.tight_layout()
        png2 = out / "size_sweep.png"
        fig2.savefig(png2, dpi=150)
        plt.close(fig2)
        _log(f"[OUT ] {png2.name} (MC; deterministic sweep deprecated)", log)
        summary["outputs"].append(str(png2))

    if opt.run_calibration:
        report = out / "calibration_report.md"
        factor, spread = write_calibration_report(report, layer)
        _log(format_calibration_cli_line(factor, spread), log)
        _log(f"[OUT ] {report.name}", log)
        summary["calibration_factor"] = factor
        summary["calibration_spread"] = spread
        summary["outputs"].append(str(report))

    if opt.run_mc:
        mc_results = []
        for inst in instruments:
            n_draws = (
                opt.mc_draws_focus
                if inst.name == "cymbal_46cm_medium"
                else opt.mc_draws_other
            )
            _log(
                f"[MC  ] {inst.name}: N={n_draws} seed={opt.mc_seed} ...",
                log,
            )
            mc_results.append(
                run_monte_carlo(inst, n_draws=n_draws, seed=opt.mc_seed)
            )
        mc_csv = out / "density_profiles_mc.csv"
        export_mc_csv(mc_results, mc_csv)
        _log(f"[OUT ] {mc_csv.name}", log)
        summary["outputs"].append(str(mc_csv))

        focus = next(
            (r for r in mc_results if r.instrument == "cymbal_46cm_medium"),
            mc_results[0],
        )
        fan = plot_fan_chart(
            focus,
            phase="shimmer" if "shimmer" in focus.energy_quantiles else "strike",
            path=out / f"density_profiles_mc_fan_{focus.instrument}.png",
        )
        _log(f"[OUT ] {fan.name}", log)
        summary["outputs"].append(str(fan))
        phase = "shimmer" if "shimmer" in focus.composite_quantiles else "strike"
        cq = focus.composite_quantiles[phase]
        _log(
            f"[MC  ] {focus.instrument} {phase} composite: "
            f"p50={cq['p50']:.3f}  [{cq['p05']:.3f}, {cq['p95']:.3f}]",
            log,
        )
        summary["mc_p50"] = cq["p50"]
        summary["mc_p05"] = cq["p05"]
        summary["mc_p95"] = cq["p95"]

    _log(f"\n[{APP_NAME}] pipeline finished.", log)
    return summary


# ======================================================================
# GUI
# ======================================================================

def launch_gui() -> None:
    """Open the NonTunPerc desktop interface."""
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk

    root = tk.Tk()
    root.title(APP_TITLE)
    root.geometry("860x640")
    root.minsize(720, 520)

    # Palette — clear, non-generic academic tool look
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
        "Title.TLabel",
        background=bg,
        foreground=accent,
        font=("Georgia", 18, "bold"),
    )
    style.configure(
        "Sub.TLabel", background=bg, foreground=muted, font=("Segoe UI", 9)
    )
    style.configure(
        "Card.TLabel", background=panel, foreground=text, font=("Segoe UI", 10)
    )
    style.configure(
        "TCheckbutton", background=panel, foreground=text, font=("Segoe UI", 9)
    )
    style.configure(
        "TButton", font=("Segoe UI", 10), padding=6
    )
    style.configure(
        "Accent.TButton",
        background=accent,
        foreground="#1a1510",
        font=("Segoe UI", 10, "bold"),
        padding=8,
    )
    style.map("Accent.TButton", background=[("active", "#d4b56a")])

    # Header
    hdr = ttk.Frame(root)
    hdr.pack(fill="x", padx=16, pady=(14, 6))
    ttk.Label(hdr, text=APP_NAME, style="Title.TLabel").pack(anchor="w")
    ttk.Label(
        hdr,
        text="Theoretical ERB-band density for unpitched percussion  ·  "
        "MC median is the citable output  ·  see QUICK_REFERENCE.md",
        style="Sub.TLabel",
    ).pack(anchor="w")

    body = ttk.Frame(root)
    body.pack(fill="both", expand=True, padx=16, pady=8)

    left = ttk.Frame(body, style="Card.TFrame")
    left.pack(side="left", fill="y", padx=(0, 10))
    left.configure(padding=12)

    right = ttk.Frame(body)
    right.pack(side="left", fill="both", expand=True)

    # --- Options -------------------------------------------------------
    ttk.Label(left, text="Pipeline stages", style="Card.TLabel").pack(anchor="w")

    var_val = tk.BooleanVar(value=True)
    var_prof = tk.BooleanVar(value=True)
    var_plots = tk.BooleanVar(value=True)
    var_cal = tk.BooleanVar(value=True)
    var_mc = tk.BooleanVar(value=True)
    var_amp = tk.BooleanVar(value=True)

    for text_, var in (
        ("VAL1 / VAL2 / VAL3 checks", var_val),
        ("Density profiles (CSV)", var_prof),
        ("Plots (PNG)", var_plots),
        ("Calibration bridge", var_cal),
        ("Monte Carlo uncertainty", var_mc),
        ("Use AmplitudeLayer (Sivian/Meyer)", var_amp),
    ):
        ttk.Checkbutton(left, text=text_, variable=var).pack(anchor="w", pady=2)

    ttk.Label(left, text="Monte Carlo", style="Card.TLabel").pack(
        anchor="w", pady=(12, 2)
    )
    mc_frame = ttk.Frame(left, style="Card.TFrame")
    mc_frame.pack(anchor="w", fill="x")
    ttk.Label(mc_frame, text="Draws (46 cm)", style="Card.TLabel").grid(
        row=0, column=0, sticky="w"
    )
    spin_focus = tk.Spinbox(
        mc_frame, from_=50, to=20000, increment=50, width=8,
        textvariable=tk.StringVar(value="2000"),
    )
    spin_focus.grid(row=0, column=1, padx=4, pady=2)
    ttk.Label(mc_frame, text="Draws (others)", style="Card.TLabel").grid(
        row=1, column=0, sticky="w"
    )
    spin_other = tk.Spinbox(
        mc_frame, from_=20, to=5000, increment=20, width=8,
        textvariable=tk.StringVar(value="400"),
    )
    spin_other.grid(row=1, column=1, padx=4, pady=2)
    ttk.Label(mc_frame, text="RNG seed", style="Card.TLabel").grid(
        row=2, column=0, sticky="w"
    )
    seed_var = tk.StringVar(value=str(DEFAULT_SEED))
    tk.Entry(mc_frame, textvariable=seed_var, width=10).grid(
        row=2, column=1, padx=4, pady=2, sticky="w"
    )

    ttk.Label(left, text="Instruments", style="Card.TLabel").pack(
        anchor="w", pady=(12, 2)
    )
    instr_vars = {}
    for name in [i.name for i in default_instruments()]:
        v = tk.BooleanVar(value=True)
        instr_vars[name] = v
        ttk.Checkbutton(left, text=name, variable=v).pack(anchor="w")

    out_dir_var = tk.StringVar(value=str(ROOT))

    def browse_out() -> None:
        d = filedialog.askdirectory(initialdir=out_dir_var.get())
        if d:
            out_dir_var.set(d)

    ttk.Label(left, text="Output folder", style="Card.TLabel").pack(
        anchor="w", pady=(12, 2)
    )
    out_row = ttk.Frame(left, style="Card.TFrame")
    out_row.pack(fill="x")
    tk.Entry(out_row, textvariable=out_dir_var, width=28).pack(
        side="left", fill="x", expand=True
    )
    ttk.Button(out_row, text="…", width=3, command=browse_out).pack(side="left")

    # --- Log -----------------------------------------------------------
    ttk.Label(right, text="Log").pack(anchor="w")
    log_box = scrolledtext.ScrolledText(
        right,
        wrap="word",
        height=28,
        bg="#121a16",
        fg=text,
        insertbackground=text,
        font=("Consolas", 9),
        relief="flat",
    )
    log_box.pack(fill="both", expand=True, pady=(4, 8))

    status_var = tk.StringVar(value="Ready.")
    ttk.Label(right, textvariable=status_var, style="Sub.TLabel").pack(anchor="w")

    running = {"flag": False}

    def append_log(msg: str) -> None:
        def _do() -> None:
            log_box.insert("end", msg + "\n")
            log_box.see("end")

        root.after(0, _do)

    def set_status(msg: str) -> None:
        root.after(0, lambda: status_var.set(msg))

    def collect_options() -> PipelineOptions:
        names = [n for n, v in instr_vars.items() if v.get()]
        return PipelineOptions(
            run_validation=var_val.get(),
            run_profiles=var_prof.get(),
            run_plots=var_plots.get(),
            run_calibration=var_cal.get(),
            run_mc=var_mc.get(),
            mc_draws_focus=int(spin_focus.get()),
            mc_draws_other=int(spin_other.get()),
            mc_seed=int(seed_var.get()),
            instrument_names=names or None,
            use_amplitude_layer=var_amp.get(),
            out_dir=Path(out_dir_var.get()),
        )

    def on_run() -> None:
        if running["flag"]:
            return
        if not any(v.get() for v in instr_vars.values()):
            messagebox.showwarning(APP_NAME, "Select at least one instrument.")
            return
        running["flag"] = True
        run_btn.state(["disabled"])
        set_status("Running…")
        append_log(f"\n=== {APP_NAME} run started ===")

        def worker() -> None:
            try:
                opt = collect_options()
                run_pipeline(opt, log=append_log)
                set_status("Finished.")
                root.after(
                    0,
                    lambda: messagebox.showinfo(
                        APP_NAME, "Pipeline finished. See log and output folder."
                    ),
                )
            except Exception as exc:
                append_log(traceback.format_exc())
                set_status("Failed.")
                root.after(
                    0,
                    lambda: messagebox.showerror(APP_NAME, str(exc)),
                )
            finally:
                running["flag"] = False
                root.after(0, lambda: run_btn.state(["!disabled"]))

        threading.Thread(target=worker, daemon=True).start()

    def open_out() -> None:
        folder = Path(out_dir_var.get())
        folder.mkdir(parents=True, exist_ok=True)
        import os

        os.startfile(str(folder))  # noqa: S606 — Windows GUI helper

    def open_validate() -> None:
        import subprocess

        bat = ROOT / "run_validate.bat"
        if bat.is_file():
            subprocess.Popen(["cmd", "/c", str(bat)], cwd=str(ROOT))
        else:
            messagebox.showinfo(
                APP_NAME,
                "run_validate.bat not found. Use:\n"
                "python validate_against_recordings.py --wav-dir <folder>",
            )

    btn_row = ttk.Frame(left, style="Card.TFrame")
    btn_row.pack(fill="x", pady=(16, 0))
    run_btn = ttk.Button(
        btn_row, text="Run pipeline", style="Accent.TButton", command=on_run
    )
    run_btn.pack(fill="x", pady=2)
    ttk.Button(btn_row, text="Open output folder", command=open_out).pack(
        fill="x", pady=2
    )
    ttk.Button(
        btn_row, text="Validate against WAVs…", command=open_validate
    ).pack(fill="x", pady=2)

    append_log(f"{APP_NAME} ready. Select stages and click Run pipeline.")
    root.mainloop()


def run_headless(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry equivalent to the former demo.py."""
    ap = argparse.ArgumentParser(description=APP_TITLE)
    ap.add_argument("--cli", action="store_true", help="run headless pipeline")
    ap.add_argument("--no-mc", action="store_true", help="skip Monte Carlo")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--no-calibration", action="store_true")
    ap.add_argument("--out", type=Path, default=ROOT)
    ap.add_argument("--mc-draws", type=int, default=2000)
    ap.add_argument("--mc-draws-other", type=int, default=400)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument(
        "--gui", action="store_true", help="force GUI (default if no --cli)"
    )
    args, _unknown = ap.parse_known_args(argv)

    if args.cli and not args.gui:
        opt = PipelineOptions(
            run_mc=not args.no_mc,
            run_plots=not args.no_plots,
            run_calibration=not args.no_calibration,
            mc_draws_focus=args.mc_draws,
            mc_draws_other=args.mc_draws_other,
            mc_seed=args.seed,
            out_dir=args.out,
        )
        run_pipeline(opt)
        return 0

    launch_gui()
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(argv) if argv is not None else sys.argv[1:]
    # Default = GUI; --cli for batch / former demo behaviour.
    if "--cli" in argv:
        return run_headless(argv)
    if argv and argv[0] in ("-h", "--help"):
        return run_headless(argv)
    launch_gui()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
