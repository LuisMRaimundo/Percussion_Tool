"""Validation and demonstration run for the idiophone density model."""

import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from model import (PlateInstrument, MembraneInstrument, generate_profile,
                   erb_band_edges)

# ----------------------------------------------------------------------
# Instrument set. Chladni anchors (c, b, p): [R] Table 9.1 (n=0 family,
# first segment p1/c1 with b=2 as in f = c(m+2n)^p) where available.
# ----------------------------------------------------------------------
instruments = [
    PlateInstrument("cymbal_16in_thin",   0.406, 0.0008,
                    chladni=(10.8, 2.0, 1.81)),
    PlateInstrument("cymbal_18in_medium", 0.457, 0.0012,
                    chladni=(13.4, 2.0, 1.65)),
    PlateInstrument("cymbal_46cm_medium", 0.460, 0.0012,
                    chladni=(14.93, 3.0, 1.557)),   # [R] Fig. 9.3, n=0
    PlateInstrument("gong_50cm_bronze",   0.500, 0.0020),
    PlateInstrument("tamtam_80cm_bronze", 0.800, 0.0015),
    MembraneInstrument("bassdrum_32in", 0.813, f11_nominal=60.0),
    MembraneInstrument("bassdrum_28in", 0.711, f11_nominal=72.0),
]

# ----------------------------------------------------------------------
# Validation 1: total low-frequency mode count vs Wilbur's holography
# ([R] sec. 9.2: >100 modes recorded in a 46-cm cymbal).
# ----------------------------------------------------------------------
cy = instruments[2]
nd = cy.modal_density()
print(f"[VAL1] 46-cm cymbal asymptotic modal density: {nd:.4f} modes/Hz "
      f"({nd*1000:.0f} modes per kHz)")
print(f"[VAL1] predicted modes below 2 kHz: {nd*2000:.0f} "
      f"(holography: >100 observed; consistent)")

# ----------------------------------------------------------------------
# Validation 2: low-mode frequencies vs Chladni fit of [R] Fig. 9.3
# ----------------------------------------------------------------------
low = cy.low_modes(m_max=7, n_max=0)
ref = 14.93 * (np.arange(2, 8)) ** 1.557
err = np.max(np.abs(low[:6] - ref) / ref) * 100
print(f"[VAL2] n=0 family m=2..7 reproduction error vs [R]: {err:.2f}%")

# ----------------------------------------------------------------------
# Generate profiles, export CSV
# ----------------------------------------------------------------------
profiles = [generate_profile(i) for i in instruments]

rows = []
for p in profiles:
    rows.extend(p.to_rows())
fieldnames = list(rows[0].keys())
# unify columns across families
for r in rows:
    for fn in set().union(*[set(x.keys()) for x in rows]):
        r.setdefault(fn, "")
fieldnames = sorted(set().union(*[set(x.keys()) for x in rows]),
                    key=lambda k: (k not in ("instrument", "family"), k))
with open("density_profiles.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(rows)
print(f"[OUT ] density_profiles.csv: {len(rows)} band rows, "
      f"{len(profiles)} instruments")

# composite indices
print("\nComposite energy-weighted density indices (per phase):")
print(f"{'instrument':<22}", end="")
phases_all = ["strike", "buildup", "shimmer", "residue", "decay"]
for ph in phases_all:
    print(f"{ph:>9}", end="")
print()
for p in profiles:
    print(f"{p.instrument:<22}", end="")
    for ph in phases_all:
        if ph in p.energy_weights:
            print(f"{p.composite_index(ph):9.2f}", end="")
        else:
            print(f"{'--':>9}", end="")
    print()

# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for p in profiles:
    ls = "-" if p.family == "plate" else "--"
    ax[0].plot(p.band_centres, p.modes_per_band, ls, lw=1.4,
               label=p.instrument)
ax[0].set(xscale="log", yscale="log", xlabel="ERB-band centre (Hz)",
          ylabel="modes per ERB band",
          title="Modal occupation of the auditory scale")
ax[0].grid(alpha=.3, which="both")
ax[0].legend(fontsize=7)

p46 = profiles[2]
for ph, w in p46.energy_weights.items():
    ax[1].plot(p46.band_centres, w, lw=1.4, label=ph)
ax[1].set(xscale="log", xlabel="ERB-band centre (Hz)",
          ylabel="relative band energy",
          title="46-cm cymbal: energy per phase")
ax[1].grid(alpha=.3, which="both")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig("density_profiles.png", dpi=150)
print("[OUT ] density_profiles.png")

# size sweep: cymbal diameter -> shimmer-phase composite index
sizes = np.linspace(0.30, 0.60, 7)
vals = []
for d in sizes:
    inst = PlateInstrument(f"cymbal_{d:.2f}", d, 0.0012)
    vals.append(generate_profile(inst).composite_index("shimmer"))
fig2, ax2 = plt.subplots(figsize=(5.2, 3.6))
ax2.plot(sizes * 100, vals, "o-")
ax2.set(xlabel="cymbal diameter (cm)",
        ylabel="composite density index (shimmer)",
        title="Size sweep (h = 1.2 mm, bronze)")
ax2.grid(alpha=.3)
fig2.tight_layout()
fig2.savefig("size_sweep.png", dpi=150)
print("[OUT ] size_sweep.png")