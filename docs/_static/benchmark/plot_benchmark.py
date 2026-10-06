# Wall time of a neighbor search, adaptive CNA and q6 for one million atoms,
# for pyscal 3.3, pyscal 4, OVITO and freud with their default settings:
# pyscal 4, OVITO and freud on all 14 threads, pyscal 3.3 (no threading) on one.
#
# Data: benchmark.json next to this script, taken from a benchmark that ran
# every code one after the other on one machine (details in the file).
# Each time is the minimum over repeated calls. OVITO has no q6 routine and
# freud no CNA in the comparison, so those bars are missing.
#
# Output: benchmark.png (next to this script)
import json
import os

import numpy as np
from pychromatic import Multiplot, Palette

HERE = os.path.dirname(os.path.abspath(__file__))
DARK = "#363636"
p = Palette("tableau10")

data = json.load(open(os.path.join(HERE, "benchmark.json")))
times = {(r["code"], r["threads"], r["operation"]): r["time_s"] for r in data["rows"]}

operations = ["neighbor search", "adaptive CNA", "q6 with neighbor search"]
labels = ["neighbor\nsearch", "adaptive\nCNA", "$q_6$ with\nneighbor search"]
# (code, threads, legend, colour)
series = [
    ("pyscal 3.3", 1, "pyscal 3.3", p.grey.hex),
    ("pyscal 4", 14, "pyscal 4", p.blue.hex),
    ("OVITO 3.16", 14, "OVITO 3.16", p.red.hex),
    ("freud 3.4", 14, "freud 3.4", p.brown.hex),
]

mp = Multiplot(width=510, ratio=0.42)
ax = mp[0, 0]
w = 0.2
for k, (code, threads, legend, colour) in enumerate(series):
    x = np.arange(len(operations)) + (k - 1.5) * w
    t = [times.get((code, threads, op), np.nan) for op in operations]
    ax.bar(x, t, width=w, color=colour, edgecolor=DARK, lw=0.8, label=legend, zorder=3)
    for xi, ti in zip(x, t):
        if np.isfinite(ti):
            text = f"{ti:.0f}" if ti >= 10 else f"{ti:.1f}" if ti >= 1 else f"{ti:.2f}"
            ax.text(xi, ti * 1.25, text, ha="center", va="bottom", fontsize=7, color=DARK,
                    rotation=90)

ax.set_yscale("log")
ax.set_ylim(0.05, 600)
ax.set_xlim(-0.55, len(operations) - 0.45)
ax.set_xticks(np.arange(len(operations)), labels, fontsize=9)
ax.set_ylabel("Wall time  (s)", fontsize=10)
ax.tick_params(axis="y", labelsize=9)
ax.text(0.99, 0.97, "1 000 188 fcc atoms, 12 neighbors each\nno bar: no such routine in OVITO or freud",
        transform=ax.transAxes, ha="right", va="top", fontsize=8, color=DARK)
mp.fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=8, ncol=4,
              loc="upper center", bbox_to_anchor=(0.5, -0.02), columnspacing=1.5)
mp.fig.text(0.5, -0.11, "14 threads for pyscal 4, OVITO and freud. pyscal 3.3 has no threading and runs on one.",
            ha="center", va="top", fontsize=8, color=DARK)
mp.fig.savefig(os.path.join(HERE, "benchmark.png"), dpi=300, bbox_inches="tight")
print("wrote benchmark.png")
