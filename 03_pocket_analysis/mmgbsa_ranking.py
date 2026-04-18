"""
CRYPTAD — MM-GBSA preliminary ranking chart
All 4 confirmed cryptic pockets (partial results, updated as jobs complete).
"""

import os, statistics
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Data collection ──────────────────────────────────────────────────────────
BASE = "/Volumes/PortableSSD/untitled folder/CRYPTAD/04_virtual_screening/mmgbsa"

POCKETS = {
    "S1_site688":  ("BIN1 BAR\nsite688",  "#1A6FA8"),   # blue (primary lead)
    "S1_site680":  ("BIN1 BAR\nsite680",  "#5BA3D9"),   # light blue
    "S1_site1813": ("BIN1 BAR\nsite1813", "#D95F3B"),   # warm red-orange
    "S3_site473":  ("PICALM ANTH\nsite473",  "#2E9E9A"),  # teal
}

results = []
for pocket_id, (label, color) in POCKETS.items():
    pocket_dir = f"{BASE}/{pocket_id}"
    for lig in os.listdir(pocket_dir):
        csv_path = f"{pocket_dir}/{lig}/mmgbsa/FINAL_RESULTS_MMPBSA.csv"
        if not os.path.exists(csv_path):
            continue
        delta_totals = []
        in_delta = False
        with open(csv_path) as f:
            for line in f:
                line = line.strip()
                if line == "Delta Energy Terms":
                    in_delta = True; continue
                if in_delta and line.startswith("Frame #"): continue
                if in_delta and line and line[0].isdigit():
                    parts = line.split(",")
                    try:
                        total = float(parts[-1])
                        if abs(total) < 5000:
                            delta_totals.append(total)
                    except:
                        pass
        if delta_totals:
            mean_dg = statistics.mean(delta_totals)
            sd = statistics.stdev(delta_totals) if len(delta_totals) > 1 else 0
            results.append({
                "dg": mean_dg, "sd": sd,
                "lig": lig, "pocket": pocket_id,
                "label": label, "color": color,
            })

results.sort(key=lambda r: r["dg"])
n = len(results)

# ── Figure setup ─────────────────────────────────────────────────────────────
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "pdf.fonttype": 42,   # editable text in Illustrator
    "svg.fonttype": "none",
})

fig, ax = plt.subplots(figsize=(14, 6))
fig.patch.set_facecolor("white")
ax.set_facecolor("#F8F8F8")

# ── Bars ─────────────────────────────────────────────────────────────────────
x = np.arange(n)
bar_colors = [r["color"] for r in results]
bar_alphas = [1.0 if i < 10 else 0.75 for i in range(n)]

bars = ax.bar(
    x,
    [r["dg"] for r in results],
    color=bar_colors,
    alpha=0.9,
    width=0.75,
    zorder=3,
    linewidth=0,
)

# Error bars
ax.errorbar(
    x,
    [r["dg"] for r in results],
    yerr=[r["sd"] for r in results],
    fmt="none",
    ecolor="black",
    elinewidth=0.7,
    capsize=2.5,
    capthick=0.7,
    zorder=4,
    alpha=0.6,
)

# ── Highlight top 3 per pocket ───────────────────────────────────────────────
top_ligs = set()
for pocket_id in POCKETS:
    top3 = sorted([r for r in results if r["pocket"] == pocket_id], key=lambda r: r["dg"])[:3]
    top_ligs.update(r["lig"] for r in top3)

for i, r in enumerate(results):
    if r["lig"] in top_ligs:
        # Bold outline on top hits
        bars[i].set_linewidth(1.2)
        bars[i].set_edgecolor("black")
        # Compound label (short ID)
        short = r["lig"].replace("CHEMBL", "")
        ax.text(
            i, r["dg"] - r["sd"] - 0.6,
            short,
            ha="center", va="top",
            fontsize=6.5, fontweight="bold",
            color="black", rotation=90,
            zorder=5,
        )

# ── Reference line ───────────────────────────────────────────────────────────
ax.axhline(y=-20, color="#888888", linewidth=0.7, linestyle="--", zorder=2, alpha=0.8)
ax.text(n - 0.3, -20.5, "−20 kcal/mol", ha="right", va="top",
        fontsize=7.5, color="#666666", style="italic")

# ── Axes styling ─────────────────────────────────────────────────────────────
ax.set_xlim(-0.8, n - 0.2)
ax.set_ylim(min(r["dg"] - r["sd"] for r in results) - 5, 2)
ax.set_ylabel("ΔG$_{bind}$ (kcal/mol, MM-GBSA)", fontsize=10, labelpad=8)
ax.set_xticks([])
ax.yaxis.grid(True, color="white", linewidth=0.8, zorder=1)
ax.set_axisbelow(True)

# Spine cleanup
for spine in ["top", "right", "bottom"]:
    ax.spines[spine].set_visible(False)
ax.spines["left"].set_color("#AAAAAA")

# ── Pocket separators + span annotations ─────────────────────────────────────
pocket_xs = {pid: [i for i, r in enumerate(results) if r["pocket"] == pid]
             for pid in POCKETS}
# Vertical separators between pockets
prev_pid = None
for pid in POCKETS:
    if prev_pid and pocket_xs[pid]:
        bx = pocket_xs[pid][0] - 0.5
        ax.axvline(x=bx, color="#CCCCCC", linewidth=1.0, linestyle=":", zorder=2)
    prev_pid = pid

# Pocket label annotations
pocket_labels = {
    "S1_site688":  "BIN1 BAR — site688 [primary]",
    "S1_site680":  "BIN1 BAR — site680",
    "S1_site1813": "BIN1 BAR — site1813",
    "S3_site473":  "PICALM ANTH — site473",
}
for pid, (_, color) in POCKETS.items():
    xs = pocket_xs[pid]
    if not xs:
        continue
    n_poc = len(xs)
    ax.annotate(
        "", xy=(max(xs) + 0.4, 0.5), xytext=(min(xs) - 0.4, 0.5),
        arrowprops=dict(arrowstyle="|-|", color=color, lw=1.5),
    )
    ax.text(np.mean(xs), 1.3,
            f"{pocket_labels[pid]}  (n={n_poc})",
            ha="center", va="bottom", fontsize=8.0,
            color=color, fontweight="bold")

# ── Title and caption ────────────────────────────────────────────────────────
ax.set_title(
    f"CRYPTAD — MM-GBSA Rescoring  ({n} / 201 compounds complete)",
    fontsize=11, fontweight="bold", pad=14, loc="left",
)
fig.text(
    0.01, 0.01,
    "GB5, AMBER99SB-ILDN, 1 ns production MD, 101 frames per compound  |  "
    "Error bars = ±1 SD across frames  |  Bold outline = top-3 per pocket  |  "
    "[primary] = primary lead pocket (site688)",
    fontsize=7, color="#888888", va="bottom",
)

# ── Legend ───────────────────────────────────────────────────────────────────
legend_patches = [
    mpatches.Patch(color=POCKETS["S1_site688"][1],  label="BIN1 BAR site688 [primary lead] (tau axis)"),
    mpatches.Patch(color=POCKETS["S1_site680"][1],  label="BIN1 BAR site680 (tau axis)"),
    mpatches.Patch(color=POCKETS["S1_site1813"][1], label="BIN1 BAR site1813 (tau axis)"),
    mpatches.Patch(color=POCKETS["S3_site473"][1],  label="PICALM ANTH site473 (Aβ axis)"),
]
ax.legend(handles=legend_patches, loc="lower right", frameon=True,
          framealpha=0.9, edgecolor="#CCCCCC", fontsize=8.5)

plt.tight_layout(rect=[0, 0.04, 1, 1])

OUTDIR = "/Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis"
plt.savefig(f"{OUTDIR}/mmgbsa_ranking.png", dpi=300, bbox_inches="tight")
plt.savefig(f"{OUTDIR}/mmgbsa_ranking.pdf", bbox_inches="tight")
print(f"Saved to {OUTDIR}/mmgbsa_ranking.png (.pdf)")
plt.show()
