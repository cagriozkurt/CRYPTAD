"""
Figure 06 — Ibufenac COM drift: 6-replica time series.

Reads ibufenac_drift.csv (rep, frame, drift_A) and plots one panel per
replica in a 2×3 grid. Each panel shows COM drift (Å) vs simulation time
(ns), with a horizontal dashed line at the engagement cutoff (15 Å).
Engaged replicas (mean drift < 15 Å) are coloured teal; partial/dissociated
are grey, matching the colour coding used in the manuscript.
% engaged (frames below cutoff) is shown in each panel subtitle.
A shared legend is placed below the panels.

Output:
  06_figures/fig_ibufenac_drift.png                       (600 dpi)
  06_figures/publication/Figure_06_ibufenac_drift.png     (600 dpi)
  06_figures/publication/Figure_06_ibufenac_drift.tif     (600 dpi, LZW)

Usage:
  python3 09_scripts/10_figures/10_ibufenac_drift_figure.py
  python3 09_scripts/10_figures/10_ibufenac_drift_figure.py --project-root /path/to/CRYPTAD
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from PIL import Image

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

FRAME_DT_NS    = 0.2          # 200 ps per saved frame
ENGAGED_CUTOFF = 15.0         # Å — COM drift threshold for engagement
DPI            = 600

# Engagement classification (from complex MD analysis)
# 4/6 replicas engaged: rep1, rep2, rep5, rep6
ENGAGED = {"rep1", "rep2", "rep5", "rep6"}

COLOR_ENGAGED     = "#2a9d8f"   # teal
COLOR_PARTIAL     = "#e9c46a"   # amber — mean drift borderline
COLOR_DISSOCIATED = "#adb5bd"   # grey

REP_LABELS = {
    "rep1": "Rep 1",
    "rep2": "Rep 2",
    "rep3": "Rep 3",
    "rep4": "Rep 4",
    "rep5": "Rep 5",
    "rep6": "Rep 6",
}

# Mean drift determines colour: engaged < 5 Å, partial 5–10 Å, else dissociated
def rep_color(mean_drift: float, rep: str) -> str:
    if rep in ENGAGED:
        return COLOR_ENGAGED
    if mean_drift < 10.0:
        return COLOR_PARTIAL
    return COLOR_DISSOCIATED


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Ibufenac COM drift (Figure 06)")
    parser.add_argument("--project-root", default=str(_default_root))
    args = parser.parse_args()

    root    = Path(args.project_root).resolve()
    csv     = root / "04_virtual_screening" / "complex_md_results" / "ibufenac_drift.csv"
    fig_dir = root / "06_figures"
    pub_dir = fig_dir / "publication"
    pub_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv)
    reps = sorted(df["rep"].unique())   # rep1 … rep6

    fig, axes = plt.subplots(2, 3, figsize=(7, 4.67), sharey=True)
    axes_flat = axes.flatten()

    plt.rcParams.update({
        "font.family":  "sans-serif",
        "font.size":    10,
        "axes.linewidth": 0.8,
    })

    for ax, rep in zip(axes_flat, reps):
        sub   = df[df["rep"] == rep].sort_values("frame")
        time  = sub["frame"].values * FRAME_DT_NS
        drift = sub["drift_A"].values
        mean_d = float(drift.mean())
        pct_engaged = float((drift < ENGAGED_CUTOFF).mean() * 100)

        color = rep_color(mean_d, rep)
        status = "engaged" if rep in ENGAGED else ("partial" if mean_d < 10 else "dissociated")

        # Color trace: teal below cutoff, grey above cutoff
        drift_below = np.ma.masked_where(drift >= ENGAGED_CUTOFF, drift)
        drift_above = np.ma.masked_where(drift <  ENGAGED_CUTOFF, drift)
        ax.plot(time, drift_below, color=COLOR_ENGAGED,     linewidth=0.8, alpha=0.85)
        ax.plot(time, drift_above, color=COLOR_DISSOCIATED, linewidth=0.8, alpha=0.85)

        cutoff_line = ax.axhline(ENGAGED_CUTOFF, color="#e76f51", linewidth=1.0,
                                 linestyle="--", alpha=0.8)

        # Shade engaged region
        ax.axhspan(0, ENGAGED_CUTOFF, color=COLOR_ENGAGED, alpha=0.07)

        # Title: rep label + % engaged subtitle
        ax.set_title(f"{REP_LABELS[rep]}\n({pct_engaged:.0f}% engaged)",
                     fontsize=10, fontweight="bold", pad=4)
        ax.set_xlabel("Time (ns)", fontsize=9)
        ax.set_ylabel("COM Drift (Å)", fontsize=9)
        ax.set_xlim(0, time.max())

        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.tick_params(which="both", direction="in", top=True, right=True,
                       labelsize=8)

        # Status label top-right
        ax.text(0.97, 0.95, status, transform=ax.transAxes,
                ha="right", va="top", fontsize=8,
                color=color if status != "dissociated" else "#6c757d",
                fontweight="bold")

    # Build legend once after loop (trace colours + cutoff line)
    import matplotlib.lines as mlines
    legend_handles = [
        mlines.Line2D([], [], color=COLOR_ENGAGED,     linewidth=1.5,
                      label=f"Drift < {ENGAGED_CUTOFF:.0f} Å (engaged)"),
        mlines.Line2D([], [], color=COLOR_DISSOCIATED, linewidth=1.5,
                      label=f"Drift ≥ {ENGAGED_CUTOFF:.0f} Å (dissociated)"),
        mlines.Line2D([], [], color="#e76f51", linewidth=1.0, linestyle="--",
                      alpha=0.8, label=f"{ENGAGED_CUTOFF:.0f} Å cutoff"),
    ]

    # Set shared y-axis range once after loop, using global data max
    axes_flat[0].set_ylim(0, df["drift_A"].max() * 1.05)

    # Legend below the panels
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.legend(handles=legend_handles, loc="lower center",
               ncol=len(legend_handles), fontsize=9,
               frameon=True, framealpha=0.9,
               bbox_to_anchor=(0.5, 0.01))

    # Save PNG
    out_png = fig_dir / "fig_ibufenac_drift.png"
    fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
    print(f"Saved: {out_png.relative_to(root)}")

    pub_png = pub_dir / "Figure_06_ibufenac_drift.png"
    fig.savefig(pub_png, dpi=DPI, bbox_inches="tight")
    print(f"Saved: {pub_png.relative_to(root)}")

    plt.close(fig)

    # Save TIFF
    pub_tif = pub_dir / "Figure_06_ibufenac_drift.tif"
    img = Image.open(pub_png).convert("RGB")
    img.save(pub_tif, dpi=(DPI, DPI), compression="tiff_lzw")
    print(f"Saved: {pub_tif.relative_to(root)}  ({img.width}×{img.height} px @ {DPI} dpi)")


if __name__ == "__main__":
    main()
