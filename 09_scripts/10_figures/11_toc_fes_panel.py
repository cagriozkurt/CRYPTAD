"""
CRYPTAD — TOC Graphic: Styled PICALM S3 FES Panel

Reads fes_200ns.dat for the PICALM ANTH (S3) system and produces a
publication-quality, aestheticised free energy surface panel intended
for use as a component in the graphical abstract.

Design choices:
  - Dark navy background (#0a1628) for cinematic TOC aesthetic
  - Custom colormap: deep indigo → teal → gold → near-white
  - Filled contours + sparse thin white contour lines
  - Global minimum marked with a white star
  - Cryptic-pocket open region (cv_site4 > 3.0 nm) annotated subtly
  - CV axes labelled minimally in clean white sans-serif
  - No title, no colourbar — panel is meant to be merged in Illustrator
  - Separate version WITH colourbar for reference

Output (06_figures/publication/):
  TOC_fes_panel.pdf / .png          — clean panel, no colourbar
  TOC_fes_panel_ref.pdf / .png      — same with colourbar for reference

Usage:
  python3 09_scripts/10_figures/11_toc_fes_panel.py
  python3 09_scripts/10_figures/11_toc_fes_panel.py --project-root /path
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

# ── constants ────────────────────────────────────────────────────────────────
KJ_TO_KCAL       = 1 / 4.184
TRUNCATE_KCAL    = 14.0          # clip FES above this value (kcal/mol)
POCKET_OPEN_NM   = 3.0           # cv_site4 threshold for open state (nm)
BG_COLOR         = "#0a1628"     # dark navy background
STAR_COLOR       = "#ffffff"
DASHED_COLOR     = "#7ecfcf"     # teal dashed threshold line

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]


def load_fes(dat_path: Path):
    """
    Parse a PLUMED 2D FES .dat file using np.loadtxt (handles blank lines
    and comment lines automatically).
    Returns (cv1, cv2, fes_kcal) each as 2-D arrays shaped (n1, n2).
    Unsampled PLUMED sentinel values (~166.87 kJ/mol) are set to NaN.
    """
    data     = np.loadtxt(dat_path, comments="#")   # shape (N, 5)
    cv1_vals = np.unique(data[:, 0])
    cv2_vals = np.unique(data[:, 1])
    n1, n2   = len(cv1_vals), len(cv2_vals)

    fes_kj = data[:, 2].reshape(n1, n2)

    # Mask PLUMED sentinel (unsampled) values — they sit at ~166.87 kJ/mol
    SENTINEL_THRESHOLD_KJ = 160.0
    fes_kj = fes_kj.astype(float)
    fes_kj[fes_kj > SENTINEL_THRESHOLD_KJ] = np.nan

    # Zero at sampled global minimum
    fes_kj = fes_kj - np.nanmin(fes_kj)

    # Convert to kcal/mol; clip at TRUNCATE_KCAL
    fes_kcal = fes_kj * KJ_TO_KCAL
    fes_kcal = np.where(np.isnan(fes_kcal), np.nan,
                        np.clip(fes_kcal, 0, TRUNCATE_KCAL))

    cv1_grid, cv2_grid = np.meshgrid(cv1_vals, cv2_vals, indexing="ij")
    return cv1_grid, cv2_grid, fes_kcal


def make_cmap():
    """
    Inverted colormap: low energy glows warm against a dark background.
      0 kcal/mol  (global minimum) = bright amber/gold
      2–5         = warm yellow → lime
      6–9         = teal → cyan-blue
      10–14       = dark slate blue → near-navy (blends into background)
    Unsampled NaN regions show as the navy background colour.
    """
    colors = [
        (1.00, 0.85, 0.20),   # bright amber/gold  (0 — minimum, glowing)
        (0.85, 0.95, 0.35),   # yellow-lime
        (0.30, 0.88, 0.72),   # bright teal
        (0.08, 0.62, 0.78),   # medium blue-teal
        (0.06, 0.35, 0.58),   # dark slate blue
        (0.04, 0.10, 0.22),   # near-navy  (14 — high energy, fades to bg)
    ]
    cmap = LinearSegmentedColormap.from_list("toc_fes", colors, N=512)
    cmap.set_bad(color=BG_COLOR)   # NaN (unsampled) → navy background
    return cmap


def plot_fes(cv1, cv2, fes, ax, show_cbar=False):
    cmap    = make_cmap()
    levels  = np.linspace(0, TRUNCATE_KCAL, 256)
    c_lines = np.arange(2, TRUNCATE_KCAL + 1, 2)   # every 2 kcal/mol

    # Filled contours
    cf = ax.contourf(cv2, cv1, fes, levels=levels, cmap=cmap, extend="max")

    # Thin white contour lines
    ax.contour(cv2, cv1, fes, levels=c_lines,
               colors="white", linewidths=0.4, alpha=0.35)

    # Dashed line: cryptic pocket open threshold
    ax.axvline(POCKET_OPEN_NM, color=DASHED_COLOR, linewidth=1.2,
               linestyle="--", alpha=0.8)

    # Global minimum star
    min_idx  = np.unravel_index(np.argmin(fes), fes.shape)
    min_cv2  = cv2[min_idx]
    min_cv1  = cv1[min_idx]
    ax.plot(min_cv2, min_cv1, "*", color=STAR_COLOR,
            markersize=11, markeredgewidth=0.6,
            markeredgecolor="#aaaaaa", zorder=10)

    # Axis styling
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors="white", labelsize=8, length=3, width=0.6)
    for spine in ax.spines.values():
        spine.set_edgecolor("#334466")
        spine.set_linewidth(0.8)

    ax.set_xlabel("site473 aperture, CV₂ (nm)", color="white",
                  fontsize=9, labelpad=5)
    ax.set_ylabel("PIP₂-loop, CV₁ (nm)", color="white",
                  fontsize=9, labelpad=5)

    # Subtle annotation for open-state region
    ax.text(POCKET_OPEN_NM + 0.06, cv1.max() * 0.95,
            "open\nstate", color=DASHED_COLOR, fontsize=7.5,
            va="top", ha="left",
            fontfamily="DejaVu Sans")

    if show_cbar:
        cbar = plt.colorbar(cf, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label("ΔG (kcal mol⁻¹)", color="white", fontsize=8)
        cbar.ax.yaxis.set_tick_params(color="white", labelsize=7)
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")
        cbar.outline.set_edgecolor("#334466")


def save_figure(fig, stem: Path):
    for ext in ("pdf", "png"):
        out = stem.with_suffix(f".{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight",
                    facecolor=BG_COLOR, edgecolor="none")
        print(f"Saved: {out.name}")


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD TOC — Styled PICALM S3 FES panel")
    parser.add_argument("--project-root", type=Path, default=_default_root)
    args = parser.parse_args()
    root = args.project_root.resolve()

    dat   = root / "02_md_simulations/PICALM/metadynamics/S3_ANTH/fes_200ns.dat"
    outdir = root / "06_figures/publication"
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading FES: {dat}")
    cv1, cv2, fes = load_fes(dat)
    print(f"  Grid: {cv1.shape}  FES range: {np.nanmin(fes):.2f}–{np.nanmax(fes):.2f} kcal/mol  NaN: {np.isnan(fes).sum()}")

    # ── clean panel (no colourbar) ────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(4.5, 3.8),
                           facecolor=BG_COLOR)
    fig.patch.set_facecolor(BG_COLOR)
    plot_fes(cv1, cv2, fes, ax, show_cbar=False)
    save_figure(fig, outdir / "TOC_fes_panel")
    plt.close()

    # ── reference panel (with colourbar) ─────────────────────────────────────
    fig, ax = plt.subplots(figsize=(5.2, 3.8),
                           facecolor=BG_COLOR)
    fig.patch.set_facecolor(BG_COLOR)
    plot_fes(cv1, cv2, fes, ax, show_cbar=True)
    save_figure(fig, outdir / "TOC_fes_panel_ref")
    plt.close()

    print("Done.")


if __name__ == "__main__":
    main()
