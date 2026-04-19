"""
CRYPTAD — TOC Graphic: Allosteric coupling scatter (PICALM site473 vs PIP2-helix)

Reproduces the r = +0.718 scatter plot from the conformational coupling analysis,
stripped down for use as a TOC graphic component:
  - No axis labels or tick labels
  - Points coloured by simulation time (dark-blue → amber gradient)
  - Clean white regression line
  - Single r-value annotation, bottom-right
  - Dark navy background to match TOC_fes_panel

Data sources:
  02_md_simulations/PICALM/metadynamics/S3_ANTH/COLVAR    (cv_site4)
  02_md_simulations/PICALM/metadynamics/S3_ANTH/meta.tpr  (topology)
  02_md_simulations/PICALM/metadynamics/S3_ANTH/meta.xtc  (trajectory)

Cache: 03_pocket_analysis/conformational/pip2_scatter_cache.npz
  Written on first run; subsequent runs load from cache (~instant).

Output (06_figures/publication/):
  TOC_scatter_pip2.pdf / .png          — clean, no labels (for Illustrator)
  TOC_scatter_pip2_ref.pdf / .png      — same with axes labels (reference)

Usage:
  python3 09_scripts/10_figures/12_toc_scatter_pip2.py
  python3 09_scripts/10_figures/12_toc_scatter_pip2.py --project-root /path
"""

import argparse
import warnings
import logging
from pathlib import Path

warnings.filterwarnings("ignore")
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

BG_COLOR   = "#0a1628"
PT_ALPHA   = 0.45
PT_SIZE    = 6
MAX_PTS    = 2000       # subsample for rendering speed


def time_cmap():
    """Blue → teal → amber: simulation time gradient for scatter points."""
    colors = [
        (0.10, 0.30, 0.65),   # early frames — deep blue
        (0.08, 0.62, 0.72),   # mid — teal
        (0.95, 0.72, 0.15),   # late frames — amber
    ]
    return LinearSegmentedColormap.from_list("toc_time", colors, N=256)


def compute_pip2_data(root: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load or compute (cv_site4, pip2_rmsd_A, time_ns).
    Saves cache to 03_pocket_analysis/conformational/pip2_scatter_cache.npz.
    """
    cache = root / "03_pocket_analysis/conformational/pip2_scatter_cache.npz"

    if cache.exists():
        print(f"Loading cache: {cache.name}")
        d = np.load(cache)
        return d["cv_site4"], d["pip2_rmsd_A"], d["time_ns"]

    print("Cache not found — computing from trajectory (this may take a few minutes)...")
    import MDAnalysis as mda

    tpr = root / "02_md_simulations/PICALM/metadynamics/S3_ANTH/meta.tpr"
    xtc = root / "02_md_simulations/PICALM/metadynamics/S3_ANTH/meta.xtc"
    col = root / "02_md_simulations/PICALM/metadynamics/S3_ANTH/COLVAR"

    # Load COLVAR: col 0 = time (ps), col 2 = cv_site4 (nm)
    colvar      = np.loadtxt(str(col), comments="#")
    u           = mda.Universe(str(tpr), str(xtc))
    n_frames    = len(u.trajectory)

    # COLVAR written every 2 ps; trajectory every 200 ps → stride 10
    cv_site4   = colvar[::10, 2][:n_frames]
    time_ns    = colvar[::10, 0][:n_frames] / 1000.0

    # PIP2-helix RMSD: align core (resids 50–200), measure helix (resids 1–22)
    anth_core  = u.select_atoms("protein and name CA and resid 50:200")
    pip2_helix = u.select_atoms("protein and name CA and resid 1:22")

    # Reference = first frame
    u.trajectory[0]
    ref_core   = anth_core.positions.copy()
    ref_pip2   = pip2_helix.positions.copy()
    c_ref      = ref_core.mean(axis=0)

    pip2_rmsd_A = []
    for _ in u.trajectory:
        core_now = anth_core.positions.copy()
        pip2_now = pip2_helix.positions.copy()
        c_now    = core_now.mean(axis=0)

        # Kabsch rotation: align core
        H   = (core_now - c_now).T @ (ref_core - c_ref)
        U_, _, Vt = np.linalg.svd(H)
        d   = np.sign(np.linalg.det(Vt.T @ U_.T))
        R   = Vt.T @ np.diag([1, 1, d]) @ U_.T

        pip2_rot = (pip2_now - c_now) @ R.T + c_ref
        rmsd_val = np.sqrt(np.mean(np.sum((pip2_rot - ref_pip2) ** 2, axis=1)))
        pip2_rmsd_A.append(rmsd_val)

    pip2_rmsd_A = np.array(pip2_rmsd_A)

    np.savez(str(cache), cv_site4=cv_site4, pip2_rmsd_A=pip2_rmsd_A, time_ns=time_ns)
    print(f"Cache saved: {cache.name}")
    return cv_site4, pip2_rmsd_A, time_ns


def make_figure(cv_site4, pip2_rmsd_A, time_ns,
                show_labels: bool, show_rval: bool = True):
    """
    Build the scatter figure.
    show_labels=False → TOC version (no axis text, no ticks).
    show_labels=True  → reference version.
    """
    r, _p = stats.pearsonr(cv_site4, pip2_rmsd_A)

    # Subsample
    step  = max(1, len(cv_site4) // MAX_PTS)
    x     = cv_site4[::step]
    y     = pip2_rmsd_A[::step]
    t     = time_ns[::step]

    # Regression line
    m, b  = np.polyfit(cv_site4, pip2_rmsd_A, 1)
    xfit  = np.linspace(cv_site4.min(), cv_site4.max(), 200)
    yfit  = m * xfit + b

    fig, ax = plt.subplots(figsize=(3.8, 3.5), facecolor=BG_COLOR)
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    cmap = time_cmap()
    sc   = ax.scatter(x, y, c=t, cmap=cmap,
                      s=PT_SIZE, alpha=PT_ALPHA, linewidths=0,
                      vmin=time_ns.min(), vmax=time_ns.max(), zorder=2)

    ax.plot(xfit, yfit, color="white", lw=1.8, alpha=0.90, zorder=3)

    if show_rval:
        ax.text(0.97, 0.05, f"r = +{r:.3f}",
                transform=ax.transAxes,
                color="white", fontsize=11, fontweight="bold",
                ha="right", va="bottom",
                fontfamily="DejaVu Sans")

    # Spine styling
    for key in ("top", "bottom", "left", "right"):
        ax.spines[key].set_edgecolor("#2a3f5f")
        ax.spines[key].set_linewidth(0.8)

    if show_labels:
        ax.set_xlabel("site473 aperture, CV₂ (nm)", color="white", fontsize=9)
        ax.set_ylabel("PIP₂-helix RMSD (Å)", color="white", fontsize=9)
        ax.tick_params(colors="white", labelsize=8, length=3, width=0.6)
        cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label("Time (ns)", color="white", fontsize=8)
        cbar.ax.yaxis.set_tick_params(color="white", labelsize=7)
        plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")
        cbar.outline.set_edgecolor("#2a3f5f")
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.tight_layout(pad=0.4)
    return fig


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD TOC — PIP2/site473 allosteric coupling scatter")
    parser.add_argument("--project-root", type=Path, default=_default_root)
    args  = parser.parse_args()
    root  = args.project_root.resolve()

    outdir = root / "06_figures/publication"
    outdir.mkdir(parents=True, exist_ok=True)

    cv_site4, pip2_rmsd_A, time_ns = compute_pip2_data(root)
    r, p = stats.pearsonr(cv_site4, pip2_rmsd_A)
    print(f"  r = {r:.3f}  p = {p:.2e}  N = {len(cv_site4)}")
    print(f"  cv_site4: {cv_site4.min():.2f}–{cv_site4.max():.2f} nm")
    print(f"  pip2_rmsd: {pip2_rmsd_A.min():.2f}–{pip2_rmsd_A.max():.2f} Å")

    # Clean TOC panel
    fig = make_figure(cv_site4, pip2_rmsd_A, time_ns,
                      show_labels=False, show_rval=True)
    for ext in ("pdf", "png"):
        out = outdir / f"TOC_scatter_pip2.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight",
                    facecolor=BG_COLOR, edgecolor="none")
        print(f"Saved: {out.name}")
    plt.close()

    # Reference panel (with labels + colourbar)
    fig = make_figure(cv_site4, pip2_rmsd_A, time_ns,
                      show_labels=True, show_rval=True)
    for ext in ("pdf", "png"):
        out = outdir / f"TOC_scatter_pip2_ref.{ext}"
        fig.savefig(out, dpi=300, bbox_inches="tight",
                    facecolor=BG_COLOR, edgecolor="none")
        print(f"Saved: {out.name}")
    plt.close()

    print("Done.")


if __name__ == "__main__":
    main()
