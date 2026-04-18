"""
CRYPTAD — Supplementary Figure: RMSF apo vs. holo (PICALM ANTH)

Computes per-Cα RMSF for:
  Apo : S3_ANTH metadynamics trajectory (200 ns, stride=10 → ~1000 frames)
  Holo: Ibufenac (CHEMBL341812) complex MD, rep1 + rep5 (both 100% engaged)

RMSF is computed after fitting each trajectory to its own mean structure.

Key residues highlighted:
  site473 pocket  (holo numbering): LEU219, TYR237, PHE240, LEU267, SER268, GLN269, ALA270, LEU278
  PIP2-helix      (apo  numbering):  1– 22
  Ubiquitin face  (apo  numbering): ~140–180  (approx C-terminal unique region)

Residue numbering:
  Apo  (meta.tpr)  : resids  1–268
  Holo (ref_nw.pdb): resids 19–286   (offset +18)
  Plot uses APO numbering (subtract 18 from holo resids for pocket labels).

Output: 06_figures/fig_rmsf_apo_holo.png  (300 dpi)

Usage:
  python3 09_scripts/10_figures/02_rmsf_apo_holo.py
  python3 09_scripts/10_figures/02_rmsf_apo_holo.py --project-root /path/to/CRYPTAD
  python3 09_scripts/10_figures/02_rmsf_apo_holo.py --no-cache
"""

import argparse
import logging
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage import uniform_filter1d
import MDAnalysis as mda

warnings.filterwarnings("ignore")
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

# ── Site473 pocket residues (holo numbering → apo = holo − 18) ───────────────
POCKET_HOLO = [219, 237, 240, 267, 268, 269, 270, 278]
POCKET_APO  = [r - 18 for r in POCKET_HOLO]

HOLO_REPS = ["rep1", "rep5"]    # 100% engaged

LABEL_OFFSETS = {
    201: (  0,  32),
    219: (-20,  28),
    222: ( 28,  38),
    249: (-42,  32),
    250: ( 42,  48),
    251: (-42,  80),
    252: ( 42,  96),
    260: ( 24,  28),
}


# ── Helper: compute RMSF ──────────────────────────────────────────────────────
def compute_rmsf(tpr: Path, xtc: Path, sel="protein and name CA", stride=1):
    """
    Returns (resids, resnames, rmsf_angstrom).
    RMSF computed after iterative Kabsch fitting to mean structure (3 passes).
    """
    print(f"  Loading: {xtc.name}  (stride={stride})", flush=True)
    u  = mda.Universe(str(tpr), str(xtc))
    ca = u.select_atoms(sel)
    res = ca.resids.copy()
    rn  = np.array([a.resname for a in ca])

    frames = []
    for ts in u.trajectory[::stride]:
        frames.append(ca.positions.copy())
    frames = np.array(frames)
    n_frames = len(frames)
    print(f"    {n_frames} frames  |  {len(res)} Cα atoms")

    mean_pos = frames.mean(axis=0)
    for _ in range(3):
        centroid_frames = frames.mean(axis=1, keepdims=True)
        centroid_mean   = mean_pos.mean(axis=0, keepdims=True)
        frames_cent     = frames - centroid_frames
        mean_cent       = mean_pos - centroid_mean
        aligned = np.zeros_like(frames)
        for i, f in enumerate(frames_cent):
            H        = f.T @ mean_cent
            U, S, Vt = np.linalg.svd(H)
            d        = np.linalg.det(Vt.T @ U.T)
            D        = np.diag([1, 1, d])
            R        = Vt.T @ D @ U.T
            aligned[i] = f @ R.T
        frames   = aligned
        mean_pos = frames.mean(axis=0)

    diff = frames - mean_pos[np.newaxis, :, :]
    rmsf = np.sqrt((diff**2).sum(axis=2).mean(axis=0))
    return res, rn, rmsf


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD RMSF apo vs. holo figure (PICALM ANTH)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    parser.add_argument("--no-cache", action="store_true",
                        help="Recompute RMSF even if cache files exist")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    outdir       = project_root / "06_figures"
    outdir.mkdir(parents=True, exist_ok=True)

    s3_tpr       = project_root / "02_md_simulations" / "PICALM" / "metadynamics" / "S3_ANTH" / "meta.tpr"
    s3_xtc       = project_root / "02_md_simulations" / "PICALM" / "metadynamics" / "S3_ANTH" / "meta.xtc"
    complex_base = project_root / "02_md_simulations" / "complex_md" / "S3_site473_CHEMBL341812"

    cache_apo  = outdir / "rmsf_cache_apo.npy"
    cache_holo = outdir / "rmsf_cache_holo.npy"
    cache_res  = outdir / "rmsf_cache_res.npy"
    cache_rn   = outdir / "rmsf_cache_rn.npy"

    have_cache = (cache_apo.exists() and cache_holo.exists()
                  and cache_res.exists() and cache_rn.exists())

    if have_cache and not args.no_cache:
        print("Loading cached RMSF arrays …")
        rmsf_apo  = np.load(cache_apo)
        rmsf_holo = np.load(cache_holo)
        res_apo   = np.load(cache_res)
        rn_apo    = np.load(cache_rn, allow_pickle=True)
        holo_res  = res_apo.copy()
        print(f"  apo: {len(rmsf_apo)} residues  |  holo: {len(rmsf_holo)} residues")
    else:
        print("=" * 60)
        print("RMSF — Apo (metadynamics, stride=10)")
        print("=" * 60)
        res_apo, rn_apo, rmsf_apo = compute_rmsf(s3_tpr, s3_xtc, stride=10)

        print("\n" + "=" * 60)
        print("RMSF — Holo (Ibufenac rep1+rep5, 100% engaged)")
        print("=" * 60)
        holo_rmsf_list = []
        holo_res = None
        for rep in HOLO_REPS:
            tpr_h = complex_base / rep / "ref_nw.pdb"
            xtc_h = complex_base / rep / "complexmd_fit_nw.xtc"
            res_h, rn_h, rmsf_h = compute_rmsf(tpr_h, xtc_h, stride=1)
            holo_rmsf_list.append(rmsf_h)
            if holo_res is None:
                holo_res = res_h

        rmsf_holo = np.mean(holo_rmsf_list, axis=0)

        np.save(cache_apo,  rmsf_apo)
        np.save(cache_holo, rmsf_holo)
        np.save(cache_res,  res_apo)
        np.save(cache_rn,   rn_apo)
        print("RMSF arrays cached.")

    assert len(res_apo) == len(holo_res), (
        f"Cα count mismatch: apo={len(res_apo)}, holo={len(holo_res)}")

    n_res   = len(res_apo)
    x_idx   = np.arange(1, n_res + 1)

    pocket_idx = [r for r in POCKET_APO if r in res_apo]
    pocket_pos = [list(res_apo).index(r) + 1 for r in pocket_idx]

    pip2_mask = (res_apo >= 1)   & (res_apo <= 22)
    ubiq_mask = (res_apo >= 140) & (res_apo <= 180)

    print(f"\nApo  RMSF: mean={rmsf_apo.mean():.2f}, max={rmsf_apo.max():.2f} Å")
    print(f"Holo RMSF: mean={rmsf_holo.mean():.2f}, max={rmsf_holo.max():.2f} Å")
    pocket_apo_rmsf  = rmsf_apo[[i for i, r in enumerate(res_apo) if r in POCKET_APO]]
    pocket_holo_rmsf = rmsf_holo[[i for i, r in enumerate(res_apo) if r in POCKET_APO]]
    print(f"\nsite473 residues:")
    print(f"  Apo  RMSF (mean at pocket): {pocket_apo_rmsf.mean():.2f} Å")
    print(f"  Holo RMSF (mean at pocket): {pocket_holo_rmsf.mean():.2f} Å")

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1.2]})
    ax  = axes[0]
    ax2 = axes[1]

    rmsf_apo_sm  = uniform_filter1d(rmsf_apo,  size=3)
    rmsf_holo_sm = uniform_filter1d(rmsf_holo, size=3)

    COLOR_APO  = "#2171B5"
    COLOR_HOLO = "#E6550D"
    COLOR_DIFF = "#31A354"

    ax.axvspan(1,   22,  alpha=0.08, color="#9ECAE1", zorder=0)
    ax.axvspan(140, 180, alpha=0.06, color="#FDAE6B", zorder=0)

    ax.fill_between(x_idx, 0, rmsf_apo_sm,  alpha=0.15, color=COLOR_APO)
    ax.fill_between(x_idx, 0, rmsf_holo_sm, alpha=0.20, color=COLOR_HOLO)
    ax.plot(x_idx, rmsf_apo_sm,  color=COLOR_APO,  lw=1.4,
            label="Apo (200 ns metadynamics)")
    ax.plot(x_idx, rmsf_holo_sm, color=COLOR_HOLO, lw=1.4,
            label="Holo (Ibufenac rep1+rep5, 100% engaged)")

    for pos, rid in zip(pocket_pos, pocket_idx):
        y_apo  = rmsf_apo[pos - 1]
        y_holo = rmsf_holo[pos - 1]
        rn     = rn_apo[pos - 1]
        ax.scatter(pos, y_apo,  color=COLOR_APO,  s=55, zorder=5)
        ax.scatter(pos, y_holo, color=COLOR_HOLO, s=55, zorder=5)
        y_anchor  = max(y_apo, y_holo)
        dx, dy    = LABEL_OFFSETS.get(rid, (0, 32))
        label_txt = f"{rn}{rid}"
        ha = "right" if dx < 0 else ("center" if dx == 0 else "left")
        ax.annotate(
            label_txt,
            xy=(pos, y_anchor),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=ha, va="bottom",
            fontsize=7.5, fontweight="bold", color="#333333",
            arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.7,
                            shrinkA=2, shrinkB=3),
            annotation_clip=False,
        )

    ax.set_ylabel("Cα RMSF (Å)", fontsize=10)
    ax.set_ylim(bottom=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9, loc="upper left", framealpha=0.9)
    ax.set_title(
        "PICALM ANTH — Cα RMSF: Apo (metadynamics) vs. Holo (Ibufenac complex MD)\n"
        "site473 pocket residues marked  |  Blue band = PIP2-helix (resids 1–22)  |"
        "  Orange band = C-terminal ubiquitin face",
        fontsize=10, fontweight="bold")

    ax.text(11,  ax.get_ylim()[1] * 0.95, "PIP2-helix",
            ha="center", va="top", fontsize=8, color="#2166AC", style="italic")
    ax.text(160, ax.get_ylim()[1] * 0.95, "Ub-face",
            ha="center", va="top", fontsize=8, color="#D94801", style="italic")

    diff = rmsf_holo_sm - rmsf_apo_sm
    ax2.axhline(0, color="black", lw=0.8, ls="-")
    ax2.fill_between(x_idx, 0, diff,
                     where=(diff < 0), alpha=0.55, color=COLOR_DIFF,
                     label="Stabilised in holo (ΔRMSF < 0)")
    ax2.fill_between(x_idx, 0, diff,
                     where=(diff > 0), alpha=0.35, color="#FB6A4A",
                     label="More flexible in holo (ΔRMSF > 0)")
    ax2.plot(x_idx, diff, color="#666666", lw=0.8, alpha=0.8)

    for pos, rid in zip(pocket_pos, pocket_idx):
        ax2.scatter(pos, diff[pos - 1], color="#333333", s=40, zorder=5)

    ax2.axvspan(1,   22,  alpha=0.08, color="#9ECAE1", zorder=0)
    ax2.axvspan(140, 180, alpha=0.06, color="#FDAE6B", zorder=0)
    ax2.set_ylabel("ΔRMSF\n(holo−apo, Å)", fontsize=9)
    ax2.set_xlabel("Residue (apo numbering, 1–268)", fontsize=10)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(fontsize=8, loc="lower right", framealpha=0.9, ncol=2)

    ax2.set_xticks(list(range(0, n_res + 1, 50)))
    pocket_str = "/".join([str(r) for r in sorted(POCKET_APO)])
    fig.text(0.5, -0.01,
             f"site473 pocket residues (apo numbering): {pocket_str}\n"
             "(holo numbering = apo + 18; offset from CHARMM-GUI construct start at residue 19)",
             ha="center", fontsize=8, color="#555555", style="italic")

    plt.tight_layout()
    out = outdir / "fig_rmsf_apo_holo.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\nFigure saved: {out.relative_to(project_root)}")

    # ── Text summary ──────────────────────────────────────────────────────────
    summary_lines = [
        "RMSF Apo vs Holo — PICALM ANTH (site473)",
        "=" * 50,
        f"Apo  (metadynamics, 200 ns, stride=10): {len(res_apo)} Cα, "
        f"mean RMSF = {rmsf_apo.mean():.2f} Å, max = {rmsf_apo.max():.2f} Å",
        f"Holo (Ibufenac rep1+rep5 avg):          {len(holo_res)} Cα, "
        f"mean RMSF = {rmsf_holo.mean():.2f} Å, max = {rmsf_holo.max():.2f} Å",
        "",
        "site473 pocket residues (holo → apo numbering):",
    ]
    for holo_r, apo_r in zip(POCKET_HOLO, POCKET_APO):
        idx  = list(res_apo).index(apo_r)
        rn   = rn_apo[idx]
        ra   = rmsf_apo[idx]
        rh   = rmsf_holo[idx]
        ddr  = rh - ra
        flag = ("↓ stabilised" if ddr < -0.05
                else ("↑ more flexible" if ddr > 0.05 else "≈ unchanged"))
        summary_lines.append(
            f"  {rn}{holo_r:>4} (apo {apo_r:>3}): apo {ra:.2f} Å  holo {rh:.2f} Å  "
            f"Δ {ddr:+.2f}  {flag}"
        )
    summary_lines += [
        "",
        f"PIP2-helix (resids 1–22):  apo {rmsf_apo[pip2_mask].mean():.2f} Å  "
        f"holo {rmsf_holo[pip2_mask].mean():.2f} Å",
        f"Ubiquitin face (~140–180): apo {rmsf_apo[ubiq_mask].mean():.2f} Å  "
        f"holo {rmsf_holo[ubiq_mask].mean():.2f} Å",
        "",
        "Caveat: Apo RMSF reflects enhanced sampling (metadynamics) — higher values expected.",
        "Meaningful comparison: RELATIVE RMSF change at site473 vs rest of domain.",
    ]
    summary = "\n".join(summary_lines)
    print("\n" + summary)

    summary_path = outdir / "rmsf_apo_holo_summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary + "\n")
    print(f"\nSummary saved: {summary_path.relative_to(project_root)}")


if __name__ == "__main__":
    main()
