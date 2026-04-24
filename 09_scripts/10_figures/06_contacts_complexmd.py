"""
CRYPTAD — Figure 4: PICALM site473 contact heatmap from 100 ns complex MD.

Two-panel figure:
  Left:  Ibufenac (CHEMBL341812) engaged replicas 1,2,5,6
         Contact occupancy computed in engaged frames (COM < 15 Å from frame-0).
  Right: Phensuximide (CHEMBL797) all replicas 1,2,3 (dissociated negative control)
         Contact occupancy computed across all frames.

Residue set: top 12 by Ibufenac rep1 mean occupancy, in holo residue numbering.

Output: 06_figures/publication/Figure_04_contacts_site473.{pdf,png}

Usage:
  python3 09_scripts/10_figures/06_contacts_complexmd.py
  python3 09_scripts/10_figures/06_contacts_complexmd.py --project-root /path/to/CRYPTAD
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import MDAnalysis as mda
from MDAnalysis.analysis import distances

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

CONTACT_CUTOFF  = 4.5   # Å
ENGAGED_CUTOFF  = 15.0  # Å — only used for Ibufenac
TOP_N_RESIDUES  = 12    # residues anchored by Ibufenac rep1

IBUFENAC_REPS   = ["rep1", "rep2", "rep5", "rep6"]   # engaged
PHENSUXIMIDE_REPS = ["rep1", "rep2", "rep3"]           # all (dissociated)


def compute_contacts(top_path, xtc_path, engaged_only=True):
    """Return (n_frames_used, dict of (resid,resname) -> occupancy fraction)."""
    u    = mda.Universe(str(top_path), str(xtc_path))
    lig  = u.select_atoms("resname MOL")
    prot = u.select_atoms("protein")

    u.trajectory[0]
    com0 = lig.center_of_geometry()

    res_counts = {}
    n_used = 0

    for ts in u.trajectory:
        com   = lig.center_of_geometry()
        drift = np.linalg.norm(com - com0)

        if engaged_only and drift >= ENGAGED_CUTOFF:
            continue

        n_used += 1
        dm = distances.distance_array(prot.positions, lig.positions, box=None)
        min_d = dm.min(axis=1)

        seen = set()
        for i, atom in enumerate(prot):
            if min_d[i] <= CONTACT_CUTOFF:
                key = (atom.resid, atom.resname)
                seen.add(key)
        for key in seen:
            res_counts[key] = res_counts.get(key, 0) + 1

    occ = {k: v / n_used for k, v in res_counts.items()} if n_used > 0 else {}
    return n_used, occ


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD Figure 4 — site473 contact heatmap (100 ns complex MD)")
    parser.add_argument("--project-root", type=Path, default=_default_root)
    args = parser.parse_args()

    root  = args.project_root.resolve()
    cmdir = root / "02_md_simulations" / "complex_md"
    outdir = root / "06_figures" / "publication"
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Ibufenac: engaged reps ────────────────────────────────────────────────
    print("Computing Ibufenac contacts (engaged frames)...")
    ibu_data = {}   # rep -> {(resid,resname): occupancy}
    ibu_frames = {}
    for rep in IBUFENAC_REPS:
        wdir = cmdir / "S3_site473_CHEMBL341812" / rep
        n, occ = compute_contacts(
            wdir / "ref_nw.pdb", wdir / "complexmd_fit_nw.xtc",
            engaged_only=True)
        ibu_data[rep]   = occ
        ibu_frames[rep] = n
        print(f"  {rep}: {n} engaged frames, {len(occ)} residues contacted")

    # ── Phensuximide: all frames ───────────────────────────────────────────────
    print("Computing Phensuximide contacts (all frames, dissociated control)...")
    phe_data = {}
    phe_frames = {}
    for rep in PHENSUXIMIDE_REPS:
        wdir = cmdir / "S3_site473_CHEMBL797" / rep
        n, occ = compute_contacts(
            wdir / "ref_nw.pdb", wdir / "complexmd_fit_nw.xtc",
            engaged_only=False)
        phe_data[rep]   = occ
        phe_frames[rep] = n
        print(f"  {rep}: {n} total frames, {len(occ)} residues contacted")

    # ── Select residues: top 12 by Ibufenac rep1 occupancy ───────────────────
    rep1_occ = ibu_data["rep1"]
    top_residues = sorted(rep1_occ, key=lambda k: -rep1_occ[k])[:TOP_N_RESIDUES]
    # Sort display order: highest rep1 occupancy at top
    res_labels = [f"{rn}{rid}" for rid, rn in top_residues]
    print(f"\nTop {TOP_N_RESIDUES} residues (by Ibufenac rep1): {res_labels}")

    # ── Build matrices ────────────────────────────────────────────────────────
    # Ibufenac: rows=residues, cols=reps
    ibu_mat = np.zeros((len(top_residues), len(IBUFENAC_REPS)))
    for j, rep in enumerate(IBUFENAC_REPS):
        for i, key in enumerate(top_residues):
            ibu_mat[i, j] = ibu_data[rep].get(key, 0.0)

    # Phensuximide: rows=residues, cols=reps
    phe_mat = np.zeros((len(top_residues), len(PHENSUXIMIDE_REPS)))
    for j, rep in enumerate(PHENSUXIMIDE_REPS):
        for i, key in enumerate(top_residues):
            phe_mat[i, j] = phe_data[rep].get(key, 0.0)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5),
                              gridspec_kw={"wspace": 0.45})
    fig.suptitle(
        "PICALM site473 — Residue Contact Occupancy from 100 ns Complex MD",
        fontweight="bold", fontsize=11, y=1.01)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})

    cmap_ibu = LinearSegmentedColormap.from_list("ibu", ["#ffffff", "#2ca02c"])
    cmap_phe = LinearSegmentedColormap.from_list("phe", ["#ffffff", "#d62728"])

    def annotate_heatmap(ax, mat):
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if v >= 0.10:
                    ax.text(j, i, f"{v*100:.0f}%", ha="center", va="center",
                            fontsize=8,
                            color="white" if v > 0.6 else "#222222")

    # Left — Ibufenac
    ax = axes[0]
    im = ax.imshow(ibu_mat, cmap=cmap_ibu, vmin=0, vmax=1,
                   aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(IBUFENAC_REPS)))
    rep_labels = [f"rep{r[-1]}\n({ibu_frames[r]:,} fr)" for r in IBUFENAC_REPS]
    ax.set_xticklabels(rep_labels, fontsize=8.5)
    ax.set_yticks(range(len(res_labels)))
    ax.set_yticklabels(res_labels, fontsize=9)
    ax.set_title("Ibufenac (CHEMBL341812)\nengaged replicas (4/6)",
                 fontweight="bold", fontsize=10)
    ax.set_xlabel("Replica", fontsize=9)
    ax.set_ylabel("Residue", fontsize=9)
    annotate_heatmap(ax, ibu_mat)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Contact occupancy")

    # Right — Phensuximide (dissociated control)
    ax = axes[1]
    im2 = ax.imshow(phe_mat, cmap=cmap_phe, vmin=0, vmax=1,
                    aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(PHENSUXIMIDE_REPS)))
    phe_labels = [f"rep{r[-1]}\n({phe_frames[r]:,} fr)" for r in PHENSUXIMIDE_REPS]
    ax.set_xticklabels(phe_labels, fontsize=8.5)
    ax.set_yticks(range(len(res_labels)))
    ax.set_yticklabels(res_labels, fontsize=9)
    ax.set_title("Phensuximide (CHEMBL797)\ndissociated control (0/3)",
                 fontweight="bold", fontsize=10)
    ax.set_xlabel("Replica", fontsize=9)
    annotate_heatmap(ax, phe_mat)
    plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.04, label="Contact occupancy")

    # Panel labels
    for ax, label in zip(axes, ["a", "b"]):
        ax.text(-0.12, 1.04, f"({label})", transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="top")

    plt.tight_layout()

    for ext in ("pdf", "png"):
        out = outdir / f"Figure_04_contacts_site473.{ext}"
        plt.savefig(out, dpi=600, bbox_inches="tight")
        print(f"Saved: {out.relative_to(root)}")

    plt.close()
    print("Done.")


if __name__ == "__main__":
    main()
