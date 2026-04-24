"""
CRYPTAD — Figure 5: BIN1 BAR contact heatmap from MM-GBSA 3 ns production MD.

Two-panel figure:
  Left:  BIN1 site688 — four final hits (Norethindrone, Segesterone Ac.,
         Zuranolone, Ethynodiol DA) — 101 frames each.
  Right: BIN1 site680 — Zuranolone only (dual-pocket binder) — 101 frames.

Residues displayed:
  Left : top 12 by mean occupancy across all four site688 compounds.
  Right: top 10 by Zuranolone site680 occupancy.

Contact cutoff: 4.5 Å heavy atom–heavy atom distance.

Output: 06_figures/publication/Figure_05_contacts_bin1.{pdf,png}

Usage:
  python3 09_scripts/10_figures/07_contacts_bin1.py
  python3 09_scripts/10_figures/07_contacts_bin1.py --project-root /path
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

CONTACT_CUTOFF_688 = 4.5   # Å — hydrophobic packing
CONTACT_CUTOFF_680 = 3.5   # Å — polar/H-bond contacts (tighter, gives occupancy gradient)

SITE688_COMPOUNDS = [
    ("CHEMBL1162",   "Norethindrone"),
    ("CHEMBL3707377","Segesterone Ac."),
    ("CHEMBL4105630","Zuranolone"),
    ("CHEMBL1200624","Ethynodiol DA"),
]
TOP_SITE688 = 12
TOP_SITE680 = 10


def compute_contacts_all_frames(gro_path, xtc_path, cutoff=CONTACT_CUTOFF_688):
    """Return (n_frames, dict (resid, resname) -> occupancy fraction, all frames)."""
    u    = mda.Universe(str(gro_path), str(xtc_path))
    lig  = u.select_atoms("resname MOL")
    prot = u.select_atoms("protein")

    res_counts = {}
    n_frames   = 0

    for ts in u.trajectory:
        n_frames += 1
        dm    = distances.distance_array(prot.positions, lig.positions, box=None)
        min_d = dm.min(axis=1)
        seen  = set()
        for i, atom in enumerate(prot):
            if min_d[i] <= cutoff:
                seen.add((atom.resid, atom.resname))
        for key in seen:
            res_counts[key] = res_counts.get(key, 0) + 1

    occ = {k: v / n_frames for k, v in res_counts.items()}
    return n_frames, occ


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD Figure 5 — BIN1 BAR contact heatmap (MM-GBSA 3 ns MD)")
    parser.add_argument("--project-root", type=Path, default=_default_root)
    args  = parser.parse_args()
    root  = args.project_root.resolve()
    base  = root / "04_virtual_screening" / "mmgbsa"
    outdir = root / "06_figures" / "publication"
    outdir.mkdir(parents=True, exist_ok=True)

    # ── site688: four compounds ───────────────────────────────────────────────
    print("Computing BIN1 site688 contacts...")
    site688_data   = {}
    site688_frames = {}
    for cid, name in SITE688_COMPOUNDS:
        wdir = base / "S1_site688" / cid
        n, occ = compute_contacts_all_frames(wdir / "prod.gro", wdir / "prod.xtc")
        site688_data[cid]   = occ
        site688_frames[cid] = n
        print(f"  {name} ({cid}): {n} frames, {len(occ)} residues contacted")

    # ── site688 residue selection: top N by mean occupancy across all 4 ───────
    all_keys = set()
    for occ in site688_data.values():
        all_keys.update(occ.keys())
    mean_occ = {k: np.mean([site688_data[cid].get(k, 0.0)
                             for cid, _ in SITE688_COMPOUNDS])
                for k in all_keys}
    top688 = sorted(mean_occ, key=lambda k: -mean_occ[k])[:TOP_SITE688]
    labels688 = [f"{rn}{rid}" for rid, rn in top688]
    print(f"Top {TOP_SITE688} site688 residues: {labels688}")

    # ── site688 matrix ─────────────────────────────────────────────────────────
    mat688 = np.zeros((len(top688), len(SITE688_COMPOUNDS)))
    for j, (cid, _) in enumerate(SITE688_COMPOUNDS):
        for i, key in enumerate(top688):
            mat688[i, j] = site688_data[cid].get(key, 0.0)

    # ── site680: Zuranolone only (3.5 Å polar cutoff for occupancy gradient) ────
    print("\nComputing BIN1 site680 contacts (Zuranolone, 3.5 Å cutoff)...")
    wdir680 = base / "S1_site680" / "CHEMBL4105630"
    n680, occ680 = compute_contacts_all_frames(
        wdir680 / "prod.gro", wdir680 / "prod.xtc",
        cutoff=CONTACT_CUTOFF_680)
    print(f"  Zuranolone site680: {n680} frames, {len(occ680)} residues contacted")

    top680 = sorted(occ680, key=lambda k: -occ680[k])[:TOP_SITE680]
    labels680 = [f"{rn}{rid}" for rid, rn in top680]
    print(f"Top {TOP_SITE680} site680 residues: {labels680}")
    mat680 = np.array([[occ680.get(k, 0.0)] for k in top680])

    # ── Figure ────────────────────────────────────────────────────────────────
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})

    fig = plt.figure(figsize=(14, 6))
    # Manual gridspec: left panel wider, right panel narrower, explicit spacing
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1.6],
                          left=0.08, right=0.95, wspace=0.55)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    axes = [ax0, ax1]

    fig.suptitle(
        "BIN1 BAR Domain — Residue Contact Occupancy (MM-GBSA 3 ns production MD)",
        fontweight="bold", fontsize=11, y=1.00)

    cmap688 = LinearSegmentedColormap.from_list("bin1", ["#ffffff", "#1f77b4"])
    cmap680 = LinearSegmentedColormap.from_list("bin1l", ["#ffffff", "#aec7e8"])

    def annotate(ax, mat, fs=8):
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if v >= 0.10:
                    ax.text(j, i, f"{v*100:.0f}%", ha="center", va="center",
                            fontsize=fs,
                            color="white" if v > 0.6 else "#222222")

    # Left — site688
    ax = axes[0]
    im = ax.imshow(mat688, cmap=cmap688, vmin=0, vmax=1,
                   aspect="auto", interpolation="nearest")
    col_labels = [f"{name}\n({site688_frames[cid]:,} fr)"
                  for cid, name in SITE688_COMPOUNDS]
    ax.set_xticks(range(len(SITE688_COMPOUNDS)))
    ax.set_xticklabels(col_labels, fontsize=8.5)
    ax.set_yticks(range(len(labels688)))
    ax.set_yticklabels(labels688, fontsize=9)
    ax.set_title("(a)  BIN1 site688 — four final hits",
                 fontweight="bold", fontsize=10, loc="left")
    ax.set_ylabel("Residue", fontsize=9)
    annotate(ax, mat688, fs=8)
    cb1 = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb1.set_label("Contact occupancy", fontsize=8)
    cb1.ax.tick_params(labelsize=8)

    # Right — site680 Zuranolone
    ax = axes[1]
    im2 = ax.imshow(mat680, cmap=cmap680, vmin=0, vmax=1,
                    aspect="auto", interpolation="nearest")
    ax.set_xticks([0])
    ax.set_xticklabels([f"Zuranolone\n({n680:,} fr)"], fontsize=8.5)
    ax.set_yticks(range(len(labels680)))
    ax.set_yticklabels(labels680, fontsize=9)
    ax.set_title("(b)  BIN1 site680\nZuranolone (3.5 Å cutoff)",
                 fontweight="bold", fontsize=10, loc="left")
    annotate(ax, mat680, fs=8)
    cb2 = fig.colorbar(im2, ax=ax, fraction=0.12, pad=0.04)
    cb2.set_label("Contact occupancy", fontsize=8)
    cb2.ax.tick_params(labelsize=8)

    for ext in ("pdf", "png"):
        out = outdir / f"Figure_05_contacts_bin1.{ext}"
        plt.savefig(out, dpi=600, bbox_inches="tight")
        print(f"Saved: {out.relative_to(root)}")

    plt.close()

    # ── Print key residues for legend writing ──────────────────────────────────
    print("\n--- Key site688 residues (mean occ ≥ 0.80) ---")
    for rid, rn in top688:
        m = mean_occ[(rid, rn)]
        if m >= 0.80:
            print(f"  {rn}{rid}: {m*100:.1f}%")

    print("\n--- Key site680 residues (Zuranolone occ ≥ 0.30) ---")
    for rid, rn in top680:
        v = occ680[(rid, rn)]
        if v >= 0.30:
            print(f"  {rn}{rid}: {v*100:.1f}%")

    print("\nDone.")


if __name__ == "__main__":
    main()
