"""
CRYPTAD — Complex MD RMSD traces + Ibufenac H-bond occupancy
=============================================================
Generates two publication outputs:

1. fig_complex_md_rmsd.png  (4-panel, one per system)
   Cα backbone RMSD vs. time for all replicas.
   Ibufenac reps colour-coded: green = engaged, orange = partial, red = dissociated.

2. fig_ibufenac_hbond.png  (heatmap)
   Distance-based H-bond occupancy at Ibufenac carboxylate (O1/O2) for engaged
   reps (rep1, rep2, rep5, rep6).  Criterion: heavy-atom donor/acceptor distance
   ≤ 3.5 Å between MOL O1 or O2 and any protein N or O atom.

3. ibufenac_hbond_table.csv
   Full per-replica, per-residue occupancy table (supplementary data).

Notes
-----
- All trajectories loaded from complexmd_fit_nw.xtc (PBC-corrected, protein-centred).
- RMSD reference = frame 0 of each replica's own trajectory (relative stability).
- H-bond analysis uses distance criterion only (element types absent in PDB topology);
  this is equivalent to a standard 3.5 Å acceptor–donor distance cutoff widely used
  in MD contact analyses.
- Ibufenac (CHEMBL341812) = p-isobutylphenylacetic acid, C12H16O2; only polar atoms
  are O1 and O2 (carboxylate).

Usage:
  python3 09_scripts/10_figures/01_complex_md_rmsd_hbond.py
  python3 09_scripts/10_figures/01_complex_md_rmsd_hbond.py --project-root /path/to/CRYPTAD
"""

import argparse
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import uniform_filter1d
import MDAnalysis as mda

warnings.filterwarnings("ignore")
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

# ── System definitions ────────────────────────────────────────────────────────
# engagement: "engaged" | "partial" | "dissociated"
SYSTEMS = [
    {
        "id":    "BIN1_site688_Zuranolone",
        "label": "BIN1 site688 · Zuranolone",
        "dir":   "S1_site688_CHEMBL4105630",
        "reps":  ["rep1", "rep2", "rep3"],
        "engage": {"rep1": "dissociated", "rep2": "dissociated", "rep3": "dissociated"},
        "rmsd_sel": "segid A and name CA",
        "hbond": False,
    },
    {
        "id":    "BIN1_site680_Zuranolone",
        "label": "BIN1 site680 · Zuranolone",
        "dir":   "S1_site680_CHEMBL4105630",
        "reps":  ["rep1", "rep2", "rep3"],
        "engage": {"rep1": "dissociated", "rep2": "dissociated", "rep3": "dissociated"},
        "rmsd_sel": "segid A and name CA",
        "hbond": False,
    },
    {
        "id":    "PICALM_site473_Ibufenac",
        "label": "PICALM site473 · Ibufenac",
        "dir":   "S3_site473_CHEMBL341812",
        "reps":  ["rep1", "rep2", "rep3", "rep4", "rep5", "rep6"],
        "engage": {
            "rep1": "engaged",     # 100%
            "rep2": "engaged",     # 84%
            "rep3": "partial",     # 35%
            "rep4": "partial",     # 56%
            "rep5": "engaged",     # 100%
            "rep6": "engaged",     # 84%
        },
        "rmsd_sel": "protein and name CA",
        "hbond": True,
        "hbond_reps": ["rep1", "rep2", "rep5", "rep6"],
    },
    {
        "id":    "PICALM_site473_Phensuximide",
        "label": "PICALM site473 · Phensuximide",
        "dir":   "S3_site473_CHEMBL797",
        "reps":  ["rep1", "rep2", "rep3"],
        "engage": {"rep1": "dissociated", "rep2": "dissociated", "rep3": "dissociated"},
        "rmsd_sel": "protein and name CA",
        "hbond": False,
    },
]

ENGAGE_COLOR = {
    "engaged":     "#2CA02C",
    "partial":     "#FF7F0E",
    "dissociated": "#D62728",
}
ENGAGE_ALPHA = {
    "engaged": 0.85, "partial": 0.65, "dissociated": 0.50,
}


# ── Helper: load and compute Cα RMSD (Kabsch superposition) ──────────────────
def compute_rmsd(pdb: Path, xtc: Path, sel="protein and name CA", stride=1):
    """Return (time_ns, rmsd_A) arrays.  Reference = frame 0."""
    u  = mda.Universe(str(pdb), str(xtc))
    ca = u.select_atoms(sel)

    u.trajectory[0]
    ref_pos  = ca.positions.copy()
    ref_com  = ref_pos.mean(axis=0)
    ref_cent = ref_pos - ref_com

    rmsds = []
    times = []
    for ts in u.trajectory[::stride]:
        pos   = ca.positions.copy()
        com   = pos.mean(axis=0)
        cent  = pos - com
        H     = cent.T @ ref_cent
        U, S, Vt = np.linalg.svd(H)
        d     = np.linalg.det(Vt.T @ U.T)
        D     = np.diag([1.0, 1.0, d])
        R     = Vt.T @ D @ U.T
        aligned = cent @ R.T
        diff  = aligned - ref_cent
        rmsd  = np.sqrt((diff ** 2).sum(axis=1).mean())
        rmsds.append(rmsd)
        times.append(ts.time / 1000.0)   # ps → ns

    return np.array(times), np.array(rmsds)


# ── Helper: distance-based H-bond contacts ───────────────────────────────────
def compute_hbond_contacts(pdb: Path, xtc: Path, dist_cutoff=3.5, stride=1):
    """
    For each frame, find protein residues whose N or O atoms are within
    dist_cutoff of the ligand O1 or O2 atoms.
    Returns: (dict {res_label: occupancy_fraction}, n_frames)
    """
    u    = mda.Universe(str(pdb), str(xtc))
    lig  = u.select_atoms("resname MOL and (name O1 or name O2)")
    prot = u.select_atoms("protein and (name N* or name O*)")

    if len(lig) == 0:
        print("    WARNING: no O1/O2 atoms found in ligand selection")
        return {}, 0

    contact_counts = {}
    n_frames = 0

    for ts in u.trajectory[::stride]:
        n_frames += 1
        lig_pos    = lig.positions
        prot_pos   = prot.positions
        prot_names = np.array([f"{a.resname}{a.resid}_{a.name}" for a in prot])

        for lp in lig_pos:
            dists = np.sqrt(((prot_pos - lp) ** 2).sum(axis=1))
            hits  = prot_names[dists <= dist_cutoff]
            for h in hits:
                contact_counts[h] = contact_counts.get(h, 0) + 1

    occupancy = {k: v / n_frames for k, v in contact_counts.items()}
    return occupancy, n_frames


def atom_to_residue(occ_dict: dict) -> dict:
    """Collapse atom-level contacts to residue level (max per residue)."""
    res_occ = {}
    for atom_key, frac in occ_dict.items():
        res_id = atom_key.split("_")[0]
        res_occ[res_id] = max(res_occ.get(res_id, 0.0), frac)
    return res_occ


def res_sort_key(r: str) -> int:
    num = "".join(c for c in r if c.isdigit())
    return int(num) if num else 0


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD complex MD RMSD + H-bond figures")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    cmbase  = project_root / "02_md_simulations" / "complex_md"
    outdir  = project_root / "06_figures"
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Compute all RMSDs ─────────────────────────────────────────────────────
    print("Computing Cα RMSD for all systems and replicas ...")
    rmsd_data = {}
    for sys in SYSTEMS:
        rmsd_data[sys["id"]] = {}
        for rep in sys["reps"]:
            d   = cmbase / sys["dir"] / rep
            pdb = d / "ref_nw.pdb"
            xtc = d / "complexmd_fit_nw.xtc"
            if not (pdb.exists() and xtc.exists()):
                print(f"  SKIP (missing): {sys['id']} {rep}")
                continue
            print(f"  {sys['id']:35s} {rep} ...", end=" ", flush=True)
            t, r = compute_rmsd(pdb, xtc,
                                sel=sys.get("rmsd_sel", "protein and name CA"),
                                stride=1)
            rmsd_data[sys["id"]][rep] = (t, r)
            print(f"{len(t)} frames, max RMSD {r.max():.2f} Å")

    # ── Figure 1: RMSD traces (4 panels) ─────────────────────────────────────
    print("\nPlotting RMSD figure ...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), sharey=False)
    axes = axes.flatten()

    for ax_i, sys in enumerate(SYSTEMS):
        ax  = axes[ax_i]
        sid = sys["id"]
        for rep in sys["reps"]:
            if rep not in rmsd_data[sid]:
                continue
            t, r = rmsd_data[sid][rep]
            eng  = sys["engage"].get(rep, "dissociated")
            col  = ENGAGE_COLOR[eng]
            alp  = ENGAGE_ALPHA[eng]
            lw   = 1.6 if eng == "engaged" else 1.1
            r_sm = uniform_filter1d(r, size=5)
            ax.plot(t, r_sm, color=col, lw=lw, alpha=alp, label=f"{rep} ({eng})")

        ax.set_title(sys["label"], fontsize=10, fontweight="bold")
        ax.set_xlabel("Time (ns)", fontsize=9)
        ax.set_ylabel("Cα RMSD (Å)", fontsize=9)
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(fontsize=7.5, loc="upper right", framealpha=0.85,
                  handlelength=1.5, labelspacing=0.3)

    legend_elements = [
        mpatches.Patch(color=ENGAGE_COLOR["engaged"],     label="Engaged (≥80% simulation)"),
        mpatches.Patch(color=ENGAGE_COLOR["partial"],     label="Partial (35–56%)"),
        mpatches.Patch(color=ENGAGE_COLOR["dissociated"], label="Dissociated"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=3,
               fontsize=9, framealpha=0.9, bbox_to_anchor=(0.5, -0.03))
    fig.suptitle(
        "Complex MD — Protein Cα RMSD (reference: frame 0)\n"
        "100 ns replicas · CHARMM36m · GROMACS 2023.3",
        fontsize=11, fontweight="bold"
    )
    plt.tight_layout()
    out_rmsd = outdir / "fig_complex_md_rmsd.png"
    plt.savefig(out_rmsd, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_rmsd.relative_to(project_root)}")

    # ── Compute H-bond contacts for Ibufenac ─────────────────────────────────
    print("\nComputing H-bond contacts for Ibufenac engaged reps ...")
    ibufenac_sys = next(s for s in SYSTEMS if s["id"] == "PICALM_site473_Ibufenac")
    hbond_reps   = ibufenac_sys["hbond_reps"]

    hbond_results = {}
    for rep in hbond_reps:
        d   = cmbase / ibufenac_sys["dir"] / rep
        pdb = d / "ref_nw.pdb"
        xtc = d / "complexmd_fit_nw.xtc"
        if not (pdb.exists() and xtc.exists()):
            print(f"  SKIP: {rep}")
            continue
        print(f"  {rep} ...", end=" ", flush=True)
        occ, nf = compute_hbond_contacts(pdb, xtc, dist_cutoff=3.5, stride=1)
        hbond_results[rep] = occ
        n_contacts = sum(1 for v in occ.values() if v >= 0.05)
        print(f"{nf} frames, {n_contacts} contacts ≥5%")

    # ── Aggregate to residue level ────────────────────────────────────────────
    res_occ_per_rep = {}
    all_residues    = set()
    for rep, occ in hbond_results.items():
        ro = atom_to_residue(occ)
        res_occ_per_rep[rep] = ro
        for res, frac in ro.items():
            if frac >= 0.10:
                all_residues.add(res)

    sorted_residues = sorted(all_residues, key=res_sort_key)
    print(f"\nH-bond residues (≥10% in ≥1 rep): {sorted_residues}")

    # ── Build DataFrame for heatmap + CSV ────────────────────────────────────
    df = pd.DataFrame(index=sorted_residues, columns=hbond_reps, dtype=float)
    for rep in hbond_reps:
        ro = res_occ_per_rep.get(rep, {})
        for res in sorted_residues:
            df.loc[res, rep] = ro.get(res, 0.0) * 100.0

    csv_path = outdir / "ibufenac_hbond_table.csv"
    df.to_csv(csv_path)
    print(f"Saved: {csv_path.relative_to(project_root)}")
    print("\nH-bond occupancy table (%):")
    print(df.to_string())

    # ── Figure 2: H-bond heatmap ──────────────────────────────────────────────
    if sorted_residues:
        n_res = len(sorted_residues)
        n_rep = len(hbond_reps)
        fig2, ax2 = plt.subplots(figsize=(6, max(4.5, n_res * 0.45 + 2.0)))

        cmap = LinearSegmentedColormap.from_list(
            "hbond", ["#FFFFFF", "#C6DBEF", "#4292C6", "#08306B"])
        mat  = df.values.astype(float)

        im = ax2.imshow(mat, cmap=cmap, vmin=0, vmax=100,
                        aspect="auto", interpolation="nearest")
        ax2.set_xticks(range(n_rep))
        ax2.set_xticklabels(hbond_reps, fontsize=9)
        ax2.set_yticks(range(n_res))
        ax2.set_yticklabels(sorted_residues, fontsize=8.5)

        for i in range(n_res):
            for j in range(n_rep):
                val = mat[i, j]
                if val >= 5:
                    txt_col = "white" if val > 55 else "#222222"
                    ax2.text(j, i, f"{val:.0f}%", ha="center", va="center",
                             fontsize=7.5, color=txt_col, fontweight="bold")

        cbar = plt.colorbar(im, ax=ax2, fraction=0.035, pad=0.02)
        cbar.set_label("H-bond occupancy (%)\n(O1/O2 – protein N/O, d ≤ 3.5 Å)",
                       fontsize=8)
        ax2.set_title(
            "PICALM site473 · Ibufenac (CHEMBL341812)\n"
            "H-bond contact occupancy — engaged replicas (rep1, rep2, rep5, rep6)",
            fontsize=10, fontweight="bold"
        )
        ax2.set_xlabel("Replica", fontsize=9)
        ax2.set_ylabel("Protein residue", fontsize=9)

        plt.tight_layout()
        out_hbond = outdir / "fig_ibufenac_hbond.png"
        plt.savefig(out_hbond, dpi=300, bbox_inches="tight")
        plt.close()
        print(f"\nSaved: {out_hbond.relative_to(project_root)}")
    else:
        print("\nNo H-bond contacts ≥10% found — skipping heatmap.")

    # ── Text summary ──────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("COMPLEX MD RMSD + H-BOND SUMMARY")
    print("=" * 60)
    for sys in SYSTEMS:
        sid = sys["id"]
        print(f"\n{sys['label']}:")
        for rep in sys["reps"]:
            if rep not in rmsd_data.get(sid, {}):
                continue
            t, r = rmsd_data[sid][rep]
            eng  = sys["engage"].get(rep, "?")
            print(f"  {rep}: mean RMSD {r.mean():.2f} Å  max {r.max():.2f} Å  [{eng}]")

    if hbond_results:
        print("\nIbufenac H-bond summary (top contacts per rep):")
        for rep in hbond_reps:
            ro  = res_occ_per_rep.get(rep, {})
            top = sorted(ro.items(), key=lambda x: -x[1])[:5]
            print(f"  {rep}: " + "  ".join(f"{r} {v*100:.0f}%" for r, v in top))


if __name__ == "__main__":
    main()
