"""
CRYPTAD — Ibufenac (CHEMBL341812) residue contact analysis at PICALM site473.

NOTE: CHEMBL341812 = Ibufenac (p-isobutylphenylacetic acid, C12H16O2, no nitrogen).
Phensuximide = CHEMBL797 (S3_site473_CHEMBL797); Phensuximide showed 0/3 binding.

For each replica, tracks:
  - Per-frame ligand COM distance from initial position (drift)
  - Per-residue contact occupancy in "engaged" frames (drift < ENGAGED_CUTOFF)
  - Time series of drift to visualise when ligand stays vs. leaves

Outputs (04_virtual_screening/complex_md_results/):
  - ibufenac_contacts.txt  : residue occupancy table per rep
  - ibufenac_drift.csv     : per-frame drift for all reps (for plotting)

Usage:
  python3 09_scripts/09_complex_md/06_ibufenac_contacts.py
  python3 09_scripts/09_complex_md/06_ibufenac_contacts.py --project-root /path/to/CRYPTAD
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

SYS  = "S3_site473_CHEMBL341812"
REPS = ["rep1", "rep2", "rep3", "rep4", "rep5", "rep6"]

CONTACT_CUTOFF = 4.5   # Å — residue counts as contacting ligand
ENGAGED_CUTOFF = 15.0  # Å — frames where ligand is considered pocket-engaged
TOP_N_RESIDUES = 20    # residues to report per rep


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD Ibufenac contact analysis (Step 5.3)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    complexmd_base = project_root / "02_md_simulations" / "complex_md"
    outdir = project_root / "04_virtual_screening" / "complex_md_results"
    outdir.mkdir(parents=True, exist_ok=True)

    drift_records   = []
    contact_results = {}

    for rep in REPS:
        wdir = complexmd_base / SYS / rep
        top  = wdir / "ref_nw.pdb"
        xtc  = wdir / "complexmd_fit_nw.xtc"

        if not top.exists() or not xtc.exists():
            print(f"\n=== {rep} — SKIP (ref_nw.pdb or complexmd_fit_nw.xtc missing) ===")
            continue

        print(f"\n=== {rep} ===")
        u = mda.Universe(str(top), str(xtc))

        lig     = u.select_atoms("resname MOL")
        protein = u.select_atoms("protein")

        # Frame-0 ligand COM reference (use geometry — MOL has unknown elements/masses)
        u.trajectory[0]
        com0 = lig.center_of_geometry()

        drift_vals  = []
        res_contacts = {}  # (resid, resname, segid) -> frame count
        n_engaged   = 0

        for ts in u.trajectory:
            com   = lig.center_of_geometry()
            drift = np.linalg.norm(com - com0)
            drift_vals.append(drift)

            if drift < ENGAGED_CUTOFF:
                n_engaged += 1
                dist_matrix      = distances.distance_array(
                    protein.positions, lig.positions, box=None
                )
                min_dist_per_atom = dist_matrix.min(axis=1)

                contacted_this_frame = set()
                for i, atom in enumerate(protein):
                    if min_dist_per_atom[i] <= CONTACT_CUTOFF:
                        rid = (atom.resid, atom.resname, atom.segid)
                        contacted_this_frame.add(rid)
                for rid in contacted_this_frame:
                    res_contacts[rid] = res_contacts.get(rid, 0) + 1

        drift_vals  = np.array(drift_vals)
        n_frames    = len(drift_vals)
        mean_drift  = float(drift_vals.mean())
        max_drift   = float(drift_vals.max())
        final_drift = float(drift_vals[-1])
        pct_engaged = 100.0 * n_engaged / n_frames

        print(f"  Frames: {n_frames} | Engaged (<{ENGAGED_CUTOFF} Å): "
              f"{n_engaged} ({pct_engaged:.1f}%)")
        print(f"  Drift — mean: {mean_drift:.1f} Å  max: {max_drift:.1f} Å  "
              f"final: {final_drift:.1f} Å")

        for f, d in enumerate(drift_vals):
            drift_records.append({"rep": rep, "frame": f, "drift_A": round(float(d), 3)})

        if n_engaged > 0:
            sorted_res = sorted(res_contacts.items(), key=lambda x: -x[1])
            contact_results[rep] = {
                "n_frames":   n_frames,
                "n_engaged":  n_engaged,
                "pct_engaged": pct_engaged,
                "mean_drift": mean_drift,
                "max_drift":  max_drift,
                "final_drift": final_drift,
                "residues":   [(rid, cnt, 100.0 * cnt / n_engaged)
                               for rid, cnt in sorted_res[:TOP_N_RESIDUES]],
            }
        else:
            contact_results[rep] = {
                "n_frames": n_frames, "n_engaged": 0, "pct_engaged": 0.0,
                "mean_drift": mean_drift, "max_drift": max_drift,
                "final_drift": final_drift, "residues": [],
            }

    # ── Write drift CSV ───────────────────────────────────────────────────────
    drift_df  = pd.DataFrame(drift_records)
    drift_csv = outdir / "ibufenac_drift.csv"
    drift_df.to_csv(drift_csv, index=False)
    print(f"\nDrift CSV written: {drift_csv.relative_to(project_root)}")

    # ── Write contacts report ─────────────────────────────────────────────────
    report_path = outdir / "ibufenac_contacts.txt"
    with open(report_path, "w") as f:
        f.write("CRYPTAD — Ibufenac (CHEMBL341812) residue contact analysis\n")
        f.write(f"System: {SYS}\n")
        f.write(f"Contact cutoff: {CONTACT_CUTOFF} Å | Engaged cutoff: {ENGAGED_CUTOFF} Å\n")
        f.write("=" * 70 + "\n\n")

        for rep in REPS:
            if rep not in contact_results:
                f.write(f"--- {rep} --- SKIPPED (files missing)\n\n")
                continue
            r = contact_results[rep]
            f.write(f"--- {rep} ---\n")
            f.write(f"  Total frames  : {r['n_frames']}\n")
            f.write(f"  Engaged frames: {r['n_engaged']} ({r['pct_engaged']:.1f}%)\n")
            f.write(f"  Drift (mean/max/final): "
                    f"{r['mean_drift']:.1f} / {r['max_drift']:.1f} / {r['final_drift']:.1f} Å\n")
            if r["residues"]:
                f.write(f"\n  Top residues by contact occupancy (in engaged frames):\n")
                f.write(f"  {'Resid':<8} {'Resname':<8} {'Chain':<6} {'Frames':<8} {'Occupancy'}\n")
                f.write(f"  {'-'*45}\n")
                for (resid, resname, segid), cnt, occ in r["residues"]:
                    f.write(f"  {resid:<8} {resname:<8} {segid:<6} {cnt:<8} {occ:.1f}%\n")
            else:
                f.write("  No engaged frames — ligand fully dissociated.\n")
            f.write("\n")

    print(f"Contact report written: {report_path.relative_to(project_root)}")
    print("\nDone.")


if __name__ == "__main__":
    main()
