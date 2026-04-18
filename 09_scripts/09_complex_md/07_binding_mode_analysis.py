"""
CRYPTAD — Step 4.6: Binding mode analysis for final hits.

Contact types computed:
  Contact     — any lig heavy atom within 4.0 Å of protein heavy atom
  Hydrophobic — lig C within 4.5 Å of protein C (non-polar)
  HBond       — MDAnalysis HydrogenBondAnalysis (vectorised, D-H···A ≤ 3.5 Å, ≥ 150°)

Per-residue occupancy over all trajectory frames (≥ 20% reported).
Representative frame = frame nearest to median ligand RMSD.
PyMOL .pml scripts for 3D visualisation.

Outputs (03_pocket_analysis/binding_mode/):
  contacts_all.csv           All contacts occ ≥ 20%
  contacts_summary.csv       Per hit/pocket summary
  {pocket}_{chembl}_rep.pdb  Representative frame PDB
  {pocket}_{chembl}.pml      PyMOL script

Usage:
  python3 09_scripts/09_complex_md/07_binding_mode_analysis.py
  python3 09_scripts/09_complex_md/07_binding_mode_analysis.py --project-root /path/to/CRYPTAD
"""

import argparse
import csv
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

warnings.filterwarnings("ignore")
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

CUT_CONTACT = 4.0   # Å
CUT_HYDRO   = 4.5   # Å
OCC_THRESH  = 0.20
TOP_N       = 3     # top hits per pocket (by dG_mean) for full analysis

# Ibufenac at site473 is the MD lead — always include even if outside top-3
FORCE_INCLUDE = [("S3_site473", "CHEMBL341812", "IBUFENAC")]


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD binding mode analysis (Step 4.6)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    mmgbsa_base  = project_root / "04_virtual_screening" / "mmgbsa"
    hits_csv     = project_root / "04_virtual_screening" / "mmgbsa_results" / "final_hits.csv"
    outdir       = project_root / "03_pocket_analysis" / "binding_mode"
    outdir.mkdir(parents=True, exist_ok=True)

    # ── Load final hits ───────────────────────────────────────────────────────
    df        = pd.read_csv(hits_csv)
    df_sorted = df.sort_values(["pocket", "dG_mean"])

    hits = []
    for pocket, grp in df_sorted.groupby("pocket"):
        top = list(zip(grp["pocket"], grp["chembl_id"], grp["name"]))[:TOP_N]
        hits.extend(top)
    present = {(p, c) for p, c, _ in hits}
    for entry in FORCE_INCLUDE:
        if (entry[0], entry[1]) not in present:
            hits.append(entry)

    print(f"Analysing {len(hits)} hit/pocket pairs.\n")

    # ── Per-trajectory analysis ───────────────────────────────────────────────
    all_contacts = []
    summaries    = []

    for pocket, chembl_id, name in hits:
        wdir = mmgbsa_base / pocket / chembl_id
        tpr  = wdir / "prod.tpr"
        xtc  = wdir / "prod_pbc.xtc"

        if not tpr.exists() or not xtc.exists():
            print(f"  [SKIP] {pocket}/{chembl_id} — trajectory missing")
            continue

        print(f"  {pocket}/{chembl_id} ({name})")

        try:
            u = mda.Universe(str(tpr), str(xtc))
        except Exception as e:
            print(f"    [ERROR] load: {e}"); continue

        lig  = u.select_atoms("resname LIG MOL")
        prot = u.select_atoms("protein")

        if len(lig) == 0:
            print(f"    [ERROR] ligand not found"); continue

        n_frames = len(u.trajectory)

        # ── H-bonds via built-in vectorised analysis ──────────────────────────
        hb = HydrogenBondAnalysis(
            universe=u,
            donors_sel="(resname MOL and name N* O* S*) or (protein and name N* O* S*)",
            acceptors_sel="(resname MOL and name N* O* S*) or (protein and name N* O* S*)",
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150.0,
            update_selections=False,
        )
        try:
            hb.run(verbose=False)
            hb_results = hb.results.hbonds
        except Exception as e:
            print(f"    [WARN] HBond analysis failed: {e}")
            hb_results = np.empty((0, 6))

        idx_to_res = {a.index: (a.resid, a.resname, a.segid) for a in u.atoms}

        hb_frames = {}
        if len(hb_results) > 0:
            lig_resnames = set(lig.residues.resnames)
            for row in hb_results:
                frame = int(row[0])
                d_idx = int(row[1])
                a_idx = int(row[3])
                d_res = idx_to_res.get(d_idx)
                a_res = idx_to_res.get(a_idx)
                d_rn  = d_res[1] if d_res else ""
                a_rn  = a_res[1] if a_res else ""
                if d_rn in lig_resnames and a_rn not in lig_resnames:
                    key = (*a_res, "HBond")
                elif a_rn in lig_resnames and d_rn not in lig_resnames:
                    key = (*d_res, "HBond")
                else:
                    continue
                hb_frames.setdefault(key, set()).add(frame)

        # ── Contact + Hydrophobic: vectorised over trajectory ─────────────────
        lig_heavy  = lig.select_atoms("not name H*")
        prot_heavy = prot.select_atoms("not name H*")
        lig_C      = lig.select_atoms("name C*")
        prot_C     = prot.select_atoms("name C*")

        contact_frames = {}

        for ts in u.trajectory:
            dm      = distances.distance_array(lig_heavy.positions, prot_heavy.positions)
            close_j = np.unique(np.where(dm < CUT_CONTACT)[1])
            seen    = set()
            for j in close_j:
                a   = prot_heavy[j]
                key = (a.resid, a.resname, a.segid, "Contact")
                if key not in seen:
                    contact_frames[key] = contact_frames.get(key, 0) + 1
                    seen.add(key)

            if len(lig_C) > 0 and len(prot_C) > 0:
                hd       = distances.distance_array(lig_C.positions, prot_C.positions)
                hclose_j = np.unique(np.where(hd < CUT_HYDRO)[1])
                for j in hclose_j:
                    a   = prot_C[j]
                    key = (a.resid, a.resname, a.segid, "Hydrophobic")
                    if key not in seen:
                        contact_frames[key] = contact_frames.get(key, 0) + 1
                        seen.add(key)

        # ── Representative frame ──────────────────────────────────────────────
        lig_heavy_pos = []
        for ts in u.trajectory:
            lig_heavy_pos.append(lig_heavy.positions.copy())
        lig_heavy_pos = np.array(lig_heavy_pos)
        ref       = lig_heavy_pos[0]
        rmsds     = np.sqrt(np.mean(np.sum((lig_heavy_pos - ref)**2, axis=-1), axis=-1))
        rep_frame = int(np.argmin(np.abs(rmsds - np.median(rmsds))))
        u.trajectory[rep_frame]
        rep_pdb = outdir / f"{pocket}_{chembl_id}_rep.pdb"
        with mda.Writer(str(rep_pdb), multiframe=False) as W:
            W.write(u.select_atoms("protein or resname LIG MOL"))

        # ── Compile contacts ──────────────────────────────────────────────────
        contacts = []
        for (resid, resname, segid, itype), count in contact_frames.items():
            occ = count / n_frames
            if occ >= OCC_THRESH:
                contacts.append({
                    "pocket": pocket, "chembl_id": chembl_id, "name": name,
                    "residue": f"{resname}{resid}.{segid}",
                    "resid": resid, "resname": resname, "chain": segid,
                    "itype": itype, "occupancy": round(occ, 3),
                })
        for (resid, resname, segid, itype), frames in hb_frames.items():
            occ = len(frames) / n_frames
            if occ >= OCC_THRESH:
                contacts.append({
                    "pocket": pocket, "chembl_id": chembl_id, "name": name,
                    "residue": f"{resname}{resid}.{segid}",
                    "resid": resid, "resname": resname, "chain": segid,
                    "itype": itype, "occupancy": round(occ, 3),
                })
        contacts.sort(key=lambda c: (-c["occupancy"], c["itype"]))
        all_contacts.extend(contacts)

        hbonds  = [c for c in contacts if c["itype"] == "HBond"]
        hydrop  = sorted([c for c in contacts if c["itype"] == "Hydrophobic"],
                         key=lambda c: -c["occupancy"])
        generic = [c for c in contacts if c["itype"] == "Contact"]

        print(f"    {n_frames} frames | {len(generic)} contacts | "
              f"{len(hbonds)} H-bonds | {len(hydrop)} hydrophobic  "
              f"[rep frame {rep_frame}, RMSD={rmsds[rep_frame]:.2f}Å]")
        for c in hbonds[:5]:
            print(f"      HB  {c['residue']:<20} {c['occupancy']*100:>5.1f}%")
        for c in hydrop[:3]:
            print(f"      Hφ  {c['residue']:<20} {c['occupancy']*100:>5.1f}%")

        # ── PyMOL script ──────────────────────────────────────────────────────
        contact_resi = "+".join(sorted(set(str(c["resid"]) for c in generic))) or "none"
        hbond_resi   = "+".join(sorted(set(str(c["resid"]) for c in hbonds))) or "none"
        pml_path     = outdir / f"{pocket}_{chembl_id}.pml"
        png_path     = outdir / f"{pocket}_{chembl_id}.png"
        with open(pml_path, "w") as f:
            f.write(f"""\
# CRYPTAD binding mode — {name} @ {pocket}
# Rep frame {rep_frame} | {n_frames} frames | contacts ≥ 20% occupancy

load {rep_pdb}, complex

bg_color white
set cartoon_fancy_helices, 1
hide everything
show cartoon, complex and polymer
show sticks, complex and resname MOL
color slate, complex and polymer
color yellow, complex and resname MOL
set stick_radius, 0.15

select contacts, complex and polymer and resi {contact_resi}
show sticks, contacts
color cyan, contacts
set stick_radius, 0.12, contacts

select hbond_res, complex and polymer and resi {hbond_resi}
color tv_orange, hbond_res

distance hbonds, complex and resname MOL, hbond_res, 3.5, mode=2
set dash_color, firebrick
set dash_width, 2.5

label hbond_res and name CA, "%s%s" % (resn, resi)
set label_size, 11
set label_color, black

center complex and resname MOL
zoom complex and resname MOL, 8
set ray_shadows, 0
set ray_opaque_background, 1
bg_color white
ray 1200, 1200
png {png_path}, dpi=300
quit
""")

        summaries.append({
            "pocket": pocket, "chembl_id": chembl_id, "name": name,
            "n_frames": n_frames, "n_contacts": len(generic),
            "n_hbonds": len(hbonds), "n_hydrophobic": len(hydrop),
            "hbond_residues": "; ".join(
                f"{c['residue']}({c['occupancy']*100:.0f}%)" for c in hbonds[:5]),
            "top_hydrophobic": "; ".join(
                f"{c['residue']}({c['occupancy']*100:.0f}%)" for c in hydrop[:3]),
            "rep_frame": rep_frame,
        })
        print()

    # ── Write outputs ─────────────────────────────────────────────────────────
    contact_fields = ["pocket", "chembl_id", "name", "residue", "resid",
                      "resname", "chain", "itype", "occupancy"]
    out_all = outdir / "contacts_all.csv"
    with open(out_all, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=contact_fields)
        w.writeheader(); w.writerows(all_contacts)
    print(f"→ contacts_all.csv  ({len(all_contacts)} rows)")

    summary_fields = ["pocket", "chembl_id", "name", "n_frames", "n_contacts",
                      "n_hbonds", "n_hydrophobic", "hbond_residues",
                      "top_hydrophobic", "rep_frame"]
    out_summary = outdir / "contacts_summary.csv"
    with open(out_summary, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields)
        w.writeheader(); w.writerows(summaries)
    print(f"→ contacts_summary.csv  ({len(summaries)} rows)")
    print(f"\nAll outputs in: {outdir.relative_to(project_root)}")


if __name__ == "__main__":
    main()
