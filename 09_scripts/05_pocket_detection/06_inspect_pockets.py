#!/usr/bin/env python3
"""
CRYPTAD — Visual inspection helper for top fpocket sites

For each system, identifies the top Tier-2 pocket sites from the cross-replica
analysis, extracts the single best open-pocket frame per site from the trajectory,
and generates a ready-to-run PyMOL inspection script.

Two modes:
  Default    : Extract the best frame from protein_traj.pdb for each site
               (requires protein_traj.pdb to be available locally — it is).
  --ref-only : Use protein_ref.pdb with centroid spheres only; no frame
               extraction needed. Useful for a quick first look.

Output in 03_pocket_analysis/inspection/{SYS}/:
  ref.pdb               — copy of protein_ref.pdb (reference structure)
  site{N}_best.pdb      — best open-pocket frame for site N (if not --ref-only)
  top_sites_{SYS}.csv   — top N sites table for this system
  inspect_{SYS}.pml     — PyMOL script: load structures + markers + environment

Usage:
  python 09_scripts/05_pocket_detection/06_inspect_pockets.py
  python 09_scripts/05_pocket_detection/06_inspect_pockets.py --sys S2_BIN1_SH3
  python 09_scripts/05_pocket_detection/06_inspect_pockets.py --top 3
  python 09_scripts/05_pocket_detection/06_inspect_pockets.py --ref-only
  python 09_scripts/05_pocket_detection/06_inspect_pockets.py --runs metad
"""

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SYSTEMS = [
    "S1_BIN1_BAR",
    "S2_BIN1_SH3",
    "S3_PICALM_ANTH",
    "S4_CD2AP_SH3-2",
    "S5_CD2AP_SH3-1",
]
RUNS = ["run1", "run2", "run3"]

# PyMOL color names for sites 1–8 (cycles if more than 8 requested)
SITE_COLORS = ["firebrick", "marine", "forest", "gold", "purple",
               "cyan", "orange", "magenta"]

# Same eps as parse_fpocket.py — must match for correct proximity matching
CLUSTER_EPS = 8.0   # Å
CROSS_DIST  = 10.0  # Å — cross-replica match radius

# Biological context notes shown in the PyMOL script header
SYSTEM_NOTES = {
    "S1_BIN1_BAR":    "Look for sites on the CONCAVE membrane-binding face "
                      "(inner curved surface); avoid the convex outer face.",
    "S2_BIN1_SH3":    "Look for sites adjacent to the PxxP-binding groove "
                      "(RT-loop, n-Src loop, distal loop region).",
    "S3_PICALM_ANTH": "Look for sites near the PI(4,5)P2-binding loops "
                      "(helix α0, β1-α1 linker) and the clathrin-box region.",
    "S4_CD2AP_SH3-2": "Look for sites in or near the RT-loop groove "
                      "(canonical SH3 binding cleft, Pro-rich peptide face).",
    "S5_CD2AP_SH3-1": "Look for sites in or near the RT-loop groove "
                      "(canonical SH3 binding cleft, Pro-rich peptide face).",
}


# ---------------------------------------------------------------------------
# Frame identification and extraction
# ---------------------------------------------------------------------------
def find_best_frame(pocket_csv, run_sites_csv, site_cx, site_cy, site_cz):
    """
    Find the trajectory frame where the pocket at (site_cx, site_cy, site_cz)
    was most druggable/open.

    Strategy:
      1. From {SYS}_{RUN}_sites.csv, find the run-specific site centroid
         nearest to the cross-replica centroid (must be within CROSS_DIST Å).
      2. In pocket_data.csv, select all rows whose centroid is within
         CLUSTER_EPS Å of that run-specific centroid.
      3. Return the row with the highest drugg_score (tie-break: volume).

    Returns a dict with frame info, or None if the site has no match in
    this run.
    """
    df = pd.read_csv(pocket_csv)
    df = df[df["cx"].abs() + df["cy"].abs() + df["cz"].abs() > 0]
    if df.empty:
        return None

    # Step 1 — nearest run-specific site centroid
    run_sites = pd.read_csv(run_sites_csv)
    if run_sites.empty:
        return None

    run_coords = run_sites[["cx", "cy", "cz"]].values
    dists_to_site = cdist([[site_cx, site_cy, site_cz]], run_coords)[0]
    best_site_idx = int(dists_to_site.argmin())
    if dists_to_site[best_site_idx] > CROSS_DIST:
        return None

    rcx = float(run_sites.iloc[best_site_idx]["cx"])
    rcy = float(run_sites.iloc[best_site_idx]["cy"])
    rcz = float(run_sites.iloc[best_site_idx]["cz"])

    # Step 2 — filter pocket_data rows that fall within the cluster radius
    obs_coords = df[["cx", "cy", "cz"]].values
    obs_dists = cdist([[rcx, rcy, rcz]], obs_coords)[0]
    matching = df[obs_dists <= CLUSTER_EPS].copy()
    if matching.empty:
        return None

    # Step 3 — best row: highest drugg_score, tie-break by volume
    best = matching.sort_values(
        ["drugg_score", "volume"], ascending=[False, False]
    ).iloc[0]

    return {
        "frame":       int(best["frame"]),
        "time_ns":     float(best["time_ns"]),
        "drugg_score": float(best["drugg_score"]),
        "volume":      float(best["volume"]),
        "cx": float(best["cx"]),
        "cy": float(best["cy"]),
        "cz": float(best["cz"]),
    }


def extract_frame(traj_pdb_path, frame_idx, output_pdb_path):
    """
    Extract a single MODEL block (0-based frame_idx) from a multi-model PDB.
    Writes ATOM/HETATM/TER lines + END to output_pdb_path.
    Returns True on success, False if the frame was not found.
    """
    current = -1
    lines = []
    found = False

    with open(traj_pdb_path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                current += 1
                if current == frame_idx:
                    found = True
            elif line.startswith("ENDMDL"):
                if found:
                    break
                found = False
            elif found:
                lines.append(line)

    if not lines:
        return False

    with open(output_pdb_path, "w") as fh:
        fh.writelines(lines)
        fh.write("END\n")
    return True


# ---------------------------------------------------------------------------
# PyMOL script generation
# ---------------------------------------------------------------------------
def _pml_pos(cx, cy, cz):
    return f"[{cx:.3f}, {cy:.3f}, {cz:.3f}]"


def write_pymol_script(sys_name, ref_pdb_name, top_sites, site_frames, outpath):
    """
    Write a ready-to-run PyMOL .pml inspection script.

    Parameters
    ----------
    sys_name     : system label (e.g. "S2_BIN1_SH3")
    ref_pdb_name : filename of the reference PDB in the inspection directory
    top_sites    : list of row dicts from selected_pockets.csv (up to 8)
    site_frames  : list of frame-info dicts (or None) parallel to top_sites
    outpath      : Path to write the .pml file
    """
    L = []
    a = L.append

    note = SYSTEM_NOTES.get(sys_name, "")
    n_sites = len(top_sites)
    has_frames = any(f is not None for f in site_frames)

    a(f"# CRYPTAD — Pocket Inspection: {sys_name}")
    a(f"# {n_sites} top cross-replica Tier-2 sites (sorted: replicas desc, persistence desc)")
    a(f"#")
    a(f"# Biological context:")
    a(f"#   {note}")
    a(f"#")
    a(f"# Color guide:")
    for i, (s, c) in enumerate(zip(top_sites, SITE_COLORS), 1):
        runs_label = ", ".join(f"run{r}" for r in s["runs"].split(","))
        a(f"#   {c:<12} — Site {i} | persist={float(s['persistence']):.1%} | "
          f"drugg={float(s['mean_drugg']):.3f} | vol={float(s['mean_volume']):.0f} A3 | "
          f"{runs_label}")
    a(f"#")
    a(f"# Usage: File > Run Script... > select this .pml")
    a(f"#")
    a(f"# Quick exploration commands:")
    a(f"#   hide all, open_frames      show only the reference + centroid markers")
    a(f"#   show sticks, pocket1_env   show stick repr for site-1 pocket residues")
    a(f"#   set sphere_scale, 1.0, centroids   smaller centroid spheres")
    a(f"#   color white, reference     make protein background white for figures")
    a(f"")
    a(f"reinitialize")
    a(f"")

    # --- Reference structure ---
    a(f"# Reference structure (protein_ref.pdb, frame 0 of production)")
    a(f"load {ref_pdb_name}, reference")
    a(f"show_as cartoon, reference")
    a(f"color gray70, reference")
    a(f"show surface, reference")
    a(f"set surface_color, gray85, reference")
    a(f"set transparency, 0.55, reference")
    a(f"set cartoon_smooth_loops, 1, reference")
    a(f"")

    # --- Per-site: best frame + centroid sphere + pocket environment ---
    for i, (site, frame_info, color) in enumerate(
            zip(top_sites, site_frames, SITE_COLORS), 1):

        cx = float(site["cx"])
        cy = float(site["cy"])
        cz = float(site["cz"])
        persist_pct = f"{float(site['persistence']):.1%}"
        drugg       = float(site["mean_drugg"])
        vol         = float(site["mean_volume"])
        runs_label  = ", ".join(f"run{r}" for r in site["runs"].split(","))

        a(f"# --- Site {i} ---  persist={persist_pct}  "
          f"drugg={drugg:.3f}  vol={vol:.0f} A3  ({runs_label})")

        # Load best open-pocket frame if available
        if frame_info is not None:
            a(f"load site{i}_best.pdb, site{i}_open")
            a(f"align site{i}_open, reference")
            a(f"show_as cartoon, site{i}_open")
            a(f"color {color}, site{i}_open")
            a(f"set cartoon_transparency, 0.40, site{i}_open")

        # Centroid pseudoatom sphere
        a(f"pseudoatom centroid_{i}, pos={_pml_pos(cx, cy, cz)}")
        a(f"show spheres, centroid_{i}")
        a(f"set sphere_scale, 2.5, centroid_{i}")
        a(f"color {color}, centroid_{i}")
        a(f'label centroid_{i}, "S{i} ({persist_pct})"')

        # Highlight residues within 8 Å on the reference surface
        a(f"select pocket{i}_env, (byres reference within 8 of centroid_{i})")
        a(f"show surface, pocket{i}_env")
        a(f"color {color}, pocket{i}_env")
        a(f"set transparency, 0.20, pocket{i}_env")
        a(f"")

    # --- Display settings ---
    a(f"# Display settings")
    a(f"set ray_shadow, 0")
    a(f"set bg_rgb, white")
    a(f"set depth_cue, 0")
    a(f"set specular, 0.3")
    a(f"set cartoon_fancy_helices, 1")
    a(f"set label_size, 14")
    a(f"set label_color, black")
    a(f"zoom all, buffer=5")
    a(f"")

    # --- Groups for easy toggling ---
    a(f"# Group objects for easy show/hide in the PyMOL object panel")
    a(f"group centroids, centroid_*")
    a(f"group pocket_envs, pocket*_env")
    if has_frames:
        a(f"group open_frames, site*_open")
    a(f"")

    # --- Suggested selections for specific systems ---
    a(f"# System-specific orientation hint:")
    a(f"# {note}")

    with open(outpath, "w") as fh:
        fh.write("\n".join(L) + "\n")


# ---------------------------------------------------------------------------
# Per-system processor
# ---------------------------------------------------------------------------
def process_system(sys_name, pocket_root, sites_dir, inspect_dir,
                   n_top, ref_only, runs):
    sel_csv = sites_dir / "selected_pockets.csv"
    if not sel_csv.exists():
        print(f"  [SKIP] selected_pockets.csv not found — "
              f"run 05_parse_fpocket.py first")
        return 0

    sel = pd.read_csv(sel_csv)
    sys_sel = sel[sel["system"] == sys_name].sort_values(
        ["n_replicas", "persistence"], ascending=[False, False]
    ).head(n_top).reset_index(drop=True)

    if sys_sel.empty:
        print(f"  No selected pockets found for {sys_name}")
        return 0

    total_available = (sel["system"] == sys_name).sum()
    print(f"  Top {len(sys_sel)} of {total_available} selected cross-replica sites:")
    for i, row in sys_sel.iterrows():
        runs_label = ", ".join(f"run{r}" for r in str(row["runs"]).split(","))
        print(f"    Site {i+1}: persist={row['persistence']:.1%}  "
              f"drugg={row['mean_drugg']:.3f}  vol={row['mean_volume']:.0f} A3  "
              f"({runs_label})  "
              f"({row['cx']:.1f}, {row['cy']:.1f}, {row['cz']:.1f})")

    # Output directory
    out = inspect_dir / sys_name
    out.mkdir(parents=True, exist_ok=True)

    # Copy protein_ref.pdb → ref.pdb
    ref_pdb_name = "ref.pdb"
    for run in runs:
        cand = pocket_root / sys_name / run / "protein_ref.pdb"
        if cand.exists():
            shutil.copy2(cand, out / ref_pdb_name)
            print(f"  Copied {cand.parent.name}/{cand.name} → ref.pdb")
            break
    else:
        print(f"  [WARN] protein_ref.pdb not found for any run of {sys_name}")

    # Save top-sites table
    sys_sel.to_csv(out / f"top_sites_{sys_name}.csv", index=False)

    # Extract best frames
    site_frames = []
    for i, row in sys_sel.iterrows():
        if ref_only:
            site_frames.append(None)
            continue

        run_indices = [int(r) - 1 for r in str(row["runs"]).split(",")]
        frame_info = None

        for run_idx in run_indices:
            if run_idx >= len(runs):
                continue
            run = runs[run_idx]
            pocket_csv    = pocket_root / sys_name / run / "pocket_data.csv"
            run_sites_csv = sites_dir / f"{sys_name}_{run}_sites.csv"
            traj_pdb      = pocket_root / sys_name / run / "protein_traj.pdb"

            if not pocket_csv.exists() or not run_sites_csv.exists():
                continue

            info = find_best_frame(
                pocket_csv, run_sites_csv,
                float(row["cx"]), float(row["cy"]), float(row["cz"])
            )
            if info is None:
                continue

            if not traj_pdb.exists():
                print(f"    Site {i+1}: best is frame {info['frame']} "
                      f"({info['time_ns']:.1f} ns, drugg={info['drugg_score']:.3f}) "
                      f"in {run} — protein_traj.pdb not found locally")
                break

            out_pdb = out / f"site{i+1}_best.pdb"
            ok = extract_frame(traj_pdb, info["frame"], out_pdb)
            if ok:
                info["run"] = run
                frame_info = info
                print(f"    Site {i+1}: frame {info['frame']} "
                      f"({info['time_ns']:.1f} ns, drugg={info['drugg_score']:.3f}) "
                      f"from {run} → site{i+1}_best.pdb")
            else:
                print(f"    Site {i+1}: frame {info['frame']} not found in "
                      f"{run}/protein_traj.pdb")
            break  # use the first run that has the site

        if frame_info is None and not ref_only:
            print(f"    Site {i+1}: no frame extracted")

        site_frames.append(frame_info)

    # Write PyMOL script
    pml_path = out / f"inspect_{sys_name}.pml"
    write_pymol_script(
        sys_name,
        ref_pdb_name,
        sys_sel.to_dict("records"),
        site_frames,
        pml_path,
    )
    print(f"  PyMOL script → {pml_path.relative_to(pml_path.parent.parent.parent)}")

    return len(sys_sel)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    _script = Path(__file__).resolve()
    _default_root = _script.parents[2]

    parser = argparse.ArgumentParser(
        description="CRYPTAD pocket inspection: extract best frames + PyMOL scripts")
    parser.add_argument("--sys",      default=None,
                        help="Single system name (default: all 5 systems)")
    parser.add_argument("--top",      type=int, default=5,
                        help="Top N sites per system (default: 5)")
    parser.add_argument("--ref-only", action="store_true",
                        help="Skip frame extraction; use protein_ref.pdb only")
    parser.add_argument("--runs",     nargs="+", default=list(RUNS),
                        help="Run names to search for protein_ref.pdb and best frames "
                             "(default: run1 run2 run3). Use --runs metad for "
                             "metadynamics data.")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    pocket_root  = project_root / "03_pocket_analysis" / "fpocket_results"
    sites_dir    = project_root / "03_pocket_analysis" / "pocket_sites"
    inspect_dir  = project_root / "03_pocket_analysis" / "inspection"
    inspect_dir.mkdir(parents=True, exist_ok=True)

    if not (sites_dir / "selected_pockets.csv").exists():
        sys.exit(f"[ERROR] {sites_dir / 'selected_pockets.csv'} not found.\n"
                 f"       Run 05_parse_fpocket.py first.")

    systems = [args.sys] if args.sys else SYSTEMS
    mode = "reference-only (centroid markers)" if args.ref_only else \
           "reference + best-frame extraction"

    print(f"CRYPTAD pocket inspection")
    print(f"  Mode  : {mode}")
    print(f"  Runs  : {args.runs}")
    print(f"  Sites : top {args.top} per system (sorted by n_replicas, persistence)")
    print(f"  Output: {inspect_dir}")

    total = 0
    for sys_name in systems:
        print(f"\n{'='*60}")
        print(f"  {sys_name}")
        print(f"{'='*60}")
        total += process_system(sys_name, pocket_root, sites_dir, inspect_dir,
                                args.top, args.ref_only, args.runs)

    print(f"\n===== DONE =====")
    print(f"Total sites prepared: {total}")
    print(f"")
    print(f"Open PyMOL, then for each system:")
    print(f"  File > Run Script > 03_pocket_analysis/inspection/{{SYS}}/inspect_{{SYS}}.pml")
    print(f"")
    print(f"What to look for:")
    print(f"  S1 BIN1 BAR   — sites on the concave membrane-binding face")
    print(f"  S2 BIN1 SH3   — sites adjacent to the PxxP-binding groove")
    print(f"  S3 PICALM ANTH — sites near PI(4,5)P2-binding loops / clathrin-box")
    print(f"  S4/S5 CD2AP   — sites in or flanking the RT-loop SH3 groove")
    print(f"")
    print(f"Record promising sites → these become metadynamics CV targets (max 2/system)")

    # Write manifest
    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "project_root": str(project_root),
        "systems": systems,
        "runs": args.runs,
        "top_n": args.top,
        "ref_only": args.ref_only,
        "total_sites_prepared": total,
    }
    manifest_path = inspect_dir / "inspect_pockets_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"  Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
