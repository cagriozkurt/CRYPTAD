"""
extract_cryptic_frames.py  —  CRYPTAD
Extract best representative frames for cryptic pocket candidates and
generate PyMOL inspection scripts.

For each cryptic candidate in cryptic_candidates.csv:
  1. Finds the frame with the highest fpocket score near the site centroid
  2. Extracts that MODEL block from protein_traj_metad.pdb → single PDB
  3. Writes a PyMOL .pml script to visualise the candidate with fpocket spheres

Output → 03_pocket_analysis/cryptic_inspection/{SYS}/
  site{N}_frame{F}.pdb      — representative protein frame
  inspect_cryptic_{SYS}.pml — PyMOL session for all candidates of that system

Usage:
    source .venv/bin/activate
    python3 09_scripts/07_conformational/04_extract_cryptic_frames.py
    python3 09_scripts/07_conformational/04_extract_cryptic_frames.py --sys S3_PICALM_ANTH
"""

import argparse
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

SITE_RADIUS = 8.0   # Å — match pocket observations to site centroid

SYSTEM_PATHS = {
    "S1_BIN1_BAR":     "BIN1/metadynamics/S1_BAR",
    "S2_BIN1_SH3":     "BIN1/metadynamics/S2_SH3",
    "S3_PICALM_ANTH":  "PICALM/metadynamics/S3_ANTH",
    "S4_CD2AP_SH3-2":  "CD2AP/metadynamics/S4_SH3-2",
    "S5_CD2AP_SH3-1":  "CD2AP/metadynamics/S5_SH3-1",
}

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; CYAN  = "\033[96m"; RED    = "\033[91m"
def log(m, c=RESET): print(f"{c}{m}{RESET}")


# ---------------------------------------------------------------------------
def extract_model(pdb_path, model_idx):
    """
    Extract the model_idx-th MODEL block (0-indexed) from a multi-model PDB.
    Returns list of ATOM/HETATM lines (without MODEL/ENDMDL wrappers).
    Also returns the time_ns parsed from the TITLE line.
    """
    atoms = []
    time_ns = None
    current_model = -1
    in_target = False

    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("MODEL"):
                current_model += 1
                in_target = (current_model == model_idx)
                atoms = []
            elif line.startswith("TITLE") and current_model == model_idx:
                m = re.search(r't=\s*([\d.eE+\-]+)', line)
                if m:
                    time_ns = float(m.group(1)) / 1000.0
            elif line.startswith("ENDMDL"):
                if in_target:
                    return atoms, time_ns
                in_target = False
            elif in_target and (line.startswith("ATOM") or line.startswith("HETATM")):
                atoms.append(line)

    return atoms, time_ns


def write_pdb(atom_lines, outpath, remark=""):
    with open(outpath, "w") as fh:
        if remark:
            fh.write(f"REMARK {remark}\n")
        fh.writelines(atom_lines)
        fh.write("END\n")


def make_pymol_script(sys_name, candidates, frame_pdbs, outpath):
    """
    Write a PyMOL .pml that loads all candidate frames and draws
    a sphere at each cryptic site centroid.
    """
    lines = [
        "# CRYPTAD — Cryptic pocket inspection",
        f"# System: {sys_name}",
        "# Colour legend: magenta sphere = cryptic site centroid",
        "reinitialize",
        "",
    ]
    for i, (row, pdb) in enumerate(zip(candidates, frame_pdbs)):
        obj = f"site{int(row['site'])}"
        lines += [
            f"load {pdb.name}, {obj}",
            f"hide everything, {obj}",
            f"show cartoon, {obj}",
            f"color slate, {obj}",
            # Draw pseudoatom at cryptic site centroid
            f"pseudoatom crypt_{obj}, pos=[{row['cx']:.2f},{row['cy']:.2f},{row['cz']:.2f}]",
            f"show spheres, crypt_{obj}",
            f"set sphere_scale, 3.0, crypt_{obj}",
            f"color magenta, crypt_{obj}",
            "",
        ]
    lines += [
        "bg_color white",
        "set ray_shadows, 0",
        "orient",
    ]
    with open(outpath, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD: extract representative frames for cryptic pocket candidates")
    parser.add_argument("--sys", default=None,
                        help="Process a single system only (default: all systems in "
                             "cryptic_candidates.csv)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    pocket_root  = project_root / "03_pocket_analysis" / "fpocket_results"
    metad_dir    = project_root / "03_pocket_analysis" / "metadynamics"
    outdir       = project_root / "03_pocket_analysis" / "cryptic_inspection"

    cands_path = metad_dir / "cryptic_candidates.csv"
    if not cands_path.exists():
        log(f"[ERROR] cryptic_candidates.csv not found at {cands_path}", RED)
        log("        Run 09_scripts/07_conformational/03_analyze_metad_pockets.py first.", RED)
        return

    cands = pd.read_csv(cands_path)
    cands = cands[cands["tier"] != "Noise"]
    if args.sys:
        cands = cands[cands["system"] == args.sys]
    if cands.empty:
        log("No cryptic candidates to process.", YELLOW)
        return

    log(f"\n{BOLD}Extracting representative frames for {len(cands)} cryptic candidates{RESET}")

    manifest_entries = []

    for sys_name, grp in cands.groupby("system"):
        log(f"\n{CYAN}▶ {sys_name}  ({len(grp)} candidates){RESET}")

        pocket_csv = pocket_root / sys_name / "metad" / "pocket_data.csv"
        traj_pdb   = pocket_root / sys_name / "metad" / "protein_traj.pdb"

        if not pocket_csv.exists():
            log(f"  [SKIP] pocket_data.csv not found", RED)
            continue
        if not traj_pdb.exists():
            log(f"  [SKIP] protein_traj.pdb not found", RED)
            continue

        df = pd.read_csv(pocket_csv)
        df = df[(df["cx"].abs() + df["cy"].abs() + df["cz"].abs()) > 0]

        out_sys = outdir / sys_name
        out_sys.mkdir(parents=True, exist_ok=True)

        frame_pdbs = []

        for _, row in grp.iterrows():
            site_id = int(row["site"])
            cx, cy, cz = row["cx"], row["cy"], row["cz"]

            # Find all pocket observations near this site centroid
            dists = np.sqrt(
                (df["cx"] - cx)**2 + (df["cy"] - cy)**2 + (df["cz"] - cz)**2
            )
            near = df[dists <= SITE_RADIUS]

            if near.empty:
                log(f"  site {site_id}: no observations within {SITE_RADIUS} Å of centroid", YELLOW)
                frame_pdbs.append(None)
                continue

            # Best frame = highest fpocket score among nearby observations
            best_row  = near.loc[near["score"].idxmax()]
            best_frame = int(best_row["frame"])
            best_time  = float(best_row["time_ns"])
            best_score = float(best_row["score"])

            # Extract that MODEL from protein_traj_metad.pdb
            # MODEL index = frame index in the PDB (0-based)
            # Note: protein_traj.pdb uses stride=2 for S1 (every other frame)
            # The 'frame' column in pocket_data.csv is the original frame index
            # (0, 2, 4, ... for stride=2), so MODEL index = frame / stride
            # But run_fpocket_frames.py tracks 'frame' as MODEL index directly
            model_idx = best_frame

            atom_lines, t_ns = extract_model(traj_pdb, model_idx)
            if not atom_lines:
                log(f"  site {site_id}: could not extract MODEL {model_idx}", YELLOW)
                frame_pdbs.append(None)
                continue

            out_pdb = out_sys / f"site{site_id}_frame{best_frame}.pdb"
            remark  = (f"CRYPTAD cryptic candidate | sys={sys_name} "
                       f"site={site_id} frame={best_frame} "
                       f"t={best_time:.1f}ns score={best_score:.4f} "
                       f"drugg={row['mean_drugg']:.3f} vol={row['mean_volume']:.0f}A3 "
                       f"centroid=({cx:.1f},{cy:.1f},{cz:.1f})")
            write_pdb(atom_lines, out_pdb, remark)
            frame_pdbs.append(out_pdb)

            log(f"  site {site_id}: best frame {best_frame} "
                f"(t={best_time:.1f} ns, score={best_score:.3f}) → {out_pdb.name}",
                GREEN)
            manifest_entries.append({
                "system":      sys_name,
                "site":        site_id,
                "frame":       best_frame,
                "time_ns":     round(best_time, 3),
                "score":       round(best_score, 4),
                "output_pdb":  str(out_pdb),
            })

        # PyMOL script
        valid = [(r, p) for r, p in zip(grp.itertuples(), frame_pdbs) if p is not None]
        if valid:
            pml_rows = [pd.Series(r._asdict()) for r, _ in valid]
            pml_pdbs = [p for _, p in valid]
            pml_out  = out_sys / f"inspect_cryptic_{sys_name}.pml"
            make_pymol_script(sys_name, pml_rows, pml_pdbs, pml_out)
            log(f"  PyMOL script → {pml_out.name}", GREEN)
            log(f"  Open with: pymol {pml_out}")

    log(f"\n{BOLD}Done. Frames in: {outdir}{RESET}")
    log("Visual inspection checklist:")
    log("  1. Is the pocket enclosed (not solvent-exposed surface groove)?")
    log("  2. Is it away from crystal contacts / artefact-prone regions?")
    log("  3. Does it have hydrophobic interior + polar rim (druggable geometry)?")
    log("  4. Unbiased persistence gate: seed 3×20 ns unbiased MD from this frame")
    log("     → pocket persists in ≥1 run → confirmed cryptic pocket for VS")

    # Write manifest
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "generated_at":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":        "09_scripts/07_conformational/04_extract_cryptic_frames.py",
        "project_root":  str(project_root),
        "site_radius_A": SITE_RADIUS,
        "n_extracted":   len(manifest_entries),
        "frames":        manifest_entries,
    }
    manifest_path = outdir / "extract_cryptic_frames_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
