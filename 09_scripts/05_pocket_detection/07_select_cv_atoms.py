#!/usr/bin/env python3
"""
select_cv_atoms.py  —  CRYPTAD
Select CV (collective variable) atom pairs for metadynamics from ref.pdb files.

For S1 BIN1 BAR:
  CV1 = DISTANCE between two concave-face Cα atoms (one per chain) — inter-monomer jaw
  CV2 = ANGLE  between tip_A — center — tip_B along BAR long axis (hinge bending)

For S2–S5 (SH3 / ANTH domains):
  One DISTANCE CV per targeted pocket: the maximum-distance Cα pair within
  ENV_RADIUS of the pocket centroid, constrained to [JAW_MIN, JAW_MAX] Å.

Output: 03_pocket_analysis/cv_atoms_recommended.txt

Usage:
  python 09_scripts/05_pocket_detection/07_select_cv_atoms.py
  python 09_scripts/05_pocket_detection/07_select_cv_atoms.py --project-root /path/to/CRYPTAD
"""

import argparse
import json
import sys
import time
import numpy as np
from pathlib import Path

# ─── Constants ────────────────────────────────────────────────────────────────
ENV_RADIUS  = 8.0   # Å — pocket environment radius for SH3/ANTH jaw search
JAW_MIN     = 8.0   # Å — minimum useful CV distance for SH3/ANTH
JAW_MAX     = 22.0  # Å — maximum useful CV distance for SH3/ANTH
BAR_JAW_MIN = 10.0  # Å — S1 BAR concave inter-chain jaw (larger domain)
BAR_JAW_MAX = 40.0  # Å — S1 BAR concave inter-chain jaw
CONCAVE_FRAC = 0.40  # fraction of CA atoms considered on concave face


# ─── PDB parsing ─────────────────────────────────────────────────────────────

def parse_pdb_ca(pdb_path):
    """
    Parse CA atoms from a PDB file (single model or first MODEL block).
    Returns list of (chain, resnum, resname, xyz_array).
    """
    atoms = []
    in_model = False

    def _parse_atom(line):
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            return None
        if line[12:16].strip() != "CA":
            return None
        chain   = line[21]
        resnum  = int(line[22:26].strip())
        resname = line[17:20].strip()
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        return (chain, resnum, resname, np.array([x, y, z]))

    with open(pdb_path) as fh:
        for raw in fh:
            rec = raw[:6].strip()
            if rec == "MODEL":
                if in_model:
                    break          # stop after first MODEL
                in_model = True
                continue
            if rec == "ENDMDL":
                break
            a = _parse_atom(raw)
            if a:
                atoms.append(a)

    # If no MODEL records, parse the whole file
    if not atoms:
        with open(pdb_path) as fh:
            for raw in fh:
                a = _parse_atom(raw)
                if a:
                    atoms.append(a)

    return atoms


# ─── CV atom selection ────────────────────────────────────────────────────────

def jaw_pair(atoms, centroid, env_radius=ENV_RADIUS,
             jaw_min=JAW_MIN, jaw_max=JAW_MAX):
    """
    Find the maximum-distance Cα pair within env_radius of centroid,
    with distance in [jaw_min, jaw_max].
    Returns (atom1, atom2, distance) or None.
    Falls back to overall max-distance pair with a printed warning.
    """
    c = np.array(centroid)
    env = [a for a in atoms if np.linalg.norm(a[3] - c) <= env_radius]
    if len(env) < 2:
        print(f"    WARNING: only {len(env)} CA within {env_radius} Å of centroid — "
              f"try increasing ENV_RADIUS")
        return None

    best, best_d = None, -1.0
    for i in range(len(env)):
        for j in range(i + 1, len(env)):
            d = np.linalg.norm(env[i][3] - env[j][3])
            if jaw_min <= d <= jaw_max and d > best_d:
                best_d = d
                best = (env[i], env[j], d)

    if best is None:
        # Relax: return the widest pair with a warning
        for i in range(len(env)):
            for j in range(i + 1, len(env)):
                d = np.linalg.norm(env[i][3] - env[j][3])
                if d > best_d:
                    best_d = d
                    best = (env[i], env[j], d)
        if best:
            print(f"    WARNING: no pair in [{jaw_min},{jaw_max}] Å — "
                  f"returning max-dist pair ({best_d:.1f} Å); consider widening range")

    return best


def concave_jaw_pair(atoms, convex_centroid,
                     jaw_min=BAR_JAW_MIN, jaw_max=BAR_JAW_MAX,
                     concave_frac=CONCAVE_FRAC):
    """
    For S1 BIN1 BAR dimer: find an inter-chain Cα jaw pair on the CONCAVE face.

    Method:
      1. Compute protein COM from all CA atoms.
      2. 'Convex direction' = unit vector from COM toward convex_centroid.
      3. Project every CA atom onto this direction.
      4. Atoms in the bottom concave_frac (lowest projections) face OPPOSITE
         the convex surface → concave / membrane-binding face.
      5. Among those concave atoms, find the inter-chain (A×B) pair with the
         largest distance in [jaw_min, jaw_max].
    """
    coords    = np.array([a[3] for a in atoms])
    com       = coords.mean(axis=0)
    conv_dir  = np.array(convex_centroid) - com
    conv_dir /= np.linalg.norm(conv_dir)

    projs     = sorted([(np.dot(a[3] - com, conv_dir), a) for a in atoms],
                       key=lambda x: x[0])

    n_keep    = max(10, int(len(projs) * concave_frac))
    concave   = [a for _, a in projs[:n_keep]]

    chains = {a[0] for a in concave}
    if len(chains) < 2:
        # Try wider window
        n_keep  = int(len(projs) * 0.60)
        concave = [a for _, a in projs[:n_keep]]
        chains  = {a[0] for a in concave}
        print(f"  NOTE: expanded concave window to 60%; chains now: {chains}")

    ch_a = [a for a in concave if a[0] == "A"]
    ch_b = [a for a in concave if a[0] == "B"]
    print(f"  Concave face: {len(ch_a)} chain-A atoms, {len(ch_b)} chain-B atoms")

    best, best_d = None, -1.0
    for ai in ch_a:
        for bi in ch_b:
            d = np.linalg.norm(ai[3] - bi[3])
            if jaw_min <= d <= jaw_max and d > best_d:
                best_d = d
                best = (ai, bi, d)

    if best is None:
        print(f"  WARNING: no inter-chain pair in [{jaw_min},{jaw_max}] Å — "
              f"returning widest inter-chain pair")
        for ai in ch_a:
            for bi in ch_b:
                d = np.linalg.norm(ai[3] - bi[3])
                if d > best_d:
                    best_d = d
                    best = (ai, bi, d)

    return best


def bar_hinge_atoms(atoms):
    """
    Identify 3 Cα atoms for the BAR hinge ANGLE CV along the long axis.
    Uses PCA on all CA coordinates:
      tip_A  = atom with most-negative projection (one banana tip)
      center = atom nearest the median projection  (bend vertex)
      tip_B  = atom with most-positive projection (other banana tip)
    Returns (tip_A, center, tip_B) as full atom tuples.
    """
    coords   = np.array([a[3] for a in atoms])
    com      = coords.mean(axis=0)
    _, _, Vt = np.linalg.svd(coords - com, full_matrices=False)
    long_ax  = Vt[0]

    projs = sorted([(np.dot(a[3] - com, long_ax), a) for a in atoms],
                   key=lambda x: x[0])

    tip_a  = projs[0][1]
    tip_b  = projs[-1][1]
    median = np.median([p for p, _ in projs])
    center = min(projs, key=lambda x: abs(x[0] - median))[1]

    return tip_a, center, tip_b


# ─── Formatting ───────────────────────────────────────────────────────────────

def fmt(atom):
    ch, rn, rname, xyz = atom
    return (f"chain {ch}  res {rn:4d} {rname}"
            f"  xyz=({xyz[0]:.2f},{xyz[1]:.2f},{xyz[2]:.2f})")


# ─── Target pocket definitions (centroids from top_sites CSVs + visual inspection)
# pdb paths are built at runtime relative to --project-root

_POCKET_DEFS = {
    # S2 BIN1 SH3 — Sites 4 + 5 (novel lateral face; user confirmed biologically relevant)
    "S2_BIN1_SH3": [
        {"name": "Site4 (lateral face)",  "centroid": (21.76, 32.28, 32.98)},
        {"name": "Site5 (lateral face)",  "centroid": (24.98, 25.30, 35.80)},
    ],
    # S3 PICALM ANTH — Sites 1 + 4 (PI(4,5)P2 loop region)
    "S3_PICALM_ANTH": [
        {"name": "Site1 (PIP2 loop face)", "centroid": (41.26, 42.15, 42.59)},
        {"name": "Site4 (PIP2 loop face)", "centroid": (49.25, 55.70, 53.12)},
    ],
    # S4 CD2AP SH3-2 — Sites 1 + 2 (RT-loop flank + groove)
    "S4_CD2AP_SH3-2": [
        {"name": "Site1 (RT-loop flank)", "centroid": (37.00, 27.26, 30.56)},
        {"name": "Site2 (groove)",         "centroid": (28.40, 20.88, 26.86)},
    ],
    # S5 CD2AP SH3-1 — Sites 1 + 4 (groove + RT-loop flank; visual rank order)
    "S5_CD2AP_SH3-1": [
        {"name": "Site1 (groove face)",   "centroid": (29.04, 28.65, 19.88)},
        {"name": "Site4 (RT-loop flank)", "centroid": (32.74, 36.65, 29.05)},
    ],
}

# S1 BIN1 BAR — geometric approach; fpocket found only convex-face artefacts
S1_CONVEX_CENTROID = (81.84, 71.23, 74.61)   # Site 1 in top_sites (convex face)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    _script = Path(__file__).resolve()
    _default_root = _script.parents[2]

    parser = argparse.ArgumentParser(
        description="CRYPTAD: select CV atom pairs from pocket inspection PDBs")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    insp         = project_root / "03_pocket_analysis" / "inspection"
    out_file     = project_root / "03_pocket_analysis" / "cv_atoms_recommended.txt"

    # Build runtime TARGETS dict with pdb paths from project_root
    targets = {
        sys_name: {
            "pdb":     insp / sys_name / "ref.pdb",
            "pockets": pockets,
        }
        for sys_name, pockets in _POCKET_DEFS.items()
    }
    s1_pdb = insp / "S1_BIN1_BAR" / "ref.pdb"

    # Verify required PDBs exist
    missing = [str(cfg["pdb"]) for cfg in targets.values()
               if not cfg["pdb"].exists()]
    if not s1_pdb.exists():
        missing.append(str(s1_pdb))
    if missing:
        print("[WARN] The following ref.pdb files are missing "
              "(run 06_inspect_pockets.py first):")
        for m in missing:
            print(f"  {m}")

    lines = []

    def h(title):
        lines.append("")
        lines.append("=" * 72)
        lines.append(title)
        lines.append("=" * 72)

    lines.append("CRYPTAD — CV Atom Recommendations")
    lines.append(f"Generated by: 09_scripts/05_pocket_detection/07_select_cv_atoms.py")
    lines.append(f"Generated at: {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    lines.append("All Cα atoms; distances in Å.")
    lines.append("GROMACS indices (1-based) must be looked up with gmx select or gmx dump.")

    # ── S1 BIN1 BAR ──────────────────────────────────────────────────────────
    h("S1_BIN1_BAR — concave-face geometric selection")
    print("\n=== S1_BIN1_BAR ===")
    if not s1_pdb.exists():
        print(f"  [SKIP] {s1_pdb} not found")
        lines.append(f"  SKIPPED: {s1_pdb} not found")
        atoms = []
    else:
        atoms = parse_pdb_ca(s1_pdb)
    print(f"  {len(atoms)} CA atoms loaded")

    # CV1 — inter-monomer concave jaw (DISTANCE)
    result = concave_jaw_pair(atoms, S1_CONVEX_CENTROID) if atoms else None
    if result:
        a1, a2, d = result
        lines.append("")
        lines.append("CV1  DISTANCE — inter-monomer concave jaw")
        lines.append(f"  Atom 1 : {fmt(a1)}")
        lines.append(f"  Atom 2 : {fmt(a2)}")
        lines.append(f"  Ref dist: {d:.2f} Å")
        lines.append("  PLUMED : DISTANCE ATOMS=IDX1,IDX2 LABEL=cv1")
        print(f"  CV1  {fmt(a1)}")
        print(f"   <-> {fmt(a2)}   d={d:.2f} Å")
    else:
        lines.append("\n  ERROR: could not find concave jaw pair")

    # CV2 — BAR hinge angle (ANGLE)
    if not atoms:
        lines.append("\n  SKIPPED CV2: no atoms loaded")
        tip_a = center = tip_b = None
    else:
        tip_a, center, tip_b = bar_hinge_atoms(atoms)
    lines.append("")
    lines.append("CV2  ANGLE — BAR long-axis hinge (tip_A — center — tip_B)")
    if tip_a is not None:
        lines.append(f"  tip_A : {fmt(tip_a)}")
        lines.append(f"  center: {fmt(center)}")
        lines.append(f"  tip_B : {fmt(tip_b)}")
        lines.append("  PLUMED : ANGLE ATOMS=IDX_A,IDX_C,IDX_B LABEL=cv2")
        print(f"  CV2  tip_A  = {fmt(tip_a)}")
        print(f"       center = {fmt(center)}")
        print(f"       tip_B  = {fmt(tip_b)}")

    # ── S2–S5 ────────────────────────────────────────────────────────────────
    for sys_name, cfg in targets.items():
        h(sys_name)
        print(f"\n=== {sys_name} ===")
        if not cfg["pdb"].exists():
            print(f"  [SKIP] {cfg['pdb']} not found")
            lines.append(f"  SKIPPED: ref.pdb not found")
            continue
        atoms = parse_pdb_ca(cfg["pdb"])
        print(f"  {len(atoms)} CA atoms loaded")

        for pk in cfg["pockets"]:
            c = pk["centroid"]
            print(f"\n  Pocket: {pk['name']}  centroid=({c[0]},{c[1]},{c[2]})")
            result = jaw_pair(atoms, c)
            lines.append("")
            lines.append(f"Pocket: {pk['name']}")
            lines.append(f"  centroid: ({c[0]}, {c[1]}, {c[2]})")
            if result:
                a1, a2, d = result
                lines.append(f"  Atom 1 : {fmt(a1)}")
                lines.append(f"  Atom 2 : {fmt(a2)}")
                lines.append(f"  Ref dist: {d:.2f} Å")
                lines.append("  PLUMED : DISTANCE ATOMS=IDX1,IDX2 LABEL=cv_...")
                print(f"    Jaw: {fmt(a1)}")
                print(f"     <-> {fmt(a2)}   d={d:.2f} Å")
            else:
                lines.append("  WARNING: no jaw pair found — check centroid / ENV_RADIUS")
                print("    WARNING: no jaw pair found")

    # ── Save ─────────────────────────────────────────────────────────────────
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(lines) + "\n")
    print(f"\nSaved → {out_file}")

    # Write manifest
    manifest = {
        "generated_at":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "project_root":  str(project_root),
        "output":        str(out_file),
        "env_radius_A":  ENV_RADIUS,
        "jaw_min_A":     JAW_MIN,
        "jaw_max_A":     JAW_MAX,
        "bar_jaw_min_A": BAR_JAW_MIN,
        "bar_jaw_max_A": BAR_JAW_MAX,
        "concave_frac":  CONCAVE_FRAC,
        "systems":       ["S1_BIN1_BAR"] + list(_POCKET_DEFS.keys()),
    }
    manifest_path = out_file.parent / "cv_atoms_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
