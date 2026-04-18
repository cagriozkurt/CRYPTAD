"""
CRYPTAD — PBC fixer for complex MD trajectories.

Uses MDAnalysis per-frame unwrap(compound='fragments') which handles each
protein chain independently.  This avoids the cumulative drift that the
transformation-pipeline approach accumulates over frames, which caused
BOND=*** (Fortran overflow) in 1-6 frames per compound and gmx_MMPBSA failure.

Usage (single trajectory — called by SLURM array):
    python3 09_scripts/09_complex_md/04_fix_pbc_mdanalysis.py <tpr> <xtc> <out_xtc>

Usage (auto mode — all replicas under 02_md_simulations/complex_md/):
    python3 09_scripts/09_complex_md/04_fix_pbc_mdanalysis.py
    python3 09_scripts/09_complex_md/04_fix_pbc_mdanalysis.py --project-root /path/to/CRYPTAD
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import MDAnalysis as mda

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]


def fix_xtc(top_path: str, xtc_path: str, out_path: str) -> int:
    """
    Unwrap PBC for one trajectory.

    Parameters
    ----------
    top_path : str  — topology file (TPR preferred; GRO accepted)
    xtc_path : str  — input trajectory
    out_path : str  — output (PBC-fixed) trajectory

    Returns number of broken residues found in frame 0 of the output (0 = success).
    """
    u = mda.Universe(top_path, xtc_path)

    protein     = u.select_atoms("protein")
    system      = u.select_atoms("all")
    not_protein = u.select_atoms("not protein")

    n = len(u.trajectory)
    print(f"  Writing {n} frames → {os.path.basename(out_path)}", flush=True)
    with mda.Writer(out_path, n_atoms=len(system)) as W:
        for ts in u.trajectory:
            # 1. Make each protein chain whole (bond-graph traversal per fragment)
            protein.atoms.unwrap(compound='fragments', reference='cog', inplace=True)
            # 2. Centre simulation box on protein centre of mass
            box_center = ts.dimensions[:3] / 2.0
            prot_com   = protein.center_of_mass()
            system.atoms.translate(box_center - prot_com)
            # 3. Wrap solvent + ligand residues back into the box
            not_protein.wrap(compound='residues')
            W.write(system)
    print(f"  Done.", flush=True)

    # Quick sanity check — count broken residues in frame 0.
    # Use resindex (not resid) to avoid homodimer false positives where chains A+B
    # share the same residue numbers (e.g. BIN1 BAR 1-260 in both chains).
    u2      = mda.Universe(top_path, out_path)
    prot2   = u2.select_atoms("protein")
    u2.trajectory[0]
    res_atoms = defaultdict(list)
    for i, atom in enumerate(prot2.atoms):
        res_atoms[atom.residue.resindex].append(i)
    coords = prot2.positions
    broken = 0
    for idxs in res_atoms.values():
        if len(idxs) < 2:
            continue
        spread = max(
            max(coords[i][j] for i in idxs) - min(coords[i][j] for i in idxs)
            for j in range(3)
        )
        if spread > 10:
            broken += 1
    print(f"  Sanity check — broken residues in frame 0: {broken}", flush=True)
    return broken


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD PBC fixer (Step 5.2)")
    parser.add_argument("top",  nargs="?", help="Topology file (TPR or GRO)")
    parser.add_argument("xtc",  nargs="?", help="Input trajectory (XTC)")
    parser.add_argument("out",  nargs="?", help="Output PBC-fixed trajectory (XTC)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    # ── Single-trajectory mode (called by SLURM array) ───────────────────────
    if args.top and args.xtc and args.out:
        fix_xtc(args.top, args.xtc, args.out)
        return

    # ── Auto mode — all replicas under 02_md_simulations/complex_md/ ─────────
    project_root  = args.project_root.resolve()
    complexmd_dir = project_root / "02_md_simulations" / "complex_md"

    if not complexmd_dir.exists():
        print(f"[ERROR] complex_md directory not found: {complexmd_dir}", file=sys.stderr)
        sys.exit(1)

    total_fixed  = 0
    total_failed = 0

    for sys_dir in sorted(complexmd_dir.iterdir()):
        if not sys_dir.is_dir():
            continue
        for rep_dir in sorted(sys_dir.iterdir()):
            if not rep_dir.is_dir():
                continue
            top = rep_dir / "complexmd.tpr"
            xtc = rep_dir / "complexmd.xtc"
            out = rep_dir / "complexmd_pbc.xtc"

            if not xtc.exists():
                continue

            label = f"{sys_dir.name}/{rep_dir.name}"
            print(f"\n{'='*60}")
            print(f"  {label}")

            if out.exists():
                print(f"  re-fixing existing complexmd_pbc.xtc")
                out.unlink()

            try:
                broken = fix_xtc(str(top), str(xtc), str(out))
                if broken == 0:
                    total_fixed += 1
                else:
                    print(f"  WARNING: {broken} residues still broken after fix!")
                    total_failed += 1
            except Exception as e:
                print(f"  ERROR: {e}")
                total_failed += 1

    print(f"\n{'='*60}")
    print(f"DONE — fixed: {total_fixed}, failed/skipped: {total_failed}")


if __name__ == "__main__":
    main()
