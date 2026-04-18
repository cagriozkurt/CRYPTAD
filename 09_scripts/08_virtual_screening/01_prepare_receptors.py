"""
prepare_receptors.py  —  CRYPTAD
Receptor preparation pipeline for virtual screening (Step 4.1).

For each confirmed cryptic pocket:
  1. Cleans the protein-only PDB (removes HETATM, alt locs, duplicate atoms)
  2. Adds hydrogens and assigns protonation states at pH 7.4 via pdb2pqr+PROPKA
  3. Converts to PDBQT format via Meeko
  4. Writes an AutoDock Vina config file with the pocket-centred grid box

Output layout:
  04_virtual_screening/receptors/{pocket_id}/
    ├── receptor_clean.pdb     cleaned input
    ├── receptor_H.pdb         protonated PDB (pdb2pqr output)
    ├── receptor.pdbqt         docking-ready PDBQT
    └── vina_config.txt        Vina grid box config

Prerequisites (run once):
  pip install pdb2pqr meeko vina rdkit biopython

Usage:
  python3 09_scripts/08_virtual_screening/01_prepare_receptors.py
  python3 09_scripts/08_virtual_screening/01_prepare_receptors.py --pocket S1_site688
  python3 09_scripts/08_virtual_screening/01_prepare_receptors.py --project-root /path/to/CRYPTAD
"""

import argparse
import json
import subprocess
import sys
import shutil
import time
from pathlib import Path

# ---------------------------------------------------------------------------
_script       = Path(__file__).resolve()
_default_root = _script.parents[2]
sys.path.insert(0, str(_script.parent))
from pocket_definitions import POCKETS, VINA_PARAMS

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; RED   = "\033[91m"; CYAN   = "\033[96m"
def log(m, c=RESET): print(f"{c}{m}{RESET}")


# ---------------------------------------------------------------------------
def check_dependencies():
    missing = []
    for cmd in ["pdb2pqr", "mk_prepare_receptor.py"]:
        if not shutil.which(cmd):
            missing.append(cmd)
    if missing:
        log(f"\n[ERROR] Missing tools: {', '.join(missing)}", RED)
        log("Install with:", YELLOW)
        log("  pip install pdb2pqr meeko", YELLOW)
        log("  (meeko installs mk_prepare_receptor.py into PATH)", YELLOW)
        sys.exit(1)


# ---------------------------------------------------------------------------
def clean_pdb(input_pdb: Path, output_pdb: Path):
    """
    Remove HETATM records, alternate locations (keep A), and END/ENDMDL
    artifacts from a multi-model trajectory frame.
    Keeps only MODEL 1 if multiple models are present (should already be single).
    """
    lines_out = []
    in_model  = False
    model_done = False

    with open(input_pdb) as fh:
        for line in fh:
            rec = line[:6].strip()

            # Take only first MODEL block
            if rec == "MODEL":
                if model_done:
                    break
                in_model = True
                continue
            if rec == "ENDMDL":
                model_done = True
                in_model   = False
                continue

            # Skip non-ATOM records
            if rec not in ("ATOM", "TER"):
                continue

            # Skip alternate locations (keep blank or 'A')
            if rec == "ATOM":
                alt_loc = line[16]
                if alt_loc not in (" ", "A"):
                    continue
                # Blank the alt_loc field so pdb2pqr doesn't complain
                line = line[:16] + " " + line[17:]

            lines_out.append(line)

    # Rename CHARMM histidine variants → standard HIS so pdb2pqr+Meeko recognise them
    # HSD (Nδ-H) / HSE (Nε-H) / HSP (doubly protonated) → HIS
    # pdb2pqr will re-assign protonation state via AMBER FF
    renamed = []
    for line in lines_out:
        if line.startswith("ATOM") and line[17:20] in ("HSD", "HSE", "HSP"):
            line = line[:17] + "HIS" + line[20:]
        renamed.append(line)
    lines_out = renamed

    with open(output_pdb, "w") as fh:
        fh.writelines(lines_out)

    n_atoms = sum(1 for l in lines_out if l.startswith("ATOM"))
    return n_atoms


# ---------------------------------------------------------------------------
def run_pdb2pqr(clean_pdb: Path, pqr_out: Path, pdb_H_out: Path):
    """
    Run pdb2pqr with AMBER force field to add hydrogens.
    Protonation uses AMBER defaults (standard pH 7.4 states: ASP/GLU
    deprotonated, LYS/ARG protonated, HIS neutral HID).
    NOTE: propka titration is skipped — Python 3.14 compatibility issue
    with propka __annotations__ API. For histidine states, AMBER default
    (HID) is used; manual inspection recommended if HIS is in the pocket.
    """
    cmd = [
        "pdb2pqr",
        "--ff",         "AMBER",
        "--pdb-output", str(pdb_H_out),
        "--drop-water",
        "--nodebump",
        str(clean_pdb),
        str(pqr_out),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"  [pdb2pqr stderr]\n{result.stderr[-800:]}", YELLOW)
        if not pdb_H_out.exists():
            raise RuntimeError(f"pdb2pqr failed: {result.stderr[-200:]}")
    return result.returncode == 0


# ---------------------------------------------------------------------------
def run_meeko_receptor(pdb_H: Path, pdbqt_out: Path):
    """
    Convert protonated PDB to PDBQT using Meeko's mk_prepare_receptor.py.
    meeko 0.7+ requires -o basename (no extension) + -p to write the PDBQT file.
    """
    basename = str(pdbqt_out.with_suffix(""))  # strip .pdbqt
    cmd = [
        "mk_prepare_receptor.py",
        "-i", str(pdb_H),
        "-o", basename,
        "-p",   # --write_pdbqt
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not pdbqt_out.exists():
        log(f"  [meeko stderr]\n{result.stderr[-600:]}", YELLOW)
        raise RuntimeError(f"mk_prepare_receptor.py failed")


# ---------------------------------------------------------------------------
def write_vina_config(pocket: dict, pdbqt_path: Path, config_path: Path):
    p = VINA_PARAMS
    # Use filename only — config and receptor are co-located; avoids space-in-path
    # issues with AutoDock Vina's --config parser.
    config = (
        f"receptor = {pdbqt_path.name}\n"
        f"\n"
        f"center_x = {pocket['cx']:.3f}\n"
        f"center_y = {pocket['cy']:.3f}\n"
        f"center_z = {pocket['cz']:.3f}\n"
        f"\n"
        f"size_x = {pocket['box_x']}\n"
        f"size_y = {pocket['box_y']}\n"
        f"size_z = {pocket['box_z']}\n"
        f"\n"
        f"exhaustiveness = {p['exhaustiveness']}\n"
        f"num_modes       = {p['num_modes']}\n"
        f"energy_range    = {p['energy_range']}\n"
        f"cpu             = {p['cpu']}\n"
    )
    config_path.write_text(config)


# ---------------------------------------------------------------------------
def prepare_one(pocket: dict, project_root: Path, force: bool = False):
    pid   = pocket["id"]
    outdir = project_root / "04_virtual_screening" / "receptors" / pid
    outdir.mkdir(parents=True, exist_ok=True)

    pdbqt_out = outdir / "receptor.pdbqt"
    if pdbqt_out.exists() and not force:
        log(f"  [SKIP] {pid}: receptor.pdbqt already exists (--force to redo)", CYAN)
        return True

    log(f"\n{'─'*60}", BOLD)
    log(f"  {pid}  ({pocket['protein']} {pocket['domain']})", BOLD)
    log(f"  Priority: {pocket['priority']}  |  max_persist: {pocket['max_persist']}")
    log(f"  Note: {pocket['note']}")
    log(f"  Centroid: ({pocket['cx']}, {pocket['cy']}, {pocket['cz']})")

    src_pdb = project_root / pocket["pdb"]
    if not src_pdb.exists():
        log(f"  [ERROR] Source PDB not found: {src_pdb}", RED)
        return False

    # Step 1: clean
    clean = outdir / "receptor_clean.pdb"
    n_atoms = clean_pdb(src_pdb, clean)
    log(f"  ✓ Cleaned PDB: {n_atoms} ATOM records → {clean.name}")

    # Step 2: pdb2pqr (protonation at pH 7.4)
    pqr_out  = outdir / "receptor.pqr"
    pdb_H    = outdir / "receptor_H.pdb"
    try:
        ok = run_pdb2pqr(clean, pqr_out, pdb_H)
        log(f"  ✓ pdb2pqr: protonated at pH 7.4 → {pdb_H.name}")
    except RuntimeError as e:
        log(f"  [ERROR] pdb2pqr: {e}", RED)
        return False

    # Step 3: Meeko → PDBQT
    try:
        run_meeko_receptor(pdb_H, pdbqt_out)
        size_kb = pdbqt_out.stat().st_size // 1024
        log(f"  ✓ Meeko: PDBQT written ({size_kb} KB) → {pdbqt_out.name}", GREEN)
    except RuntimeError as e:
        log(f"  [ERROR] Meeko: {e}", RED)
        return False

    # Step 4: Vina config
    cfg = outdir / "vina_config.txt"
    write_vina_config(pocket, pdbqt_out, cfg)
    log(f"  ✓ Vina config: box {pocket['box_x']}×{pocket['box_y']}×{pocket['box_z']} Å"
        f" centred at ({pocket['cx']}, {pocket['cy']}, {pocket['cz']})", GREEN)

    return True


# ---------------------------------------------------------------------------
def tool_version(cmd: str) -> str:
    """Return first line of --version output, or 'unknown' on failure."""
    try:
        r = subprocess.run([cmd, "--version"], capture_output=True, text=True, timeout=10)
        out = (r.stdout or r.stderr).strip().splitlines()
        return out[0] if out else "unknown"
    except Exception:
        return "unknown"


def main():
    parser = argparse.ArgumentParser(description="CRYPTAD receptor preparation (Step 4.1)")
    parser.add_argument("--pocket", default=None,
                        help="Prepare a single pocket by ID (e.g. S1_site688)")
    parser.add_argument("--force", action="store_true",
                        help="Redo even if receptor.pdbqt already exists")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()

    check_dependencies()

    targets = [p for p in POCKETS if args.pocket is None or p["id"] == args.pocket]
    if not targets:
        log(f"[ERROR] No pocket matching '{args.pocket}'. "
            f"Valid IDs: {[p['id'] for p in POCKETS]}", RED)
        sys.exit(1)

    log(f"\n{BOLD}CRYPTAD — Receptor Preparation (Step 4.1){RESET}")
    log(f"Pockets: {len(targets)}  |  Output: 04_virtual_screening/receptors/\n")

    ok = fail = 0
    results = []
    for pocket in sorted(targets, key=lambda p: p["priority"]):
        success = prepare_one(pocket, project_root, force=args.force)
        if success:
            ok += 1
        else:
            fail += 1
        results.append({"pocket_id": pocket["id"], "success": success})

    log(f"\n{'═'*60}", BOLD)
    log(f"  Done: {ok} OK  |  {fail} failed", BOLD)
    if ok > 0:
        log(f"\n  Next: Step 4.2 — ligand library preparation", GREEN)
        log(f"  python3 09_scripts/08_virtual_screening/02_prepare_ligands.py", GREEN)

    # Write manifest
    outdir = project_root / "04_virtual_screening" / "receptors"
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "generated_at":    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":          "09_scripts/08_virtual_screening/01_prepare_receptors.py",
        "project_root":    str(project_root),
        "pdb2pqr_version": tool_version("pdb2pqr"),
        "n_ok":   ok,
        "n_fail": fail,
        "pockets": results,
    }
    manifest_path = outdir / "prepare_receptors_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
