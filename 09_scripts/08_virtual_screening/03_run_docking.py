"""
run_docking.py  —  CRYPTAD
AutoDock Vina docking campaign (Step 4.3).

Docks CNS-penetrant approved ligands against cryptic-pocket receptors
using the Vina configs prepared in Step 4.1.

Local mode  (--local):
  Runs docking sequentially on this machine.  Use for testing only.

TRUBA mode  (default):
  Writes a SLURM array job script to 09_scripts/08_virtual_screening/run_docking_truba.sh
  and transfers all required files to TRUBA via rsync.

  Each SLURM task docks one batch of ligands against all receptors.
  Array size = ceil(n_ligands / BATCH_SIZE).

Output layout (local):
  04_virtual_screening/docking_results/{pocket_id}/{ligand_id}_out.pdbqt
  04_virtual_screening/docking_results/{pocket_id}/{ligand_id}_vina.txt

Output layout (TRUBA):
  /arf/scratch/mozkurt/CRYPTAD/04_virtual_screening/docking_results/
  (same tree — rsync back after completion)

Usage:
  # Write SLURM script + rsync to TRUBA (default):
  python3 09_scripts/08_virtual_screening/03_run_docking.py

  # Test locally on a single ligand:
  python3 09_scripts/08_virtual_screening/03_run_docking.py --local --ligand CHEMBL3

  # Test locally on first N ligands:
  python3 09_scripts/08_virtual_screening/03_run_docking.py --local --n 10

  # Custom project root:
  python3 09_scripts/08_virtual_screening/03_run_docking.py --project-root /path/to/CRYPTAD
"""

import argparse
import csv
import json
import math
import subprocess
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

TRUBA_USER    = "mozkurt"
TRUBA_HOST    = "hamsi.truba.gov.tr"
TRUBA_PROJECT = "/arf/scratch/mozkurt/CRYPTAD"

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; RED   = "\033[91m"; CYAN   = "\033[96m"
def log(m, c=RESET): print(f"{c}{m}{RESET}", flush=True)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_manifest(project_root: Path, limit=None, single=None):
    """Return list of (ligand_id, pdbqt_path) for PASS compounds."""
    manifest = project_root / "04_virtual_screening" / "ligands" / "library_manifest.csv"
    rows = []
    with open(manifest) as f:
        for row in csv.DictReader(f):
            if row["status"] != "PASS":
                continue
            if single and row["id"] != single:
                continue
            pdbqt = project_root / row["pdbqt"]
            if pdbqt.exists():
                rows.append((row["id"], pdbqt))
            if limit and len(rows) >= limit:
                break
    return rows


def get_receptors(recdir: Path):
    """Return list of (pocket_id, receptor_pdbqt, vina_config) tuples."""
    receptors = []
    for pocket_dir in sorted(recdir.iterdir()):
        pdbqt  = pocket_dir / "receptor.pdbqt"
        config = pocket_dir / "vina_config.txt"
        if pdbqt.exists() and config.exists():
            receptors.append((pocket_dir.name, pdbqt, config))
    return receptors


def parse_vina_config(config_path: Path) -> dict:
    """Parse a Vina config file into a {key: value} dict."""
    cfg = {}
    for line in config_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            k, _, v = line.partition("=")
            cfg[k.strip()] = v.strip()
    return cfg


# ---------------------------------------------------------------------------
# Local docking (single ligand against all receptors)
# ---------------------------------------------------------------------------
def dock_one(lig_id: str, lig_pdbqt: Path, receptors: list, dockdir: Path):
    """Dock lig_pdbqt against all receptors. Returns list of (pocket_id, score).

    Passes all Vina params as direct CLI flags (not --config) to avoid
    failures when paths contain spaces.
    """
    results = []
    for pocket_id, rec_pdbqt, config in receptors:
        out_dir  = dockdir / pocket_id
        out_dir.mkdir(parents=True, exist_ok=True)
        out_pdbqt = out_dir / f"{lig_id}_out.pdbqt"
        out_log   = out_dir / f"{lig_id}_vina.txt"

        cfg = parse_vina_config(config)

        cmd = [
            "vina",
            "--receptor",  str(rec_pdbqt),
            "--ligand",    str(lig_pdbqt),
            "--out",       str(out_pdbqt),
            "--log",       str(out_log),
            "--center_x",  cfg.get("center_x", "0"),
            "--center_y",  cfg.get("center_y", "0"),
            "--center_z",  cfg.get("center_z", "0"),
            "--size_x",    cfg.get("size_x", "20"),
            "--size_y",    cfg.get("size_y", "20"),
            "--size_z",    cfg.get("size_z", "20"),
            "--exhaustiveness", cfg.get("exhaustiveness", "16"),
            "--num_modes",      cfg.get("num_modes", "9"),
            "--energy_range",   cfg.get("energy_range", "3"),
            "--cpu",            cfg.get("cpu", "4"),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0 or not out_pdbqt.exists():
            log(f"  [WARN] {lig_id} @ {pocket_id}: {result.stderr[-200:]}", YELLOW)
            results.append((pocket_id, None))
            continue

        # Parse best score from log
        score = None
        for line in out_log.read_text().splitlines():
            if line.strip().startswith("1 "):
                try:
                    score = float(line.split()[1])
                except (IndexError, ValueError):
                    pass
                break
        results.append((pocket_id, score))

    return results


def run_local(ligands, receptors, project_root: Path):
    """Run docking locally (sequential). For testing only."""
    dockdir = project_root / "04_virtual_screening" / "docking_results"
    dockdir.mkdir(parents=True, exist_ok=True)

    log(f"\n{BOLD}Local docking — {len(ligands)} ligands × {len(receptors)} receptors{RESET}")
    log(f"Output: {dockdir.relative_to(project_root)}\n")

    results_csv = dockdir / "docking_scores.csv"
    with open(results_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand_id", "pocket_id", "vina_score_kcal_mol"])
        for i, (lig_id, lig_pdbqt) in enumerate(ligands, 1):
            log(f"  [{i}/{len(ligands)}] {lig_id}", CYAN)
            scores = dock_one(lig_id, lig_pdbqt, receptors, dockdir)
            for pocket_id, score in scores:
                writer.writerow([lig_id, pocket_id, score if score is not None else "NA"])

    log(f"\n  Results → {results_csv.relative_to(project_root)}", GREEN)


# ---------------------------------------------------------------------------
# TRUBA SLURM script generator
# ---------------------------------------------------------------------------
# hamsi: 56 cores/node (preferred — fairshare not depressed).
# Strategy: 56 cores per task, batch BATCH_SIZE ligands per task,
# N_PARALLEL parallel Vina processes × VINA_CPU CPUs each = 56 cores fully utilised.
BATCH_SIZE  = 50   # ligands per SLURM task
VINA_CPU    = 4    # cores per Vina call
N_PARALLEL  = 56 // VINA_CPU  # = 14 parallel Vina calls

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=cryptad_vina
#SBATCH --account=mozkurt
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-{n_tasks}

# ── Environment ──────────────────────────────────────────────────────────────
export PATH="$HOME/apps/vina:$PATH"

# ── Paths ────────────────────────────────────────────────────────────────────
PROJ="{truba_project}"
LIGDIR="$PROJ/04_virtual_screening/ligands/pdbqt"
RECDIR="$PROJ/04_virtual_screening/receptors"
DOCKDIR="$PROJ/04_virtual_screening/docking_results"
MANIFEST="$PROJ/04_virtual_screening/ligands/ligand_ids.txt"

BATCH_SIZE={batch_size}
N_PARALLEL={n_parallel}
VINA_CPU={vina_cpu}

mkdir -p "$DOCKDIR/logs"

# ── Select ligand batch for this task ────────────────────────────────────────
TOTAL=$(wc -l < "$MANIFEST")
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
[[ $END -gt $TOTAL ]] && END=$TOTAL

LIGANDS=$(awk "NR>=${{START}} && NR<=${{END}}" "$MANIFEST")

# ── Per-ligand docking function (called by GNU parallel) ─────────────────────
dock_ligand() {{
    LIG_ID="$1"
    LIG_PDBQT="$LIGDIR/${{LIG_ID}}.pdbqt"
    [[ ! -f "$LIG_PDBQT" ]] && {{ echo "SKIP $LIG_ID (pdbqt missing)"; return; }}

    for POCKET_DIR in "$RECDIR"/*/; do
        POCKET_ID=$(basename "$POCKET_DIR")
        CONFIG="$POCKET_DIR/vina_config.txt"
        RECEPTOR="$POCKET_DIR/receptor.pdbqt"
        [[ ! -f "$CONFIG"   ]] && continue
        [[ ! -f "$RECEPTOR" ]] && continue

        OUTDIR="$DOCKDIR/$POCKET_ID"
        mkdir -p "$OUTDIR"
        OUT_PDBQT="$OUTDIR/${{LIG_ID}}_out.pdbqt"
        OUT_LOG="$OUTDIR/${{LIG_ID}}_vina.txt"

        [[ -f "$OUT_PDBQT" ]] && continue   # already done — safe resume

        # Extract box params from config to avoid --config with space-in-path issues
        CX=$(awk -F'[[:space:]]*=[[:space:]]*' '/^center_x/{{print $2}}' "$CONFIG")
        CY=$(awk -F'[[:space:]]*=[[:space:]]*' '/^center_y/{{print $2}}' "$CONFIG")
        CZ=$(awk -F'[[:space:]]*=[[:space:]]*' '/^center_z/{{print $2}}' "$CONFIG")
        SX=$(awk -F'[[:space:]]*=[[:space:]]*' '/^size_x/{{print $2}}' "$CONFIG")
        SY=$(awk -F'[[:space:]]*=[[:space:]]*' '/^size_y/{{print $2}}' "$CONFIG")
        SZ=$(awk -F'[[:space:]]*=[[:space:]]*' '/^size_z/{{print $2}}' "$CONFIG")
        EX=$(awk -F'[[:space:]]*=[[:space:]]*' '/^exhaustiveness/{{print $2}}' "$CONFIG")

        vina \\
            --receptor "$RECEPTOR" \\
            --ligand   "$LIG_PDBQT" \\
            --out      "$OUT_PDBQT" \\
            --log      "$OUT_LOG" \\
            --center_x "$CX" --center_y "$CY" --center_z "$CZ" \\
            --size_x   "$SX" --size_y   "$SY" --size_z   "$SZ" \\
            --exhaustiveness "$EX" \\
            --cpu      $VINA_CPU
    done

    SCORE=$(awk '/^REMARK VINA RESULT/{{print $4; exit}}' "$DOCKDIR"/S1_site688/"${{LIG_ID}}_out.pdbqt" 2>/dev/null || echo "NA")
    echo "DONE $LIG_ID score=$SCORE" \\
        >> "$DOCKDIR/logs/done_${{SLURM_ARRAY_JOB_ID}}_${{SLURM_ARRAY_TASK_ID}}.txt"
}}
export -f dock_ligand
export LIGDIR RECDIR DOCKDIR VINA_CPU

# ── Run in parallel (14 ligands at a time × 4 CPUs = 56 cores) ───────────────
echo "$LIGANDS" | /usr/bin/parallel --jobs $N_PARALLEL dock_ligand {{}}
"""


def write_slurm_script(n_ligands: int, script_dir: Path):
    n_tasks = math.ceil(n_ligands / BATCH_SIZE)
    script_path = script_dir / "run_docking_truba.sh"
    content = SLURM_TEMPLATE.format(
        truba_project=TRUBA_PROJECT,
        n_tasks=n_tasks,
        batch_size=BATCH_SIZE,
        n_parallel=N_PARALLEL,
        vina_cpu=VINA_CPU,
    )
    script_path.write_text(content)
    return script_path, n_tasks


# ---------------------------------------------------------------------------
# Ligand ID list for TRUBA (plain text, one ID per line)
# ---------------------------------------------------------------------------
def write_ligand_id_list(project_root: Path):
    id_file = project_root / "04_virtual_screening" / "ligands" / "ligand_ids.txt"
    ligands = load_manifest(project_root)
    id_file.write_text("\n".join(lid for lid, _ in ligands) + "\n")
    return id_file, len(ligands)


# ---------------------------------------------------------------------------
# rsync helper
# ---------------------------------------------------------------------------
def rsync_to_truba(paths: list, dest_prefix: str, project_root: Path):
    """rsync a list of local paths to TRUBA."""
    for local_path in paths:
        dest = f"{TRUBA_USER}@{TRUBA_HOST}:{dest_prefix}/{local_path.relative_to(project_root)}"
        log(f"  rsync {local_path.relative_to(project_root)} → TRUBA …", CYAN)
        result = subprocess.run(
            ["rsync", "-avz", "--mkpath", str(local_path), dest],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            log(f"  [rsync ERROR]\n{result.stderr[-400:]}", RED)
        else:
            log(f"  ✓ transferred", GREEN)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="CRYPTAD docking campaign (Step 4.3)")
    parser.add_argument("--local", action="store_true",
                        help="Run docking locally (testing mode)")
    parser.add_argument("--ligand", default=None,
                        help="Single ligand ChEMBL ID to dock (local mode)")
    parser.add_argument("--n", type=int, default=None,
                        help="Limit to first N ligands (local mode)")
    parser.add_argument("--no-rsync", action="store_true",
                        help="Write SLURM script but skip rsync to TRUBA")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    recdir  = project_root / "04_virtual_screening" / "receptors"

    receptors = get_receptors(recdir)
    if not receptors:
        log(f"[ERROR] No receptor PDBQT files found in "
            f"{recdir.relative_to(project_root)}. Run Step 4.1 first.", RED)
        sys.exit(1)

    log(f"\n{BOLD}CRYPTAD — Docking Campaign (Step 4.3){RESET}")
    log(f"Receptors : {len(receptors)} pockets — "
        + ", ".join(r[0] for r in receptors))

    # ── Local mode ────────────────────────────────────────────────────────────
    if args.local:
        if not __import__("shutil").which("vina"):
            log("[ERROR] 'vina' not found in PATH. "
                "Install with: conda install -c conda-forge vina", RED)
            sys.exit(1)
        ligands = load_manifest(project_root, limit=args.n, single=args.ligand)
        manifest_path = project_root / "04_virtual_screening" / "ligands" / "library_manifest.csv"
        if not ligands:
            log(f"[ERROR] No ligands found (manifest: {manifest_path})", RED)
            sys.exit(1)
        log(f"Ligands   : {len(ligands)}\n")
        run_local(ligands, receptors, project_root)
        # Write manifest
        manifest = {
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "script":       "09_scripts/08_virtual_screening/03_run_docking.py",
            "project_root": str(project_root),
            "mode":         "local",
            "n_ligands":    len(ligands),
            "n_receptors":  len(receptors),
            "receptors":    [r[0] for r in receptors],
        }
        dockdir = project_root / "04_virtual_screening" / "docking_results"
        dockdir.mkdir(parents=True, exist_ok=True)
        with open(dockdir / "run_docking_manifest.json", "w") as fh:
            json.dump(manifest, fh, indent=2)
        return

    # ── TRUBA mode ────────────────────────────────────────────────────────────
    id_file, n_ligands = write_ligand_id_list(project_root)
    log(f"Ligands   : {n_ligands}  |  batch size: {BATCH_SIZE}  |  "
        f"{N_PARALLEL} parallel × {VINA_CPU} CPUs = 56 cores/node\n")

    slurm_script, n_tasks = write_slurm_script(n_ligands, _script.parent)
    log(f"  ✓ SLURM script: {slurm_script.relative_to(project_root)}", GREEN)
    log(f"  ✓ Ligand ID list: {id_file.relative_to(project_root)}", GREEN)

    total_runs = n_ligands * len(receptors)
    log(f"\n  Array: {n_tasks} tasks × {BATCH_SIZE} ligands × {len(receptors)} receptors "
        f"= {total_runs:,} docking runs", BOLD)
    log(f"  Estimated wall time: ~3–4 h/task (56 cores, exhaustiveness=16)", CYAN)

    if args.no_rsync:
        log("\n  [--no-rsync] Skipping transfer to TRUBA.", YELLOW)
    else:
        log("\n  Transferring files to TRUBA …", BOLD)
        rsync_to_truba([
            project_root / "04_virtual_screening" / "receptors",
            project_root / "04_virtual_screening" / "ligands",
            slurm_script,
        ], TRUBA_PROJECT, project_root)

    log(f"\n{'═'*60}", BOLD)
    log(f"  Next steps on TRUBA:", BOLD)
    log(f"  1. ssh {TRUBA_USER}@{TRUBA_HOST}", GREEN)
    log(f"  2. cd {TRUBA_PROJECT}", GREEN)
    log(f"  3. sbatch 09_scripts/08_virtual_screening/run_docking_truba.sh", GREEN)
    log(f"  4. squeue -u {TRUBA_USER}   # monitor", GREEN)
    log(f"\n  After completion:", GREEN)
    log(f"  5. python3 09_scripts/08_virtual_screening/04_parse_docking_results.py", GREEN)

    # Write manifest
    manifest = {
        "generated_at":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":        "09_scripts/08_virtual_screening/03_run_docking.py",
        "project_root":  str(project_root),
        "mode":          "truba",
        "n_ligands":     n_ligands,
        "n_receptors":   len(receptors),
        "receptors":     [r[0] for r in receptors],
        "n_tasks":       n_tasks,
        "batch_size":    BATCH_SIZE,
        "vina_cpu":      VINA_CPU,
        "slurm_script":  str(slurm_script),
        "truba_project": TRUBA_PROJECT,
    }
    outdir = project_root / "04_virtual_screening"
    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "run_docking_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"Manifest → {outdir / 'run_docking_manifest.json'}")


if __name__ == "__main__":
    main()
