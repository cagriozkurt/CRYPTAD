"""
setup_mmgbsa.py  —  CRYPTAD
MM-GBSA rescoring setup (Step 4.4).

Workflow:
  LOCAL (this script):
    1. Read top50 CSVs from docking results (all 4 pockets)
    2. Deduplicate: build unique (pocket, ligand) pair list
    3. Extract best docked pose from each PDBQT → SDF (meeko)
    4. Write per-complex directory structure
    5. Write SLURM array script for TRUBA
    6. rsync to TRUBA

  TRUBA (SLURM job — one task per complex):
    1. acpype: GAFF2 topology + AM1-BCC charges for ligand
    2. pdb2gmx: CHARMM36m topology for protein
    3. Build combined topology + merge GRO files
    4. Solvate + add ions
    5. Energy minimise → NVT 100 ps → NPT 100 ps → 3 ns production
    6. gmx_MMPBSA: compute ΔG_bind with GB model

Output layout:
  04_virtual_screening/mmgbsa/
    pairs.txt                   one "pocket_id ligand_id" per line
    {pocket_id}/{ligand_id}/
      ligand_pose.sdf           best docked pose (from Vina PDBQT)
      ligand_pose.pdb           same, PDB format
      [TRUBA adds:]
      ligand_GMX.itp / .gro     GAFF2 topology (acpype)
      complex.tpr               solvated complex topology
      md.xtc                    3 ns production trajectory
      mmgbsa/FINAL_RESULTS.csv  per-frame ΔG_bind

Usage:
  python3 09_scripts/08_virtual_screening/05_setup_mmgbsa.py             # prepare + rsync
  python3 09_scripts/08_virtual_screening/05_setup_mmgbsa.py --no-rsync  # local prep only
  python3 09_scripts/08_virtual_screening/05_setup_mmgbsa.py --pocket S1_site688
  python3 09_scripts/08_virtual_screening/05_setup_mmgbsa.py --project-root /path/to/CRYPTAD
"""

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate

# ---------------------------------------------------------------------------
_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

TRUBA_USER    = "mozkurt"
TRUBA_HOST    = "hamsi.truba.gov.tr"
TRUBA_PROJECT = "/arf/scratch/mozkurt/CRYPTAD"

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; RED   = "\033[91m"; CYAN   = "\033[96m"
def log(m, c=RESET): print(f"{c}{m}{RESET}", flush=True)

TOP_N = 50


# ---------------------------------------------------------------------------
# Load top-N hits per pocket from docking CSVs
# ---------------------------------------------------------------------------
def load_top_pairs(pockets: list, top_n: int, dockdir: Path) -> list:
    """Return sorted unique (pocket, ligand_id) pairs from top-N CSVs."""
    pairs = []
    seen  = set()
    for pocket in pockets:
        csv_path = dockdir / f"top{top_n}_{pocket}.csv"
        if not csv_path.exists():
            log(f"  [WARN] {csv_path.name} not found — run 04_parse_docking_results.py first",
                YELLOW)
            continue
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                lig_id = row["id"] if "id" in row else list(row.values())[0]
                key = (pocket, lig_id)
                if key not in seen:
                    pairs.append(key)
                    seen.add(key)
        log(f"  {pocket}: loaded top {top_n} → {sum(1 for p, _ in pairs if p == pocket)} pairs")
    return pairs


# ---------------------------------------------------------------------------
# Extract best docked pose from PDBQT → SDF + PDB
# ---------------------------------------------------------------------------
def extract_pose(pdbqt_path: Path, out_sdf: Path, out_pdb: Path) -> bool:
    """
    Extract best (rank 1) docked pose from Vina PDBQT output.
    Returns True on success.
    """
    try:
        pdbqt_mol = PDBQTMolecule.from_file(str(pdbqt_path), poses_to_read=1)
        sdf_str, failures = RDKitMolCreate.write_sd_string(pdbqt_mol)
        if failures or not sdf_str.strip():
            return False

        out_sdf.write_text(sdf_str)

        # SDF → PDB via RDKit
        suppl = Chem.SDMolSupplier(str(out_sdf), removeHs=False)
        mol   = next((m for m in suppl if m is not None), None)
        if mol is None:
            return False
        Chem.MolToPDBFile(mol, str(out_pdb))
        return True
    except Exception as e:
        log(f"    [pose extract error] {e}", YELLOW)
        return False


# ---------------------------------------------------------------------------
# Build per-complex directory with pose files
# ---------------------------------------------------------------------------
def prepare_local(pairs: list, project_root: Path) -> list:
    """
    For each (pocket, ligand) pair:
      - Create mmgbsa/{pocket}/{ligand}/
      - Extract pose SDF + PDB
      - Copy receptor_H.pdb
    Returns list of successfully prepared pairs.
    """
    mmgbsadir = project_root / "04_virtual_screening" / "mmgbsa"
    recdir    = project_root / "04_virtual_screening" / "receptors"
    dockdir   = project_root / "04_virtual_screening" / "docking_results"

    ok_pairs = []
    n_fail   = 0

    for pocket, lig_id in pairs:
        out_dir = mmgbsadir / pocket / lig_id
        out_dir.mkdir(parents=True, exist_ok=True)

        # Receptor
        rec_src = recdir / pocket / "receptor_H.pdb"
        rec_dst = out_dir / "receptor_H.pdb"
        if rec_src.exists() and not rec_dst.exists():
            shutil.copy2(rec_src, rec_dst)

        # Ligand pose
        pdbqt_src = dockdir / pocket / f"{lig_id}_out.pdbqt"
        out_sdf   = out_dir / "ligand_pose.sdf"
        out_pdb   = out_dir / "ligand_pose.pdb"

        if out_sdf.exists() and out_pdb.exists():
            ok_pairs.append((pocket, lig_id))
            continue

        if not pdbqt_src.exists():
            log(f"  [SKIP] {pocket}/{lig_id}: PDBQT not found", YELLOW)
            n_fail += 1
            continue

        if extract_pose(pdbqt_src, out_sdf, out_pdb):
            ok_pairs.append((pocket, lig_id))
        else:
            log(f"  [FAIL] pose extraction: {pocket}/{lig_id}", YELLOW)
            n_fail += 1

    log(f"\n  Prepared: {len(ok_pairs)} pairs  |  failed: {n_fail}", GREEN)
    return ok_pairs


# ---------------------------------------------------------------------------
# Write pairs list for SLURM
# ---------------------------------------------------------------------------
def write_pairs_file(pairs: list, mmgbsadir: Path) -> Path:
    pairs_file = mmgbsadir / "pairs.txt"
    pairs_file.write_text("\n".join(f"{p} {l}" for p, l in pairs) + "\n")
    return pairs_file


# ---------------------------------------------------------------------------
# MDP files
# ---------------------------------------------------------------------------
MDP_MINIM = """\
; Energy minimisation
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 10000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
"""

MDP_NVT = """\
; NVT equilibration 100 ps
define      = -DPOSRES
integrator  = md
nsteps      = 25000
dt          = 0.004
nstxout     = 0
nstvout     = 0
nstenergy   = 500
nstlog      = 500
continuation = no
gen_vel     = yes
gen_temp    = 310
gen_seed    = -1
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 20
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype = PME
pme_order   = 4
ewald_rtol  = 1e-5
tcoupl      = V-rescale
tc-grps     = Protein_LIG Water_and_ions
tau_t       = 0.1  0.1
ref_t       = 310  310
pcoupl      = no
pbc         = xyz
DispCorr    = EnerPres
constraints = h-bonds
constraint_algorithm = LINCS
"""

MDP_NPT = """\
; NPT equilibration 100 ps
define      = -DPOSRES
integrator  = md
nsteps      = 25000
dt          = 0.004
nstxout     = 0
nstenergy   = 500
nstlog      = 500
continuation = yes
gen_vel     = no
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 20
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype = PME
pme_order   = 4
ewald_rtol  = 1e-5
tcoupl      = V-rescale
tc-grps     = Protein_LIG Water_and_ions
tau_t       = 0.1  0.1
ref_t       = 310  310
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
DispCorr    = EnerPres
constraints = h-bonds
constraint_algorithm = LINCS
"""

MDP_PROD = """\
; Production MD 3 ns
integrator  = md
nsteps      = 750000
dt          = 0.004
nstxout     = 0
nstvout     = 0
nstfout     = 0
nstenergy   = 5000
nstlog      = 5000
nstxout-compressed = 5000
compressed-x-grps  = System
continuation = yes
gen_vel     = no
cutoff-scheme = Verlet
ns_type     = grid
nstlist     = 20
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype = PME
pme_order   = 4
ewald_rtol  = 1e-5
tcoupl      = V-rescale
tc-grps     = Protein_LIG Water_and_ions
tau_t       = 0.1  0.1
ref_t       = 310  310
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
pbc         = xyz
DispCorr    = EnerPres
constraints = h-bonds
constraint_algorithm = LINCS
"""

MMGBSA_IN = """\
&general
  sys_name = "CRYPTAD_complex",
  startframe = 1,
  endframe = 999999,
  interval = 5,
  verbose = 1,
  use_sander = 0,
/
&gb
  igb = 5,
  saltcon = 0.150,
/
"""


def write_mdp_files(mdp_dir: Path):
    mdp_dir.mkdir(parents=True, exist_ok=True)
    (mdp_dir / "minim.mdp").write_text(MDP_MINIM)
    (mdp_dir / "nvt.mdp").write_text(MDP_NVT)
    (mdp_dir / "npt.mdp").write_text(MDP_NPT)
    (mdp_dir / "prod.mdp").write_text(MDP_PROD)
    (mdp_dir / "mmgbsa.in").write_text(MMGBSA_IN)


# ---------------------------------------------------------------------------
# SLURM script
# ---------------------------------------------------------------------------
SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH --job-name=cryptad_mmgbsa
#SBATCH --account=mozkurt
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=2
#SBATCH --time=08:00:00
#SBATCH --array=1-{n_pairs}
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -uo pipefail

# ── Environment ──────────────────────────────────────────────────────────────
module load apps/gromacs/2024.1-oneapi2024
source $HOME/.venv/bin/activate 2>/dev/null || true
# gmx_MMPBSA — install once: pip install gmx_MMPBSA
# acpype     — install once: pip install acpype

export OMP_NUM_THREADS=${{SLURM_CPUS_PER_TASK:-2}}
export GMX="gmx_mpi"
export MPI="mpirun -np ${{SLURM_NTASKS:-28}}"

# ── Paths ────────────────────────────────────────────────────────────────────
PROJ="${{CRYPTAD:-{truba_project}}}"
MMGBSADIR="$PROJ/04_virtual_screening/mmgbsa"
MDPDIR="$MMGBSADIR/mdp"
PAIRS="$MMGBSADIR/pairs.txt"

# ── Select pair for this task ─────────────────────────────────────────────────
LINE=$(awk "NR==${{SLURM_ARRAY_TASK_ID}}" "$PAIRS")
POCKET=$(echo "$LINE" | awk '{{print $1}}')
LIG_ID=$(echo "$LINE" | awk '{{print $2}}')
WORKDIR="$MMGBSADIR/$POCKET/$LIG_ID"
cd "$WORKDIR" || {{ echo "ERROR: cannot cd to $WORKDIR"; exit 1; }}

echo "=== Task $SLURM_ARRAY_TASK_ID: $POCKET / $LIG_ID ==="

# ── 1. GAFF2 ligand topology (acpype) ────────────────────────────────────────
if [[ ! -f "LIG_GMX.itp" ]]; then
    echo "--- acpype ---"
    acpype -i ligand_pose.sdf -c bcc -a gaff2 -n 0 -o gmx 2>&1
    # acpype creates a directory like ligand_pose.acpype/
    ACPYPE_DIRS=( *.acpype )
    ACPYPE_DIR="${{ACPYPE_DIRS[0]:-}}"
    if [[ -z "$ACPYPE_DIR" || ! -d "$ACPYPE_DIR" ]]; then
        echo "ERROR: acpype failed for $LIG_ID" >&2; exit 1
    fi
    cp "$ACPYPE_DIR/LIG_GMX.gro"  LIG_GMX.gro
    cp "$ACPYPE_DIR/LIG_GMX.itp"  LIG_GMX.itp
    cp "$ACPYPE_DIR/LIG_GMX.top"  LIG_GMX.top
    cp "$ACPYPE_DIR/posre_LIG.itp" posre_LIG.itp 2>/dev/null || true
fi

# ── 2. Protein topology (CHARMM36m) ──────────────────────────────────────────
if [[ ! -f "topol.top" ]]; then
    echo "--- pdb2gmx ---"
    $GMX pdb2gmx \
        -f receptor_H.pdb \
        -o protein.gro \
        -p topol.top \
        -i posre.itp \
        -ff charmm36m \
        -water tip3p \
        -ignh \
        -nobackup

    # Inject ligand ITP into topol.top (before first [ molecules ] section)
    sed -i '/; Include water topology/i #include "LIG_GMX.itp"\n#ifdef POSRES\n#include "posre_LIG.itp"\n#endif' topol.top
    # Add LIG to [ molecules ] section
    echo "LIG  1" >> topol.top
fi

# ── 3. Build complex GRO ─────────────────────────────────────────────────────
if [[ ! -f "complex.gro" ]]; then
    echo "--- building complex GRO ---"
    N_PROT=$(awk 'NR==2{{print; exit}}' protein.gro)
    N_LIG=$(awk 'NR==2{{print; exit}}' LIG_GMX.gro)
    N_TOTAL=$(( N_PROT + N_LIG ))
    BOX=$(tail -1 protein.gro)
    {{
        echo "Complex"
        echo "  $N_TOTAL"
        head -n $(( N_PROT + 2 )) protein.gro | tail -n $N_PROT
        head -n $(( N_LIG + 2 )) LIG_GMX.gro  | tail -n $N_LIG
        echo "$BOX"
    }} > complex.gro
fi

# ── 4. Solvate + ions ────────────────────────────────────────────────────────
if [[ ! -f "solvated.gro" ]]; then
    echo "--- solvate ---"
    $GMX editconf -f complex.gro -o complex_box.gro -bt dodecahedron -d 1.2 -nobackup
    $GMX solvate  -cp complex_box.gro -cs spc216.gro -o solvated.gro -p topol.top -nobackup
    $GMX grompp   -f "$MDPDIR/minim.mdp" -c solvated.gro -p topol.top \
                  -o ions.tpr -maxwarn 1 -nobackup
    echo "SOL" | $GMX genion -s ions.tpr -o solvated_ions.gro -p topol.top \
                  -pname NA -nname CL -neutral -conc 0.15 -nobackup
    mv solvated_ions.gro solvated.gro
fi

# ── 5. Make index ─────────────────────────────────────────────────────────────
if [[ ! -f "index.ndx" ]]; then
    echo "--- make_ndx ---"
    printf 'q\\n' | $GMX make_ndx -f solvated.gro -o index.ndx -nobackup
    printf 'name 20 Protein_LIG\\n1 | 13\\nq\\n' | \
        $GMX make_ndx -f solvated.gro -n index.ndx -o index.ndx -nobackup || true
fi

# ── 6. Energy minimisation ───────────────────────────────────────────────────
if [[ ! -f "minim.gro" ]]; then
    echo "--- energy minimisation ---"
    $GMX grompp -f "$MDPDIR/minim.mdp" -c solvated.gro -p topol.top \
                -o minim.tpr -nobackup
    $MPI $GMX mdrun -v -deffnm minim -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
fi

# ── 7. NVT equilibration ─────────────────────────────────────────────────────
if [[ ! -f "nvt.gro" ]]; then
    echo "--- NVT ---"
    $GMX grompp -f "$MDPDIR/nvt.mdp" -c minim.gro -r minim.gro \
                -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1 -nobackup
    $MPI $GMX mdrun -deffnm nvt -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
fi

# ── 8. NPT equilibration ─────────────────────────────────────────────────────
if [[ ! -f "npt.gro" ]]; then
    echo "--- NPT ---"
    $GMX grompp -f "$MDPDIR/npt.mdp" -c nvt.gro -r nvt.gro \
                -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 1 -nobackup
    $MPI $GMX mdrun -deffnm npt -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
fi

# ── 9. Production MD 3 ns ────────────────────────────────────────────────────
if [[ ! -f "prod.gro" ]]; then
    echo "--- production MD ---"
    $GMX grompp -f "$MDPDIR/prod.mdp" -c npt.gro -t npt.cpt \
                -p topol.top -n index.ndx -o prod.tpr -nobackup
    $MPI $GMX mdrun -deffnm prod -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
fi

# ── 10. gmx_MMPBSA ───────────────────────────────────────────────────────────
if [[ ! -f "mmgbsa/FINAL_RESULTS_MMPBSA.csv" ]]; then
    echo "--- gmx_MMPBSA ---"
    mkdir -p mmgbsa
    cd mmgbsa
    gmx_MMPBSA \
        -O \
        -i "$MDPDIR/mmgbsa.in" \
        -cs ../prod.tpr \
        -ct ../prod.xtc \
        -ci ../index.ndx \
        -cg 1 13 \
        -cp ../topol.top \
        -nogui \
        -o FINAL_RESULTS_MMPBSA.dat \
        -eo FINAL_RESULTS_MMPBSA.csv
    cd ..
fi

echo "=== DONE: $POCKET / $LIG_ID ==="
"""


def write_slurm_script(n_pairs: int, script_dir: Path) -> Path:
    script_path = script_dir / "run_mmgbsa_truba.sh"
    script_path.write_text(SLURM_TEMPLATE.format(
        truba_project=TRUBA_PROJECT,
        n_pairs=n_pairs,
    ))
    return script_path


# ---------------------------------------------------------------------------
# rsync
# ---------------------------------------------------------------------------
def rsync_to_truba(paths: list, project_root: Path):
    for local in paths:
        dest = f"{TRUBA_USER}@{TRUBA_HOST}:{TRUBA_PROJECT}/{local.relative_to(project_root)}"
        log(f"  rsync {local.relative_to(project_root)} → TRUBA …", CYAN)
        r = subprocess.run(
            ["rsync", "-avz", "--mkpath", str(local), dest],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            log(f"  [rsync ERROR]\n{r.stderr[-300:]}", RED)
        else:
            log(f"  ✓", GREEN)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="CRYPTAD MM-GBSA setup (Step 4.4)")
    parser.add_argument("--pocket", default=None,
                        help="Limit to one pocket (e.g. S1_site688)")
    parser.add_argument("--top-n", type=int, default=TOP_N,
                        help=f"Top N compounds per pocket (default: {TOP_N})")
    parser.add_argument("--no-rsync", action="store_true",
                        help="Prepare locally but skip rsync to TRUBA")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    dockdir   = project_root / "04_virtual_screening" / "docking_results"
    recdir    = project_root / "04_virtual_screening" / "receptors"
    mmgbsadir = project_root / "04_virtual_screening" / "mmgbsa"
    mmgbsadir.mkdir(parents=True, exist_ok=True)

    log(f"\n{BOLD}CRYPTAD — MM-GBSA Setup (Step 4.4){RESET}")

    pockets = sorted(p.name for p in recdir.iterdir() if p.is_dir())
    if args.pocket:
        pockets = [p for p in pockets if p == args.pocket]
        if not pockets:
            log(f"[ERROR] pocket '{args.pocket}' not found", RED); sys.exit(1)

    log(f"Pockets : {pockets}")
    log(f"Top N   : {args.top_n} per pocket\n")

    # 1. Collect pairs
    log("Loading top-hit pairs …")
    pairs = load_top_pairs(pockets, args.top_n, dockdir)
    log(f"\n  Total (pocket, ligand) pairs: {len(pairs)}", BOLD)

    # 2. Prepare local files
    log("\nExtracting docked poses …")
    pairs = prepare_local(pairs, project_root)

    # 3. Write MDP files
    mdp_dir = mmgbsadir / "mdp"
    write_mdp_files(mdp_dir)
    log(f"  ✓ MDP files → {mdp_dir.relative_to(project_root)}", GREEN)

    # 4. Write pairs file
    pairs_file = write_pairs_file(pairs, mmgbsadir)
    log(f"  ✓ Pairs list: {len(pairs)} pairs → {pairs_file.relative_to(project_root)}", GREEN)

    # 5. Write SLURM script
    slurm_script = write_slurm_script(len(pairs), _script.parent)
    log(f"  ✓ SLURM script → {slurm_script.relative_to(project_root)}", GREEN)

    log(f"\n  Array: {len(pairs)} tasks  |  ~6–8 h/task (56 cores)", BOLD)
    log(f"  Est. total compute: ~{len(pairs) * 56 * 7 / 60:.0f} CPU-hours", CYAN)

    # 6. rsync
    if args.no_rsync:
        log("\n  [--no-rsync] Skipping transfer to TRUBA.", YELLOW)
    else:
        log("\nTransferring to TRUBA …")
        rsync_to_truba([mmgbsadir, slurm_script], project_root)

    log(f"\n{'═'*60}", BOLD)
    log(f"  On TRUBA — install dependencies (once):", BOLD)
    log(f"  pip install gmx_MMPBSA acpype", GREEN)
    log(f"\n  Then submit:", BOLD)
    log(f"  cd {TRUBA_PROJECT}", GREEN)
    log(f"  sbatch 09_scripts/08_virtual_screening/run_mmgbsa_truba.sh", GREEN)
    log(f"  squeue -u {TRUBA_USER}", GREEN)
    log(f"\n  After completion:", GREEN)
    log(f"  python3 09_scripts/08_virtual_screening/06_parse_mmgbsa_results.py", GREEN)

    # Write manifest
    manifest = {
        "generated_at":  time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":        "09_scripts/08_virtual_screening/05_setup_mmgbsa.py",
        "project_root":  str(project_root),
        "pockets":       pockets,
        "top_n":         args.top_n,
        "n_pairs":       len(pairs),
        "slurm_script":  str(slurm_script),
        "pairs_file":    str(pairs_file),
        "truba_project": TRUBA_PROJECT,
    }
    manifest_path = mmgbsadir / "setup_mmgbsa_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
