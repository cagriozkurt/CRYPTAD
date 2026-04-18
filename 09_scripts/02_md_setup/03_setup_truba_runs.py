"""
CRYPTAD — TRUBA Run Setup Script
Creates run directories, production MDPs, NPT equilibration MDPs,
and SLURM job scripts for all 5 MD systems.

Reproducibility notes
---------------------
- nsteps is computed using integer arithmetic (ps→fs then integer division)
  to avoid floating-point rounding from the 4 fs time step.
- SLURM scripts target the hamsi partition (56 cores/node, preferred over
  orfoz — see project MEMORY for fairshare rationale).
  hamsi is CPU-only; standard GROMACS 2024.1 (no GPU) is used here.
  Metadynamics runs use a separate PLUMED-patched build.
- maxwarn is set to 1 (the single known CHARMM-GUI CMAP note).
  Raise only if grompp explicitly reports a second known-harmless warning.
- A setup_manifest.json is written on completion for audit purposes.

Usage:
    python 09_scripts/02_md_setup/03_setup_truba_runs.py [--project-root PATH]

Options:
    --project-root    Explicit path to CRYPTAD project root.
                      Default: inferred as three directories above this script.
"""

import argparse
import json
import os
import stat
import sys
from datetime import datetime, timezone
from pathlib import Path

# ── Path resolution ────────────────────────────────────────────────────────────
# Script lives at <root>/09_scripts/02_md_setup/03_setup_truba_runs.py
# parents[0] = 02_md_setup/
# parents[1] = 09_scripts/
# parents[2] = project root
_SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_ROOT = _SCRIPT_PATH.parents[2]

# ── System definitions ────────────────────────────────────────────────────────
#   ns_prod  : target production simulation length in nanoseconds
#   nreps    : number of independent replicas
SYSTEMS = {
    "S1_BIN1_BAR":    {"ns_prod": 500, "label": "BIN1BAR",   "nreps": 3, "chain": 30},
    "S2_BIN1_SH3":    {"ns_prod": 500, "label": "BIN1SH3",   "nreps": 3, "chain":  4},
    "S3_PICALM_ANTH": {"ns_prod": 500, "label": "PICALMANT", "nreps": 3, "chain":  8},
    "S4_CD2AP_SH3-2": {"ns_prod": 200, "label": "CD2APSH32", "nreps": 3, "chain":  3},
    "S5_CD2AP_SH3-1": {"ns_prod": 200, "label": "CD2APSH31", "nreps": 3, "chain":  3},
}

# Time steps in femtoseconds (integers to avoid floating-point arithmetic downstream)
DT_PROD_FS  = 4    # 4 fs — HMR enabled
DT_EQUIL_FS = 1    # 1 fs — NVT equilibration from CHARMM-GUI
DT_NPT_FS   = 4    # 4 fs — NPT equilibration

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; CYAN  = "\033[96m"; RED = "\033[91m"
def log(m, c=RESET): print(f"{c}{m}{RESET}", flush=True)


# ── Helper: find the extracted charmm-gui-XXXXXXXXXX directory ────────────────

def get_charmm_dir(system_dir: Path) -> Path | None:
    for entry in sorted(system_dir.iterdir()):
        if entry.name.startswith("charmm-gui-") and entry.is_dir():
            return entry
    return None


# ── MDP generators ────────────────────────────────────────────────────────────

NPT_EQUIL_TEMPLATE = """\
; NPT equilibration — 1 ns, weak position restraints, 4 fs step (HMR)
; Restraints: BB 40 kJ/mol/nm², SC 10 kJ/mol/nm² (10x weaker than NVT)
define                  = -DPOSRES -DPOSRES_FC_BB=40.0 -DPOSRES_FC_SC=10.0
integrator              = md
dt                      = 0.004
nsteps                  = 250000
nstxout-compressed      = 25000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
rlist                   = 1.2
rcoulomb                = 1.2
coulombtype             = PME
;
tcoupl                  = v-rescale
tc_grps                 = SOLU SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU SOLV
"""

PROD_TEMPLATE = """\
; Production MD — {ns} ns  ({nsteps} steps at 4 fs, HMR)
; NPT ensemble, CHARMM36m, TIP3P, 310.15 K, 1 bar, KCl 0.15 M
integrator              = md
dt                      = 0.004
nsteps                  = {nsteps}
nstxout-compressed      = 25000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
rlist                   = 1.2
rcoulomb                = 1.2
coulombtype             = PME
;
tcoupl                  = v-rescale
tc_grps                 = SOLU SOLV
tau_t                   = 1.0 1.0
ref_t                   = 310.15 310.15
;
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU SOLV
"""


# ── SLURM script template ─────────────────────────────────────────────────────
# Targets hamsi partition: 56 cores/node, CPU-only (no GPU on hamsi).
# hamsi convention: --ntasks=28 --cpus-per-task=2, OMP_NUM_THREADS=2, -npme 4.
# Standard GROMACS 2024.1 (no PLUMED patch) — metadynamics runs use a separate
# PLUMED-patched build loaded via $HOME/load_gmx_plumed.sh.

SLURM_TEMPLATE = """\
#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A mozkurt
#SBATCH -J {label}_r{rep}
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

module purge
module load apps/gromacs/2024.1-oneapi2024

export OMP_NUM_THREADS=2

GMX=gmx_mpi

SETUP={cgdir}/gromacs
GRO=$SETUP/step3_input.gro
TOP=$SETUP/topol.top
NDX=$SETUP/index.ndx
EM_MDP=$SETUP/step4.0_minimization.mdp
NVT_MDP=$SETUP/step4.1_equilibration.mdp
NPT_MDP=../npt_equil.mdp
PROD_MDP=../prod.mdp

echo "===== {label} replica {rep} ====="
echo "SETUP : $SETUP"
echo "Start : $(date)"

# ── Stage 1: Energy minimisation ─────────────────────────────────────────────
if [ ! -f em.gro ]; then
    echo "--- Stage 1: Energy minimisation ---"
    $GMX grompp -f $EM_MDP -c $GRO -r $GRO -p $TOP -n $NDX \\
                -o em.tpr -maxwarn 1
    mpirun -np $SLURM_NTASKS $GMX mdrun -v -deffnm em -ntomp $OMP_NUM_THREADS
fi

# ── Stage 2: NVT equilibration (125 ps, 1 fs, heavy restraints) ──────────────
if [ ! -f nvt.gro ]; then
    echo "--- Stage 2: NVT equilibration ---"
    $GMX grompp -f $NVT_MDP -c em.gro -r $GRO -p $TOP -n $NDX \\
                -o nvt.tpr -maxwarn 1
    mpirun -np $SLURM_NTASKS $GMX mdrun -v -deffnm nvt -ntomp $OMP_NUM_THREADS
fi

# ── Stage 3: NPT equilibration (1 ns, 4 fs, light restraints) ────────────────
if [ ! -f npt.gro ]; then
    echo "--- Stage 3: NPT equilibration ---"
    $GMX grompp -f $NPT_MDP -c nvt.gro -r $GRO -p $TOP -n $NDX \\
                -o npt.tpr -maxwarn 1
    mpirun -np $SLURM_NTASKS $GMX mdrun -v -deffnm npt -ntomp $OMP_NUM_THREADS
fi

# ── Stage 4: Production MD (checkpoint-restartable) ──────────────────────────
# Re-submitting the same job continues from prod.cpt automatically.
echo "--- Stage 4: Production MD ---"
if [ ! -f prod.tpr ]; then
    $GMX grompp -f $PROD_MDP -c npt.gro -p $TOP -n $NDX \\
                -o prod.tpr -maxwarn 1
fi
# Only pass -cpi if a checkpoint exists (guards against mdrun error on first run)
CPT_FLAG=""
[ -f prod.cpt ] && CPT_FLAG="-cpi prod.cpt"
mpirun -np $SLURM_NTASKS $GMX mdrun -v -deffnm prod \\
       $CPT_FLAG -npme 4 -ntomp $OMP_NUM_THREADS -maxh 11.9

echo "Done : $(date)"
"""


# ── Per-system setup ──────────────────────────────────────────────────────────

def setup_system(sname: str, cfg: dict, md_dir: Path) -> dict | None:
    sdir  = md_dir / sname
    cgdir = get_charmm_dir(sdir)

    if cgdir is None:
        log(f"  [SKIP] {sname}: no charmm-gui-* directory found in {sdir}", RED)
        return None

    log(f"\n{CYAN}▶ {sname}{RESET}")

    # Write npt_equil.mdp at system level (shared by all replicas)
    npt_path = sdir / "npt_equil.mdp"
    npt_path.write_text(NPT_EQUIL_TEMPLATE)
    log(f"  npt_equil.mdp  → 1 ns NPT equilibration", GREEN)

    # nsteps computed with integer arithmetic to avoid floating-point rounding.
    # ns_prod × 1_000_000 converts ns → fs; dividing by DT_PROD_FS (integer)
    # is exact for all production lengths used in this project.
    nsteps = cfg["ns_prod"] * 1_000_000 // DT_PROD_FS
    prod_mdp = PROD_TEMPLATE.format(ns=cfg["ns_prod"], nsteps=nsteps)
    prod_path = sdir / "prod.mdp"
    prod_path.write_text(prod_mdp)
    log(f"  prod.mdp       → {cfg['ns_prod']} ns  ({nsteps:,} steps)", GREEN)

    # Create run directories and SLURM scripts
    run_dirs = []
    for rep in range(1, cfg["nreps"] + 1):
        run_dir = sdir / f"run{rep}"
        run_dir.mkdir(exist_ok=True)

        cgdir_rel = os.path.relpath(cgdir, run_dir)
        script = SLURM_TEMPLATE.format(
            label  = cfg["label"],
            rep    = rep,
            cgdir  = cgdir_rel,
        )
        script_path = run_dir / "run.sh"
        script_path.write_text(script)
        script_path.chmod(script_path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP)
        run_dirs.append(str(run_dir))

    label_dirs = "  ".join(f"run{r}/" for r in range(1, cfg["nreps"] + 1))
    log(f"  {label_dirs} created with run.sh", GREEN)

    return {
        "system":    sname,
        "cgdir":     str(cgdir),
        "ns_prod":   cfg["ns_prod"],
        "nsteps":    nsteps,
        "nreps":     cfg["nreps"],
        "run_dirs":  run_dirs,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Set up TRUBA MD run directories for CRYPTAD.")
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT,
                        help=f"Path to CRYPTAD project root (default: {DEFAULT_ROOT})")
    args = parser.parse_args()

    root = args.project_root.resolve()
    md_dir = root / "02_md_inputs"

    if not (root / "01_structures").is_dir():
        log(f"[ERROR] Project root does not contain 01_structures/: {root}", RED)
        log( "        Pass the correct path with --project-root.", RED)
        sys.exit(1)

    if not md_dir.is_dir():
        log(f"[ERROR] MD inputs directory not found: {md_dir}", RED)
        log( "        This directory (~73 GB) is not tracked in git.", RED)
        log( "        Retrieve from the Zenodo data deposit or TRUBA scratch.", RED)
        sys.exit(1)

    log(f"\n{'═'*60}", BOLD)
    log("  CRYPTAD — Setting up TRUBA run directories", BOLD)
    log(f"{'═'*60}", BOLD)
    log(f"  Project root : {root}")
    log(f"  MD inputs    : {md_dir}\n")

    manifest_entries = []
    for sname, cfg in SYSTEMS.items():
        entry = setup_system(sname, cfg, md_dir)
        if entry:
            manifest_entries.append(entry)

    # Write manifest
    if manifest_entries:
        manifest = {
            "generated_at":  datetime.now(timezone.utc).isoformat(),
            "dt_prod_fs":    DT_PROD_FS,
            "dt_npt_fs":     DT_NPT_FS,
            "dt_equil_fs":   DT_EQUIL_FS,
            "slurm_partition": "hamsi",
            "gromacs_module":  "apps/gromacs/2024.1-oneapi2024",
            "systems":       manifest_entries,
        }
        manifest_path = md_dir / "setup_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2))
        log(f"\n  Manifest written → {manifest_path.relative_to(root)}", CYAN)

    log(f"\n{'═'*60}", BOLD)
    log("  CRYPTAD — TRUBA run setup complete", BOLD)
    log(f"{'═'*60}\n", BOLD)
    log(f"""\
  Transfer to TRUBA (from Mac):
    rsync -avz --progress \\
      "02_md_inputs/" \\
      mozkurt@levrek.ulakbim.gov.tr:/arf/scratch/mozkurt/CRYPTAD/02_md_inputs/

  Submit on TRUBA (per system, per replica):
    cd /arf/scratch/mozkurt/CRYPTAD/02_md_inputs/S1_BIN1_BAR/run1 && sbatch run.sh
    (repeat for run2, run3, and S2–S5)

  Expected walltime (hamsi, 28×2 CPU):
    S1  ~12–17 days/replica  (checkpoint-restartable, 12 h job slices)
    S2  ~1 day/replica
    S3  ~3–4 days/replica
    S4  <0.5 days/replica
    S5  <0.5 days/replica
""", YELLOW)


if __name__ == "__main__":
    main()
