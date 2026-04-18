#!/bin/bash
# 01_run_metad_S3_PICALM_ANTH.sh  —  CRYPTAD
# WTMetaD for S3 PICALM ANTH domain
# Partition: orfoz (56 cores per node, PLUMED-patched GROMACS 2023.3)
#
# Strategy: 72 h self-terminating jobs chained via SLURM dependency.
#   mdrun -maxh 71.9 writes a checkpoint and exits cleanly before the 3-day
#   SLURM limit.  Re-submit with --dependency=afterany:$PREV_JOBID to extend.
#
# Submit (from CRYPTAD project root on TRUBA):
#   sbatch 09_scripts/03_metadynamics/01_run_metad_S3_PICALM_ANTH.sh
#
# Re-submit to extend (checkpoint restart):
#   sbatch --dependency=afterany:$PREV_JOBID \
#          09_scripts/03_metadynamics/01_run_metad_S3_PICALM_ANTH.sh
#
#SBATCH --job-name=meta_S3_ANTH
#SBATCH -p orfoz
#SBATCH -A mozkurt
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err

set -euo pipefail

# ── Environment (PLUMED-patched GROMACS 2023.3) ───────────────────────────────
source "$HOME/load_gmx_plumed.sh"
export OMP_NUM_THREADS=1
GMX=gmx_mpi

# ── Paths ─────────────────────────────────────────────────────────────────────
CRYPTAD=/arf/scratch/mozkurt/CRYPTAD
SCRIPTS=$CRYPTAD/09_scripts/03_metadynamics
RUNDIR=$CRYPTAD/02_md_simulations/PICALM/metadynamics/S3_ANTH

RUN1=$CRYPTAD/02_md_inputs/S3_PICALM_ANTH/run1
GRO=$RUN1/npt.gro
CPT=$RUN1/npt.cpt

# Locate topology and index inside the charmm-gui-XXXXXXXXXX directory.
# Use bash glob array instead of ls|head-1 — avoids SIGPIPE under set -o pipefail.
TOP_CANDIDATES=($CRYPTAD/02_md_inputs/S3_PICALM_ANTH/charmm-gui-*/gromacs/topol.top)
NDX_CANDIDATES=($CRYPTAD/02_md_inputs/S3_PICALM_ANTH/charmm-gui-*/gromacs/index.ndx)
TOP=${TOP_CANDIDATES[0]:-}
NDX=${NDX_CANDIDATES[0]:-}
if [[ -z "$TOP" || ! -f "$TOP" || -z "$NDX" || ! -f "$NDX" ]]; then
    echo "ERROR: Cannot find topol.top or index.ndx in charmm-gui-* directory."
    echo "  Searched: $CRYPTAD/02_md_inputs/S3_PICALM_ANTH/charmm-gui-*/gromacs/"
    exit 1
fi

MDP=$SCRIPTS/meta.mdp
PLUMED_DAT=$SCRIPTS/plumed_S3_PICALM_ANTH.dat

# ── Verify PLUMED dat file exists ─────────────────────────────────────────────
# plumed.dat defines the CVs, Gaussian hills width/height, and bias factor —
# it entirely determines the simulation.  Abort if not present.
if [[ ! -f "$PLUMED_DAT" ]]; then
    echo "ERROR: PLUMED input not found: $PLUMED_DAT"
    exit 1
fi

# ── Run directory ─────────────────────────────────────────────────────────────
mkdir -p "$RUNDIR"
cd "$RUNDIR"

echo "=== S3 PICALM ANTH — WTMetaD ==="
echo "RUNDIR  : $RUNDIR"
echo "TOP     : $TOP"
echo "GRO     : $GRO"
echo "PLUMED  : $PLUMED_DAT"
echo "SHA256  : $(sha256sum "$PLUMED_DAT" | awk '{print $1}')"
echo "Start   : $(date)"
echo ""

# Copy plumed.dat to run directory (PLUMED reads it from CWD).
cp "$PLUMED_DAT" ./plumed.dat

# ── Step 1: grompp (only if meta.tpr does not exist) ─────────────────────────
if [[ ! -f meta.tpr ]]; then
    echo "--- grompp: creating meta.tpr ---"
    $GMX grompp \
        -f "$MDP" \
        -c "$GRO" \
        -p "$TOP" \
        -n "$NDX" \
        -t "$CPT" \
        -o meta.tpr \
        -maxwarn 1
fi

# ── Step 2: mdrun with PLUMED (checkpoint-restartable) ───────────────────────
# -maxh 71.9 stops mdrun cleanly after 71.9 h so SLURM can finish before
# the 3-day wall limit.  On re-submission, -cpi meta.cpt resumes the run.
# Guard: only pass -cpi if a checkpoint exists (first run has none).
echo "--- mdrun: WTMetaD ---"
CPT_FLAG=""
[[ -f meta.cpt ]] && CPT_FLAG="-cpi meta.cpt"
mpirun -np "$SLURM_NTASKS" $GMX mdrun \
    -v \
    -deffnm meta \
    -s meta.tpr \
    -plumed plumed.dat \
    $CPT_FLAG \
    -maxh 71.9

echo ""
echo "Done : $(date)"
echo "HILLS lines : $(awk 'END{print NR}' HILLS  2>/dev/null || echo 'HILLS not found')"
echo "COLVAR lines: $(awk 'END{print NR}' COLVAR 2>/dev/null || echo 'COLVAR not found')"
