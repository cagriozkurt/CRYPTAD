#!/bin/bash
#SBATCH --job-name=cryptad_cplxmd
#SBATCH --account=mozkurt
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=2
#SBATCH --time=3-00:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-12
#SBATCH --signal=B:USR1@120

# ── CRYPTAD Phase 4 — Complex MD 100 ns (Step 5.1) ───────────────────────────
# 12 tasks: 4 systems × 3 replicas
# Each task: grompp + mdrun 100 ns with checkpoint restart
#
# pairs.txt layout (no header):
#   S1_site688_CHEMBL4105630 rep1
#   S1_site688_CHEMBL4105630 rep2
#   ...  (12 lines total)

set -uo pipefail

CRYPTAD="${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}"
PAIRS="$CRYPTAD/09_scripts/09_complex_md/complex_md_pairs.txt"
COMPLEXMD_BASE="$CRYPTAD/02_md_simulations/complex_md"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}

# ── Load GROMACS (standard, no PLUMED needed) ─────────────────────────────────
module load apps/gromacs/2024.1-oneapi2024
GMX="gmx_mpi"

# ── Resolve system + replica from task ID ────────────────────────────────────
LINE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" "$PAIRS")
SYS=$(echo "$LINE" | awk '{print $1}')
REP=$(echo "$LINE" | awk '{print $2}')
WDIR="$COMPLEXMD_BASE/$SYS/$REP"

echo "=== Task ${SLURM_ARRAY_TASK_ID}: $SYS / $REP ==="
echo "    Working dir: $WDIR"
cd "$WDIR"

# ── Checkpoint restart handler ────────────────────────────────────────────────
_checkpoint() {
    echo "[$(date '+%H:%M:%S')] SIGTERM/USR1 received — saving checkpoint"
    if [ -f complexmd.cpt ]; then
        cp complexmd.cpt complexmd_prev.cpt
    fi
}
trap '_checkpoint' USR1 TERM

# ── Step 1: Berendsen re-equilibration (200 ps, gen_vel=yes) ─────────────────
if [ ! -f equil.tpr ]; then
    echo "[$(date '+%H:%M:%S')] Running grompp (equil)..."
    $GMX grompp \
        -f equil.mdp \
        -c npt.gro \
        -t npt.cpt \
        -p topol_flat.top \
        -n index.ndx \
        -o equil.tpr \
        -maxwarn 1 \
        2>&1 | tee grompp_equil.log
fi

if [ ! -f equil.gro ]; then
    echo "[$(date '+%H:%M:%S')] Running mdrun (equil)..."
    mpirun -np "$SLURM_NTASKS" $GMX mdrun \
        -v -deffnm equil \
        -ntomp "$OMP_NUM_THREADS" -npme 4 \
        -nobackup \
        2>&1 | tee mdrun_equil.log
    echo "[$(date '+%H:%M:%S')] Equil done."
else
    echo "[$(date '+%H:%M:%S')] equil.gro exists — skipping equil."
fi

# ── Step 2: Production (100 ns, continuation=yes) ────────────────────────────
if [ ! -f complexmd.tpr ]; then
    echo "[$(date '+%H:%M:%S')] Running grompp (production)..."
    $GMX grompp \
        -f complexmd.mdp \
        -c equil.gro \
        -t equil.cpt \
        -p topol_flat.top \
        -n index.ndx \
        -o complexmd.tpr \
        -maxwarn 1 \
        2>&1 | tee grompp_prod.log
    echo "[$(date '+%H:%M:%S')] grompp done."
else
    echo "[$(date '+%H:%M:%S')] complexmd.tpr exists — skipping grompp."
fi

CPT_FLAG=""
[[ -f complexmd.cpt ]] && CPT_FLAG="-cpi complexmd.cpt"

if [ -n "$CPT_FLAG" ]; then
    echo "[$(date '+%H:%M:%S')] Resuming from checkpoint..."
else
    echo "[$(date '+%H:%M:%S')] Starting fresh mdrun..."
fi

mpirun -np "$SLURM_NTASKS" $GMX mdrun \
    -v -deffnm complexmd \
    -maxh 71.9 \
    -ntomp "$OMP_NUM_THREADS" -npme 4 \
    $CPT_FLAG \
    -nobackup \
    2>&1 | tee mdrun_prod.log

echo "=== DONE: $SYS / $REP ==="
