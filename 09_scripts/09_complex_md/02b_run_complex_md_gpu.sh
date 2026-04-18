#!/bin/bash
#SBATCH --job-name=cryptad_cplxmd_gpu
#SBATCH --account=mozkurt
#SBATCH --partition=barbun-cuda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --gres=gpu:1
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --array=1-6
#SBATCH --signal=B:USR1@120

# ── CRYPTAD Phase 4 — Complex MD (BIN1 tasks, GPU acceleration) ──────────────
# Container: $HOME/containers/gromacs-2024.sif (GROMACS 2023.2, CUDA 12.1, thread_mpi)
# Binary:    /usr/local/gromacs/avx2_256/bin/gmx  (thread_mpi — NO mpirun)
# GPU:       barbun-cuda V100 (sm_70); 1 GPU + 20 CPU threads per task
#
# Checkpoint logic:
#   equil.cpt absent  → fresh start: run equil, then wipe complexmd.cpt → fresh production
#   equil.cpt present → resume:      skip equil, keep complexmd.cpt   → resume production
#
# Before first submit: delete equil.cpt + complexmd.cpt in each replica dir on TRUBA.

set -uo pipefail

CRYPTAD="${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}"
PAIRS="${PAIRS:-$CRYPTAD/09_scripts/09_complex_md/complex_md_pairs.txt}"
COMPLEXMD_BASE="$CRYPTAD/02_md_simulations/complex_md"
CONTAINER="$HOME/containers/gromacs-2024.sif"
GMX_BIN="/usr/local/gromacs/avx2_256/bin/gmx"
GMX="apptainer exec --nv $CONTAINER $GMX_BIN"

# ── Resolve system + replica ──────────────────────────────────────────────────
LINE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" "$PAIRS")
SYS=$(echo "$LINE" | awk '{print $1}')
REP=$(echo "$LINE" | awk '{print $2}')
WDIR="$COMPLEXMD_BASE/$SYS/$REP"

echo "=== Task ${SLURM_ARRAY_TASK_ID}: $SYS / $REP ==="
echo "    Working dir: $WDIR"
cd "$WDIR"

# ── Checkpoint handler (SIGTERM / USR1 at 2 min before wall time) ─────────────
_checkpoint() {
    echo "[$(date '+%H:%M:%S')] SIGTERM/USR1 — saving checkpoint"
    [ -f complexmd.cpt ] && cp complexmd.cpt complexmd_prev.cpt
}
trap '_checkpoint' USR1 TERM

# ── Always regenerate TPRs (ensures correct MDP parameters) ───────────────────
rm -f equil.tpr complexmd.tpr

# ── Step 1: Berendsen re-equilibration (200 ps, dt=0.002, gen_vel=yes) ────────
if [ ! -f equil.cpt ]; then
    echo "[$(date '+%H:%M:%S')] Fresh start — running equil..."
    $GMX grompp \
        -f equil.mdp -c npt.gro -t npt.cpt \
        -p topol_flat.top -n index.ndx \
        -o equil.tpr -maxwarn 1 \
        2>&1 | tee grompp_equil.log
    $GMX mdrun -v -deffnm equil -maxh 1.0 \
        -ntmpi 1 -ntomp 20 -nb gpu -pme gpu -bonded gpu \
        -nobackup \
        2>&1 | tee mdrun_equil.log
    echo "[$(date '+%H:%M:%S')] equil done."
    # Discard any stale production checkpoint — it may be from a different
    # dt or rank count and will cause GROMACS to report 0 s remaining.
    rm -f complexmd.cpt
else
    echo "[$(date '+%H:%M:%S')] equil.cpt exists — resuming session, skipping equil."
fi

# ── Step 2: C-rescale production (100 ns, dt=0.002) ───────────────────────────
echo "[$(date '+%H:%M:%S')] grompp — production..."
$GMX grompp \
    -f complexmd.mdp -c equil.gro -t equil.cpt \
    -p topol_flat.top -n index.ndx \
    -o complexmd.tpr -maxwarn 1 \
    2>&1 | tee grompp_prod.log
echo "[$(date '+%H:%M:%S')] grompp done."

CPT_FLAG=""
[[ -f complexmd.cpt ]] && CPT_FLAG="-cpi complexmd.cpt"

if [ -n "$CPT_FLAG" ]; then
    echo "[$(date '+%H:%M:%S')] Resuming production from checkpoint..."
else
    echo "[$(date '+%H:%M:%S')] Starting fresh production mdrun..."
fi

$GMX mdrun -v -deffnm complexmd -maxh 23.9 \
    -ntmpi 1 -ntomp 20 -nb gpu -pme gpu -bonded gpu \
    $CPT_FLAG \
    -nobackup \
    2>&1 | tee mdrun_prod.log

echo "=== DONE: $SYS / $REP ==="
