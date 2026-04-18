#!/bin/bash
# run_persist_S1.sh  —  CRYPTAD
# SLURM array job: 3 × 20 ns unbiased persistence gate for S1 BIN1 BAR site688
# Each array task is an independent replica with a unique velocity seed.
#
# Prerequisites (run first on login node):
#   bash 09_scripts/06_persistence_gate/01_extract_persist_start_S1.sh
#
# Submit:
#   sbatch 09_scripts/06_persistence_gate/02_run_persist_S1.sh
#
# Output:
#   02_md_simulations/BIN1/persistence_gate/rep{1,2,3}/persist.xtc
#
#SBATCH --job-name=persist_S1
#SBATCH -p hamsi
#SBATCH -A mozkurt
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=2
#SBATCH --array=1-3
#SBATCH --time=0-12:00:00
#SBATCH --output=slurm_persist_S1_%A_%a.out
#SBATCH --error=slurm_persist_S1_%A_%a.err

set -euo pipefail

# ── Environment (standard GROMACS — no PLUMED needed) ────────────────────────
module load apps/gromacs/2024.1-oneapi2024
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
GMX=gmx_mpi

# ── Paths ─────────────────────────────────────────────────────────────────────
CRYPTAD=${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}
SCRIPTS=$CRYPTAD/09_scripts/06_persistence_gate
GATEDIR=$CRYPTAD/02_md_simulations/BIN1/persistence_gate
START=$GATEDIR/site688_start.gro
TOP_CANDIDATES=($CRYPTAD/02_md_inputs/S1_BIN1_BAR/charmm-gui-*/gromacs/topol.top)
NDX_CANDIDATES=($CRYPTAD/02_md_inputs/S1_BIN1_BAR/charmm-gui-*/gromacs/index.ndx)
TOP=${TOP_CANDIDATES[0]:-}
NDX=${NDX_CANDIDATES[0]:-}

if [[ ! -f "$START" ]]; then
    echo "ERROR: $START not found."
    echo "       Run 01_extract_persist_start_S1.sh on login node first."
    exit 1
fi
if [[ -z "$TOP" || -z "$NDX" ]]; then
    echo "ERROR: topol.top or index.ndx not found in charmm-gui-* directory."
    exit 1
fi

# ── Replica setup ─────────────────────────────────────────────────────────────
REP=$SLURM_ARRAY_TASK_ID
REPDIR=$GATEDIR/rep${REP}
mkdir -p "$REPDIR"
cd "$REPDIR"

# Unique velocity seed per replica (deterministic, reproducible)
SEED=$((REP * 314159 + 271828))

echo "=== S1 BIN1 BAR — Persistence gate replica ${REP} ==="
echo "REPDIR : $REPDIR"
echo "START  : $START"
echo "SEED   : $SEED"
echo "Start  : $(date)"
echo ""

# ── Generate replica MDP (substitute GEN_SEED) ───────────────────────────────
sed "s/GEN_SEED/${SEED}/" "$SCRIPTS/persist_template.mdp" > persist.mdp

# ── grompp (only if persist.tpr not yet built) ───────────────────────────────
if [[ ! -f persist.tpr ]]; then
    echo "--- grompp ---"
    $GMX grompp \
        -f persist.mdp \
        -c "$START" \
        -p "$TOP" \
        -n "$NDX" \
        -o persist.tpr \
        -maxwarn 1
fi

# ── mdrun (checkpoint-restartable) ───────────────────────────────────────────
CPT_FLAG=""
[[ -f persist.cpt ]] && CPT_FLAG="-cpi persist.cpt"

echo "--- mdrun: 20 ns unbiased MD (replica ${REP}) ---"
mpirun -np $SLURM_NTASKS $GMX mdrun \
    -v \
    -deffnm persist \
    -s persist.tpr \
    $CPT_FLAG \
    -npme 4 \
    -ntomp $OMP_NUM_THREADS \
    -maxh 11.5

echo ""
echo "Done : $(date)"
NFRAMES=$(python3 -c "
import subprocess, re
out = subprocess.check_output(
    ['$GMX', 'check', '-f', 'persist.xtc'], stderr=subprocess.STDOUT
).decode()
m = re.search(r'Last frame\s+(\d+)', out)
print(m.group(1) if m else '?')
" 2>/dev/null || echo "?")
echo "Frames written: $NFRAMES"

# Write per-replica manifest
cat > "$REPDIR/persist_manifest_rep${REP}.json" <<EOF
{
  "generated_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "script": "09_scripts/06_persistence_gate/02_run_persist_S1.sh",
  "system": "S1_BIN1_BAR",
  "site": "site688",
  "replica": $REP,
  "seed": $SEED,
  "start_gro": "$START",
  "top": "$TOP",
  "ndx": "$NDX",
  "repdir": "$REPDIR",
  "slurm_job_id": "${SLURM_JOB_ID:-unknown}",
  "slurm_array_task_id": "${SLURM_ARRAY_TASK_ID:-unknown}",
  "n_frames": "$NFRAMES"
}
EOF
echo "Manifest → $REPDIR/persist_manifest_rep${REP}.json"
