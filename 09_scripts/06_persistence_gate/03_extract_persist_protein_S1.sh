#!/bin/bash
# extract_persist_protein_S1.sh  —  CRYPTAD
# SLURM job: extract protein-only frames from S1 persistence gate replicas.
# Run AFTER run_persist_S1.sh completes (all 3 replicas done).
#
# For each rep{1,2,3}/persist.xtc → protein_traj.pdb
# Also copies to fpocket_results/S1_BIN1_BAR/persist/rep{N}/ for run_fpocket_frames.py
#
# Submit (after persistence gate replicas finish):
#   sbatch 09_scripts/06_persistence_gate/03_extract_persist_protein_S1.sh
#
#SBATCH --job-name=prot_persist_S1
#SBATCH -p hamsi
#SBATCH -A mozkurt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --output=slurm_prot_persist_S1_%j.out
#SBATCH --error=slurm_prot_persist_S1_%j.err

set -euo pipefail
module load apps/gromacs/2024.1-oneapi2024
GMX=gmx_mpi

CRYPTAD=${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}
RUNDIR=$CRYPTAD/02_md_simulations/BIN1/metadynamics/S1_BAR
GATEDIR=$CRYPTAD/02_md_simulations/BIN1/persistence_gate
FPDIR=$CRYPTAD/03_pocket_analysis/fpocket_results/S1_BIN1_BAR/persist
TPR=$RUNDIR/meta.tpr    # same topology as metadynamics run
DT_PS=100               # every 100 ps stride = 200 frames per 20 ns run
GRP_PROTEIN=1

echo "=== S1 persistence gate — protein frame extraction ==="
echo "Start: $(date)"

for REP in 1 2 3; do
    REPDIR=$GATEDIR/rep${REP}
    XTC=$REPDIR/persist.xtc
    OUT=$REPDIR/protein_traj.pdb

    if [[ ! -f "$XTC" ]]; then
        echo "[SKIP] rep${REP}: persist.xtc not found"
        continue
    fi

    echo ""
    echo "--- rep${REP} ---"

    NOJUMP=$REPDIR/nojump.xtc
    CLUSTER=$REPDIR/cluster.xtc
    CENTER=$REPDIR/center.xtc
    rm -f "$NOJUMP" "$CLUSTER" "$CENTER" "$OUT"

    # Step 1: nojump + stride
    printf "%d\n" "0" | $GMX trjconv \
        -s "$TPR" -f "$XTC" -o "$NOJUMP" -pbc nojump -dt "$DT_PS" -quiet

    # Step 2: cluster chains into same image (dimer)
    printf "%d\n%d\n" "$GRP_PROTEIN" "0" | $GMX trjconv \
        -s "$TPR" -f "$NOJUMP" -o "$CLUSTER" -pbc cluster -quiet

    # Step 3: center + compact
    printf "%d\n%d\n" "$GRP_PROTEIN" "0" | $GMX trjconv \
        -s "$TPR" -f "$CLUSTER" -o "$CENTER" -center -pbc mol -ur compact -quiet

    # Step 4: protein-only → multi-model PDB
    printf "%d\n" "$GRP_PROTEIN" | $GMX trjconv \
        -s "$TPR" -f "$CENTER" -o "$OUT" -quiet

    rm -f "$NOJUMP" "$CLUSTER" "$CENTER"

    NFRAMES=$(awk '/^MODEL/{n++} END{print n+0}' "$OUT")
    SIZE=$(du -h "$OUT" | awk '{print $1}' || echo "?")
    echo "  rep${REP}: $NFRAMES frames, $SIZE → $OUT"

    # Copy to fpocket_results directory (expected by run_fpocket_frames.py)
    FPREPDIR=$FPDIR/rep${REP}
    mkdir -p "$FPREPDIR"
    cp "$OUT" "$FPREPDIR/protein_traj.pdb"
    echo "  Copied to $FPREPDIR/protein_traj.pdb"
done

echo ""
echo "Done: $(date)"

# Write manifest
cat > "$GATEDIR/extract_persist_protein_S1_manifest.json" <<EOF
{
  "generated_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "script": "09_scripts/06_persistence_gate/03_extract_persist_protein_S1.sh",
  "system": "S1_BIN1_BAR",
  "tpr": "$TPR",
  "dt_ps": $DT_PS,
  "gatedir": "$GATEDIR",
  "fpdir": "$FPDIR",
  "slurm_job_id": "${SLURM_JOB_ID:-unknown}"
}
EOF
echo "Manifest → $GATEDIR/extract_persist_protein_S1_manifest.json"

echo ""
echo "Next steps:"
echo "  1. rsync $FPDIR to local 03_pocket_analysis/fpocket_results/S1_BIN1_BAR/persist/"
echo "  2. Locally: python3 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --sys S1_BIN1_BAR --run persist/rep1"
echo "              python3 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --sys S1_BIN1_BAR --run persist/rep2"
echo "              python3 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --sys S1_BIN1_BAR --run persist/rep3"
echo "  3. Locally: python3 09_scripts/06_persistence_gate/04_analyze_persistence_gate.py"
