#!/bin/bash
# extract_persist_protein_S3.sh  —  CRYPTAD
# SLURM job: extract protein-only frames from S3 persistence gate replicas.
# Run AFTER run_persist_S3.sh completes (all 3 replicas done).
#
# S3 PICALM ANTH is a monomer — 2-pass PBC (nojump → center+mol), no cluster step.
#
# Submit:
#   sbatch 09_scripts/06_persistence_gate/03_extract_persist_protein_S3.sh
#
#SBATCH --job-name=prot_persist_S3
#SBATCH -p hamsi
#SBATCH -A mozkurt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:20:00
#SBATCH --output=slurm_prot_persist_S3_%j.out
#SBATCH --error=slurm_prot_persist_S3_%j.err

set -euo pipefail
module load apps/gromacs/2024.1-oneapi2024
GMX=gmx_mpi

CRYPTAD=${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}
RUNDIR=$CRYPTAD/02_md_simulations/PICALM/metadynamics/S3_ANTH
GATEDIR=$CRYPTAD/02_md_simulations/PICALM/persistence_gate
FPDIR=$CRYPTAD/03_pocket_analysis/fpocket_results/S3_PICALM_ANTH/persist
TPR=$RUNDIR/meta.tpr
DT_PS=100
GRP_PROTEIN=1

echo "=== S3 PICALM ANTH — Persistence gate protein frame extraction ==="
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
    CENTER=$REPDIR/center.xtc
    rm -f "$NOJUMP" "$CENTER" "$OUT"

    # Step 1: nojump + stride
    printf "%d\n" "0" | $GMX trjconv \
        -s "$TPR" -f "$XTC" -o "$NOJUMP" -pbc nojump -dt "$DT_PS" -quiet

    # Step 2: center + compact → protein PDB (monomer: no cluster step)
    printf "%d\n%d\n" "$GRP_PROTEIN" "0" | $GMX trjconv \
        -s "$TPR" -f "$NOJUMP" -o "$CENTER" -center -pbc mol -ur compact -quiet

    # Step 3: protein-only multi-model PDB
    printf "%d\n" "$GRP_PROTEIN" | $GMX trjconv \
        -s "$TPR" -f "$CENTER" -o "$OUT" -quiet

    rm -f "$NOJUMP" "$CENTER"

    NFRAMES=$(awk '/^MODEL/{n++} END{print n+0}' "$OUT")
    SIZE=$(du -h "$OUT" | awk '{print $1}' || echo "?")
    echo "  rep${REP}: $NFRAMES frames, $SIZE → $OUT"

    FPREPDIR=$FPDIR/rep${REP}
    mkdir -p "$FPREPDIR"
    cp "$OUT" "$FPREPDIR/protein_traj.pdb"
    echo "  Copied to $FPREPDIR/protein_traj.pdb"
done

echo ""
echo "Done: $(date)"

# Write manifest
cat > "$GATEDIR/extract_persist_protein_S3_manifest.json" <<EOF
{
  "generated_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "script": "09_scripts/06_persistence_gate/03_extract_persist_protein_S3.sh",
  "system": "S3_PICALM_ANTH",
  "tpr": "$TPR",
  "dt_ps": $DT_PS,
  "gatedir": "$GATEDIR",
  "fpdir": "$FPDIR",
  "slurm_job_id": "${SLURM_JOB_ID:-unknown}"
}
EOF
echo "Manifest → $GATEDIR/extract_persist_protein_S3_manifest.json"

echo ""
echo "Next steps:"
echo "  1. rsync $FPDIR to local 03_pocket_analysis/fpocket_results/S3_PICALM_ANTH/persist/"
echo "  2. Locally: for rep in 1 2 3; do"
echo "       python3 09_scripts/05_pocket_detection/04_run_fpocket_frames.py \\"
echo "           --sys S3_PICALM_ANTH --run persist/rep\${rep}"
echo "     done"
echo "  3. Locally: python3 09_scripts/06_persistence_gate/04_analyze_persistence_gate.py --sys S3_PICALM_ANTH"
