#!/bin/bash
# 02b_extract_metad_frames_S1.sh  —  CRYPTAD
# SLURM job to extract S1 BIN1 BAR protein frames from meta.xtc.
# Running trjconv on the login node is blocked by CPU time limit (379k atom system).
# Needs ~10 min on a compute node.
#
# Reproducibility notes
# ---------------------
# - Single task allocated: trjconv is serial; --ntasks=56 was wasteful.
# - No grep filter downstream of trjconv — grep exiting 1 (no matches) under
#   set -o pipefail aborts the job on a successful run.
# - Intermediate file sizes reported via du (not grep on binary XTC).
#
# Submit (from CRYPTAD project root on TRUBA):
#   sbatch 09_scripts/05_pocket_detection/02b_extract_metad_frames_S1.sh
#
#SBATCH --job-name=extract_S1
#SBATCH -p hamsi
#SBATCH -A mozkurt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --output=slurm_extract_S1_%j.out
#SBATCH --error=slurm_extract_S1_%j.err

set -euo pipefail
source "$HOME/load_gmx_plumed.sh"
GMX=gmx_mpi

CRYPTAD=/arf/scratch/mozkurt/CRYPTAD
RUNDIR=$CRYPTAD/02_md_simulations/BIN1/metadynamics/S1_BAR
XTC=$RUNDIR/meta.xtc
TPR=$RUNDIR/meta.tpr
OUT=$RUNDIR/protein_traj_metad.pdb

NOJUMP=$RUNDIR/meta_nojump.xtc
CLUSTER=$RUNDIR/meta_cluster.xtc
CENTER=$RUNDIR/meta_center.xtc

DT_PS=500
GRP_SYSTEM=0
GRP_PROTEIN=1

echo "=== S1 BIN1 BAR — metad frame extraction ==="
echo "Start : $(date)"
echo "RUNDIR: $RUNDIR"

# Clean up any partial output from previous login-node attempt
rm -f "$NOJUMP" "$CLUSTER" "$CENTER" "$OUT"

# ── Step 1: nojump + stride (reads all frames, writes every 500 ps) ──────────
echo "--- Step 1: nojump (dt=${DT_PS} ps) ---"
printf "%d\n" "$GRP_SYSTEM" | \
    "$GMX" trjconv \
        -s "$TPR" -f "$XTC" -o "$NOJUMP" \
        -pbc nojump \
        -dt "$DT_PS" \
        -quiet || {
    echo "ERROR Step 1 (nojump) failed"; exit 1
}
echo "Done: $(du -h "$NOJUMP" 2>/dev/null | awk '{print $1}' || echo '?')"

# ── Step 2: cluster dimer chains into same periodic image ────────────────────
echo "--- Step 2: pbc cluster ---"
printf "%d\n%d\n" "$GRP_PROTEIN" "$GRP_SYSTEM" | \
    "$GMX" trjconv \
        -s "$TPR" -f "$NOJUMP" -o "$CLUSTER" \
        -pbc cluster \
        -quiet || {
    echo "ERROR Step 2 (cluster) failed"; exit 1
}
echo "Done: $(du -h "$CLUSTER" 2>/dev/null | awk '{print $1}' || echo '?')"

# ── Step 3: center + compact wrapping ────────────────────────────────────────
echo "--- Step 3: center + pbc mol ---"
printf "%d\n%d\n" "$GRP_PROTEIN" "$GRP_SYSTEM" | \
    "$GMX" trjconv \
        -s "$TPR" -f "$CLUSTER" -o "$CENTER" \
        -center -pbc mol -ur compact \
        -quiet || {
    echo "ERROR Step 3 (center) failed"; exit 1
}
echo "Done: $(du -h "$CENTER" 2>/dev/null | awk '{print $1}' || echo '?')"

# ── Step 4: extract protein-only → multi-model PDB ───────────────────────────
echo "--- Step 4: extract Protein → PDB ---"
printf "%d\n" "$GRP_PROTEIN" | \
    "$GMX" trjconv \
        -s "$TPR" -f "$CENTER" -o "$OUT" \
        -quiet || {
    echo "ERROR Step 4 (protein extract) failed"; exit 1
}

NFRAMES=$(grep -c "^MODEL" "$OUT" 2>/dev/null || echo "?")
SIZE=$(du -h "$OUT" 2>/dev/null | awk '{print $1}' || echo "?")
echo "Done: $NFRAMES frames, $SIZE"

# ── Cleanup intermediates ─────────────────────────────────────────────────────
rm -f "$NOJUMP" "$CLUSTER" "$CENTER"

# ── Manifest ─────────────────────────────────────────────────────────────────
MANIFEST="$RUNDIR/extract_S1_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "gromacs_binary": "%s",\n' "$(command -v "$GMX")"
    printf '  "dt_ps": %d,\n' "$DT_PS"
    printf '  "nframes": "%s",\n' "$NFRAMES"
    printf '  "output_size": "%s",\n' "$SIZE"
    printf '  "output": "%s"\n' "$OUT"
    printf '}\n'
} > "$MANIFEST"
echo "Manifest → $MANIFEST"

echo "Finished: $(date)"
