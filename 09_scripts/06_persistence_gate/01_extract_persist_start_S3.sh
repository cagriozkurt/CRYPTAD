#!/bin/bash
# extract_persist_start_S3.sh  —  CRYPTAD
# Extract full-system GRO at t=192 ns (frame 384) from S3 meta.xtc.
# Starting structure for 3-replica unbiased persistence gate of S3 PICALM ANTH site473.
#
# S3 is a monomer — single chain, no cluster step needed.
# Smaller system than S1 — safe to run on login node.
#
# Usage (on TRUBA login node):
#   source $HOME/load_gmx_plumed.sh
#   bash 09_scripts/06_persistence_gate/01_extract_persist_start_S3.sh

set -euo pipefail
GMX=gmx_mpi

CRYPTAD=${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}
RUNDIR=$CRYPTAD/02_md_simulations/PICALM/metadynamics/S3_ANTH
OUTDIR=$CRYPTAD/02_md_simulations/PICALM/persistence_gate
OUTGRO=$OUTDIR/site473_start.gro

# Frame 384 at dt=500 ps stride → t = 384 × 500 = 192000 ps = 192 ns
DUMP_PS=192000

mkdir -p "$OUTDIR"

echo "=== S3 PICALM ANTH — Persistence gate starting structure ==="
echo "Source  : $RUNDIR/meta.xtc"
echo "Dumping : t = ${DUMP_PS} ps (frame 384, site473 best-score frame)"
echo "Output  : $OUTGRO"
echo ""

# Group 0 = System (protein + water + ions)
echo "0" | $GMX trjconv \
    -s  "$RUNDIR/meta.tpr" \
    -f  "$RUNDIR/meta.xtc" \
    -o  "$OUTGRO" \
    -dump "$DUMP_PS" \
    -quiet

if [[ -f "$OUTGRO" ]]; then
    NATOMS=$(awk 'NR==2{print; exit}' "$OUTGRO" | tr -d ' ')
    echo ""
    echo "Done: $OUTGRO  ($NATOMS atoms)"

    # Write manifest
    cat > "$OUTDIR/extract_persist_start_S3_manifest.json" <<EOF
{
  "generated_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "script": "09_scripts/06_persistence_gate/01_extract_persist_start_S3.sh",
  "source_xtc": "$RUNDIR/meta.xtc",
  "source_tpr": "$RUNDIR/meta.tpr",
  "dump_ps": $DUMP_PS,
  "output_gro": "$OUTGRO",
  "n_atoms": "$NATOMS"
}
EOF
    echo "Manifest → $OUTDIR/extract_persist_start_S3_manifest.json"
    echo "Next: sbatch 09_scripts/06_persistence_gate/02_run_persist_S3.sh"
else
    echo "ERROR: output GRO not created — check that meta.xtc contains t=${DUMP_PS} ps"
    exit 1
fi
