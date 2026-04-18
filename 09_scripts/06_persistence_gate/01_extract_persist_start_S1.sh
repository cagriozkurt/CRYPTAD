#!/bin/bash
# extract_persist_start_S1.sh  —  CRYPTAD
# Extract full-system GRO at t=151 ns (frame 302) from S1 meta.xtc.
# This is the starting structure for the 3-replica unbiased persistence gate.
#
# Single-frame dump is fast — run directly on login node (no SLURM needed).
#
# Usage (on TRUBA login node):
#   source $HOME/load_gmx_plumed.sh   # or module load standard GROMACS
#   bash 09_scripts/06_persistence_gate/01_extract_persist_start_S1.sh

set -euo pipefail
GMX=gmx_mpi

CRYPTAD=${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}
RUNDIR=$CRYPTAD/02_md_simulations/BIN1/metadynamics/S1_BAR
OUTDIR=$CRYPTAD/02_md_simulations/BIN1/persistence_gate
OUTGRO=$OUTDIR/site688_start.gro

# Frame 302 at dt=500 ps stride → t = 302 × 500 = 151000 ps = 151 ns
DUMP_PS=151000

mkdir -p "$OUTDIR"

echo "=== S1 BIN1 BAR — Persistence gate starting structure ==="
echo "Source  : $RUNDIR/meta.xtc"
echo "Dumping : t = ${DUMP_PS} ps (frame 302, site688 best-score frame)"
echo "Output  : $OUTGRO"
echo ""

# Group 0 = System (protein + water + ions) — full system needed for unbiased MD
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
    cat > "$OUTDIR/extract_persist_start_S1_manifest.json" <<EOF
{
  "generated_at": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "script": "09_scripts/06_persistence_gate/01_extract_persist_start_S1.sh",
  "source_xtc": "$RUNDIR/meta.xtc",
  "source_tpr": "$RUNDIR/meta.tpr",
  "dump_ps": $DUMP_PS,
  "output_gro": "$OUTGRO",
  "n_atoms": "$NATOMS"
}
EOF
    echo "Manifest → $OUTDIR/extract_persist_start_S1_manifest.json"
    echo "Next: sbatch 09_scripts/06_persistence_gate/02_run_persist_S1.sh"
else
    echo "ERROR: output GRO not created — check that meta.xtc contains t=${DUMP_PS} ps"
    exit 1
fi
