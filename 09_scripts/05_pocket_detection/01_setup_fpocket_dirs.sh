#!/usr/bin/env bash
# =============================================================================
# CRYPTAD — Local setup: copy metad protein PDBs into fpocket_results tree
#
# Run this LOCALLY after rsyncing protein_traj_metad.pdb from TRUBA.
# Creates 03_pocket_analysis/fpocket_results/{SYS}/metad/protein_traj.pdb
# so that run_fpocket_frames.py --run metad can find the files.
#
# Usage (from CRYPTAD project root):
#   bash 09_scripts/05_pocket_detection/01_setup_fpocket_dirs.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SIMDIR="$PROJECT_ROOT/02_md_simulations"
POCKETDIR="$PROJECT_ROOT/03_pocket_analysis/fpocket_results"

# Format: "SYS_NAME:relative/path/to/metad/rundir"
SYSTEMS=(
    "S1_BIN1_BAR:BIN1/metadynamics/S1_BAR"
    "S2_BIN1_SH3:BIN1/metadynamics/S2_SH3"
    "S3_PICALM_ANTH:PICALM/metadynamics/S3_ANTH"
    "S4_CD2AP_SH3-2:CD2AP/metadynamics/S4_SH3-2"
    "S5_CD2AP_SH3-1:CD2AP/metadynamics/S5_SH3-1"
)

OK=0; MISSING=0

for ENTRY in "${SYSTEMS[@]}"; do
    SYS="${ENTRY%%:*}"
    RELPATH="${ENTRY#*:}"
    SRC="$SIMDIR/$RELPATH/protein_traj_metad.pdb"
    DESTDIR="$POCKETDIR/$SYS/metad"
    DEST="$DESTDIR/protein_traj.pdb"

    if [[ ! -f "$SRC" ]]; then
        echo "[MISSING] $SYS: $SRC not found — rsync from TRUBA first"
        ((MISSING++)) || true
        continue
    fi

    mkdir -p "$DESTDIR"

    if [[ -f "$DEST" ]]; then
        echo "[SKIP]    $SYS: protein_traj.pdb already in fpocket_results"
    else
        cp "$SRC" "$DEST"
        NFRAMES=$(grep -c "^MODEL" "$DEST" 2>/dev/null || echo "?")
        echo "[OK]      $SYS → $DESTDIR/ ($NFRAMES frames)"
    fi
    ((OK++)) || true
done

echo ""
if [[ "$MISSING" -gt 0 ]]; then
    echo "$MISSING system(s) missing. Rsync command:"
    echo "  rsync -avz --include='protein_traj_metad.pdb' --include='*/' --exclude='*' \\"
    echo "    mozkurt@arf-ui.yonetim.truba.gov.tr:/arf/scratch/mozkurt/CRYPTAD/02_md_simulations/ \\"
    echo "    \"$SIMDIR/\""
else
    echo "All $OK systems set up. Run fpocket:"
    echo "  source .venv/bin/activate"
    for ENTRY in "${SYSTEMS[@]}"; do
        SYS="${ENTRY%%:*}"
        STRIDE=1; [[ "$SYS" == "S1_BIN1_BAR" ]] && STRIDE=2
        echo "  python3 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --sys $SYS --run metad --stride $STRIDE"
    done
fi

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
MANIFEST="$POCKETDIR/setup_fpocket_dirs_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "ok": %d,\n' "$OK"
    printf '  "missing": %d\n' "$MISSING"
    printf '}\n'
} > "$MANIFEST"
echo "Manifest → $MANIFEST"
