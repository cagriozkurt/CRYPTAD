#!/usr/bin/env bash
# =============================================================================
# CRYPTAD — Extract protein frames from WTMetaD trajectories for fpocket
#
# Generates a multi-model protein-only PDB from meta.xtc for all 5 systems.
# Mirrors the production PBC pipeline (01_preprocess_trajectories.sh +
# 02_fix_pbc_dimer_s1.sh) adapted for metadynamics output.
#
# PBC treatment:
#   All systems  : nojump → center + pbc mol + ur compact → extract Protein
#   S1_BIN1_BAR  : nojump → pbc cluster (keeps dimer together) → center → Protein
#
# Stride: one frame every 500 ps  →  ~400 frames per 200 ns system
#
# Output per system:
#   $RUNDIR/protein_traj_metad.pdb   (multi-model PDB, protein only)
#
# Reproducibility notes
# ---------------------
# - No grep filter downstream of gmx trjconv — grep exiting 1 (no matches)
#   under set -o pipefail falsely triggers error handlers on a successful run.
# - A manifest is written on completion.
#
# Usage (TRUBA login node, interactive):
#   source "$HOME/load_gmx_plumed.sh"   # or: module load apps/gromacs/2024.1-oneapi2024
#   bash 09_scripts/05_pocket_detection/02_extract_metad_frames.sh
#
# After completion, rsync back to local and run:
#   bash 09_scripts/05_pocket_detection/01_setup_fpocket_dirs.sh   # local
#   python3 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --run metad --stride 1
# =============================================================================

set -euo pipefail

source "$HOME/load_gmx_plumed.sh"
GMX=gmx_mpi

CRYPTAD=/arf/scratch/mozkurt/CRYPTAD

# Frame interval: 500 ps.
# meta.xtc is written every 20 ps (nstxout-compressed=5000, dt=0.004 ps)
# → -dt 500 selects every 25th frame → ~400 frames / 200 ns
DT_PS=500

# GROMACS index groups (auto-generated from topology)
GRP_SYSTEM=0    # System  (all atoms)
GRP_PROTEIN=1   # Protein (backbone + side chain, no water/ions)

# System → run directory mapping (colon-split; preserves order, bash 3.2-safe)
SYSTEMS=(
    "S1_BIN1_BAR:$CRYPTAD/02_md_simulations/BIN1/metadynamics/S1_BAR"
    "S2_BIN1_SH3:$CRYPTAD/02_md_simulations/BIN1/metadynamics/S2_SH3"
    "S3_PICALM_ANTH:$CRYPTAD/02_md_simulations/PICALM/metadynamics/S3_ANTH"
    "S4_CD2AP_SH3-2:$CRYPTAD/02_md_simulations/CD2AP/metadynamics/S4_SH3-2"
    "S5_CD2AP_SH3-1:$CRYPTAD/02_md_simulations/CD2AP/metadynamics/S5_SH3-1"
)

log() { echo "[$(date '+%H:%M:%S')] $*"; }

FAILED=0
COMPLETED=()

for ENTRY in "${SYSTEMS[@]}"; do
    SYS="${ENTRY%%:*}"
    RUNDIR="${ENTRY#*:}"
    XTC="$RUNDIR/meta.xtc"
    TPR="$RUNDIR/meta.tpr"
    OUT="$RUNDIR/protein_traj_metad.pdb"

    NOJUMP="$RUNDIR/meta_nojump.xtc"
    CLUSTER="$RUNDIR/meta_cluster.xtc"   # S1 only
    CENTER="$RUNDIR/meta_center.xtc"

    if [[ ! -f "$XTC" || ! -f "$TPR" ]]; then
        log "SKIP  $SYS: meta.xtc or meta.tpr not found at $RUNDIR"
        continue
    fi

    if [[ -f "$OUT" ]]; then
        NFRAMES=$(grep -c "^MODEL" "$OUT" 2>/dev/null || echo "?")
        log "SKIP  $SYS: protein_traj_metad.pdb already exists ($NFRAMES frames)"
        continue
    fi

    log "INFO  ===== $SYS ====="
    log "INFO  XTC: $XTC"

    # ── Step 1: nojump — remove PBC discontinuities; apply stride ────────────
    log "RUN   Step 1: nojump (dt=${DT_PS} ps) → meta_nojump.xtc"
    printf "%d\n" "$GRP_SYSTEM" | \
        "$GMX" trjconv \
            -s "$TPR" -f "$XTC" -o "$NOJUMP" \
            -pbc nojump \
            -dt "$DT_PS" \
            -quiet || {
        log "ERROR Step 1 failed for $SYS"
        ((FAILED++)) || true
        continue
    }
    log "DONE  meta_nojump.xtc"

    # ── Step 2a (S1 only): cluster dimer chains into same periodic image ──────
    if [[ "$SYS" == "S1_BIN1_BAR" ]]; then
        log "RUN   Step 2a [S1 dimer]: pbc cluster → meta_cluster.xtc"
        printf "%d\n%d\n" "$GRP_PROTEIN" "$GRP_SYSTEM" | \
            "$GMX" trjconv \
                -s "$TPR" -f "$NOJUMP" -o "$CLUSTER" \
                -pbc cluster \
                -quiet || {
            log "ERROR Step 2a failed for $SYS"
            ((FAILED++)) || true
            continue
        }
        log "DONE  meta_cluster.xtc"
        INPUT_FOR_CENTER="$CLUSTER"
    else
        INPUT_FOR_CENTER="$NOJUMP"
    fi

    # ── Step 2b: center protein + compact wrapping → full system XTC ──────────
    log "RUN   Step 2b: center + pbc mol → meta_center.xtc"
    printf "%d\n%d\n" "$GRP_PROTEIN" "$GRP_SYSTEM" | \
        "$GMX" trjconv \
            -s "$TPR" -f "$INPUT_FOR_CENTER" -o "$CENTER" \
            -center -pbc mol -ur compact \
            -quiet || {
        log "ERROR Step 2b failed for $SYS"
        ((FAILED++)) || true
        continue
    }
    log "DONE  meta_center.xtc"

    # ── Step 3: extract protein-only → multi-model PDB ───────────────────────
    log "RUN   Step 3: extract Protein → protein_traj_metad.pdb"
    printf "%d\n" "$GRP_PROTEIN" | \
        "$GMX" trjconv \
            -s "$TPR" -f "$CENTER" -o "$OUT" \
            -quiet || {
        log "ERROR Step 3 failed for $SYS"
        ((FAILED++)) || true
        continue
    }

    NFRAMES=$(grep -c "^MODEL" "$OUT" 2>/dev/null || echo "?")
    SIZE=$(du -sh "$OUT" 2>/dev/null | awk '{print $1}' || echo "?")
    log "DONE  protein_traj_metad.pdb — $NFRAMES frames, $SIZE"

    # ── Cleanup intermediate XTCs ─────────────────────────────────────────────
    rm -f "$NOJUMP" "$CENTER"
    [[ "$SYS" == "S1_BIN1_BAR" ]] && rm -f "$CLUSTER"
    log "INFO  Intermediates cleaned up"
    COMPLETED+=("$SYS")
done

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
MANIFEST="$CRYPTAD/02_md_simulations/extract_metad_frames_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "gromacs_binary": "%s",\n' "$(command -v "$GMX")"
    printf '  "dt_ps": %d,\n' "$DT_PS"
    printf '  "grp_system": %d,\n' "$GRP_SYSTEM"
    printf '  "grp_protein": %d,\n' "$GRP_PROTEIN"
    printf '  "completed": [\n'
    for i in "${!COMPLETED[@]}"; do
        COMMA=","
        [[ $i -eq $((${#COMPLETED[@]}-1)) ]] && COMMA=""
        printf '    "%s"%s\n' "${COMPLETED[$i]}" "$COMMA"
    done
    printf '  ],\n'
    printf '  "failed": %d\n' "$FAILED"
    printf '}\n'
} > "$MANIFEST"
log "INFO  Manifest written → $MANIFEST"

echo ""
log "===== EXTRACTION COMPLETE ====="
if [[ "$FAILED" -gt 0 ]]; then
    log "  $FAILED system(s) FAILED — check output above"
else
    log "  All systems OK"
fi

echo ""
echo "Rsync command (run locally from CRYPTAD project root):"
echo "  rsync -avz --include='protein_traj_metad.pdb' --include='*/' --exclude='*' \\"
echo "    mozkurt@arf-ui.yonetim.truba.gov.tr:/arf/scratch/mozkurt/CRYPTAD/02_md_simulations/ \\"
echo "    ./02_md_simulations/"
