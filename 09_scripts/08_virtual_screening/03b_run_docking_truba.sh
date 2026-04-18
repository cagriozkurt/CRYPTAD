#!/bin/bash
#SBATCH --job-name=cryptad_vina
#SBATCH --account=mozkurt
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=04:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
# NOTE: update --array upper bound to ceil(n_ligands / BATCH_SIZE) before submitting.
# Run 03_run_docking.py --no-rsync to regenerate with the correct value.
#SBATCH --array=1-25

set -uo pipefail

# ── Environment ──────────────────────────────────────────────────────────────
export PATH="$HOME/apps/vina:$PATH"

# ── Paths ────────────────────────────────────────────────────────────────────
PROJ="${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}"
LIGDIR="$PROJ/04_virtual_screening/ligands/pdbqt"
RECDIR="$PROJ/04_virtual_screening/receptors"
DOCKDIR="$PROJ/04_virtual_screening/docking_results"
MANIFEST="$PROJ/04_virtual_screening/ligands/ligand_ids.txt"

BATCH_SIZE=50
N_PARALLEL=14
VINA_CPU=4

mkdir -p "$DOCKDIR/logs"

# ── Select ligand batch for this task ────────────────────────────────────────
TOTAL=$(wc -l < "$MANIFEST")
START=$(( (SLURM_ARRAY_TASK_ID - 1) * BATCH_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * BATCH_SIZE ))
[[ $END -gt $TOTAL ]] && END=$TOTAL

LIGANDS=$(awk "NR>=$START && NR<=$END" "$MANIFEST")

# ── Per-ligand docking function ───────────────────────────────────────────────
dock_ligand() {
    LIG_ID="$1"
    LIG_PDBQT="$LIGDIR/${LIG_ID}.pdbqt"
    [[ ! -f "$LIG_PDBQT" ]] && { echo "SKIP $LIG_ID (pdbqt missing)"; return; }

    for POCKET_DIR in "$RECDIR"/*/; do
        POCKET_ID=$(basename "$POCKET_DIR")
        CONFIG="$POCKET_DIR/vina_config.txt"
        RECEPTOR="$POCKET_DIR/receptor.pdbqt"
        [[ ! -f "$CONFIG" || ! -f "$RECEPTOR" ]] && continue

        OUTDIR="$DOCKDIR/$POCKET_ID"
        mkdir -p "$OUTDIR"
        OUT_PDBQT="$OUTDIR/${LIG_ID}_out.pdbqt"
        OUT_LOG="$OUTDIR/${LIG_ID}_vina.txt"

        [[ -f "$OUT_PDBQT" ]] && continue   # already done — safe resume

        # Parse grid params from config directly to avoid any path-with-spaces issues
        CX=$(awk -F'[[:space:]]*=[[:space:]]*' '/^center_x/{print $2}'      "$CONFIG")
        CY=$(awk -F'[[:space:]]*=[[:space:]]*' '/^center_y/{print $2}'      "$CONFIG")
        CZ=$(awk -F'[[:space:]]*=[[:space:]]*' '/^center_z/{print $2}'      "$CONFIG")
        SX=$(awk -F'[[:space:]]*=[[:space:]]*' '/^size_x/{print $2}'        "$CONFIG")
        SY=$(awk -F'[[:space:]]*=[[:space:]]*' '/^size_y/{print $2}'        "$CONFIG")
        SZ=$(awk -F'[[:space:]]*=[[:space:]]*' '/^size_z/{print $2}'        "$CONFIG")
        EX=$(awk -F'[[:space:]]*=[[:space:]]*' '/^exhaustiveness/{print $2}' "$CONFIG")
        NM=$(awk -F'[[:space:]]*=[[:space:]]*' '/^num_modes/{print $2}'     "$CONFIG")
        ER=$(awk -F'[[:space:]]*=[[:space:]]*' '/^energy_range/{print $2}'  "$CONFIG")

        vina \
            --receptor     "$RECEPTOR" \
            --ligand       "$LIG_PDBQT" \
            --out          "$OUT_PDBQT" \
            --log          "$OUT_LOG" \
            --center_x     "$CX" --center_y "$CY" --center_z "$CZ" \
            --size_x       "$SX" --size_y   "$SY" --size_z   "$SZ" \
            --exhaustiveness "$EX" \
            --num_modes    "$NM" \
            --energy_range "$ER" \
            --cpu          "$VINA_CPU"
    done

    # Use first available pocket output for score summary
    FIRST_PDBQT=$(ls "$DOCKDIR"/*/"{$LIG_ID}_out.pdbqt" 2>/dev/null | head -1 || true)
    SCORE=$(awk '/^REMARK VINA RESULT/{print $4; exit}' "${FIRST_PDBQT:-/dev/null}" 2>/dev/null || echo "NA")
    echo "DONE $LIG_ID score=$SCORE" \
        >> "$DOCKDIR/logs/done_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
}

# ── Run N_PARALLEL ligands at a time (pure bash — no GNU parallel needed) ────
PIDS=()
while IFS= read -r LIG_ID; do
    [[ -z "$LIG_ID" ]] && continue
    dock_ligand "$LIG_ID" &
    PIDS+=($!)
    while (( ${#PIDS[@]} >= N_PARALLEL )); do
        NEW_PIDS=()
        for pid in "${PIDS[@]}"; do
            kill -0 "$pid" 2>/dev/null && NEW_PIDS+=("$pid")
        done
        PIDS=("${NEW_PIDS[@]}")
        (( ${#PIDS[@]} >= N_PARALLEL )) && sleep 1
    done
done <<< "$LIGANDS"
wait
