#!/usr/bin/env bash
# =============================================================================
# CRYPTAD — mdpocket batch script for all 15 MD trajectories
#
# Runs mdpocket (trajectory mode of fpocket) on every system/replica.
# Requires prod_protein.xtc from 01_preprocess_trajectories.sh.
#
# Per-run steps:
#   1. Extract protein reference PDB (frame 0 of production run)
#   2. Convert prod_protein.xtc → prod_protein_traj.pdb (multi-model PDB)
#      Reason: some fpocket builds have a heap-corruption bug in the XTC
#      reader. Multi-model PDB uses simpler I/O and is universally supported.
#   3. Run mdpocket with --trajectory_format pdb
#
# Output: 03_pocket_analysis/fpocket_results/{SYS}/{RUN}/
#   protein_ref.pdb             — reference structure used
#   protein_ref_mdpout_dens.dx  — overall pocket density (visualise in VMD)
#   protein_ref_mdpout_freq_pockets.dx — pocket frequency per grid point
#   protein_ref_mdpout_freq_atm.dx     — atom-level frequency
#   protein_ref_out/            — per-frame pocket info (pqr, info files)
#
# Reproducibility notes
# ---------------------
# - PROJECT_ROOT resolved as two directories above SCRIPT_DIR.
# - No grep/tail filters downstream of gmx or mdpocket — early pipe close under
#   set -o pipefail causes SIGPIPE and misroutes exit codes.
# - A manifest is written on completion.
#
# Usage:
#   bash 09_scripts/05_pocket_detection/03_run_mdpocket.sh
#   GMX=gmx_mpi bash 09_scripts/05_pocket_detection/03_run_mdpocket.sh
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Load GROMACS if needed
# ---------------------------------------------------------------------------
GMX="${GMX:-gmx}"
if ! command -v "${GMX}" &>/dev/null; then
    if command -v module &>/dev/null; then
        echo "[INFO] Loading GROMACS module"
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load apps/gromacs/2024.1-oneapi2024
    fi
    command -v "${GMX}" &>/dev/null || { echo "[ERROR] ${GMX} not found"; exit 1; }
fi

# ---------------------------------------------------------------------------
# Check mdpocket is installed
# ---------------------------------------------------------------------------
if ! command -v mdpocket &>/dev/null; then
    echo "[ERROR] mdpocket not found in PATH."
    echo ""
    echo "  Install fpocket (includes mdpocket):"
    echo "    conda install -c conda-forge fpocket"
    echo "  or build from source:"
    echo "    git clone https://github.com/Discngine/fpocket.git"
    echo "    cd fpocket && make && sudo make install"
    echo ""
    echo "  On TRUBA, try: module avail fpocket"
    exit 1
fi

# awk avoids head-1 early pipe close (SIGPIPE) under set -o pipefail
FPOCKET_VERSION=$(fpocket --version 2>&1 | awk 'NR==1{print; exit}' || echo "unknown")
echo "[INFO] Using: mdpocket (${FPOCKET_VERSION})"

# ---------------------------------------------------------------------------
# Path resolution
# Script lives at <root>/09_scripts/05_pocket_detection/03_run_mdpocket.sh
# SCRIPT_DIR   = <root>/09_scripts/05_pocket_detection
# ../..         = <root>
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MD_ROOT="${PROJECT_ROOT}/02_md_inputs"
POCKET_ROOT="${PROJECT_ROOT}/03_pocket_analysis/fpocket_results"

SYSTEMS=(
    S1_BIN1_BAR
    S2_BIN1_SH3
    S3_PICALM_ANTH
    "S4_CD2AP_SH3-2"
    "S5_CD2AP_SH3-1"
)
RUNS=(run1 run2 run3)

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
TOTAL=0; SKIPPED=0; FAILED=0
COMPLETED=()

for SYS in "${SYSTEMS[@]}"; do
    for RUN in "${RUNS[@]}"; do
        RUNDIR="${MD_ROOT}/${SYS}/${RUN}"
        OUTDIR="${POCKET_ROOT}/${SYS}/${RUN}"

        TPR="${RUNDIR}/prod.tpr"
        CENTER="${RUNDIR}/prod_center.xtc"
        PROTEIN_XTC="${RUNDIR}/prod_protein.xtc"
        REF_PDB="${OUTDIR}/protein_ref.pdb"
        TRAJ_PDB="${OUTDIR}/protein_traj.pdb"

        if [[ ! -f "${PROTEIN_XTC}" ]]; then
            log "WARN  ${SYS}/${RUN}: prod_protein.xtc not found — run 01_preprocess_trajectories.sh first"
            ((FAILED++)) || true
            continue
        fi

        mkdir -p "${OUTDIR}"
        log "INFO  ===== ${SYS} / ${RUN} ====="
        TOTAL=$((TOTAL + 1))

        # ------------------------------------------------------------------
        # Step 1 — Extract protein reference PDB (frame 0 of production)
        # ------------------------------------------------------------------
        if [[ -f "${REF_PDB}" ]]; then
            log "SKIP  protein_ref.pdb already exists"
        else
            log "RUN   Extracting reference PDB (frame 0)"
            printf "1\n" | \
                "${GMX}" trjconv \
                    -s  "${TPR}" \
                    -f  "${CENTER}" \
                    -o  "${REF_PDB}" \
                    -dump 0 \
                    -quiet || {
                log "ERROR Reference PDB extraction failed for ${SYS}/${RUN}"
                ((FAILED++)) || true
                continue
            }
            log "DONE  protein_ref.pdb"
        fi

        # ------------------------------------------------------------------
        # Step 2 — Convert XTC → multi-model PDB
        # Bypasses mdpocket's XTC reader which has a heap-corruption bug
        # in some builds. PDB format uses simpler I/O, universally supported.
        # Note: file will be large for long trajectories (S1: ~1–2 GB).
        # ------------------------------------------------------------------
        if [[ -f "${TRAJ_PDB}" ]]; then
            log "SKIP  protein_traj.pdb already exists"
        else
            log "RUN   Converting XTC → multi-model PDB (may take a few minutes)"
            printf "1\n" | \
                "${GMX}" trjconv \
                    -s  "${TPR}" \
                    -f  "${PROTEIN_XTC}" \
                    -o  "${TRAJ_PDB}" \
                    -quiet || {
                log "ERROR XTC→PDB conversion failed for ${SYS}/${RUN}"
                ((FAILED++)) || true
                continue
            }
            SIZE=$(du -h "${TRAJ_PDB}" 2>/dev/null | awk '{print $1}' || echo "?")
            log "DONE  protein_traj.pdb (${SIZE})"
        fi

        # ------------------------------------------------------------------
        # Step 3 — Run mdpocket with multi-model PDB trajectory
        # ------------------------------------------------------------------
        DONE_SENTINEL="${OUTDIR}/protein_ref_mdpout_freq_pockets.dx"
        if [[ -f "${DONE_SENTINEL}" ]]; then
            log "SKIP  mdpocket output already exists"
            ((SKIPPED++)) || true
            continue
        fi

        log "RUN   mdpocket — this may take several minutes for long trajectories"
        (
            cd "${OUTDIR}"
            mdpocket \
                -f  "${REF_PDB}" \
                --trajectory_file   "${TRAJ_PDB}" \
                --trajectory_format pdb
        ) || {
            log "ERROR mdpocket failed for ${SYS}/${RUN}"
            ((FAILED++)) || true
            continue
        }

        if [[ -f "${DONE_SENTINEL}" ]]; then
            log "DONE  ${SYS}/${RUN}"
            COMPLETED+=("${SYS}/${RUN}")
        else
            log "WARN  mdpocket finished but expected output not found — check ${OUTDIR}"
            ((FAILED++)) || true
        fi
    done
done

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
mkdir -p "${POCKET_ROOT}"
MANIFEST="${POCKET_ROOT}/mdpocket_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "gromacs_binary": "%s",\n' "$(command -v "${GMX}")"
    printf '  "gromacs_version": "%s",\n' "$("${GMX}" --version 2>&1 | awk '/GROMACS version/{print $NF}')"
    printf '  "fpocket_version": "%s",\n' "${FPOCKET_VERSION}"
    printf '  "completed_runs": [\n'
    for i in "${!COMPLETED[@]}"; do
        COMMA=","
        [[ $i -eq $((${#COMPLETED[@]}-1)) ]] && COMMA=""
        printf '    "%s"%s\n' "${COMPLETED[$i]}" "$COMMA"
    done
    printf '  ],\n'
    printf '  "total": %d,\n' "${TOTAL}"
    printf '  "skipped": %d,\n' "${SKIPPED}"
    printf '  "failed": %d\n' "${FAILED}"
    printf '}\n'
} > "${MANIFEST}"
log "INFO  Manifest written → ${MANIFEST}"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
log "===== MDPOCKET SUMMARY ====="
log "  Runs attempted : ${TOTAL}"
log "  Skipped        : ${SKIPPED}"
log "  Failed         : ${FAILED}"
if [[ "${FAILED}" -eq 0 ]]; then
    log "  Status: ALL OK"
    log ""
    log "  Output in: ${POCKET_ROOT}"
    log "  Key files per run:"
    log "    *_mdpout_freq_pockets.dx — pocket frequency (load in VMD with volmap)"
    log "    *_mdpout_dens.dx         — density grid"
    log "    *_out/                   — per-frame pocket descriptors"
    log ""
    log "  Next: rsync 03_pocket_analysis/ back, then run:"
    log "    python 09_scripts/05_pocket_detection/05_parse_fpocket.py"
else
    log "  Status: ${FAILED} failure(s) — check output above"
fi
