#!/usr/bin/env bash
# =============================================================================
# CRYPTAD — PBC preprocessing for all 15 MD trajectories
#
# Three steps per run:
#   Step 1: Fix periodic boundary jumps         → prod_nojump.xtc  (System)
#   Step 2: Center protein, compact box         → prod_center.xtc  (System)
#   Step 3: Extract protein-only frames         → prod_protein.xtc (Protein)
#
# Reproducibility notes
# ---------------------
# - MD_ROOT is resolved as two directories above SCRIPT_DIR so the script
#   works regardless of the working directory it is called from.
# - trjconv pipelines use no downstream grep — grep exiting 1 (no matches)
#   under set -o pipefail would falsely trigger error handlers even on
#   a successful run.  trjconv output goes to stdout/stderr directly.
# - A manifest (preprocess_manifest.json) is written on completion.
#
# Usage:
#   bash 09_scripts/04_trajectory_qc/01_preprocess_trajectories.sh
#   GMX=gmx_mpi bash 09_scripts/04_trajectory_qc/01_preprocess_trajectories.sh
#
# Already-existing output files are skipped automatically.
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Load GROMACS module if needed (TRUBA: apps/gromacs/2024.1-oneapi2024)
# ---------------------------------------------------------------------------
GMX="${GMX:-gmx}"

if ! command -v "${GMX}" &>/dev/null; then
    if command -v module &>/dev/null; then
        echo "[INFO] '${GMX}' not found — attempting: module load apps/gromacs/2024.1-oneapi2024"
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load apps/gromacs/2024.1-oneapi2024
    fi
    if ! command -v "${GMX}" &>/dev/null; then
        echo "[ERROR] '${GMX}' still not found after module load."
        echo "        Run:  module load apps/gromacs/2024.1-oneapi2024"
        echo "        Then: GMX=gmx_mpi bash 09_scripts/04_trajectory_qc/01_preprocess_trajectories.sh"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Path resolution
# Script lives at <root>/09_scripts/04_trajectory_qc/01_preprocess_trajectories.sh
# SCRIPT_DIR  = <root>/09_scripts/04_trajectory_qc
# ../..        = <root>
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MD_ROOT="${PROJECT_ROOT}/02_md_inputs"

if [[ ! -d "${MD_ROOT}" ]]; then
    echo "[ERROR] MD inputs directory not found: ${MD_ROOT}"
    echo "        This directory is not tracked in git."
    echo "        Retrieve from the Zenodo data deposit or TRUBA scratch."
    exit 1
fi

# ---------------------------------------------------------------------------
# Systems and runs
# ---------------------------------------------------------------------------
SYSTEMS=(
    S1_BIN1_BAR
    S2_BIN1_SH3
    S3_PICALM_ANTH
    "S4_CD2AP_SH3-2"
    "S5_CD2AP_SH3-1"
)
RUNS=(run1 run2 run3)

# Group indices (standard CHARMM36m / GROMACS auto-generated):
#   0 = System   1 = Protein
GRP_SYSTEM=0
GRP_PROTEIN=1

# ---------------------------------------------------------------------------
# Helper: log with timestamp
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
TOTAL=0
SKIPPED=0
FAILED=0
MANIFEST_ENTRIES=()

for SYS in "${SYSTEMS[@]}"; do
    for RUN in "${RUNS[@]}"; do
        RUNDIR="${MD_ROOT}/${SYS}/${RUN}"

        if [[ ! -d "${RUNDIR}" ]]; then
            log "WARN  Directory not found, skipping: ${RUNDIR}"
            continue
        fi

        TPR="${RUNDIR}/prod.tpr"
        XTC="${RUNDIR}/prod.xtc"
        NOJUMP="${RUNDIR}/prod_nojump.xtc"
        CENTER="${RUNDIR}/prod_center.xtc"
        PROTEIN="${RUNDIR}/prod_protein.xtc"

        if [[ ! -f "${TPR}" || ! -f "${XTC}" ]]; then
            log "WARN  Missing prod.tpr or prod.xtc in ${SYS}/${RUN}, skipping"
            ((FAILED++)) || true
            continue
        fi

        log "INFO  ===== ${SYS} / ${RUN} ====="
        TOTAL=$((TOTAL + 1))

        # ------------------------------------------------------------------
        # Step 1 — Fix PBC jumps (nojump), output = full system
        # ------------------------------------------------------------------
        if [[ -f "${NOJUMP}" ]]; then
            log "SKIP  Step 1 already done: prod_nojump.xtc"
        else
            log "RUN   Step 1: pbc nojump → prod_nojump.xtc"
            echo "${GRP_SYSTEM}" | \
                "${GMX}" trjconv \
                    -s  "${TPR}" \
                    -f  "${XTC}" \
                    -o  "${NOJUMP}" \
                    -pbc nojump \
                    -quiet || {
                log "ERROR Step 1 failed for ${SYS}/${RUN}"
                ((FAILED++)) || true
                continue
            }
            log "DONE  prod_nojump.xtc"
        fi

        # ------------------------------------------------------------------
        # Step 2 — Center protein, wrap all molecules compactly
        #          center-group = Protein (1), output-group = System (0)
        # ------------------------------------------------------------------
        if [[ -f "${CENTER}" ]]; then
            log "SKIP  Step 2 already done: prod_center.xtc"
        else
            log "RUN   Step 2: center + pbc mol → prod_center.xtc"
            printf "%d\n%d\n" "${GRP_PROTEIN}" "${GRP_SYSTEM}" | \
                "${GMX}" trjconv \
                    -s  "${TPR}" \
                    -f  "${NOJUMP}" \
                    -o  "${CENTER}" \
                    -center \
                    -pbc mol \
                    -ur compact \
                    -quiet || {
                log "ERROR Step 2 failed for ${SYS}/${RUN}"
                ((FAILED++)) || true
                continue
            }
            log "DONE  prod_center.xtc"
        fi

        # ------------------------------------------------------------------
        # Step 3 — Extract protein-only frames from the centered trajectory
        #          (used for RMSD, RMSF, fpocket, DCCM)
        # ------------------------------------------------------------------
        if [[ -f "${PROTEIN}" ]]; then
            log "SKIP  Step 3 already done: prod_protein.xtc"
            ((SKIPPED++)) || true
        else
            log "RUN   Step 3: extract Protein → prod_protein.xtc"
            echo "${GRP_PROTEIN}" | \
                "${GMX}" trjconv \
                    -s  "${TPR}" \
                    -f  "${CENTER}" \
                    -o  "${PROTEIN}" \
                    -quiet || {
                log "ERROR Step 3 failed for ${SYS}/${RUN}"
                ((FAILED++)) || true
                continue
            }
            log "DONE  prod_protein.xtc"
        fi

        log "INFO  ${SYS}/${RUN} complete"
        MANIFEST_ENTRIES+=("${SYS}/${RUN}")
    done
done

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
MANIFEST="${MD_ROOT}/preprocess_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "gromacs_binary": "%s",\n' "$(command -v "${GMX}")"
    printf '  "gromacs_version": "%s",\n' "$("${GMX}" --version 2>&1 | awk '/GROMACS version/{print $NF}')"
    printf '  "grp_system": %d,\n' "${GRP_SYSTEM}"
    printf '  "grp_protein": %d,\n' "${GRP_PROTEIN}"
    printf '  "processed_runs": [\n'
    for i in "${!MANIFEST_ENTRIES[@]}"; do
        COMMA=","
        [[ $i -eq $((${#MANIFEST_ENTRIES[@]}-1)) ]] && COMMA=""
        printf '    "%s"%s\n' "${MANIFEST_ENTRIES[$i]}" "$COMMA"
    done
    printf '  ],\n'
    printf '  "total": %d,\n' "${TOTAL}"
    printf '  "failed": %d\n' "${FAILED}"
    printf '}\n'
} > "${MANIFEST}"
log "INFO  Manifest written → ${MANIFEST}"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
log "===== PREPROCESSING SUMMARY ====="
log "  Systems processed : ${TOTAL}"
log "  Skipped (done)    : ${SKIPPED}"
log "  Failed            : ${FAILED}"
if [[ "${FAILED}" -eq 0 ]]; then
    log "  Status: ALL OK"
else
    log "  Status: CHECK FAILURES ABOVE"
fi
