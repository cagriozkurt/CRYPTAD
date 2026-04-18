#!/usr/bin/env bash
# =============================================================================
# CRYPTAD — Fix PBC dimer artifact for S1 BIN1 BAR (3 replicas)
#
# Problem: pbc mol treated the two BAR domain chains independently,
# placing them in different periodic images → 20–37 Å RMSD artifact.
#
# Fix: insert a 'pbc cluster' step between nojump and centering.
# pbc cluster gathers all molecules around the largest one (protein),
# keeping both chains in the same periodic image before centering.
#
# Steps per run:
#   Step 1: prod_nojump.xtc already exists — KEEP (re-use from 01_preprocess)
#   Step 2: prod_nojump.xtc  → prod_cluster.xtc   (pbc cluster, center=Protein)
#   Step 3: prod_cluster.xtc → prod_center.xtc    (center + pbc mol + ur compact)
#   Step 4: prod_center.xtc  → prod_protein.xtc   (extract Protein)
#
# Reproducibility notes
# ---------------------
# - MD_ROOT resolved as two directories above SCRIPT_DIR (script is two levels
#   inside 09_scripts/), not one.
# - No grep filter downstream of trjconv — grep exiting 1 (no matches) under
#   set -o pipefail would falsely trigger error handlers on a successful run.
# - Steps 3 and 4 are skip-if-exists (idempotent), consistent with step 2.
# - A manifest is written on completion.
#
# Usage:
#   GMX=gmx_mpi bash 09_scripts/04_trajectory_qc/02_fix_pbc_dimer_s1.sh
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Load GROMACS module if needed
# ---------------------------------------------------------------------------
GMX="${GMX:-gmx}"

if ! command -v "${GMX}" &>/dev/null; then
    if command -v module &>/dev/null; then
        echo "[INFO] '${GMX}' not found — loading apps/gromacs/2024.1-oneapi2024"
        # shellcheck disable=SC1091
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load apps/gromacs/2024.1-oneapi2024
    fi
    if ! command -v "${GMX}" &>/dev/null; then
        echo "[ERROR] '${GMX}' not found. Load the GROMACS module first."
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Path resolution
# Script lives at <root>/09_scripts/04_trajectory_qc/02_fix_pbc_dimer_s1.sh
# SCRIPT_DIR  = <root>/09_scripts/04_trajectory_qc
# ../..        = <root>
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MD_ROOT="${PROJECT_ROOT}/02_md_inputs/S1_BIN1_BAR"

if [[ ! -d "${MD_ROOT}" ]]; then
    echo "[ERROR] S1_BIN1_BAR inputs directory not found: ${MD_ROOT}"
    echo "        Run 01_preprocess_trajectories.sh first."
    exit 1
fi

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
RUNS=(run1 run2 run3)
GRP_PROTEIN=1
GRP_SYSTEM=0

log() { echo "[$(date '+%H:%M:%S')] $*"; }

FAILED=0
COMPLETED=()

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for RUN in "${RUNS[@]}"; do
    RUNDIR="${MD_ROOT}/${RUN}"
    TPR="${RUNDIR}/prod.tpr"
    NOJUMP="${RUNDIR}/prod_nojump.xtc"
    CLUSTER="${RUNDIR}/prod_cluster.xtc"
    CENTER="${RUNDIR}/prod_center.xtc"
    PROTEIN="${RUNDIR}/prod_protein.xtc"

    log "INFO  ===== S1_BIN1_BAR / ${RUN} ====="

    if [[ ! -f "${NOJUMP}" ]]; then
        log "ERROR prod_nojump.xtc missing — run 01_preprocess_trajectories.sh first"
        ((FAILED++)) || true
        continue
    fi

    # Remove old (wrong) center and protein trajectories so they are rebuilt
    [[ -f "${CENTER}"  ]] && { log "RM    old prod_center.xtc";  rm "${CENTER}";  }
    [[ -f "${PROTEIN}" ]] && { log "RM    old prod_protein.xtc"; rm "${PROTEIN}"; }

    # ------------------------------------------------------------------
    # Step 2 — Cluster dimer chains together (pbc cluster)
    #          center-group = Protein, output = System
    # ------------------------------------------------------------------
    if [[ -f "${CLUSTER}" ]]; then
        log "SKIP  Step 2 already done: prod_cluster.xtc"
    else
        log "RUN   Step 2: pbc cluster → prod_cluster.xtc"
        printf "%d\n%d\n" "${GRP_PROTEIN}" "${GRP_SYSTEM}" | \
            "${GMX}" trjconv \
                -s  "${TPR}" \
                -f  "${NOJUMP}" \
                -o  "${CLUSTER}" \
                -pbc cluster \
                -quiet || {
            log "ERROR Step 2 failed for ${RUN}"
            ((FAILED++)) || true
            continue
        }
        log "DONE  prod_cluster.xtc"
    fi

    # ------------------------------------------------------------------
    # Step 3 — Center protein, wrap molecules compactly
    # ------------------------------------------------------------------
    if [[ -f "${CENTER}" ]]; then
        log "SKIP  Step 3 already done: prod_center.xtc"
    else
        log "RUN   Step 3: center + pbc mol → prod_center.xtc"
        printf "%d\n%d\n" "${GRP_PROTEIN}" "${GRP_SYSTEM}" | \
            "${GMX}" trjconv \
                -s  "${TPR}" \
                -f  "${CLUSTER}" \
                -o  "${CENTER}" \
                -center \
                -pbc mol \
                -ur compact \
                -quiet || {
            log "ERROR Step 3 failed for ${RUN}"
            ((FAILED++)) || true
            continue
        }
        log "DONE  prod_center.xtc"
    fi

    # ------------------------------------------------------------------
    # Step 4 — Extract protein-only
    # ------------------------------------------------------------------
    if [[ -f "${PROTEIN}" ]]; then
        log "SKIP  Step 4 already done: prod_protein.xtc"
    else
        log "RUN   Step 4: extract Protein → prod_protein.xtc"
        printf "%d\n" "${GRP_PROTEIN}" | \
            "${GMX}" trjconv \
                -s  "${TPR}" \
                -f  "${CENTER}" \
                -o  "${PROTEIN}" \
                -quiet || {
            log "ERROR Step 4 failed for ${RUN}"
            ((FAILED++)) || true
            continue
        }
        log "DONE  prod_protein.xtc"
    fi

    log "INFO  ${RUN} complete"
    COMPLETED+=("S1_BIN1_BAR/${RUN}")
done

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
MANIFEST="${PROJECT_ROOT}/02_md_inputs/S1_BIN1_BAR/fix_pbc_dimer_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "gromacs_binary": "%s",\n' "$(command -v "${GMX}")"
    printf '  "gromacs_version": "%s",\n' "$("${GMX}" --version 2>&1 | awk '/GROMACS version/{print $NF}')"
    printf '  "fix": "pbc cluster inserted between nojump and centering for homodimer",\n'
    printf '  "grp_protein": %d,\n' "${GRP_PROTEIN}"
    printf '  "grp_system": %d,\n' "${GRP_SYSTEM}"
    printf '  "completed_runs": [\n'
    for i in "${!COMPLETED[@]}"; do
        COMMA=","
        [[ $i -eq $((${#COMPLETED[@]}-1)) ]] && COMMA=""
        printf '    "%s"%s\n' "${COMPLETED[$i]}" "$COMMA"
    done
    printf '  ],\n'
    printf '  "failed": %d\n' "${FAILED}"
    printf '}\n'
} > "${MANIFEST}"
log "INFO  Manifest written → ${MANIFEST}"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
echo ""
log "===== SUMMARY ====="
if [[ "${FAILED}" -eq 0 ]]; then
    log "  Status: ALL OK"
    log "  Next: re-run 05_plot_qc.py for S1_BIN1_BAR only, then rsync back"
else
    log "  Status: ${FAILED} run(s) FAILED — check output above"
fi
