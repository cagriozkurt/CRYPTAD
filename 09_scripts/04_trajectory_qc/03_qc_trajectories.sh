#!/usr/bin/env bash
# =============================================================================
# CRYPTAD — Basic trajectory QC for all 15 MD runs
#
# Generates per-run .xvg files in 03_trajectory_qc/{SYS}/{RUN}/:
#   rmsd_backbone.xvg  — Backbone RMSD vs frame 0 (ns)
#   rg.xvg             — Protein radius of gyration (ns)
#   temperature.xvg    — Simulation temperature (K)
#   pressure.xvg       — Simulation pressure (bar)
#   potential.xvg      — Potential energy (kJ/mol)
#
# Reproducibility notes
# ---------------------
# - PROJECT_ROOT resolved as two directories above SCRIPT_DIR (script is
#   two levels inside 09_scripts/), not one.
# - No grep filter downstream of gmx calls — grep exiting 1 under
#   set -o pipefail falsely triggers error handlers on successful runs.
# - A manifest is written on completion.
#
# Usage:
#   bash 09_scripts/04_trajectory_qc/03_qc_trajectories.sh
#   GMX=gmx_mpi bash 09_scripts/04_trajectory_qc/03_qc_trajectories.sh
#
# Run after 01_preprocess_trajectories.sh — requires prod_center.xtc.
# Already-existing output files are skipped automatically.
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
# Script lives at <root>/09_scripts/04_trajectory_qc/03_qc_trajectories.sh
# SCRIPT_DIR   = <root>/09_scripts/04_trajectory_qc
# ../..         = <root>
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
MD_ROOT="${PROJECT_ROOT}/02_md_inputs"
QC_ROOT="${PROJECT_ROOT}/03_trajectory_qc"

SYSTEMS=(
    S1_BIN1_BAR
    S2_BIN1_SH3
    S3_PICALM_ANTH
    "S4_CD2AP_SH3-2"
    "S5_CD2AP_SH3-1"
)
RUNS=(run1 run2 run3)

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
TOTAL=0; FAILED=0
COMPLETED=()

for SYS in "${SYSTEMS[@]}"; do
    for RUN in "${RUNS[@]}"; do
        RUNDIR="${MD_ROOT}/${SYS}/${RUN}"
        OUTDIR="${QC_ROOT}/${SYS}/${RUN}"

        TPR="${RUNDIR}/prod.tpr"
        CENTER="${RUNDIR}/prod_center.xtc"
        EDR="${RUNDIR}/prod.edr"

        if [[ ! -f "${CENTER}" ]]; then
            log "WARN  ${SYS}/${RUN}: prod_center.xtc not found — run 01_preprocess_trajectories.sh first"
            ((FAILED++)) || true
            continue
        fi
        if [[ ! -f "${TPR}" || ! -f "${EDR}" ]]; then
            log "WARN  ${SYS}/${RUN}: missing prod.tpr or prod.edr"
            ((FAILED++)) || true
            continue
        fi

        mkdir -p "${OUTDIR}"
        log "INFO  ===== ${SYS} / ${RUN} ====="
        TOTAL=$((TOTAL + 1))

        # ------------------------------------------------------------------
        # RMSD — Backbone fit and calculation vs frame 0
        # ------------------------------------------------------------------
        RMSD_OUT="${OUTDIR}/rmsd_backbone.xvg"
        if [[ -f "${RMSD_OUT}" ]]; then
            log "SKIP  rmsd_backbone.xvg"
        else
            log "RUN   RMSD (Backbone vs frame 0)"
            printf "Backbone\nBackbone\n" | \
                "${GMX}" rms \
                    -s "${TPR}" \
                    -f "${CENTER}" \
                    -o "${RMSD_OUT}" \
                    -tu ns \
                    -quiet || {
                log "ERROR RMSD failed for ${SYS}/${RUN}"
                ((FAILED++)) || true
            }
        fi

        # ------------------------------------------------------------------
        # Radius of gyration — Protein
        # ------------------------------------------------------------------
        RG_OUT="${OUTDIR}/rg.xvg"
        if [[ -f "${RG_OUT}" ]]; then
            log "SKIP  rg.xvg"
        else
            log "RUN   Radius of gyration (Protein)"
            printf "Protein\n" | \
                "${GMX}" gyrate \
                    -s "${TPR}" \
                    -f "${CENTER}" \
                    -o "${RG_OUT}" \
                    -tu ns \
                    -quiet || {
                log "ERROR Rg failed for ${SYS}/${RUN}"
                ((FAILED++)) || true
            }
        fi

        # ------------------------------------------------------------------
        # Energy terms from .edr (Temperature, Pressure, Potential)
        # ------------------------------------------------------------------
        for TERM in Temperature Pressure Potential; do
            TERM_OUT="${OUTDIR}/${TERM,,}.xvg"
            if [[ -f "${TERM_OUT}" ]]; then
                log "SKIP  ${TERM,,}.xvg"
            else
                log "RUN   ${TERM}"
                printf "%s\n\n" "${TERM}" | \
                    "${GMX}" energy \
                        -f "${EDR}" \
                        -o "${TERM_OUT}" \
                        -quiet || {
                    log "ERROR ${TERM} extraction failed for ${SYS}/${RUN}"
                    ((FAILED++)) || true
                }
            fi
        done

        log "INFO  ${SYS}/${RUN} done"
        COMPLETED+=("${SYS}/${RUN}")
    done
done

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
mkdir -p "${QC_ROOT}"
MANIFEST="${QC_ROOT}/qc_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "gromacs_binary": "%s",\n' "$(command -v "${GMX}")"
    printf '  "gromacs_version": "%s",\n' "$("${GMX}" --version 2>&1 | awk '/GROMACS version/{print $NF}')"
    printf '  "qc_root": "%s",\n' "${QC_ROOT}"
    printf '  "completed_runs": [\n'
    for i in "${!COMPLETED[@]}"; do
        COMMA=","
        [[ $i -eq $((${#COMPLETED[@]}-1)) ]] && COMMA=""
        printf '    "%s"%s\n' "${COMPLETED[$i]}" "$COMMA"
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
log "===== QC SUMMARY ====="
log "  Runs processed : ${TOTAL}"
log "  Failures       : ${FAILED}"
if [[ "${FAILED}" -eq 0 ]]; then
    log "  Status: ALL OK"
    echo ""
    log "  XVG files written to: ${QC_ROOT}"
    log "  Rsync back with:"
    log "    rsync -avz truba:PATH/03_trajectory_qc/ ./03_trajectory_qc/"
    log "  Then plot locally:"
    log "    python 09_scripts/04_trajectory_qc/05_plot_qc.py"
else
    log "  Status: CHECK FAILURES ABOVE"
fi
