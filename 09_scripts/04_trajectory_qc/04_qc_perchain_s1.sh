#!/usr/bin/env bash
# =============================================================================
# CRYPTAD — Per-chain backbone RMSD for S1 BIN1 BAR dimer
#
# CHARMM-GUI puts protein atoms first in the topology.
# From make_ndx: 6402 protein atoms total → 3201 per chain.
#   Chain A atoms: 1    – 3201  (system atom indices, 1-based)
#   Chain B atoms: 3202 – 6402
#
# make_ndx group numbering after default 16 groups (0–15):
#   a 1-3201     → group 16   (chain A all atoms)  ← ATOM indices, always safe
#   16 & 4       → group 17   (chain A backbone)   ← Backbone = group 4
#   a 3202-6402  → group 18   (chain B all atoms)
#   18 & 4       → group 19   (chain B backbone)
#
# Reproducibility notes
# ---------------------
# - set -euo pipefail ensures any command or pipe failure aborts the script.
# - PROJECT_ROOT resolved as two directories above SCRIPT_DIR.
# - No tail/grep filters downstream of gmx — avoids SIGPIPE and lossy output.
# - NDX and RMSD outputs are skipped if already present (idempotent).
# - A manifest is written on completion.
#
# Usage:
#   GMX=gmx_mpi bash 09_scripts/04_trajectory_qc/04_qc_perchain_s1.sh
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
# Script lives at <root>/09_scripts/04_trajectory_qc/04_qc_perchain_s1.sh
# SCRIPT_DIR   = <root>/09_scripts/04_trajectory_qc
# ../..         = <root>
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
S1_MD="${PROJECT_ROOT}/02_md_inputs/S1_BIN1_BAR"
S1_QC="${PROJECT_ROOT}/03_trajectory_qc/S1_BIN1_BAR"

# ---------------------------------------------------------------------------
# Atom index boundaries (1-based, from make_ndx protein atom count: 6402 total)
# ---------------------------------------------------------------------------
CHAIN_A_LAST=3201
CHAIN_B_FIRST=3202
CHAIN_B_LAST=6402

# NDX group numbers (0–15 = defaults; additions start at 16)
GRP_A_ALL=16    # a 1-3201
GRP_A_BB=17     # 16 & 4  (chain A backbone)
GRP_B_ALL=18    # a 3202-6402
GRP_B_BB=19     # 18 & 4  (chain B backbone)

log() { echo "[$(date '+%H:%M:%S')] $*"; }

FAILED=0
COMPLETED=()

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
for RUN in run1 run2 run3; do
    RUNDIR="${S1_MD}/${RUN}"
    QCDIR="${S1_QC}/${RUN}"
    TPR="${RUNDIR}/prod.tpr"
    CENTER="${RUNDIR}/prod_center.xtc"
    NDX="${RUNDIR}/chain_backbone.ndx"
    OUTFILE_A="${QCDIR}/rmsd_chain_A.xvg"
    OUTFILE_B="${QCDIR}/rmsd_chain_B.xvg"

    log "===== ${RUN} ====="

    if [[ ! -f "${CENTER}" ]]; then
        log "ERROR prod_center.xtc missing — run 01_preprocess_trajectories.sh and 02_fix_pbc_dimer_s1.sh first"
        ((FAILED++)) || true
        continue
    fi

    mkdir -p "${QCDIR}"

    # ------------------------------------------------------------------
    # Build NDX with per-chain backbone groups (skip if already exists)
    # ------------------------------------------------------------------
    if [[ -f "${NDX}" ]]; then
        log "SKIP  chain_backbone.ndx already exists"
    else
        log "RUN   make_ndx → per-chain backbone groups"
        printf \
            "a 1-%d\n%d & 4\na %d-%d\n%d & 4\nq\n" \
            "${CHAIN_A_LAST}" \
            "${GRP_A_ALL}" \
            "${CHAIN_B_FIRST}" "${CHAIN_B_LAST}" \
            "${GRP_B_ALL}" | \
            "${GMX}" make_ndx -f "${TPR}" -o "${NDX}" -quiet || {
            log "ERROR make_ndx failed for ${RUN} — skipping"
            ((FAILED++)) || true
            continue
        }
        if [[ ! -f "${NDX}" ]]; then
            log "ERROR make_ndx produced no output for ${RUN} — skipping"
            ((FAILED++)) || true
            continue
        fi
    fi

    # Log group names without a multi-command pipe (awk, single process)
    log "INFO  Groups in NDX:"
    awk '/^\[/{printf "          %d: %s\n", ++n, $0}' "${NDX}"

    # ------------------------------------------------------------------
    # Chain A backbone RMSD
    # ------------------------------------------------------------------
    if [[ -f "${OUTFILE_A}" ]]; then
        log "SKIP  rmsd_chain_A.xvg"
    else
        log "RUN   gmx rms chain A (group ${GRP_A_BB})"
        printf "%d\n%d\n" "${GRP_A_BB}" "${GRP_A_BB}" | \
            "${GMX}" rms \
                -s  "${TPR}" \
                -f  "${CENTER}" \
                -n  "${NDX}" \
                -o  "${OUTFILE_A}" \
                -tu ns \
                -quiet || {
            log "ERROR chain A rms failed for ${RUN}"
            ((FAILED++)) || true
            continue
        }
        log "DONE  rmsd_chain_A.xvg"
    fi

    # ------------------------------------------------------------------
    # Chain B backbone RMSD
    # ------------------------------------------------------------------
    if [[ -f "${OUTFILE_B}" ]]; then
        log "SKIP  rmsd_chain_B.xvg"
    else
        log "RUN   gmx rms chain B (group ${GRP_B_BB})"
        printf "%d\n%d\n" "${GRP_B_BB}" "${GRP_B_BB}" | \
            "${GMX}" rms \
                -s  "${TPR}" \
                -f  "${CENTER}" \
                -n  "${NDX}" \
                -o  "${OUTFILE_B}" \
                -tu ns \
                -quiet || {
            log "ERROR chain B rms failed for ${RUN}"
            ((FAILED++)) || true
            continue
        }
        log "DONE  rmsd_chain_B.xvg"
    fi

    log "INFO  ${RUN} done"
    COMPLETED+=("S1_BIN1_BAR/${RUN}")
    echo ""
done

# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
mkdir -p "${S1_QC}"
MANIFEST="${S1_QC}/qc_perchain_manifest.json"
{
    printf '{\n'
    printf '  "generated_at": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf '  "gromacs_binary": "%s",\n' "$(command -v "${GMX}")"
    printf '  "gromacs_version": "%s",\n' "$("${GMX}" --version 2>&1 | awk '/GROMACS version/{print $NF}')"
    printf '  "chain_a_atoms": "1-%d",\n' "${CHAIN_A_LAST}"
    printf '  "chain_b_atoms": "%d-%d",\n' "${CHAIN_B_FIRST}" "${CHAIN_B_LAST}"
    printf '  "grp_a_backbone": %d,\n' "${GRP_A_BB}"
    printf '  "grp_b_backbone": %d,\n' "${GRP_B_BB}"
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
log "===== ALL DONE ====="
if [[ "${FAILED}" -eq 0 ]]; then
    log "  Status: ALL OK"
    log "  Rsync: rsync -avz truba:/arf/scratch/mozkurt/CRYPTAD/03_trajectory_qc/S1_BIN1_BAR/ ./03_trajectory_qc/S1_BIN1_BAR/"
    log "  Then:  python 09_scripts/04_trajectory_qc/06_plot_qc_perchain.py"
else
    log "  Status: ${FAILED} run(s) FAILED — check output above"
fi
