#!/bin/bash
# 02_run_sum_hills.sh  —  CRYPTAD
# FES block convergence analysis using plumed sum_hills
# Run on TRUBA (plumed available via load_gmx_plumed.sh)
#
# For each system: compute FES at t=67, 133, 200 ns (thirds of simulation)
# Output: fes_67ns.dat, fes_133ns.dat, fes_200ns.dat in each run directory
#
# Fix notes:
#   - PLUMED 2.9 sum_hills has no --endhills flag; use awk truncation instead
#   - Gaussian count thresholds are derived from the actual deposition interval
#     read from each HILLS file (col 1 time difference between first two data
#     lines) — no hardcoded magic numbers that would silently misfire if
#     meta.mdp PACE or dt were ever changed
#   - Grid bounds auto-detected from HILLS CV columns (with 10% padding)
#   - All counts use awk (no grep|wc pipes) to avoid SIGPIPE under set -o pipefail
#
# Usage (submit as interactive or short batch job):
#   source "$HOME/load_gmx_plumed.sh"
#   bash 09_scripts/03_metadynamics/02_run_sum_hills.sh
#
# After completion: rsync to local and run:
#   python3 09_scripts/07_conformational/01_analyze_fes_convergence.py

set -euo pipefail
source "$HOME/load_gmx_plumed.sh"

CRYPTAD=/arf/scratch/mozkurt/CRYPTAD

NBINS=100
NBINS2D="${NBINS},${NBINS}"  # plumed sum_hills needs vector for 2D: "100,100"

# System → run directory mapping (colon-split; bash 3.2-safe pattern not needed
# here since this runs on TRUBA bash 4+, but kept explicit for clarity)
SYSTEMS=(
    "S1_BIN1_BAR:$CRYPTAD/02_md_simulations/BIN1/metadynamics/S1_BAR"
    "S2_BIN1_SH3:$CRYPTAD/02_md_simulations/BIN1/metadynamics/S2_SH3"
    "S3_PICALM_ANTH:$CRYPTAD/02_md_simulations/PICALM/metadynamics/S3_ANTH"
    "S4_CD2AP_SH3-2:$CRYPTAD/02_md_simulations/CD2AP/metadynamics/S4_SH3-2"
    "S5_CD2AP_SH3-1:$CRYPTAD/02_md_simulations/CD2AP/metadynamics/S5_SH3-1"
)

# ── Helper: detect Gaussian deposition interval in ps from HILLS file ─────────
# Reads the time difference between the first two data lines (col 1).
# This is authoritative: it reflects the actual PACE × dt used by PLUMED,
# so truncation counts are always correct regardless of meta.mdp settings.
detect_deposition_interval_ps() {
    local hills_file="$1"
    awk '!/^#/ { if (prev!="") { printf "%.6f\n", $1-prev; exit } prev=$1 }' "$hills_file"
}

# ── Helper: compute grid bounds (min,max for cv1 and cv2) from HILLS file ─────
# Reads cv1 (col 2) and cv2 (col 3) data lines, adds 10% padding.
# Output: four space-separated values: MIN_CV1 MAX_CV1 MIN_CV2 MAX_CV2
get_bounds() {
    local hills_file="$1"
    awk '
    /^#/ { next }
    first == "" { min1=$2; max1=$2; min2=$3; max2=$3; first=1; next }
    {
        if ($2 < min1) min1=$2; if ($2 > max1) max1=$2
        if ($3 < min2) min2=$3; if ($3 > max2) max2=$3
    }
    END {
        pad1 = (max1-min1)*0.10 + 0.05
        pad2 = (max2-min2)*0.10 + 0.05
        printf "%.4f %.4f %.4f %.4f\n", min1-pad1, max1+pad1, min2-pad2, max2+pad2
    }
    ' "$hills_file"
}

# ── Helper: create a truncated HILLS file with header + first N data lines ────
# Uses awk (single process, no pipe) to avoid SIGPIPE with set -o pipefail.
truncate_hills() {
    local hills_file="$1"
    local n_data="$2"
    local outfile="$3"
    awk -v n="$n_data" '
        /^#/ { print; next }
        { if (count < n) { print; count++ } else { exit } }
    ' "$hills_file" > "$outfile"
}

for ENTRY in "${SYSTEMS[@]}"; do
    SYS="${ENTRY%%:*}"
    DIR="${ENTRY#*:}"
    HILLS="$DIR/HILLS"

    if [[ ! -f "$HILLS" ]]; then
        echo "[SKIP] $SYS: HILLS not found at $HILLS"
        continue
    fi

    echo ""
    echo "=== $SYS ==="

    # Count data lines without a pipe (awk reads all before exiting — no SIGPIPE)
    TOTAL_DATA=$(awk '!/^#/{count++} END{print count+0}' "$HILLS")
    TOTAL_LINES=$(awk 'END{print NR}' "$HILLS")
    echo "HILLS: $TOTAL_DATA data lines ($TOTAL_LINES total)"

    # ── Derive Gaussian counts from the actual deposition interval ─────────────
    DEPOSITION_PS=$(detect_deposition_interval_ps "$HILLS")
    echo "  Deposition interval: ${DEPOSITION_PS} ps"
    HILLS_67NS=$(awk  "BEGIN{printf \"%d\", int(67000  / $DEPOSITION_PS)}")
    HILLS_133NS=$(awk "BEGIN{printf \"%d\", int(133000 / $DEPOSITION_PS)}")
    echo "  Thresholds: 67 ns=$HILLS_67NS  133 ns=$HILLS_133NS  200 ns=all"

    cd "$DIR"

    # ── Compute grid bounds from full HILLS ────────────────────────────────────
    echo "  Computing grid bounds..."
    read MIN_CV1 MAX_CV1 MIN_CV2 MAX_CV2 < <(get_bounds "$HILLS")
    echo "  Grid: CV1=[${MIN_CV1}, ${MAX_CV1}]  CV2=[${MIN_CV2}, ${MAX_CV2}]"

    MINARG="${MIN_CV1},${MIN_CV2}"
    MAXARG="${MAX_CV1},${MAX_CV2}"

    # ── FES at t=67 ns (first third) ──────────────────────────────────────────
    echo "  sum_hills → fes_67ns.dat  (first ${HILLS_67NS} Gaussians)"
    truncate_hills "$HILLS" "$HILLS_67NS" HILLS_67ns.tmp
    plumed sum_hills \
        --hills HILLS_67ns.tmp \
        --mintozero \
        --bin "$NBINS2D" \
        --min "$MINARG" \
        --max "$MAXARG" \
        --outfile fes_67ns.dat
    rm -f HILLS_67ns.tmp

    # ── FES at t=133 ns (first two thirds) ────────────────────────────────────
    echo "  sum_hills → fes_133ns.dat  (first ${HILLS_133NS} Gaussians)"
    truncate_hills "$HILLS" "$HILLS_133NS" HILLS_133ns.tmp
    plumed sum_hills \
        --hills HILLS_133ns.tmp \
        --mintozero \
        --bin "$NBINS2D" \
        --min "$MINARG" \
        --max "$MAXARG" \
        --outfile fes_133ns.dat
    rm -f HILLS_133ns.tmp

    # ── FES at t=200 ns (full run) ─────────────────────────────────────────────
    echo "  sum_hills → fes_200ns.dat  (full HILLS)"
    plumed sum_hills \
        --hills HILLS \
        --mintozero \
        --bin "$NBINS2D" \
        --min "$MINARG" \
        --max "$MAXARG" \
        --outfile fes_200ns.dat

    echo "  Done: $(awk 'END{printf "%s %d lines\n", FILENAME, NR}' fes_67ns.dat)  \
$(awk 'END{printf "%s %d lines\n", FILENAME, NR}' fes_133ns.dat)  \
$(awk 'END{printf "%s %d lines\n", FILENAME, NR}' fes_200ns.dat)"
done

echo ""
echo "=== sum_hills complete ==="
echo "Rsync to local and run:"
echo "  python3 09_scripts/07_conformational/01_analyze_fes_convergence.py"
