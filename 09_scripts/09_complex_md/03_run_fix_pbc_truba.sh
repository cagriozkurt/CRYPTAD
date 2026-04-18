#!/usr/bin/env bash
# ============================================================
# CRYPTAD — PBC fix (MDAnalysis make_whole) on TRUBA
# Fixes complex MD prod.xtc → complexmd_pbc.xtc
#
# Usage:
#   # First generate the fix pairs list (run once):
#   bash 09_scripts/09_complex_md/03_run_fix_pbc_truba.sh --make-pairs
#
#   # Then submit:
#   N=$(wc -l < "$PROJ/02_md_simulations/complex_md/fix_pbc_pairs.txt")
#   sbatch --array=1-${N} 09_scripts/09_complex_md/03_run_fix_pbc_truba.sh
# ============================================================
#SBATCH --job-name=cryptad_pbc_fix
#SBATCH --account=mozkurt
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-02:00:00
#SBATCH --output=pbc_fix_%A_%a.out
#SBATCH --error=pbc_fix_%A_%a.err

set -uo pipefail

# ── Environment ───────────────────────────────────────────────
# AmberTools24 conda Python 3.12; MDAnalysis installed via pip --user
CONDA_PYTHON="/arf/sw/comp/python/miniconda3/envs/AmberTools24/bin/python3"
export PATH="$HOME/.local/bin:$PATH"
export PYTHONPATH="$HOME/.local/lib/python3.12/site-packages:$PYTHONPATH"

CRYPTAD="${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}"
COMPLEXMD_BASE="$CRYPTAD/02_md_simulations/complex_md"
PAIRS="$COMPLEXMD_BASE/fix_pbc_pairs.txt"
SCRIPT="$CRYPTAD/09_scripts/09_complex_md/04_fix_pbc_mdanalysis.py"

# ── Generate pairs file (run with --make-pairs, not via SLURM) ─
if [[ "${1:-}" == "--make-pairs" ]]; then
    python3 - <<PYEOF
import os
base = os.environ.get("CRYPTAD", "/arf/scratch/mozkurt/CRYPTAD")
complexmd_base = os.path.join(base, "02_md_simulations", "complex_md")
pairs_file = os.path.join(complexmd_base, "fix_pbc_pairs.txt")
pairs = []
for sys_dir in sorted(os.listdir(complexmd_base)):
    sys_path = os.path.join(complexmd_base, sys_dir)
    if not os.path.isdir(sys_path):
        continue
    for rep in sorted(os.listdir(sys_path)):
        rep_path = os.path.join(sys_path, rep)
        if not os.path.isdir(rep_path):
            continue
        xtc = os.path.join(rep_path, "complexmd.xtc")
        if not os.path.exists(xtc):
            continue
        pairs.append(f"{sys_dir} {rep}")
os.makedirs(complexmd_base, exist_ok=True)
with open(pairs_file, "w") as f:
    f.write("\n".join(pairs) + "\n")
print(f"Wrote {len(pairs)} entries to {pairs_file}")
PYEOF
    exit 0
fi

# ── Parse system + replica from pairs file ────────────────────
LINE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" "$PAIRS")
SYS=$(echo "$LINE" | awk '{print $1}')
REP=$(echo "$LINE" | awk '{print $2}')

if [[ -z "$SYS" || -z "$REP" ]]; then
    echo "ERROR: empty line at index $SLURM_ARRAY_TASK_ID" >&2; exit 1
fi

WORKDIR="$COMPLEXMD_BASE/$SYS/$REP"
TPR="$WORKDIR/complexmd.tpr"
XTC="$WORKDIR/complexmd.xtc"
OUT="$WORKDIR/complexmd_pbc.xtc"

echo "=== Task $SLURM_ARRAY_TASK_ID: $SYS / $REP ==="

# Remove bad existing file (same size as complexmd.xtc = old trjconv copy)
if [[ -f "$OUT" ]]; then
    sz_xtc=$(stat -c%s "$XTC" 2>/dev/null || stat -f%z "$XTC")
    sz_out=$(stat -c%s "$OUT" 2>/dev/null || stat -f%z "$OUT")
    if [[ "$sz_out" -eq "$sz_xtc" ]]; then
        echo "  Removing bad complexmd_pbc.xtc (identical to complexmd.xtc)"
        rm -f "$OUT"
    fi
fi

# Run PBC fix
$CONDA_PYTHON "$SCRIPT" "$TPR" "$XTC" "$OUT"

echo "=== DONE: $SYS / $REP ==="
