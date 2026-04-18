#!/bin/bash
#SBATCH --job-name=cryptad_mmgbsa
#SBATCH --account=mozkurt
#SBATCH --partition=hamsi
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=2
#SBATCH --time=3-00:00:00
#SBATCH --array=1-99
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

set -uo pipefail

# ── Environment ──────────────────────────────────────────────────────────────
module load apps/gromacs/2024.1-oneapi2024
# packages installed to ~/.local via pip (user install — no venv on TRUBA)
export PATH="$HOME/.local/bin:$PATH"
# AmberTools24 uses Python 3.12; MDAnalysis installed via pip --user
export PYTHONPATH="$HOME/.local/lib/python3.12/site-packages:${PYTHONPATH:-}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-2}
export GMX="gmx_mpi"
export MPI="mpirun -np ${SLURM_NTASKS:-28}"
# AmberTools24 conda env provides tleap, acpype, and MDAnalysis
export AMBERHOME="/arf/sw/comp/python/miniconda3/envs/AmberTools24"
export PATH="$AMBERHOME/bin:$PATH"
CONDA_PYTHON="$AMBERHOME/bin/python3"

# ── Paths ────────────────────────────────────────────────────────────────────
PROJ="${CRYPTAD:-/arf/scratch/mozkurt/CRYPTAD}"
MMGBSADIR="$PROJ/04_virtual_screening/mmgbsa"
MDPDIR="$MMGBSADIR/mdp"
PAIRS="$MMGBSADIR/pairs.txt"

# ── Select pair for this task ─────────────────────────────────────────────────
# OFFSET allows reuse across batches: batch2 --export=ALL,OFFSET=99; batch3 --export=ALL,OFFSET=198
OFFSET=${OFFSET:-0}
LINE=$(awk "NR==$((SLURM_ARRAY_TASK_ID + OFFSET))" "$PAIRS")
POCKET=$(echo "$LINE" | awk '{print $1}')
LIG_ID=$(echo "$LINE" | awk '{print $2}')
WORKDIR="$MMGBSADIR/$POCKET/$LIG_ID"
cd "$WORKDIR" || { echo "ERROR: cannot cd to $WORKDIR"; exit 1; }

echo "=== Task $SLURM_ARRAY_TASK_ID: $POCKET / $LIG_ID ==="

# ── 1. GAFF2 ligand topology (acpype) ────────────────────────────────────────
if [[ ! -f "LIG_GMX.itp" ]]; then
    echo "--- acpype ---"
    # acpype names output files after the molecule name in the SDF (line 1)
    # Rename to LIG so output files are LIG_GMX.gro / LIG_GMX.itp as expected
    sed -i "1s/.*/LIG/" ligand_pose.sdf
    acpype -i ligand_pose.sdf -c bcc -a gaff2 -n 0 -o gmx 2>&1
    # acpype creates a directory like ligand_pose.acpype/
    ACPYPE_DIRS=( *.acpype )
    ACPYPE_DIR="${ACPYPE_DIRS[0]:-}"
    if [[ -z "$ACPYPE_DIR" || ! -d "$ACPYPE_DIR" ]]; then
        echo "ERROR: acpype failed for $LIG_ID" >&2; exit 1
    fi
    # acpype names files after input filename (ligand_pose), not molecule name
    cp "$ACPYPE_DIR/ligand_pose_GMX.gro"  LIG_GMX.gro
    cp "$ACPYPE_DIR/ligand_pose_GMX.itp"  LIG_GMX.itp
    cp "$ACPYPE_DIR/ligand_pose_GMX.top"  LIG_GMX.top
    cp "$ACPYPE_DIR/posre_ligand_pose.itp" posre_LIG.itp 2>/dev/null || true
    if [[ ! -f "LIG_GMX.itp" ]]; then
        echo "SKIP: acpype produced no GMX files for $LIG_ID (likely charged molecule — check sqm.out)" >&2
        exit 0
    fi
    # acpype sets moleculetype = input file basename ("ligand_pose"); rename to LIG
    # so topol.top "[ molecules ]" entry "LIG  1" matches
    sed -i 's/^ ligand_pose\b/LIG/' LIG_GMX.itp
    # Also fix residue name in GRO (acpype uses MOL) so make_ndx produces a LIG group
    sed -i 's/MOL/LIG/g' LIG_GMX.gro
fi

# ── 2. Protein topology (AMBER99SB-ILDN) ─────────────────────────────────────
if [[ ! -f "topol.top" ]]; then
    echo "--- pdb2gmx ---"
    # Force consistent His protonation across all receptor frames (site1813 reference pattern):
    # Chain A: HIS86=HIE(1), HIS125=HIE(1), HIS165=HID(0), HIS166=HID(0), HIS232=HIE(1)
    # Chain B: HIS86=HIE(1), HIS125=HID(0), HIS165=HIE(1), HIS166=HIE(1), HIS232=HIE(1)
    printf "1\n1\n0\n0\n1\n1\n0\n1\n1\n1\n" | $GMX pdb2gmx \
        -f receptor_H.pdb -o protein.gro -p topol.top \
        -i posre.itp -ff amber99sb-ildn -water tip3p \
        -ignh -his -nobackup

    # Extract [ atomtypes ] from LIG_GMX.itp into separate file.
    # GROMACS requires all [ atomtypes ] before any [ moleculetype ] directive.
    python3 - << 'PYEOF'
import re
with open('LIG_GMX.itp') as f:
    itp = f.read()
m = re.search(r'(\[ atomtypes \][^\[]+)', itp, re.DOTALL)
if m:
    with open('lig_atomtypes.itp', 'w') as f:
        f.write(m.group(1).rstrip('\n') + '\n')
    with open('LIG_GMX.itp', 'w') as f:
        f.write(itp[:m.start()] + itp[m.end():])
    print("atomtypes extracted OK")
else:
    print("WARNING: no [ atomtypes ] found in LIG_GMX.itp")
PYEOF
    # Include GAFF2 atomtypes right after FF include (before protein moleculetype)
    sed -i '/^#include.*forcefield\.itp/a #include "lig_atomtypes.itp"' topol.top
    # Include LIG_GMX.itp (without atomtypes) and position restraints before water
    sed -i '/; Include water topology/i #include "LIG_GMX.itp"\n#ifdef POSRES\n#include "posre_LIG.itp"\n#endif' topol.top
    # Add LIG to [ molecules ] section
    echo "LIG  1" >> topol.top
fi

# ── 3. Build complex GRO ─────────────────────────────────────────────────────
if [[ ! -f "complex.gro" ]]; then
    echo "--- building complex GRO ---"
    N_PROT=$(awk 'NR==2{print; exit}' protein.gro)
    N_LIG=$(awk 'NR==2{print; exit}' LIG_GMX.gro)
    N_TOTAL=$(( N_PROT + N_LIG ))
    BOX=$(tail -1 protein.gro)
    {
        echo "Complex"
        echo "  $N_TOTAL"
        head -n $(( N_PROT + 2 )) protein.gro | tail -n $N_PROT
        head -n $(( N_LIG + 2 )) LIG_GMX.gro  | tail -n $N_LIG
        echo "$BOX"
    } > complex.gro
fi

# ── 4. Solvate + ions ────────────────────────────────────────────────────────
if [[ ! -f "solvated.gro" ]]; then
    echo "--- solvate ---"
    $GMX editconf -f complex.gro -o complex_box.gro -bt dodecahedron -d 1.2 -nobackup
    $GMX solvate  -cp complex_box.gro -cs spc216.gro -o solvated.gro -p topol.top -nobackup
    $GMX grompp   -f "$MDPDIR/minim.mdp" -c solvated.gro -p topol.top \
                  -o ions.tpr -maxwarn 1 -nobackup
    echo "SOL" | $GMX genion -s ions.tpr -o solvated_ions.gro -p topol.top \
                  -pname NA -nname CL -neutral -conc 0.15 -nobackup
    mv solvated_ions.gro solvated.gro
    # GROMACS rewrites residue names during solvate/genion; rename MOL back to LIG
    # so make_ndx auto-generates a LIG group (required for gmx_MMPBSA -cg)
    sed -i 's/MOL/LIG/g' solvated.gro
fi

# ── 4b. Ensure solvated.gro has LIG and index.ndx has Protein_LIG ────────────
if [[ -f "solvated.gro" ]] && grep -qm1 'MOL' solvated.gro; then
    echo "Fixing MOL→LIG in solvated.gro"
    sed -i 's/MOL/LIG/g' solvated.gro
fi
if [[ -f "index.ndx" ]] && ! grep -q '\[ Protein_LIG \]' index.ndx; then
    echo "Stale index.ndx (no Protein_LIG group) — deleting for regeneration"
    rm -f index.ndx
fi

# ── 5. Make index ─────────────────────────────────────────────────────────────
if [[ ! -f "index.ndx" ]]; then
    echo "--- make_ndx ---"
    printf 'q\n' | $GMX make_ndx -f solvated.gro -o index.ndx -nobackup

    # Parse LIG group number directly from index.ndx (0-indexed by [ ] section order)
    LIG_GRP=$(awk '/^\[/{grp++} /^\[ LIG \]/{print grp-1; exit}' index.ndx)
    echo "LIG group = ${LIG_GRP:-NOT FOUND}"

    if [[ -n "$LIG_GRP" ]]; then
        N_GRP=$(awk '/^\[/{grp++} END{print grp}' index.ndx)
        printf "1 | %s\nname %s Protein_LIG\nq\n" "$LIG_GRP" "$N_GRP" \
            | $GMX make_ndx -f solvated.gro -n index.ndx -o index.ndx -nobackup || true
    else
        echo "ERROR: LIG group not found in index.ndx — check solvated.gro residue names" >&2
        exit 1
    fi
fi

# ── 6. Energy minimisation ───────────────────────────────────────────────────
if [[ ! -f "minim.gro" ]]; then
    echo "--- energy minimisation ---"
    $GMX grompp -f "$MDPDIR/minim.mdp" -c solvated.gro -p topol.top \
                -o minim.tpr -nobackup
    $MPI $GMX mdrun -v -deffnm minim -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
fi

# ── 7. NVT equilibration ─────────────────────────────────────────────────────
if [[ ! -f "nvt.gro" ]]; then
    echo "--- NVT ---"
    NVT_CPT_FLAG=""
    if [[ ! -f "nvt.tpr" ]]; then
        $GMX grompp -f "$MDPDIR/nvt.mdp" -c minim.gro -r minim.gro \
            -p topol.top -n index.ndx -o nvt.tpr -maxwarn 1 -nobackup
    elif [[ -f "nvt.cpt" ]]; then
        NVT_CPT_FLAG="-cpi nvt.cpt"
    fi
    $MPI $GMX mdrun -deffnm nvt $NVT_CPT_FLAG -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
    if [[ ! -f "nvt.gro" ]]; then
        echo "WARNING: NVT failed — removing checkpoint so next run restarts clean" >&2
        rm -f nvt.cpt nvt.tpr; exit 1
    fi
fi

# ── 8. NPT equilibration ─────────────────────────────────────────────────────
if [[ ! -f "npt.gro" ]]; then
    echo "--- NPT ---"
    NPT_CPT_FLAG=""
    if [[ ! -f "npt.tpr" ]]; then
        $GMX grompp -f "$MDPDIR/npt.mdp" -c nvt.gro -r nvt.gro \
            -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 1 -nobackup
    elif [[ -f "npt.cpt" ]]; then
        NPT_CPT_FLAG="-cpi npt.cpt"
    fi
    $MPI $GMX mdrun -deffnm npt $NPT_CPT_FLAG -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
    if [[ ! -f "npt.gro" ]]; then
        echo "WARNING: NPT failed — removing checkpoint so next run restarts clean" >&2
        rm -f npt.cpt npt.tpr; exit 1
    fi
fi

# ── 9. Production MD (1 ns) ──────────────────────────────────────────────────
if [[ ! -f "prod.gro" ]]; then
    echo "--- production MD ---"
    PROD_CPT_FLAG=""
    if [[ ! -f "prod.tpr" ]]; then
        $GMX grompp -f "$MDPDIR/prod.mdp" -c npt.gro -t npt.cpt \
            -p topol.top -n index.ndx -o prod.tpr -maxwarn 1 -nobackup
    elif [[ -f "prod.cpt" ]]; then
        PROD_CPT_FLAG="-cpi prod.cpt"
    fi
    $MPI $GMX mdrun -deffnm prod $PROD_CPT_FLAG -ntomp $OMP_NUM_THREADS -npme 4 -nobackup
    if [[ ! -f "prod.gro" ]]; then
        echo "WARNING: Production MD failed — removing checkpoint so next run restarts clean" >&2
        rm -f prod.cpt prod_prev.cpt prod.tpr prod.xtc prod.edr prod.log; exit 1
    fi
fi

# ── 9b. PBC correction (MDAnalysis make_whole — bond-by-bond unwrapping) ─────
# GROMACS trjconv -pbc mol/cluster uses naive centre-of-mass which fails for
# homodimers whose chains are split across box boundaries.  MDAnalysis unwrap()
# traverses the bond graph from TPR topology and correctly reassembles each
# fragment before re-centring and wrapping.
# IMPORTANT: always re-run unless MM-GBSA is already complete. A prod_pbc.xtc
# from a prior trjconv run is silently broken for homodimers — do NOT skip based
# on file existence alone.
if [[ ! -f "mmgbsa/FINAL_RESULTS_MMPBSA.csv" ]]; then
    echo "--- PBC correction (MDAnalysis) ---"
    FIX_SCRIPT="$PROJ/09_scripts/08_virtual_screening/fix_pbc_mdanalysis.py"
    $CONDA_PYTHON "$FIX_SCRIPT" prod.tpr prod.xtc prod_pbc.xtc
    if [[ ! -f "prod_pbc.xtc" ]]; then
        echo "ERROR: MDAnalysis PBC correction failed" >&2; exit 1
    fi
fi

# ── 10. gmx_MMPBSA ───────────────────────────────────────────────────────────
if [[ ! -f "mmgbsa/FINAL_RESULTS_MMPBSA.csv" ]]; then
    echo "--- gmx_MMPBSA ---"
    LIG_GRP=$(awk '/^\[/{grp++} /^\[ LIG \]/{print grp-1; exit}' index.ndx)
    if [[ -z "$LIG_GRP" ]]; then
        echo "ERROR: LIG group not found in index.ndx — cannot run gmx_MMPBSA" >&2; exit 1
    fi
    echo "Using LIG group $LIG_GRP for gmx_MMPBSA -cg"
    # parmed ignores GMXLIB and cannot resolve #include directives in topol.top.
    # Generate a preprocessed flat topology with all includes inlined.
    if [[ ! -f "topol_flat.top" ]]; then
        $GMX grompp -f "$MDPDIR/prod.mdp" -c npt.gro -t npt.cpt \
            -p topol.top -n index.ndx -o _flat_tmp.tpr -pp topol_flat.top \
            -maxwarn 1 -nobackup 2>/dev/null
        rm -f _flat_tmp.tpr
    fi
    [[ ! -f "topol_flat.top" ]] && { echo "ERROR: could not generate topol_flat.top" >&2; exit 1; }
    mkdir -p mmgbsa
    cd mmgbsa
    gmx_MMPBSA \
        -O \
        -i "$MDPDIR/mmgbsa.in" \
        -cs ../prod.tpr \
        -ct ../prod_pbc.xtc \
        -ci ../index.ndx \
        -cg 1 "$LIG_GRP" \
        -cp ../topol_flat.top \
        -nogui \
        -o FINAL_RESULTS_MMPBSA.dat \
        -eo FINAL_RESULTS_MMPBSA.csv
    cd ..
    if [[ ! -f "mmgbsa/FINAL_RESULTS_MMPBSA.csv" ]]; then
        echo "WARNING: gmx_MMPBSA failed for $POCKET / $LIG_ID" >&2; exit 1
    fi
fi

echo "=== DONE: $POCKET / $LIG_ID ==="
