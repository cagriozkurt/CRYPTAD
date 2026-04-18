# CRYPTAD Pipeline — Script Reference

All scripts accept `--project-root /path/to/CRYPTAD` (default: auto-detected from
the script's own location). Run everything from the **project root** or pass the
flag explicitly.

---

## 01 · Structure Preparation

```bash
python3 09_scripts/01_structure_prep/01_download_structures.py
python3 09_scripts/01_structure_prep/02_process_af3_results.py
python3 09_scripts/01_structure_prep/03_prepare_structures.py
```

---

## 02 · MD Setup

Build GROMACS+PLUMED on the HPC cluster (run once):

```bash
bash 09_scripts/02_md_setup/01_build_gromacs_plumed.sh   # with PLUMED
bash 09_scripts/02_md_setup/02_build_gromacs_only.sh      # standard GROMACS
```

Generate SLURM run scripts and rsync commands:

```bash
python3 09_scripts/02_md_setup/03_setup_truba_runs.py
```

---

## 03 · Metadynamics

One script per system. Submit as SLURM array on the HPC cluster:

```bash
sbatch 09_scripts/03_metadynamics/01_run_metad_S1_BIN1_BAR.sh
sbatch 09_scripts/03_metadynamics/01_run_metad_S2_BIN1_SH3.sh
sbatch 09_scripts/03_metadynamics/01_run_metad_S3_PICALM_ANTH.sh
sbatch 09_scripts/03_metadynamics/01_run_metad_S4_CD2AP_SH3-2.sh
sbatch 09_scripts/03_metadynamics/01_run_metad_S5_CD2AP_SH3-1.sh
```

Sum hills (FES reconstruction) after metadynamics converges:

```bash
sbatch 09_scripts/03_metadynamics/02_run_sum_hills.sh
```

---

## 04 · Trajectory QC

```bash
# On HPC cluster:
sbatch 09_scripts/04_trajectory_qc/01_preprocess_trajectories.sh
sbatch 09_scripts/04_trajectory_qc/02_fix_pbc_dimer_s1.sh
sbatch 09_scripts/04_trajectory_qc/03_qc_trajectories.sh
sbatch 09_scripts/04_trajectory_qc/04_qc_perchain_s1.sh

# Locally (after rsync):
python3 09_scripts/04_trajectory_qc/05_plot_qc.py
python3 09_scripts/04_trajectory_qc/06_plot_qc_perchain.py
```

---

## 05 · Pocket Detection (fpocket / mdpocket)

```bash
# On HPC cluster:
sbatch 09_scripts/05_pocket_detection/01_setup_fpocket_dirs.sh
sbatch 09_scripts/05_pocket_detection/02_extract_metad_frames.sh
sbatch 09_scripts/05_pocket_detection/02b_extract_metad_frames_S1.sh   # BIN1 BAR dimer
sbatch 09_scripts/05_pocket_detection/03_run_mdpocket.sh

# Locally:
python3 09_scripts/05_pocket_detection/04_run_fpocket_frames.py
python3 09_scripts/05_pocket_detection/05_parse_fpocket.py
python3 09_scripts/05_pocket_detection/06_inspect_pockets.py
python3 09_scripts/05_pocket_detection/07_select_cv_atoms.py
```

---

## 06 · Persistence Gate

```bash
# On HPC cluster:
sbatch 09_scripts/06_persistence_gate/01_extract_persist_start_S1.sh
sbatch 09_scripts/06_persistence_gate/01_extract_persist_start_S3.sh
sbatch 09_scripts/06_persistence_gate/02_run_persist_S1.sh
sbatch 09_scripts/06_persistence_gate/02_run_persist_S3.sh
sbatch 09_scripts/06_persistence_gate/03_extract_persist_protein_S1.sh
sbatch 09_scripts/06_persistence_gate/03_extract_persist_protein_S3.sh

# Locally:
python3 09_scripts/06_persistence_gate/04_analyze_persistence_gate.py
```

---

## 07 · Conformational Analysis

```bash
python3 09_scripts/07_conformational/01_analyze_fes_convergence.py
python3 09_scripts/07_conformational/02_analyze_metad_colvar.py
python3 09_scripts/07_conformational/03_analyze_metad_pockets.py
python3 09_scripts/07_conformational/04_extract_cryptic_frames.py
python3 09_scripts/07_conformational/05_conformational_coupling.py
```

---

## 08 · Virtual Screening

### Step 1 — Prepare receptors (local)

```bash
python3 09_scripts/08_virtual_screening/01_prepare_receptors.py
```

Requires `pdb2pqr` and `mk_prepare_receptor.py` (meeko) in PATH.

### Step 2 — Prepare ligand library (local)

```bash
python3 09_scripts/08_virtual_screening/02_prepare_ligands.py
```

Requires RDKit and `mk_prepare_ligand.py` (meeko) in PATH.

### Step 3 — Docking

Generate SLURM script and rsync commands, then submit on the HPC cluster:

```bash
# Local — generates 03b_run_docking_truba.sh:
python3 09_scripts/08_virtual_screening/03_run_docking.py --no-rsync

# HPC cluster:
sbatch 09_scripts/08_virtual_screening/03b_run_docking_truba.sh
```

> **Note:** Vina requires the `vina_env` conda environment:
> `conda activate vina_env`

### Step 4 — Parse docking results (local)

```bash
python3 09_scripts/08_virtual_screening/04_parse_docking_results.py
```

Outputs `docking_scores_all.csv`, per-pocket top-50 CSVs, consensus hits, and
a score distribution plot.

### Step 5 — Set up MM-GBSA (local → HPC cluster)

```bash
# Local — generates per-compound SLURM directories:
python3 09_scripts/08_virtual_screening/05_setup_mmgbsa.py

# HPC cluster:
sbatch 09_scripts/08_virtual_screening/05b_run_mmgbsa_truba.sh
```

### Step 6 — Parse MM-GBSA results (local)

```bash
python3 09_scripts/08_virtual_screening/06_parse_mmgbsa_results.py
# Optional flags:
#   --le-threshold 0.3   (default: 0.3)
#   --top-n 20           (default: 20)
```

### Step 7 — ADMET filtering (local)

```bash
python3 09_scripts/08_virtual_screening/07_admet_filter.py
# Optional flags:
#   --dg-max -20.0       (default: -20.0 kcal/mol)
#   --cns-mpo-min 3.0    (default: 3.0)
```

### Step 8 — Selectivity docking (local, vina_env)

```bash
conda run -n vina_env \
  python3 09_scripts/08_virtual_screening/08_selectivity_docking.py
```

---

## 09 · Complex MD

### Step 1 — Set up replica directories (local)

```bash
python3 09_scripts/09_complex_md/01_setup_complex_md.py
```

Writes MDP files and copies topology into
`02_md_simulations/complex_md/{system}/rep{1..N}/`.
Prints rsync and sbatch commands to run next.

### Step 2 — Run complex MD on HPC cluster

```bash
# CPU partition (12 tasks: 4 systems × 3 replicas):
sbatch 09_scripts/09_complex_md/02_run_complex_md_truba.sh

# GPU partition (BIN1 tasks):
sbatch 09_scripts/09_complex_md/02b_run_complex_md_gpu.sh
```

### Step 3 — Fix PBC (HPC cluster)

```bash
# Generate pairs list (run once):
bash 09_scripts/09_complex_md/03_run_fix_pbc_truba.sh --make-pairs

# Submit array:
N=$(wc -l < 02_md_simulations/complex_md/fix_pbc_pairs.txt)
sbatch --array=1-${N} 09_scripts/09_complex_md/03_run_fix_pbc_truba.sh
```

Or run locally on a single trajectory:

```bash
python3 09_scripts/09_complex_md/04_fix_pbc_mdanalysis.py <tpr> <xtc> <out.xtc>
# Auto-mode (all replicas):
python3 09_scripts/09_complex_md/04_fix_pbc_mdanalysis.py
```

### Step 4 — Drift analysis (VMD)

```bash
vmd -dispdev none -e 09_scripts/09_complex_md/05_drift_analysis.tcl
# Override paths:
CRYPTAD_ROOT=/path/to/CRYPTAD vmd -dispdev none \
  -e 09_scripts/09_complex_md/05_drift_analysis.tcl
```

Output: `04_virtual_screening/complex_md_results/drift_results.txt`

### Step 5 — Ibufenac contact analysis (local)

```bash
python3 09_scripts/09_complex_md/06_ibufenac_contacts.py
```

Output: `04_virtual_screening/complex_md_results/ibufenac_contacts.txt` and
`ibufenac_drift.csv`

### Step 6 — Binding mode analysis (local)

```bash
python3 09_scripts/09_complex_md/07_binding_mode_analysis.py
```

Output: `03_pocket_analysis/binding_mode/contacts_all.csv`, per-hit PDBs and
PyMOL scripts.

---

## 10 · Figures

```bash
# Complex MD RMSD + Ibufenac H-bond heatmap:
python3 09_scripts/10_figures/01_complex_md_rmsd_hbond.py

# RMSF apo vs. holo (PICALM ANTH):
python3 09_scripts/10_figures/02_rmsf_apo_holo.py
# Use --no-cache to recompute (skips .npy cache by default):
python3 09_scripts/10_figures/02_rmsf_apo_holo.py --no-cache

# Conformational composite figure (Figure 9):
python3 09_scripts/10_figures/03_assemble_conformational.py

# STRING PPI network:
python3 09_scripts/10_figures/04_string_network.py

# Virtual screening figures (Figures 1–5):
python3 09_scripts/10_figures/05_generate_figures.py
```

All figures are saved to `06_figures/` at 300 dpi.

---

## Environment Setup

All local Python scripts (analysis, figures, virtual screening) run inside the
`cryptad` conda environment:

```bash
conda env create -f environment.yml
conda activate cryptad
```

This installs Python 3.11, RDKit, MDAnalysis, AutoDock Vina 1.2.6, meeko,
gmx_MMPBSA, acpype, and all other Python dependencies.

To generate a fully pinned lockfile after creating the environment:

```bash
conda env export --no-builds > environment.lock.yml
```

### External tools (install separately)

| Tool | Purpose | Install |
|------|---------|---------|
| fpocket | Pocket detection (steps 05, 07) | Compile from [github.com/Discngine/fpocket](https://github.com/Discngine/fpocket); macOS ARM: `make ARCH=MACOSXARM64` |
| VMD | Drift analysis TCL script (step 09.4) | [ks.uiuc.edu/Research/vmd](https://www.ks.uiuc.edu/Research/vmd/) — free, requires registration |
| GROMACS + PLUMED | Metadynamics and complex MD (steps 03, 09) | Build with `singularity/gromacs-plumed.def` (see below) |

### GROMACS + PLUMED container

Build the Apptainer image once on the HPC cluster (requires `--fakeroot` or root; ~30–60 min):

```bash
apptainer build singularity/gromacs-plumed.sif singularity/gromacs-plumed.def
```

The SLURM scripts call `gmx_mpi` directly (assuming it is in PATH after sourcing the
site module or the build helper). To use the container instead, replace `gmx_mpi` calls with:

```bash
apptainer exec --nv singularity/gromacs-plumed.sif gmx_mpi ...
# Multi-node (hybrid MPI):
mpirun -np $SLURM_NTASKS \
  apptainer exec --nv --bind $SCRATCH singularity/gromacs-plumed.sif \
  gmx_mpi mdrun -ntomp $OMP_NUM_THREADS -npme 4 -deffnm md
```

### HPC cluster notes

All SLURM scripts target the CPU partition (`--ntasks=28 --cpus-per-task=2`)
unless noted otherwise. The GPU script targets a CUDA-enabled GPU partition.
Override the project root on the HPC cluster with:

```bash
export CRYPTAD=/path/to/CRYPTAD
```
