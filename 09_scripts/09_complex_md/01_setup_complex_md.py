"""
CRYPTAD — Phase 4 complex MD setup (Step 5.1).

Creates 3 independent 100 ns replica directories for each of the 4 priority
receptor-ligand systems, starting from the MM-GBSA NPT equilibration endpoint.

Priority systems:
  S1_site688 / CHEMBL4105630  — Zuranolone @ BIN1 primary pocket
  S1_site680 / CHEMBL4105630  — Zuranolone @ BIN1 site680
  S3_site473 / CHEMBL341812   — Ibufenac   @ PICALM
  S3_site473 / CHEMBL797      — Phensuximide @ PICALM

Directory layout:
  02_md_simulations/complex_md/{system}/rep{1,2,3}/
    equil.mdp      — 200 ps Berendsen re-equilibration MDP (gen_vel=yes, unique seed)
    complexmd.mdp  — 100 ns production MDP (continuation=yes)
    topol_flat.top — preprocessed topology (copied from MM-GBSA dir)
    *.itp          — force-field include files
    index.ndx      — GROMACS index
    npt.gro        — starting coordinates (end of NPT)
    npt.cpt        — checkpoint (passed to grompp -t)

Usage:
    python3 09_scripts/09_complex_md/01_setup_complex_md.py
    python3 09_scripts/09_complex_md/01_setup_complex_md.py --project-root /path/to/CRYPTAD

Then rsync to TRUBA and sbatch the SLURM array.
"""

import argparse
import json
import shutil
import time
from pathlib import Path

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

SYSTEMS = [
    ("S1_site688", "CHEMBL4105630", "Zuranolone",   "site688"),
    ("S1_site680", "CHEMBL4105630", "Zuranolone",   "site680"),
    ("S3_site473", "CHEMBL341812",  "Ibufenac",     "site473"),
    ("S3_site473", "CHEMBL797",     "Phensuximide", "site473"),
]

N_REPLICAS = 3
SEEDS      = [42, 1337, 2718]      # one per replica

# ── 200 ps Berendsen re-equilibration MDP ────────────────────────────────────
# Assigns new velocities (gen_vel=yes) then equilibrates pressure with Berendsen.
# dt=0.002 ps — no HMR in MM-GBSA topologies; GAFF2 parameterised for 2 fs.
# 100,000 × 0.002 ps = 200 ps

EQUIL_MDP_TEMPLATE = """\
; CRYPTAD complex MD — 200 ps Berendsen re-equilibration (Phase 4, Step 5.1)
; Replica: {rep}  |  System: {sysname}  |  Seed: {seed}
integrator              = md
nsteps                  = 100000
dt                      = 0.002

; Minimal output (checkpoint only)
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 100000
nstlog                  = 100000
nstxout-compressed      = 0

; Assign new velocities for independent replica
continuation            = no
gen_vel                 = yes
gen_temp                = 310
gen_seed                = {seed}

; Neighbour search
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rcoulomb                = 1.2
rvdw                    = 1.2

; Electrostatics
coulombtype             = PME
pme_order               = 4
ewald_rtol              = 1e-5

; Temperature (V-rescale, 310 K)
tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1   0.1
ref_t                   = 310   310

; Pressure (Berendsen — stable for startup after velocity reassignment)
pcoupl                  = Berendsen
pcoupltype              = isotropic
tau_p                   = 0.5
ref_p                   = 1.0
compressibility         = 4.5e-5

; Constraints
constraints             = h-bonds
constraint_algorithm    = LINCS
pbc                     = xyz
DispCorr                = EnerPres
"""

# ── 100 ns production MDP ─────────────────────────────────────────────────────
# dt = 0.002 ps — no HMR in MM-GBSA topologies; GAFF2 parameterised for 2 fs.
# 50,000,000 × 0.002 = 100,000 ps = 100 ns
# Output every 100 ps (nstxout-compressed = 50000) → 1001 frames
# continuation=yes / gen_vel=no — velocities were set in equil step

MDP_TEMPLATE = """\
; CRYPTAD complex MD — 100 ns production (Phase 4, Step 5.1)
; Replica: {rep}  |  System: {sysname}  |  Seed: {seed}
integrator              = md
nsteps                  = 50000000
dt                      = 0.002

; Output (every 100 ps)
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 50000
nstlog                  = 50000
nstxout-compressed      = 50000
compressed-x-grps       = System

; Continue from Berendsen equil (velocities already set)
continuation            = yes
gen_vel                 = no

; Neighbour search
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rcoulomb                = 1.2
rvdw                    = 1.2

; Electrostatics
coulombtype             = PME
pme_order               = 4
ewald_rtol              = 1e-5

; Temperature (V-rescale, 310 K)
tcoupl                  = V-rescale
tc-grps                 = Protein_LIG Water_and_ions
tau_t                   = 0.1   0.1
ref_t                   = 310   310

; Pressure (C-rescale, 1 bar)
; C-rescale produces correct NPT ensemble, no startup instability.
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

; Constraints
constraints             = h-bonds
constraint_algorithm    = LINCS
pbc                     = xyz
DispCorr                = EnerPres
"""

# ── Copy files needed per replica ────────────────────────────────────────────
ITP_PATTERNS = ["*.itp", "topol_flat.top", "topol.top", "index.ndx",
                "npt.gro", "npt.cpt"]


def copy_topology(src_dir: Path, dst_dir: Path):
    """Copy all topology / coordinate files from MM-GBSA dir to replica dir."""
    for pattern in ITP_PATTERNS:
        for f in src_dir.glob(pattern):
            shutil.copy2(f, dst_dir / f.name)
    for required in ("npt.gro", "npt.cpt", "topol_flat.top", "index.ndx"):
        if not (dst_dir / required).exists():
            raise FileNotFoundError(f"Missing {required} in {dst_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD complex MD setup (Step 5.1)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    mmgbsa_dir   = project_root / "04_virtual_screening" / "mmgbsa"
    outbase      = project_root / "02_md_simulations" / "complex_md"
    scripts_dir  = project_root / "09_scripts" / "09_complex_md"
    truba_base   = "/arf/scratch/mozkurt/CRYPTAD"

    outbase.mkdir(parents=True, exist_ok=True)

    pairs_lines = []
    created     = []

    for pocket, chembl_id, compound, site in SYSTEMS:
        src_dir  = mmgbsa_dir / pocket / chembl_id
        sys_name = f"{pocket}_{chembl_id}"
        sys_dir  = outbase / sys_name

        if not src_dir.exists():
            print(f"  [SKIP] {sys_name} — MM-GBSA dir not found"); continue

        print(f"\n{sys_name}  ({compound} @ {site})")

        for rep in range(1, N_REPLICAS + 1):
            seed    = SEEDS[rep - 1]
            rep_dir = sys_dir / f"rep{rep}"
            rep_dir.mkdir(parents=True, exist_ok=True)

            try:
                copy_topology(src_dir, rep_dir)
            except FileNotFoundError as e:
                print(f"  [ERROR] rep{rep}: {e}"); continue

            equil_mdp = EQUIL_MDP_TEMPLATE.format(rep=rep, sysname=sys_name, seed=seed)
            (rep_dir / "equil.mdp").write_text(equil_mdp)

            mdp = MDP_TEMPLATE.format(rep=rep, sysname=sys_name, seed=seed)
            (rep_dir / "complexmd.mdp").write_text(mdp)

            print(f"  rep{rep}  seed={seed}  → {rep_dir.relative_to(project_root)}")
            pairs_lines.append(f"{sys_name} rep{rep}")
            created.append({"system": sys_name, "rep": rep, "seed": seed})

    # Write pairs.txt for SLURM array
    pairs_path = scripts_dir / "complex_md_pairs.txt"
    pairs_path.write_text("\n".join(pairs_lines) + "\n")
    print(f"\nPairs file: {pairs_path.relative_to(project_root)}  ({len(pairs_lines)} tasks)")

    print(f"""
─────────────────────────────────────────────────────────────
NEXT: rsync + submit on TRUBA

  rsync -av --progress \\
      "{outbase}/" \\
      mozkurt@levrek1.ulakbim.gov.tr:{truba_base}/02_md_simulations/complex_md/

  rsync -av \\
      "{scripts_dir}/02_run_complex_md_truba.sh" \\
      "{pairs_path}" \\
      mozkurt@levrek1.ulakbim.gov.tr:{truba_base}/09_scripts/09_complex_md/

  ssh truba 'cd {truba_base} && sbatch 09_scripts/09_complex_md/02_run_complex_md_truba.sh'
─────────────────────────────────────────────────────────────
""")

    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":       "09_scripts/09_complex_md/01_setup_complex_md.py",
        "project_root": str(project_root),
        "n_systems":    len(SYSTEMS),
        "n_replicas":   N_REPLICAS,
        "seeds":        SEEDS,
        "pairs_file":   str(pairs_path),
        "replicas_created": created,
    }
    manifest_path = outbase / "setup_complex_md_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"Manifest → {manifest_path.relative_to(project_root)}")


if __name__ == "__main__":
    main()
