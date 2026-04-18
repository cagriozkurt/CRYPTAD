"""
prepare_ligands.py  —  CRYPTAD
Ligand library preparation for virtual screening (Step 4.2).

Workflow:
  1. Fetch approved small molecules from ChEMBL (max_phase=4) OR read a local SDF
  2. Standardise: strip salts, neutralise charges
  3. Apply CNS/BBB filter (MW, logP, TPSA, HBD, HBA, RotBonds)
  4. Discard PAINS compounds
  5. Generate 3-D conformers (ETKDG v3 + MMFF94)
  6. Convert to PDBQT via Meeko for AutoDock Vina

Output layout:
  04_virtual_screening/ligands/
    ├── library_filtered.sdf        filtered 3-D SDF
    ├── library_manifest.csv        per-compound properties + pass/fail reason
    └── pdbqt/{chembl_id}.pdbqt     one file per docking-ready ligand

CNS filter (BBB-penetrant small molecules):
  MW      150–450 Da
  LogP   -0.5 – 5.0
  TPSA    ≤ 90 Å²
  HBD     ≤ 3
  HBA     ≤ 7
  RotBonds ≤ 8

Usage:
  python3 09_scripts/08_virtual_screening/02_prepare_ligands.py              # fetch ChEMBL
  python3 09_scripts/08_virtual_screening/02_prepare_ligands.py --sdf drugs.sdf
  python3 09_scripts/08_virtual_screening/02_prepare_ligands.py --max 500    # limit compounds
  python3 09_scripts/08_virtual_screening/02_prepare_ligands.py --project-root /path/to/CRYPTAD
"""

import argparse
import csv
import json
import sys
import time
import urllib.request
import urllib.parse
import urllib.error
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

from meeko import MoleculePreparation, PDBQTWriterLegacy

# ---------------------------------------------------------------------------
_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; RED   = "\033[91m"; CYAN   = "\033[96m"
def log(m, c=RESET): print(f"{c}{m}{RESET}", flush=True)

# CNS/BBB physicochemical filter thresholds
CNS = dict(mw_min=150, mw_max=450, logp_min=-0.5, logp_max=5.0,
           tpsa_max=90, hbd_max=3, hba_max=7, rotbonds_max=8)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data/molecule"


# ---------------------------------------------------------------------------
# Standardisation helpers
# ---------------------------------------------------------------------------
_lfc       = rdMolStandardize.LargestFragmentChooser()
_uncharger = rdMolStandardize.Uncharger()

def standardise(mol):
    """Strip salts (largest fragment) and neutralise formal charges."""
    mol = _lfc.choose(mol)
    mol = _uncharger.uncharge(mol)
    return mol


# ---------------------------------------------------------------------------
# CNS/BBB filter
# ---------------------------------------------------------------------------
def cns_props(mol):
    mw   = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
    return dict(mw=round(mw,2), logp=round(logp,2), tpsa=round(tpsa,1),
                hbd=hbd, hba=hba, rotbonds=rotb)

def passes_cns(props):
    return (CNS["mw_min"]      <= props["mw"]       <= CNS["mw_max"]   and
            CNS["logp_min"]    <= props["logp"]      <= CNS["logp_max"] and
            props["tpsa"]      <= CNS["tpsa_max"]                       and
            props["hbd"]       <= CNS["hbd_max"]                        and
            props["hba"]       <= CNS["hba_max"]                        and
            props["rotbonds"]  <= CNS["rotbonds_max"])

def cns_fail_reason(props):
    reasons = []
    if not (CNS["mw_min"] <= props["mw"] <= CNS["mw_max"]):
        reasons.append(f"MW={props['mw']}")
    if not (CNS["logp_min"] <= props["logp"] <= CNS["logp_max"]):
        reasons.append(f"LogP={props['logp']}")
    if props["tpsa"] > CNS["tpsa_max"]:
        reasons.append(f"TPSA={props['tpsa']}")
    if props["hbd"] > CNS["hbd_max"]:
        reasons.append(f"HBD={props['hbd']}")
    if props["hba"] > CNS["hba_max"]:
        reasons.append(f"HBA={props['hba']}")
    if props["rotbonds"] > CNS["rotbonds_max"]:
        reasons.append(f"RotB={props['rotbonds']}")
    return "; ".join(reasons)


# ---------------------------------------------------------------------------
# PAINS filter
# ---------------------------------------------------------------------------
def build_pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)

def passes_pains(mol, catalog):
    return catalog.GetFirstMatch(mol) is None


# ---------------------------------------------------------------------------
# 3-D conformer generation
# ---------------------------------------------------------------------------
def generate_3d(mol):
    """
    Add hydrogens, embed with ETKDGv3, minimise with MMFF94.
    Returns the mol with 3-D coords, or None on failure.
    """
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    params.enforceChirality = True
    if AllChem.EmbedMolecule(mol, params) == -1:
        return None  # embedding failed
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:
        pass  # keep unoptimised conformer
    return mol


# ---------------------------------------------------------------------------
# PDBQT conversion via Meeko
# ---------------------------------------------------------------------------
_preparator = MoleculePreparation()

def mol_to_pdbqt(mol, out_path: Path) -> bool:
    """Convert an RDKit mol (with 3-D coords) to PDBQT. Returns True on success."""
    try:
        mol_setups = _preparator.prepare(mol)
        for setup in mol_setups:
            pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setup)
            if is_ok:
                out_path.write_text(pdbqt_str)
                return True
            else:
                log(f"    [meeko] {err}", YELLOW)
    except Exception as e:
        log(f"    [meeko exception] {e}", YELLOW)
    return False


# ---------------------------------------------------------------------------
# ChEMBL fetcher
# ---------------------------------------------------------------------------
def fetch_chembl_smiles(max_compounds: int):
    """
    Yield (chembl_id, smiles, name) tuples for max_phase=4 small molecules
    from the ChEMBL REST API (paginated, JSON).
    """
    params = urllib.parse.urlencode({
        "max_phase":     4,
        "molecule_type": "Small molecule",
        "format":        "json",
        "limit":         500,
        "offset":        0,
    })
    url     = f"{CHEMBL_BASE}?{params}"
    fetched = 0
    page    = 0
    while url and fetched < max_compounds:
        page += 1
        log(f"  [ChEMBL] fetching page {page} …", CYAN)
        try:
            with urllib.request.urlopen(url, timeout=30) as resp:
                data = json.loads(resp.read())
        except Exception as e:
            log(f"  [ChEMBL fetch error] {e}", RED)
            break

        for mol in data.get("molecules", []):
            if fetched >= max_compounds:
                return
            cid     = mol.get("molecule_chembl_id", "")
            structs = mol.get("molecule_structures") or {}
            smi     = structs.get("canonical_smiles", "")
            name    = mol.get("pref_name") or cid
            if smi:
                yield cid, smi, name
                fetched += 1

        next_path = (data.get("page_meta") or {}).get("next")
        url = (f"https://www.ebi.ac.uk{next_path}" if next_path else None)
        time.sleep(0.3)   # polite delay

    log(f"  [ChEMBL] fetched {fetched} molecules total", CYAN)


# ---------------------------------------------------------------------------
# Core processing loop
# ---------------------------------------------------------------------------
def process_compounds(source_iter, outdir: Path, project_root: Path, label: str):
    """
    source_iter yields (compound_id, smiles_or_mol, name) tuples.
    Runs full pipeline: standardise → CNS filter → PAINS → 3-D → PDBQT.
    Writes library_filtered.sdf, library_manifest.csv, pdbqt/*.pdbqt.
    """
    pdbqt_dir = outdir / "pdbqt"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = outdir / "library_manifest.csv"
    sdf_path      = outdir / "library_filtered.sdf"

    writer  = Chem.SDWriter(str(sdf_path))
    catalog = build_pains_catalog()

    manifest_fields = ["id", "name", "smiles", "mw", "logp", "tpsa",
                       "hbd", "hba", "rotbonds", "status", "reason", "pdbqt"]

    n_total = n_pass = n_fail_parse = n_fail_cns = n_fail_pains = 0
    n_fail_3d = n_fail_pdbqt = 0

    with open(manifest_path, "w", newline="") as mf:
        mw = csv.DictWriter(mf, fieldnames=manifest_fields)
        mw.writeheader()

        for cid, smi_or_mol, name in source_iter:
            n_total += 1
            row = dict(id=cid, name=name, smiles="", mw="", logp="", tpsa="",
                       hbd="", hba="", rotbonds="", status="FAIL", reason="", pdbqt="")

            # 1. Parse / standardise
            if isinstance(smi_or_mol, str):
                mol = Chem.MolFromSmiles(smi_or_mol)
                row["smiles"] = smi_or_mol
            else:
                mol = smi_or_mol
                row["smiles"] = Chem.MolToSmiles(mol) if mol else ""

            if mol is None:
                n_fail_parse += 1
                row["reason"] = "parse_fail"
                mw.writerow(row)
                continue

            try:
                mol = standardise(mol)
            except Exception:
                pass
            row["smiles"] = Chem.MolToSmiles(mol)

            # 2. CNS filter
            props = cns_props(mol)
            row.update(props)
            if not passes_cns(props):
                n_fail_cns += 1
                row["reason"] = cns_fail_reason(props)
                mw.writerow(row)
                continue

            # 3. PAINS filter
            if not passes_pains(mol, catalog):
                n_fail_pains += 1
                row["reason"] = "PAINS"
                mw.writerow(row)
                continue

            # 4. 3-D conformer
            mol3d = generate_3d(mol)
            if mol3d is None:
                n_fail_3d += 1
                row["reason"] = "conformer_fail"
                mw.writerow(row)
                continue

            mol3d.SetProp("_Name", cid)
            mol3d.SetProp("ChEMBL_ID", cid)
            mol3d.SetProp("Preferred_Name", name)

            # 5. PDBQT
            pdbqt_path = pdbqt_dir / f"{cid}.pdbqt"
            if not mol_to_pdbqt(mol3d, pdbqt_path):
                n_fail_pdbqt += 1
                row["reason"] = "pdbqt_fail"
                mw.writerow(row)
                continue

            # All passed
            writer.write(mol3d)
            row["status"] = "PASS"
            row["pdbqt"]  = str(pdbqt_path.relative_to(project_root))
            mw.writerow(row)
            n_pass += 1

            if n_pass % 100 == 0:
                log(f"  … {n_pass} compounds ready  ({n_total} processed)", CYAN)

    writer.close()

    log(f"\n{'═'*60}", BOLD)
    log(f"  {label}: {n_total} input compounds", BOLD)
    log(f"  ✓ PASS (docking-ready): {n_pass}", GREEN)
    log(f"  ✗ Failed parse:         {n_fail_parse}", YELLOW)
    log(f"  ✗ Failed CNS filter:    {n_fail_cns}", YELLOW)
    log(f"  ✗ Failed PAINS:         {n_fail_pains}", YELLOW)
    log(f"  ✗ Failed 3-D embed:     {n_fail_3d}", YELLOW)
    log(f"  ✗ Failed PDBQT:         {n_fail_pdbqt}", YELLOW)
    log(f"\n  SDF  → {sdf_path.relative_to(project_root)}", GREEN)
    log(f"  CSV  → {manifest_path.relative_to(project_root)}", GREEN)
    log(f"  PDBQT dir → {pdbqt_dir.relative_to(project_root)}", GREEN)

    return n_pass, {
        "n_total":       n_total,
        "n_pass":        n_pass,
        "n_fail_parse":  n_fail_parse,
        "n_fail_cns":    n_fail_cns,
        "n_fail_pains":  n_fail_pains,
        "n_fail_3d":     n_fail_3d,
        "n_fail_pdbqt":  n_fail_pdbqt,
        "sdf_path":      str(sdf_path),
        "manifest_path": str(manifest_path),
        "pdbqt_dir":     str(pdbqt_dir),
    }


# ---------------------------------------------------------------------------
# Source: local SDF
# ---------------------------------------------------------------------------
def iter_sdf(sdf_path: Path):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=True)
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        cid  = mol.GetProp("_Name") if mol.HasProp("_Name") else f"cpd_{i:06d}"
        name = (mol.GetProp("GENERIC_NAME") if mol.HasProp("GENERIC_NAME")
                else mol.GetProp("PREFERRED_NAME") if mol.HasProp("PREFERRED_NAME")
                else cid)
        yield cid, mol, name


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="CRYPTAD ligand library prep (Step 4.2)")
    src = parser.add_mutually_exclusive_group()
    src.add_argument("--sdf",  type=Path, default=None,
                     help="Local SDF file (DrugBank approved_structures.sdf, etc.)")
    src.add_argument("--fetch-chembl", action="store_true",
                     help="Download approved drugs from ChEMBL REST API (default)")
    parser.add_argument("--max", type=int, default=10_000,
                        help="Maximum compounds to process (default: 10000)")
    parser.add_argument("--outdir", type=Path, default=None,
                        help="Output directory (default: <project-root>/04_virtual_screening/ligands)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    outdir = args.outdir if args.outdir else project_root / "04_virtual_screening" / "ligands"
    outdir.mkdir(parents=True, exist_ok=True)

    log(f"\n{BOLD}CRYPTAD — Ligand Library Preparation (Step 4.2){RESET}")
    log(f"CNS filter: MW {CNS['mw_min']}–{CNS['mw_max']} Da  |  "
        f"LogP {CNS['logp_min']}–{CNS['logp_max']}  |  "
        f"TPSA ≤{CNS['tpsa_max']} Å²  |  "
        f"HBD ≤{CNS['hbd_max']}  HBA ≤{CNS['hba_max']}  RotB ≤{CNS['rotbonds_max']}")
    log(f"Output: {outdir.relative_to(project_root)}\n")

    if args.sdf:
        if not args.sdf.exists():
            log(f"[ERROR] SDF not found: {args.sdf}", RED)
            sys.exit(1)
        log(f"Source: local SDF → {args.sdf}", CYAN)
        source = iter_sdf(args.sdf)
        label  = args.sdf.stem
    else:
        log(f"Source: ChEMBL REST API (max_phase=4, max={args.max})", CYAN)
        source = ((cid, smi, name)
                  for cid, smi, name in fetch_chembl_smiles(args.max))
        label  = "ChEMBL approved"

    n_ready, stats = process_compounds(source, outdir, project_root, label)

    if n_ready > 0:
        log(f"\n  Next: Step 4.3 — docking on TRUBA", GREEN)
        log(f"  python3 09_scripts/08_virtual_screening/03_run_docking.py", GREEN)
    else:
        log(f"\n  [WARNING] No compounds passed all filters. Check input library.", YELLOW)

    # Write manifest
    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":       "09_scripts/08_virtual_screening/02_prepare_ligands.py",
        "project_root": str(project_root),
        "source":       str(args.sdf) if args.sdf else "ChEMBL REST API (max_phase=4)",
        "max_compounds": args.max,
        "cns_filter":   CNS,
        **stats,
    }
    manifest_path = outdir / "prepare_ligands_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
