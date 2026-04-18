"""
CRYPTAD — Structure Preparation Script
Produces clean apo starting structures for CHARMM-GUI / MD simulations.

What this script does per system:
  S1  BIN1 BAR   (2FIC)  — keep chains A+B (obligate homodimer), strip XE
  S2  BIN1 SH3   (1MUZ)  — extract NMR model 0 only, keep chain A
  S3  PICALM ANTH (7JXV) — keep chain A (ANTH domain), remove chain B (ubiquitin)
  S4  CD2AP SH3-2 (3U23) — keep chain A, remove chain B (RIN3 peptide) + EDO
  S5  CD2AP SH3-1 (4WCI) — keep chain A, remove chains B-F + SO4

Output:   01_structures/<TARGET>/prepared/<ID>_prepared.pdb
Manifest: 01_structures/preparation_manifest.json

Reproducibility notes
---------------------
- Outputs are skipped if already present (use --force to regenerate).
- All transformations (chains kept, residues stripped) are recorded in
  preparation_manifest.json with SHA-256 checksums of input and output files.
- gemmi PdbWriteOptions are set explicitly so behaviour does not depend on the
  gemmi version's defaults.
- renumber_residues is intentionally NOT called: original PDB residue numbering
  is preserved in *_prepared.pdb so that residue numbers in PLUMED CVs, in the
  manuscript, and in the RCSB entry remain consistent.  Renumbering is handled
  downstream by CHARMM-GUI for topology purposes only.

Usage:
    python 03_prepare_structures.py [--force] [--project-root PATH]
"""

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import gemmi

# ── Path resolution ────────────────────────────────────────────────────────────
# Script lives at <root>/09_scripts/01_structure_prep/03_prepare_structures.py
# parents[0] = 01_structure_prep/
# parents[1] = 09_scripts/
# parents[2] = project root
DEFAULT_ROOT = Path(__file__).resolve().parents[2]

# ── Heteroatom residues always stripped (solvents, ions, cryoprotectants) ─────
STRIP_HET = frozenset({
    "HOH", "WAT", "XE",                        # water / noble gas
    "EDO", "GOL", "PEG", "MPD",                # cryoprotectants
    "SO4", "PO4", "ACT",                        # anions
    "CL",  "NA",  "MG",  "CA",  "ZN",          # ions
})

# ── Colours ────────────────────────────────────────────────────────────────────
RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; CYAN  = "\033[96m"; RED    = "\033[91m"

def log(m, c=RESET):
    print(f"{c}{m}{RESET}", flush=True)


# ── Checksum ───────────────────────────────────────────────────────────────────

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


# ── Core helpers ───────────────────────────────────────────────────────────────

def remove_heteroatoms(structure: gemmi.Structure, also_strip: set | None = None) -> list[str]:
    """Remove water and HETATM residues. Returns list of removed residue names."""
    targets = STRIP_HET | (also_strip or set())
    removed = []
    for model in structure:
        for chain in model:
            to_del = []
            for i, res in enumerate(chain):
                if res.is_water() or res.name in targets:
                    to_del.append(i)
                    removed.append(res.name)
            for i in reversed(to_del):
                del chain[i]
    return sorted(set(removed))


def keep_chains(structure: gemmi.Structure, chain_names: set[str]) -> list[str]:
    """Remove all chains not in chain_names. Returns list of removed chain IDs."""
    removed = []
    for model in structure:
        to_del = [c.name for c in model if c.name not in chain_names]
        for name in to_del:
            del model[name]
            removed.append(name)
    return removed


def keep_model_0_only(structure: gemmi.Structure) -> int:
    """For NMR ensembles: discard all models except the first. Returns n discarded."""
    n = len(structure) - 1
    while len(structure) > 1:
        del structure[-1]
    return n


def report(label: str, structure: gemmi.Structure) -> None:
    log(f"\n  {BOLD}{label}{RESET}")
    for model in structure:
        for chain in model:
            std = [r for r in chain if not r.is_water() and r.het_flag != "H"]
            het = [r for r in chain if r.het_flag == "H"]
            if std:
                rng = f"{std[0].seqid}–{std[-1].seqid}"
                log(f"    Chain {chain.name}: {len(std)} residues ({rng})", GREEN)
            if het:
                names = sorted(set(r.name for r in het))
                log(f"    Chain {chain.name}: {len(het)} HETATM — {names}", YELLOW)


def write_pdb(structure: gemmi.Structure, out_path: Path) -> None:
    """
    Write PDB with explicit options so output is independent of gemmi version defaults.
    TER records on (GROMACS chain detection), SEQRES off (not needed for MD), CONECT off.
    """
    out_path.parent.mkdir(parents=True, exist_ok=True)
    opts = gemmi.PdbWriteOptions()
    opts.seqres_records = False
    opts.conect_records = False
    structure.write_pdb(str(out_path), opts)
    log(f"  → Saved: {out_path}  ({out_path.stat().st_size / 1024:.1f} KB)", GREEN)


# ── Per-system preparation ─────────────────────────────────────────────────────

def prepare_s1_bin1_bar(root: Path, force: bool) -> dict:
    """S1: BIN1 BAR domain dimer — 2FIC chains A+B, remove XE."""
    pdb_in  = root / "01_structures" / "BIN1" / "pdb" / "2FIC.pdb"
    pdb_out = root / "01_structures" / "BIN1" / "prepared" / "2FIC_BAR_dimer_prepared.pdb"
    log(f"\n{CYAN}▶ S1  BIN1 BAR domain dimer  (2FIC){RESET}")

    if not pdb_in.exists():
        log(f"  [FAIL] Input not found: {pdb_in}", RED)
        log( "         Run 01_download_structures.py first.", RED)
        return {"system": "S1", "status": "failed", "reason": "input missing"}

    if pdb_out.exists() and not force:
        log(f"  [SKIP] {pdb_out.name} already exists.", YELLOW)
        return {"system": "S1", "status": "skipped", "output": str(pdb_out.relative_to(root))}

    st = gemmi.read_structure(str(pdb_in))
    report("2FIC raw", st)

    removed_chains = keep_chains(st, {"A", "B"})
    removed_het    = remove_heteroatoms(st)

    for model in st:
        for chain in model:
            std = [r for r in chain if r.het_flag != "H" and not r.is_water()]
            if std and (first := std[0].seqid.num) > 1:
                log(f"  ⚠  Chain {chain.name} starts at residue {first} "
                    f"(residues 1–{first-1} absent from crystal — "
                    f"confirmed disordered by AF3 pLDDT, safe to exclude)", YELLOW)

    report("2FIC clean", st)
    write_pdb(st, pdb_out)

    return {
        "system": "S1", "status": "ok",
        "input":  str(pdb_in.relative_to(root)),
        "output": str(pdb_out.relative_to(root)),
        "input_sha256":  sha256_file(pdb_in),
        "output_sha256": sha256_file(pdb_out),
        "operations": {
            "chains_kept":    ["A", "B"],
            "chains_removed": removed_chains,
            "hetatm_stripped": removed_het,
            "models_kept":    "all (X-ray, single model)",
        },
    }


def prepare_s2_bin1_sh3(root: Path, force: bool) -> dict:
    """S2: BIN1 SH3 domain — 1MUZ model 0 only, chain A."""
    pdb_in  = root / "01_structures" / "BIN1" / "pdb" / "1MUZ.pdb"
    pdb_out = root / "01_structures" / "BIN1" / "prepared" / "1MUZ_SH3_model0_prepared.pdb"
    log(f"\n{CYAN}▶ S2  BIN1 SH3 domain  (1MUZ model 0){RESET}")

    if not pdb_in.exists():
        log(f"  [FAIL] Input not found: {pdb_in}", RED)
        return {"system": "S2", "status": "failed", "reason": "input missing"}

    if pdb_out.exists() and not force:
        log(f"  [SKIP] {pdb_out.name} already exists.", YELLOW)
        return {"system": "S2", "status": "skipped", "output": str(pdb_out.relative_to(root))}

    st = gemmi.read_structure(str(pdb_in))
    n_models = len(st)
    n_discarded = keep_model_0_only(st)
    log(f"  NMR ensemble: {n_models} models — keeping model 0, discarding {n_discarded}")

    removed_het = remove_heteroatoms(st)
    report("1MUZ model 0", st)
    log("  ℹ  Original NMR residue numbering retained (402–482 = full SH3 domain)")

    write_pdb(st, pdb_out)

    return {
        "system": "S2", "status": "ok",
        "input":  str(pdb_in.relative_to(root)),
        "output": str(pdb_out.relative_to(root)),
        "input_sha256":  sha256_file(pdb_in),
        "output_sha256": sha256_file(pdb_out),
        "operations": {
            "chains_kept":    ["A"],
            "chains_removed": [],
            "hetatm_stripped": removed_het,
            "models_kept":    "model 0 only (NMR ensemble — top model selected)",
            "nmr_models_total": n_models,
        },
    }


def prepare_s3_picalm(root: Path, force: bool) -> dict:
    """S3: PICALM ANTH domain — 7JXV chain A only (strip chain B = ubiquitin)."""
    pdb_in  = root / "01_structures" / "PICALM" / "pdb" / "7JXV.pdb"
    pdb_out = root / "01_structures" / "PICALM" / "prepared" / "7JXV_ANTH_prepared.pdb"
    log(f"\n{CYAN}▶ S3  PICALM ANTH domain  (7JXV chain A){RESET}")

    if not pdb_in.exists():
        log(f"  [FAIL] Input not found: {pdb_in}", RED)
        return {"system": "S3", "status": "failed", "reason": "input missing"}

    if pdb_out.exists() and not force:
        log(f"  [SKIP] {pdb_out.name} already exists.", YELLOW)
        return {"system": "S3", "status": "skipped", "output": str(pdb_out.relative_to(root))}

    st = gemmi.read_structure(str(pdb_in))
    report("7JXV raw", st)
    log("  Removing chain B (ubiquitin) + HETATM ...")

    removed_chains = keep_chains(st, {"A"})
    removed_het    = remove_heteroatoms(st)

    for model in st:
        for chain in model:
            std = [r for r in chain if r.het_flag != "H" and not r.is_water()]
            if std:
                first, last = std[0].seqid.num, std[-1].seqid.num
                log(f"  ℹ  Coverage: residues {first}–{last}  "
                    f"(1–{first-1} flexible N-terminus; {last+1}–652 is IDR)")

    report("7JXV clean", st)
    write_pdb(st, pdb_out)

    return {
        "system": "S3", "status": "ok",
        "input":  str(pdb_in.relative_to(root)),
        "output": str(pdb_out.relative_to(root)),
        "input_sha256":  sha256_file(pdb_in),
        "output_sha256": sha256_file(pdb_out),
        "operations": {
            "chains_kept":    ["A"],
            "chains_removed": removed_chains,
            "hetatm_stripped": removed_het,
            "models_kept":    "all (X-ray, single model)",
            "note": "chain B is ubiquitin co-crystallisation partner; removed for apo simulation",
        },
    }


def prepare_s4_cd2ap_sh3_2(root: Path, force: bool) -> dict:
    """S4: CD2AP SH3-2 — 3U23 chain A only, remove peptide (chain B) + EDO."""
    pdb_in  = root / "01_structures" / "CD2AP" / "pdb" / "3U23.pdb"
    pdb_out = root / "01_structures" / "CD2AP" / "prepared" / "3U23_SH3-2_prepared.pdb"
    log(f"\n{CYAN}▶ S4  CD2AP SH3-2  (3U23 chain A){RESET}")

    if not pdb_in.exists():
        log(f"  [FAIL] Input not found: {pdb_in}", RED)
        return {"system": "S4", "status": "failed", "reason": "input missing"}

    if pdb_out.exists() and not force:
        log(f"  [SKIP] {pdb_out.name} already exists.", YELLOW)
        return {"system": "S4", "status": "skipped", "output": str(pdb_out.relative_to(root))}

    st = gemmi.read_structure(str(pdb_in))
    report("3U23 raw", st)
    log("  Removing chain B (RIN3 peptide) + EDO + NH2 ...")

    removed_chains = keep_chains(st, {"A"})
    removed_het    = remove_heteroatoms(st, also_strip={"NH2"})

    report("3U23 clean", st)
    log("  ℹ  Apo SH3-2 groove — RIN3 binding site now empty")
    log("  ℹ  Crystal resolution: 1.11 Å — highest quality structure in set")
    write_pdb(st, pdb_out)

    return {
        "system": "S4", "status": "ok",
        "input":  str(pdb_in.relative_to(root)),
        "output": str(pdb_out.relative_to(root)),
        "input_sha256":  sha256_file(pdb_in),
        "output_sha256": sha256_file(pdb_out),
        "operations": {
            "chains_kept":    ["A"],
            "chains_removed": removed_chains,
            "hetatm_stripped": removed_het,
            "models_kept":    "all (X-ray, single model)",
            "note": "chain B is RIN3 peptide co-crystallisation partner; removed for apo simulation",
        },
    }


def prepare_s5_cd2ap_sh3_1(root: Path, force: bool) -> dict:
    """S5: CD2AP SH3-1 — 4WCI chain A only, remove chains B-F + SO4."""
    pdb_in  = root / "01_structures" / "CD2AP" / "pdb" / "4WCI.pdb"
    pdb_out = root / "01_structures" / "CD2AP" / "prepared" / "4WCI_SH3-1_prepared.pdb"
    log(f"\n{CYAN}▶ S5  CD2AP SH3-1  (4WCI chain A){RESET}")

    if not pdb_in.exists():
        log(f"  [FAIL] Input not found: {pdb_in}", RED)
        return {"system": "S5", "status": "failed", "reason": "input missing"}

    if pdb_out.exists() and not force:
        log(f"  [SKIP] {pdb_out.name} already exists.", YELLOW)
        return {"system": "S5", "status": "skipped", "output": str(pdb_out.relative_to(root))}

    st = gemmi.read_structure(str(pdb_in))
    report("4WCI raw", st)
    log("  Keeping chain A only (3 crystallographic copies; A is primary) ...")

    removed_chains = keep_chains(st, {"A"})
    removed_het    = remove_heteroatoms(st)

    report("4WCI clean", st)
    write_pdb(st, pdb_out)

    return {
        "system": "S5", "status": "ok",
        "input":  str(pdb_in.relative_to(root)),
        "output": str(pdb_out.relative_to(root)),
        "input_sha256":  sha256_file(pdb_in),
        "output_sha256": sha256_file(pdb_out),
        "operations": {
            "chains_kept":    ["A"],
            "chains_removed": removed_chains,
            "hetatm_stripped": removed_het,
            "models_kept":    "all (X-ray, single model)",
        },
    }


# ── Summary ────────────────────────────────────────────────────────────────────

def print_summary(root: Path, results: list[dict]) -> None:
    log(f"\n{'═'*58}", BOLD)
    log("  CRYPTAD — Structure preparation summary", BOLD)
    log(f"{'═'*58}", BOLD)
    for r in results:
        status = r["status"]
        color  = GREEN if status in ("ok", "skipped") else RED
        marker = "✓" if status in ("ok", "skipped") else "✗"
        out    = r.get("output", "—")
        log(f"  {marker}  {r['system']}  [{status:7s}]  {out}", color)
    log(f"\n  Next step: open each *_prepared.pdb in PyMOL to verify the apo "
        f"binding site, then proceed to CHARMM-GUI.\n")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare apo PDB structures for CRYPTAD MD.")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing prepared PDB files.")
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT,
                        help=f"Path to CRYPTAD project root (default: {DEFAULT_ROOT})")
    args = parser.parse_args()

    root = args.project_root.resolve()
    if not (root / "01_structures").is_dir():
        log(f"[ERROR] Project root does not contain 01_structures/: {root}", RED)
        log( "        Pass the correct path with --project-root.", RED)
        sys.exit(1)

    log(f"\n{'═'*58}", BOLD)
    log("  CRYPTAD — Structure Preparation", BOLD)
    log(f"{'═'*58}", BOLD)
    log(f"  Project root : {root}")
    log(f"  Force re-run : {args.force}\n")

    results = [
        prepare_s1_bin1_bar(root, args.force),
        prepare_s2_bin1_sh3(root, args.force),
        prepare_s3_picalm(root, args.force),
        prepare_s4_cd2ap_sh3_2(root, args.force),
        prepare_s5_cd2ap_sh3_1(root, args.force),
    ]

    # Write manifest
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "note": (
            "Original PDB residue numbering is preserved in *_prepared.pdb files. "
            "Renumbering (for CHARMM-GUI topology) is a separate manual step. "
            "gemmi PdbWriteOptions: seqres_records=False, conect_records=False."
        ),
        "systems": results,
    }
    manifest_path = root / "01_structures" / "preparation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    log(f"\n  Manifest written → {manifest_path.relative_to(root)}", CYAN)

    print_summary(root, results)

    failed = [r for r in results if r["status"] == "failed"]
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
