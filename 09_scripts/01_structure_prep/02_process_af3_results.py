"""
CRYPTAD — AlphaFold3 Results Processor
Extracts AF3 zip files, converts best model CIF → PDB, plots pLDDT per residue,
and records model confidence metadata.

Reproducibility notes
---------------------
- Only model_0 (AF3's top-ranked prediction) is converted and used downstream.
  Model selection is explicit by filename pattern, not alphabetical sort order.
- AF3 confidence metadata (iptm, ptm, ranking_score, fraction_disordered, seed)
  is extracted from summary_confidences_0.json and written to af3_manifest.json
  for audit purposes.
- pLDDT figures are written to 06_figures/publication/ at 300 Dpi to match all
  other publication figures in this project.
- Zip extraction is skipped if model CIF files already exist (--force to override).

Usage:
    python 02_process_af3_results.py [--force] [--project-root PATH]

Options:
    --force           Re-extract zips and re-convert even if outputs already exist.
    --project-root    Explicit path to CRYPTAD project root.
                      Default: inferred as three directories above this script.
"""

import argparse
import glob
import json
import os
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import gemmi
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

# ── Path resolution ────────────────────────────────────────────────────────────
# Script lives at <root>/09_scripts/01_structure_prep/02_process_af3_results.py
# parents[0] = 01_structure_prep/
# parents[1] = 09_scripts/
# parents[2] = project root
_SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_ROOT = _SCRIPT_PATH.parents[2]

# ── Targets ────────────────────────────────────────────────────────────────────
TARGETS = ["BIN1", "PICALM", "CD2AP"]

# ── Colours ────────────────────────────────────────────────────────────────────
RESET = "\033[0m"; GREEN = "\033[92m"; RED = "\033[91m"
YELLOW = "\033[93m"; BOLD = "\033[1m"; CYAN = "\033[96m"

def log(msg, c=RESET):
    print(f"{c}{msg}{RESET}", flush=True)


# ── Step 1: Extract zip ────────────────────────────────────────────────────────

def extract_zip(af_dir: Path, force: bool) -> bool:
    """Extract the AF3 zip into af_dir. Skip if model CIFs already present."""
    model_cifs = sorted(af_dir.glob("*model_0*.cif"))
    if model_cifs and not force:
        log(f"  [SKIP] CIF files already extracted ({len(list(af_dir.glob('*model*.cif')))} models found).", YELLOW)
        return True

    zips = sorted(af_dir.glob("*.zip"))
    if not zips:
        log(f"  [FAIL] No zip file found in {af_dir}", RED)
        log(f"         Place the AF3 download zip here and re-run.", RED)
        return False
    if len(zips) > 1:
        log(f"  [WARN] Multiple zips found; using {zips[0].name}", YELLOW)

    z = zips[0]
    log(f"  Extracting {z.name} ...", CYAN)
    with zipfile.ZipFile(z, "r") as zf:
        zf.extractall(af_dir)
    log(f"  Extracted {len(zf.namelist())} files.", GREEN)
    return True


# ── Step 2: CIF → PDB (model_0 only) ──────────────────────────────────────────

def cif_to_pdb(target: str, af_dir: Path, force: bool) -> tuple[Path, Path] | None:
    """
    Convert model_0 CIF → PDB. Returns (pdb_path, cif_path) or None on failure.
    Explicitly selects model_0 by filename to avoid sort-order ambiguity.
    """
    out_pdb = af_dir / f"{target}_af3_model.pdb"

    # Find model_0 CIF explicitly
    candidates = sorted(af_dir.glob("*model_0*.cif"))
    if not candidates:
        # Fallback: any model CIF, warn loudly
        candidates = sorted(af_dir.glob("*model*.cif"))
        if not candidates:
            log(f"  [FAIL] No model CIF found for {target}", RED)
            return None
        log(f"  [WARN] model_0 CIF not found; falling back to {candidates[0].name}", YELLOW)

    chosen = candidates[0]

    if out_pdb.exists() and not force:
        log(f"  [SKIP] {out_pdb.name} already exists.", YELLOW)
        return out_pdb, chosen

    log(f"  Converting {chosen.name} → {out_pdb.name} ...", CYAN)
    structure = gemmi.read_structure(str(chosen))
    structure.name = target

    # Write only the first model (index 0) — AF3 CIFs contain a single model
    opts = gemmi.PdbWriteOptions()
    structure.write_pdb(str(out_pdb), opts)
    log(f"  Saved: {out_pdb.name}", GREEN)
    return out_pdb, chosen


# ── Step 3: pLDDT extraction ───────────────────────────────────────────────────

def extract_plddt_from_cif(cif_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    AF3 stores per-residue confidence (pLDDT) in the B-factor column.
    Returns (residue_numbers, plddt_values).
    """
    st = gemmi.read_structure(str(cif_path))
    residues, plddt = [], []
    for model in st:
        for chain in model:
            for res in chain:
                if res.is_water():
                    continue
                for atom in res:
                    if atom.name == "CA":
                        residues.append(res.seqid.num)
                        plddt.append(atom.b_iso)
                        break
    return np.array(residues), np.array(plddt)


# ── Step 4: AF3 confidence metadata ───────────────────────────────────────────

def load_af3_metadata(af_dir: Path) -> dict:
    """
    Extract confidence scores from summary_confidences_0.json and
    random seed from job_request.json.  Returns a flat metadata dict.
    """
    meta = {}

    conf_files = sorted(af_dir.glob("*summary_confidences_0*.json"))
    if conf_files:
        d = json.loads(conf_files[0].read_text())
        meta["iptm"]               = d.get("iptm")
        meta["ptm"]                = d.get("ptm")
        meta["ranking_score"]      = d.get("ranking_score")
        meta["fraction_disordered"]= d.get("fraction_disordered")
        meta["has_clash"]          = d.get("has_clash")
        meta["num_recycles"]       = d.get("num_recycles")
    else:
        log(f"  [WARN] summary_confidences_0.json not found — metadata incomplete.", YELLOW)

    req_files = sorted(af_dir.glob("*job_request*.json"))
    if req_files:
        d = json.loads(req_files[0].read_text())
        if isinstance(d, list) and d:
            d = d[0]
        meta["job_name"]     = d.get("name")
        meta["model_seeds"]  = d.get("modelSeeds")
    else:
        log(f"  [WARN] job_request.json not found — seed not recorded.", YELLOW)

    return meta


# ── Step 5: pLDDT plot ────────────────────────────────────────────────────────

def plot_plddt(target: str, residues: np.ndarray, plddt: np.ndarray,
               out_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(12, 3.5))

    # AF3 confidence colour bands (matching the AlphaFold Protein Structure Database viewer)
    ax.axhspan(90, 100, alpha=0.15, color="#0053D6", label="Very high (>90)")
    ax.axhspan(70,  90, alpha=0.15, color="#65CBF3", label="Confident (70–90)")
    ax.axhspan(50,  70, alpha=0.15, color="#FFDB13", label="Low (50–70)")
    ax.axhspan( 0,  50, alpha=0.15, color="#FF7D45", label="Very low (<50)")

    ax.plot(residues, plddt, color="#333333", lw=1.2, zorder=3)
    ax.axhline(70, color="steelblue", lw=0.8, ls="--", alpha=0.6)
    ax.axhline(50, color="orange",    lw=0.8, ls="--", alpha=0.6)

    ax.set_xlim(residues[0], residues[-1])
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue number", fontsize=11)
    ax.set_ylabel("pLDDT", fontsize=11)
    ax.set_title(f"{target} — AlphaFold3 per-residue confidence (pLDDT)",
                 fontsize=12, fontweight="bold")

    # Shade contiguous stretches with pLDDT < 50
    low_mask = plddt < 50
    if low_mask.any():
        runs, in_run, start = [], False, 0
        for i, v in enumerate(low_mask):
            if v and not in_run:
                in_run, start = True, residues[i]
            elif not v and in_run:
                in_run = False
                runs.append((start, residues[i - 1]))
        if in_run:
            runs.append((start, residues[-1]))
        for s, e in runs:
            ax.axvspan(s, e, alpha=0.25, color="red", zorder=2)
            ax.text((s + e) / 2, 5, f"{s}–{e}", ha="center", fontsize=7, color="darkred")

    legend_patches = [
        mpatches.Patch(color="#0053D6", alpha=0.4, label="Very high >90"),
        mpatches.Patch(color="#65CBF3", alpha=0.4, label="Confident 70–90"),
        mpatches.Patch(color="#FFDB13", alpha=0.4, label="Low 50–70"),
        mpatches.Patch(color="#FF7D45", alpha=0.4, label="Very low <50"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=8, framealpha=0.9)

    plt.tight_layout()
    out_path = out_dir / f"{target}_pLDDT.png"
    plt.savefig(out_path, dpi=300)   # 300 DPI — consistent with all publication figures
    plt.close()
    log(f"  pLDDT plot saved: {out_path.name} (300 dpi)", GREEN)
    return out_path


def summarise_plddt(target: str, residues: np.ndarray, plddt: np.ndarray) -> dict:
    log(f"\n  {BOLD}pLDDT summary for {target}{RESET}")
    log(f"  Mean pLDDT   : {plddt.mean():.1f}")
    log(f"  Residues >90 : {(plddt > 90).sum()} ({100*(plddt > 90).mean():.0f}%)")
    log(f"  Residues >70 : {(plddt > 70).sum()} ({100*(plddt > 70).mean():.0f}%)")
    log(f"  Residues <50 : {(plddt < 50).sum()} ({100*(plddt < 50).mean():.0f}%)")

    runs, in_run, start = [], False, 0
    low_mask = plddt < 50
    for i, v in enumerate(low_mask):
        if v and not in_run:
            in_run, start = True, i
        elif not v and in_run:
            in_run = False
            if i - start >= 10:
                runs.append((int(residues[start]), int(residues[i - 1])))
    if in_run and (len(low_mask) - start) >= 10:
        runs.append((int(residues[start]), int(residues[-1])))

    if runs:
        log(f"  ⚠  Disordered stretches (≥10 residues pLDDT<50):", YELLOW)
        for s, e in runs:
            log(f"     → residues {s}–{e}", YELLOW)
    else:
        log(f"  ✓  No long disordered stretches.", GREEN)

    return {
        "mean_plddt":        float(plddt.mean()),
        "n_residues":        int(len(plddt)),
        "frac_above_90":     float((plddt > 90).mean()),
        "frac_above_70":     float((plddt > 70).mean()),
        "frac_below_50":     float((plddt < 50).mean()),
        "disordered_runs":   runs,
    }


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Process AlphaFold3 results for CRYPTAD.")
    parser.add_argument("--force", action="store_true",
                        help="Re-extract and re-convert even if outputs already exist.")
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT,
                        help=f"Path to CRYPTAD project root (default: {DEFAULT_ROOT})")
    args = parser.parse_args()

    root = args.project_root.resolve()
    if not (root / "01_structures").is_dir():
        log(f"[ERROR] Project root does not contain 01_structures/: {root}", RED)
        log( "        Pass the correct path with --project-root.", RED)
        sys.exit(1)

    figures_dir = root / "06_figures" / "publication"
    figures_dir.mkdir(parents=True, exist_ok=True)

    log(f"\n{'═'*60}", BOLD)
    log("  CRYPTAD — AF3 Results Processor", BOLD)
    log(f"{'═'*60}", BOLD)
    log(f"  Project root : {root}")
    log(f"  Figures out  : {figures_dir.relative_to(root)}")
    log(f"  Force re-run : {args.force}\n")

    manifest_entries = []
    all_ok = True

    for target in TARGETS:
        log(f"\n{CYAN}▶ {target}{RESET}")
        af_dir = root / "01_structures" / target / "alphafold"

        if not af_dir.is_dir():
            log(f"  [FAIL] Directory not found: {af_dir}", RED)
            all_ok = False
            continue

        if not extract_zip(af_dir, args.force):
            all_ok = False
            continue

        result = cif_to_pdb(target, af_dir, args.force)
        if result is None:
            all_ok = False
            continue

        pdb_path, cif_path = result
        residues, plddt   = extract_plddt_from_cif(cif_path)
        plddt_stats        = summarise_plddt(target, residues, plddt)
        af3_meta           = load_af3_metadata(af_dir)
        plot_plddt(target, residues, plddt, figures_dir)

        manifest_entries.append({
            "target":       target,
            "cif_source":   cif_path.name,
            "pdb_output":   pdb_path.name,
            "af3_metadata": af3_meta,
            "plddt_stats":  plddt_stats,
        })

    # Write manifest
    if manifest_entries:
        manifest = {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "note": (
                "AF3 model_0 selected for each target. "
                "Confidence scores from summary_confidences_0.json. "
                "Seeds from job_request.json."
            ),
            "targets": manifest_entries,
        }
        manifest_path = root / "01_structures" / "af3_manifest.json"
        manifest_path.write_text(json.dumps(manifest, indent=2))
        log(f"\n  Manifest written → {manifest_path.relative_to(root)}", CYAN)

    log(f"\n{'─'*60}")
    if all_ok:
        log("  All AF3 models processed.", GREEN)
        log(f"  pLDDT plots → {figures_dir.relative_to(root)}/", GREEN)
        log(f"  PDB files   → 01_structures/<TARGET>/alphafold/<TARGET>_af3_model.pdb\n", GREEN)
    else:
        log("  Some targets need attention (see above).\n", YELLOW)
        sys.exit(1)


if __name__ == "__main__":
    main()
