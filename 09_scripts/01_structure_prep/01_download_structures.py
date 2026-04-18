"""
CRYPTAD — PDB Structure Downloader
Downloads selected PDB files into the correct project directories.

Reproducibility notes
---------------------
SHA-256 checksums are hardcoded for the exact depositions used in this study.
If RCSB updates an entry the checksum will fail and you will be told explicitly,
rather than silently using a different structure.  A download manifest
(download_manifest.json) is written on success for audit purposes.

Usage:
    python 01_download_structures.py [--force] [--project-root PATH]

Options:
    --force           Re-download even if file already exists (re-verifies checksum).
    --project-root    Explicit path to CRYPTAD project root.
                      Default: inferred as three directories above this script.
"""

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# ── Path resolution ────────────────────────────────────────────────────────────
# Script lives at <root>/09_scripts/01_structure_prep/01_download_structures.py
# parents[0] = 01_structure_prep/
# parents[1] = 09_scripts/
# parents[2] = project root
_SCRIPT_DIR = Path(__file__).resolve()
DEFAULT_ROOT = _SCRIPT_DIR.parents[2]

# ── RCSB URL ───────────────────────────────────────────────────────────────────
RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
TIMEOUT_SEC = 60

# ── Target structures with SHA-256 checksums ───────────────────────────────────
# Checksums computed from the exact depositions used in this study (Apr 2026).
# If RCSB re-refines an entry the checksum will differ — download the original
# from the RCSB versioned archive or from the Zenodo data deposit for this paper.
TARGETS = {
    "BIN1": {
        "subdir": Path("01_structures") / "BIN1" / "pdb",
        "structures": {
            "2FIC": {
                "description": "BAR domain (X-ray, 1.99 Å) — PRIMARY for MD",
                "sha256": "f0d6aa7cc681b79077475a63570158815be080c9987d588946950280933e2e7d",
                "size_bytes": 336717,
            },
            "1MUZ": {
                "description": "SH3 domain (NMR) — reference/validation",
                "sha256": "cfd27c2b9e96e146ff9a81be145b1774a259d4ab41aed62f8e0c1a9c02882261",
                "size_bytes": 2108349,
            },
            "1MV0": {
                "description": "SH3 domain + c-Myc peptide (NMR) — binding groove reference",
                "sha256": "3c55373d56275f52e0d95b43fe5bdf0d906eff1543af507654da02342d752dce",
                "size_bytes": 2459565,
            },
        },
    },
    "PICALM": {
        "subdir": Path("01_structures") / "PICALM" / "pdb",
        "structures": {
            "7JXV": {
                "description": "ANTH domain + Ubiquitin (X-ray, 2.35 Å) — PRIMARY for MD",
                "sha256": "83c3414aea33ffca806e0dfbf007a593931644b7c772cd40ecabdff4d5534c12",
                "size_bytes": 490374,
            },
        },
    },
    "CD2AP": {
        "subdir": Path("01_structures") / "CD2AP" / "pdb",
        "structures": {
            "3U23": {
                "description": "SH3-2 + RIN3 peptide (X-ray, 1.11 Å) — PRIMARY for MD",
                "sha256": "ef97df63cafbd844dce6a6a312cf77e3b98ca132d43ede75ae3a882986af5165",
                "size_bytes": 145557,
            },
            "4WCI": {
                "description": "SH3-1 + RIN3 peptide (X-ray, 1.65 Å) — secondary",
                "sha256": "7ca0c7e2adf29b2259fef81f5c81eb72eabb89ecc8560881c1878d1539beb61b",
                "size_bytes": 363528,
            },
            "2J6F": {
                "description": "SH3-1 + CBL-B peptide (X-ray, 1.70 Å) — reference",
                "sha256": "0a3d274fe9ed5ad144b7e07116b72f779e3e39b7e704e081757f069c287c5684",
                "size_bytes": 86589,
            },
        },
    },
}

# ── Colours ────────────────────────────────────────────────────────────────────
RESET  = "\033[0m"
GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
BOLD   = "\033[1m"
CYAN   = "\033[96m"

def log(msg, color=RESET):
    print(f"{color}{msg}{RESET}", flush=True)

# ── Checksum ───────────────────────────────────────────────────────────────────

def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def verify(path: Path, expected_sha256: str, expected_size: int) -> tuple[bool, str]:
    """Return (ok, reason). reason is empty string on success."""
    if not path.exists():
        return False, "file does not exist"
    actual_size = path.stat().st_size
    if actual_size != expected_size:
        return False, f"size mismatch: {actual_size} B vs expected {expected_size} B"
    actual_hash = sha256_file(path)
    if actual_hash != expected_sha256:
        return False, f"SHA-256 mismatch:\n    got:      {actual_hash}\n    expected: {expected_sha256}"
    return True, ""

# ── Download ───────────────────────────────────────────────────────────────────

def download_pdb(pdb_id: str, dest_dir: Path, meta: dict, force: bool, retries: int = 3) -> bool:
    dest_dir.mkdir(parents=True, exist_ok=True)
    url      = RCSB_URL.format(pdb_id=pdb_id.upper())
    out_path = dest_dir / f"{pdb_id.upper()}.pdb"

    # If file exists and not forced, verify before skipping
    if out_path.exists() and not force:
        ok, reason = verify(out_path, meta["sha256"], meta["size_bytes"])
        if ok:
            log(f"  [SKIP] {pdb_id}.pdb — already present, checksum verified.", YELLOW)
            return True
        else:
            log(f"  [WARN] {pdb_id}.pdb exists but failed verification ({reason}). Re-downloading.", YELLOW)
            out_path.unlink()

    req = urllib.request.Request(url, headers={"User-Agent": "CRYPTAD/1.0 (research pipeline)"})

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT_SEC) as response:
                data = response.read()
            out_path.write_bytes(data)

            ok, reason = verify(out_path, meta["sha256"], meta["size_bytes"])
            if not ok:
                out_path.unlink()
                log(f"  [FAIL] {pdb_id}.pdb downloaded but failed verification: {reason}", RED)
                log(f"         The RCSB entry may have been updated since this study was run.", RED)
                log(f"         Retrieve the original from the Zenodo data deposit for this paper.", RED)
                return False

            size_kb = len(data) / 1024
            log(f"  [OK]   {pdb_id}.pdb  ({size_kb:.1f} KB)  — {meta['description']}", GREEN)
            return True

        except urllib.error.URLError as e:
            if attempt < retries:
                log(f"  [WARN] Attempt {attempt} failed for {pdb_id}: {e}. Retrying in 2 s...", YELLOW)
                time.sleep(2)
            else:
                log(f"  [FAIL] Could not download {pdb_id} after {retries} attempts: {e}", RED)
                if out_path.exists():
                    out_path.unlink()
                return False

# ── Manifest ───────────────────────────────────────────────────────────────────

def write_manifest(project_root: Path, results: list[dict]) -> None:
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "rcsb_url_template": RCSB_URL,
        "note": (
            "SHA-256 checksums match the depositions used in the CRYPTAD study (Apr 2026). "
            "If a checksum fails, retrieve the original file from the Zenodo data deposit."
        ),
        "files": results,
    }
    out = project_root / "01_structures" / "download_manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    log(f"\n  Manifest written → {out.relative_to(project_root)}", CYAN)

# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Download PDB structures for CRYPTAD.")
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if file already exists.")
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT,
                        help=f"Path to CRYPTAD project root (default: {DEFAULT_ROOT})")
    args = parser.parse_args()

    root = args.project_root.resolve()
    if not (root / "01_structures").is_dir():
        log(f"[ERROR] Project root does not contain 01_structures/: {root}", RED)
        log( "        Pass the correct path with --project-root.", RED)
        sys.exit(1)

    log("\n" + "═" * 60, BOLD)
    log("  CRYPTAD — PDB Structure Download", BOLD)
    log("═" * 60, BOLD)
    log(f"  Project root : {root}")
    log(f"  Force re-dl  : {args.force}\n")

    results = []
    total, succeeded, failed_ids = 0, 0, []

    for target, config in TARGETS.items():
        log(f"\n{CYAN}▶ {target}{RESET}")
        dest_dir = root / config["subdir"]
        log(f"  → {dest_dir}")

        for pdb_id, meta in config["structures"].items():
            total += 1
            ok = download_pdb(pdb_id, dest_dir, meta, force=args.force)
            path = dest_dir / f"{pdb_id.upper()}.pdb"
            results.append({
                "pdb_id":       pdb_id,
                "target":       target,
                "description":  meta["description"],
                "expected_sha256": meta["sha256"],
                "expected_size_bytes": meta["size_bytes"],
                "status":       "ok" if ok else "failed",
                "path":         str(path.relative_to(root)),
            })
            if ok:
                succeeded += 1
            else:
                failed_ids.append(f"{target}/{pdb_id}")

    log("\n" + "─" * 60)
    log(f"  Done: {succeeded}/{total} verified successfully.", BOLD)

    if succeeded == total:
        write_manifest(root, results)
        log("  All structures downloaded and checksum-verified.\n", GREEN)
    else:
        log(f"\n  Failed: {', '.join(failed_ids)}", RED)
        log( "  Retrieve originals from the Zenodo data deposit for this paper.", RED)
        sys.exit(1)

if __name__ == "__main__":
    main()
