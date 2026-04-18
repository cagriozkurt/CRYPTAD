"""
zenodo_upload.py — CRYPTAD trajectory archive uploader
=======================================================
Uploads MD trajectories to a Zenodo deposit using the bucket (streaming) API,
which handles files of any size without memory issues.

Usage:
    # Dry run — list all files without uploading:
    python3 zenodo_upload.py --dry-run

    # Create a new draft deposit and upload:
    python3 zenodo_upload.py --token YOUR_TOKEN

    # Add files to an existing draft deposit:
    python3 zenodo_upload.py --token YOUR_TOKEN --deposit-id 12345678

    # Test against Zenodo sandbox first (recommended):
    python3 zenodo_upload.py --token SANDBOX_TOKEN --sandbox

Token:
    Generate at https://zenodo.org/account/settings/applications/tokens/new/
    Scopes needed: deposit:write, deposit:actions
    Or set env var: export ZENODO_TOKEN=...

Notes:
    - Already-uploaded files are skipped automatically (checked by filename).
    - Large files (e.g. 13 GB meta.xtc) are streamed in 32 MB chunks.
    - On failure each file is retried up to 3 times before aborting.
    - A JSON manifest is written to zenodo_upload_manifest.json on completion.
"""

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import requests

_script       = Path(__file__).resolve()
_default_root = _script.parent

CHUNK_SIZE  = 32 * 1024 * 1024   # 32 MB streaming chunks
MAX_RETRIES = 3
RETRY_WAIT  = 10                  # seconds between retries

API_BASE      = "https://zenodo.org/api"
SANDBOX_BASE  = "https://sandbox.zenodo.org/api"


# ── File manifest ─────────────────────────────────────────────────────────────
# (local_relative_path, remote_filename)
# local_relative_path is relative to project root.
# remote_filename is the flat name stored on Zenodo.

def build_file_list(root: Path) -> list[tuple[Path, str]]:
    files = []

    # ── 1. Metadynamics trajectories and inputs ───────────────────────────────
    metad_systems = {
        "BIN1":   ["S1_BAR", "S2_SH3"],
        "PICALM": ["S3_ANTH"],
        "CD2AP":  ["S4_SH3-2", "S5_SH3-1"],
    }
    metad_files = ["meta.xtc", "meta.tpr", "plumed.dat", "HILLS", "COLVAR"]

    for protein, systems in metad_systems.items():
        for sys_name in systems:
            base = root / "02_md_simulations" / protein / "metadynamics" / sys_name
            for fname in metad_files:
                p = base / fname
                if p.exists():
                    remote = f"metad_{sys_name}_{fname}"
                    files.append((p, remote))

    # ── 2. Complex MD — PBC-corrected trajectories + topology ─────────────────
    complex_keep = [
        "complexmd_fit_nw.xtc",
        "complexmd.tpr",
        "complexmd.gro",
        "ref_nw.pdb",
        "topol.top",
        "topol_flat.top",
        "LIG_GMX.itp",
        "lig_atomtypes.itp",
        "posre_LIG.itp",
        "index.ndx",
    ]
    complex_base = root / "02_md_simulations" / "complex_md"
    if complex_base.exists():
        for sys_dir in sorted(complex_base.iterdir()):
            if not sys_dir.is_dir():
                continue
            for rep_dir in sorted(sys_dir.iterdir()):
                if not rep_dir.is_dir() or not rep_dir.name.startswith("rep"):
                    continue
                for fname in complex_keep:
                    p = rep_dir / fname
                    if p.exists():
                        remote = f"complex_{sys_dir.name}_{rep_dir.name}_{fname}"
                        files.append((p, remote))

    # ── 3. Persistence gate trajectories ─────────────────────────────────────
    persist_proteins = ["BIN1", "PICALM"]
    persist_files    = ["persist.xtc", "persist.tpr"]
    for protein in persist_proteins:
        pg_base = root / "02_md_simulations" / protein / "persistence_gate"
        if not pg_base.exists():
            continue
        for rep_dir in sorted(pg_base.iterdir()):
            if not rep_dir.is_dir() or not rep_dir.name.startswith("rep"):
                continue
            for fname in persist_files:
                p = rep_dir / fname
                if p.exists():
                    remote = f"persist_{protein}_{rep_dir.name}_{fname}"
                    files.append((p, remote))

    return files


# ── Helpers ───────────────────────────────────────────────────────────────────

def md5(path: Path, chunk: int = CHUNK_SIZE) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            buf = f.read(chunk)
            if not buf:
                break
            h.update(buf)
    return h.hexdigest()


def human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


class ProgressStream:
    """Wraps a file object to print upload progress."""
    def __init__(self, path: Path):
        self._f    = open(path, "rb")
        self._size = path.stat().st_size
        self._done = 0
        self._last = 0

    def read(self, n: int = -1) -> bytes:
        buf = self._f.read(n)
        self._done += len(buf)
        pct = self._done / self._size * 100 if self._size else 100
        if pct - self._last >= 5 or not buf:
            print(f"    {pct:5.1f}%  {human_size(self._done)} / {human_size(self._size)}",
                  end="\r", flush=True)
            self._last = pct
        return buf

    def __len__(self) -> int:
        return self._size

    def close(self):
        self._f.close()


# ── Zenodo API calls ──────────────────────────────────────────────────────────

def api_headers(token: str) -> dict:
    return {"Authorization": f"Bearer {token}"}


def create_deposit(base: str, token: str) -> dict:
    r = requests.post(
        f"{base}/deposit/depositions",
        headers={**api_headers(token), "Content-Type": "application/json"},
        json={},
        timeout=30,
    )
    r.raise_for_status()
    return r.json()


def get_deposit(base: str, token: str, deposit_id: int) -> dict:
    r = requests.get(
        f"{base}/deposit/depositions/{deposit_id}",
        headers=api_headers(token),
        timeout=30,
    )
    r.raise_for_status()
    return r.json()


def list_uploaded(base: str, token: str, deposit_id: int) -> set[str]:
    r = requests.get(
        f"{base}/deposit/depositions/{deposit_id}/files",
        headers=api_headers(token),
        timeout=30,
    )
    r.raise_for_status()
    return {f["filename"] for f in r.json()}


def upload_file(bucket_url: str, token: str, local: Path, remote: str) -> dict:
    url = f"{bucket_url}/{remote}"
    stream = ProgressStream(local)
    try:
        r = requests.put(
            url,
            headers=api_headers(token),
            data=stream,
            timeout=None,
        )
        r.raise_for_status()
        print()   # newline after progress bar
        return r.json()
    finally:
        stream.close()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Upload CRYPTAD trajectories to Zenodo")
    parser.add_argument("--token",      default=os.environ.get("ZENODO_TOKEN"),
                        help="Zenodo API token (or set ZENODO_TOKEN env var)")
    parser.add_argument("--deposit-id", type=int, default=None,
                        help="Existing draft deposit ID to add files to")
    parser.add_argument("--sandbox",    action="store_true",
                        help="Use sandbox.zenodo.org instead of zenodo.org")
    parser.add_argument("--dry-run",    action="store_true",
                        help="List files to upload without uploading")
    parser.add_argument("--resume",     action="store_true",
                        help="Read manifest and skip already-uploaded files; "
                             "use --deposit-id to target a new deposit for remainder")
    parser.add_argument("--skip-also-in", type=int, default=None, metavar="DEPOSIT_ID",
                        help="Also skip files already present in this deposit ID "
                             "(use when splitting across two deposits)")
    parser.add_argument("--project-root", type=Path, default=_default_root)
    args = parser.parse_args()

    root = args.project_root.resolve()
    base = SANDBOX_BASE if args.sandbox else API_BASE

    # ── Build file list ───────────────────────────────────────────────────────
    files = build_file_list(root)
    if not files:
        print("No files found. Check --project-root.")
        sys.exit(1)

    total_bytes = sum(p.stat().st_size for p, _ in files)
    print(f"Files to upload: {len(files)}  |  Total: {human_size(total_bytes)}")
    print()

    if args.dry_run:
        for local, remote in files:
            size = human_size(local.stat().st_size)
            print(f"  {size:>10}  {remote}")
            print(f"            <- {local.relative_to(root)}")
        return

    if not args.token:
        print("ERROR: --token required (or set ZENODO_TOKEN env var)")
        sys.exit(1)

    # ── Load previous manifest (--resume) ─────────────────────────────────────
    manifest_path    = root / "zenodo_upload_manifest.json"
    previously_done: set[str] = set()
    manifest: list[dict]      = []

    if args.resume and manifest_path.exists():
        with open(manifest_path) as f:
            prev = json.load(f)
        for entry in prev.get("files", []):
            if entry["status"] in ("uploaded", "skipped"):
                previously_done.add(entry["remote"])
                manifest.append(entry)
        print(f"Resume: {len(previously_done)} file(s) already done per manifest — will skip.")
        print()

    # ── Get or create deposit ─────────────────────────────────────────────────
    if args.deposit_id:
        print(f"Using existing deposit {args.deposit_id} ...")
        deposit = get_deposit(base, args.token, args.deposit_id)
    else:
        print("Creating new draft deposit ...")
        deposit = create_deposit(base, args.token)

    deposit_id = deposit["id"]
    bucket_url = deposit["links"]["bucket"]
    edit_url   = deposit["links"]["html"]
    print(f"  Deposit ID : {deposit_id}")
    print(f"  Edit URL   : {edit_url}")
    print(f"  Bucket     : {bucket_url}")
    print()

    # Also skip files already present in this specific deposit
    in_deposit = list_uploaded(base, args.token, deposit_id)
    if in_deposit:
        print(f"  {len(in_deposit)} file(s) already in this deposit — will skip.")

    # Optionally cross-check a previous deposit (--skip-also-in)
    if args.skip_also_in:
        prev_deposit = list_uploaded(base, args.token, args.skip_also_in)
        in_deposit = in_deposit | prev_deposit
        print(f"  {len(prev_deposit)} file(s) in deposit {args.skip_also_in} — will also skip.")
    print()

    n_uploaded = 0
    n_skipped  = 0
    n_failed   = 0

    def _save_manifest():
        with open(manifest_path, "w") as f:
            json.dump({"deposits": [deposit_id], "sandbox": args.sandbox,
                       "files": manifest}, f, indent=2)

    # ── Upload loop ───────────────────────────────────────────────────────────
    for i, (local, remote) in enumerate(files, 1):
        size = human_size(local.stat().st_size)
        print(f"[{i}/{len(files)}] {remote}  ({size})")

        if remote in previously_done or remote in in_deposit:
            print("    Already uploaded — skipping.")
            n_skipped += 1
            if remote not in {e["remote"] for e in manifest}:
                manifest.append({"remote": remote,
                                  "local": str(local.relative_to(root)),
                                  "status": "skipped"})
            continue

        success = False
        for attempt in range(1, MAX_RETRIES + 1):
            try:
                result = upload_file(bucket_url, args.token, local, remote)
                checksum = result.get("checksum", "")
                print(f"    OK  checksum={checksum}")
                manifest.append({"remote": remote,
                                  "local": str(local.relative_to(root)),
                                  "size_bytes": local.stat().st_size,
                                  "deposit_id": deposit_id,
                                  "zenodo_checksum": checksum,
                                  "status": "uploaded"})
                n_uploaded += 1
                success = True
                _save_manifest()   # write after every successful upload
                break
            except Exception as e:
                print(f"    Attempt {attempt}/{MAX_RETRIES} failed: {e}")
                if attempt < MAX_RETRIES:
                    print(f"    Retrying in {RETRY_WAIT}s ...")
                    time.sleep(RETRY_WAIT)

        if not success:
            print(f"    FAILED after {MAX_RETRIES} attempts — continuing.")
            manifest.append({"remote": remote, "local": str(local.relative_to(root)),
                              "status": "failed"})
            n_failed += 1
            _save_manifest()

    _save_manifest()

    print()
    print("=" * 60)
    print(f"Uploaded : {n_uploaded}")
    print(f"Skipped  : {n_skipped}  (already done)")
    print(f"Failed   : {n_failed}")
    print(f"Manifest : {manifest_path.relative_to(root)}")
    print(f"Edit URL : {edit_url}")
    print()
    print("Next steps:")
    print("  1. Open the edit URL above and fill in title, description,")
    print("     authors, keywords, and related identifiers.")
    print("  2. Click Publish when ready.")
    if n_failed:
        print(f"  !! {n_failed} file(s) failed.")
        print(f"     Re-run with: --resume --deposit-id <new-deposit-id>")


if __name__ == "__main__":
    main()
