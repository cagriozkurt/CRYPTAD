#!/usr/bin/env python3
"""
CRYPTAD — Per-frame fpocket trajectory analysis
Workaround for fpocket 4.1 mdpocket trajectory-reader crash (malloc corruption).

Reads protein_traj.pdb (multi-model PDB already on disk), calls fpocket on
each frame individually, aggregates pocket data into a CSV.

Reproducibility notes
---------------------
- fpocket binary path and version are recorded in fpocket_frames_manifest.json.
- stride is recorded in both the manifest and pocket_data.csv header comment,
  since stride directly determines which frames are analysed.

Output per run (in 03_pocket_analysis/fpocket_results/{SYS}/{RUN}/):
    pocket_data.csv   — frame, time_ns, n_pockets, rank, score, drugg_score,
                        volume, cx, cy, cz

Usage:
    python 09_scripts/05_pocket_detection/04_run_fpocket_frames.py
    python 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --sys S4_CD2AP_SH3-2 --run run1
    python 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --stride 5
    python 09_scripts/05_pocket_detection/04_run_fpocket_frames.py --fpocket ~/fpocket/bin/fpocket

Options:
    --project-root  Path to CRYPTAD project root.
                    Default: inferred as three directories above this script.
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
SYSTEMS = [
    "S1_BIN1_BAR",
    "S2_BIN1_SH3",
    "S3_PICALM_ANTH",
    "S4_CD2AP_SH3-2",
    "S5_CD2AP_SH3-1",
]
RUNS = ["run1", "run2", "run3"]

CSV_HEADER = ["frame", "time_ns", "n_pockets",
              "rank", "score", "drugg_score", "volume",
              "cx", "cy", "cz"]


# ---------------------------------------------------------------------------
# Multi-model PDB parser
# ---------------------------------------------------------------------------
def iter_frames(pdb_path, stride=1):
    """Yield (frame_idx, time_ns, atom_lines) for every `stride`-th MODEL block."""
    model_count = 0
    time_ns = 0.0
    current_atoms = []
    in_model = False

    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("TITLE"):
                # GROMACS: "TITLE     t=   0.000 step=      0"
                m = re.search(r't=\s*([\d.eE+\-]+)', line)
                if m:
                    time_ns = float(m.group(1)) / 1000.0   # ps → ns

            elif line.startswith("MODEL"):
                in_model = True
                current_atoms = []

            elif line.startswith("ENDMDL"):
                if in_model and model_count % stride == 0:
                    yield model_count, time_ns, current_atoms
                model_count += 1
                in_model = False
                current_atoms = []
                time_ns = 0.0

            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                current_atoms.append(line)


# ---------------------------------------------------------------------------
# fpocket runner
# ---------------------------------------------------------------------------
def pocket_centroid(out_dir, rank):
    """Compute centroid of alpha spheres from pocket{rank}_atm.pdb. Returns (cx,cy,cz)."""
    atm_pdb = os.path.join(out_dir, "pockets", f"pocket{rank}_atm.pdb")
    xs, ys, zs = [], [], []
    try:
        with open(atm_pdb) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
                    except (ValueError, IndexError):
                        pass
    except FileNotFoundError:
        pass
    if xs:
        n = len(xs)
        return sum(xs)/n, sum(ys)/n, sum(zs)/n
    return 0.0, 0.0, 0.0


def run_fpocket_on_lines(atom_lines, tmpdir, fpocket_bin):
    """Write atom_lines to temp PDB, run fpocket, parse info file. Returns list of dicts."""
    tmp_pdb = os.path.join(tmpdir, "frame.pdb")
    out_dir  = os.path.join(tmpdir, "frame_out")

    # Clean up from previous frame
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    if os.path.exists(tmp_pdb):
        os.remove(tmp_pdb)

    with open(tmp_pdb, "w") as fh:
        fh.writelines(atom_lines)
        fh.write("END\n")

    result = subprocess.run(
        [fpocket_bin, "-f", tmp_pdb],
        capture_output=True, text=True, cwd=tmpdir
    )

    if result.returncode != 0:
        return []

    # fpocket names the output based on the input filename: frame_out/frame_info.txt
    info_candidates = [
        os.path.join(out_dir, "frame_info.txt"),
        os.path.join(out_dir, "fpocket.log"),
        os.path.join(out_dir, "frame_out.log"),
    ]
    for cand in info_candidates:
        if os.path.exists(cand):
            pockets = parse_info_file(cand)
            # Compute centroids from alpha-sphere PDB files (not in info file in fpocket 4.0)
            for p in pockets:
                cx, cy, cz = pocket_centroid(out_dir, p["rank"])
                p["cx"], p["cy"], p["cz"] = cx, cy, cz
            return pockets

    return []


def parse_info_file(info_path):
    """Parse fpocket _info.txt. Returns list of pocket dicts sorted by rank."""
    pockets = []
    cur = {}

    with open(info_path) as fh:
        for line in fh:
            line = line.strip()

            if re.match(r'Pocket\s+\d+\s*:', line):
                if cur:
                    pockets.append(cur)
                pocket_id = int(re.search(r'\d+', line).group())
                cur = {"rank": pocket_id}

            # "Score :" but NOT "Druggability Score :"
            elif re.match(r'Score\s*:', line):
                cur["score"] = _last_float(line)

            elif "Druggability Score" in line:
                cur["drugg_score"] = _last_float(line)

            elif re.match(r'Volume\s*:', line):
                cur["volume"] = _last_float(line)

    if cur:
        pockets.append(cur)
    return pockets


def _last_float(line):
    """Return the last float-like token in a line."""
    tokens = re.findall(r'[-\d.eE+]+', line)
    try:
        return float(tokens[-1])
    except (IndexError, ValueError):
        return 0.0


# ---------------------------------------------------------------------------
# Per-run processor
# ---------------------------------------------------------------------------
def process_run(sys_name, run_name, outdir, fpocket_bin, stride, force=False):
    traj_pdb = outdir / "protein_traj.pdb"
    csv_out  = outdir / "pocket_data.csv"

    if not traj_pdb.exists():
        print(f"  [SKIP] {traj_pdb} not found — run 03_run_mdpocket.sh first (XTC→PDB step)")
        return False

    if csv_out.exists() and not force:
        print(f"  [SKIP] pocket_data.csv already exists for {sys_name}/{run_name}")
        return True
    if csv_out.exists() and force:
        csv_out.unlink()

    with open(traj_pdb) as _f:
        total_frames = sum(1 for line in _f if line.startswith("MODEL"))
    sampled = (total_frames + stride - 1) // stride
    print(f"\n  {sys_name}/{run_name}: {total_frames} frames, "
          f"stride={stride} → {sampled} frames to process")

    t0 = time.time()
    rows_written = 0

    with tempfile.TemporaryDirectory(prefix="fpocket_") as tmpdir:
        with open(csv_out, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(CSV_HEADER)

            for frame_idx, time_ns, atom_lines in iter_frames(traj_pdb, stride=stride):
                if not atom_lines:
                    continue

                pockets = run_fpocket_on_lines(atom_lines, tmpdir, fpocket_bin)
                n_pockets = len(pockets)

                if n_pockets == 0:
                    writer.writerow([frame_idx, f"{time_ns:.3f}", 0,
                                     0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
                    rows_written += 1
                else:
                    for p in pockets:
                        writer.writerow([
                            frame_idx, f"{time_ns:.3f}", n_pockets,
                            p.get("rank", 0),
                            round(p.get("score", 0.0), 4),
                            round(p.get("drugg_score", 0.0), 4),
                            round(p.get("volume", 0.0), 2),
                            round(p.get("cx", 0.0), 3),
                            round(p.get("cy", 0.0), 3),
                            round(p.get("cz", 0.0), 3),
                        ])
                        rows_written += 1

                # Progress every 200 sampled frames
                sample_num = frame_idx // stride
                if sample_num % 200 == 0 and sample_num > 0:
                    elapsed = time.time() - t0
                    rate = sample_num / elapsed
                    remaining = (sampled - sample_num) / rate if rate > 0 else 0
                    print(f"    frame {frame_idx:5d} ({time_ns:7.1f} ns) "
                          f"— {n_pockets} pockets | "
                          f"{rate:.1f} frames/s | ~{remaining/60:.0f} min left")

    elapsed = time.time() - t0
    print(f"  Done in {elapsed/60:.1f} min → {csv_out}  ({rows_written} rows)")
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    _script = Path(__file__).resolve()
    _default_root = _script.parents[2]

    parser = argparse.ArgumentParser(
        description="Per-frame fpocket analysis (mdpocket crash workaround)")
    parser.add_argument("--sys",    default=None, help="Single system name (e.g. S4_CD2AP_SH3-2)")
    parser.add_argument("--run",    default=None, help="Single run (e.g. run1)")
    parser.add_argument("--stride", type=int, default=1,
                        help="Process every Nth frame (default: 1 = all frames). "
                             "Use 5-10 for large systems to save time.")
    parser.add_argument("--fpocket", default="fpocket",
                        help="Path to fpocket binary (default: fpocket in PATH)")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing pocket_data.csv (default: skip if exists)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help="Path to CRYPTAD project root "
                             f"(default: {_default_root})")
    args = parser.parse_args()

    # Verify fpocket
    fpocket_bin = shutil.which(args.fpocket) or args.fpocket
    if not os.path.isfile(fpocket_bin):
        sys.exit(f"[ERROR] fpocket not found: {fpocket_bin}\n"
                 "       Build from source: cd ~/fpocket && make ARCH=MACOSXARM64\n"
                 "       Or: export PATH=$HOME/fpocket/bin:$PATH")

    # Record fpocket version
    ver_result = subprocess.run([fpocket_bin, "--version"],
                                capture_output=True, text=True)
    fpocket_version = (ver_result.stdout or ver_result.stderr).strip().splitlines()[0] \
                      if (ver_result.stdout or ver_result.stderr) else "unknown"
    print(f"Using fpocket: {fpocket_bin}  ({fpocket_version})")

    # Locate project root
    project_root = args.project_root.resolve()
    pocket_root  = project_root / "03_pocket_analysis" / "fpocket_results"

    if not pocket_root.exists():
        sys.exit(f"[ERROR] {pocket_root} not found")

    # Build job list
    if args.sys and args.run:
        jobs = [(args.sys, args.run)]
    elif args.sys:
        jobs = [(args.sys, r) for r in RUNS]
    else:
        jobs = [(s, r) for s in SYSTEMS for r in RUNS]

    print(f"Jobs: {len(jobs)}  |  stride: {args.stride}")

    ok = fail = 0
    for sys_name, run_name in jobs:
        outdir = pocket_root / sys_name / run_name
        outdir.mkdir(parents=True, exist_ok=True)
        result = process_run(sys_name, run_name, outdir, fpocket_bin, args.stride,
                             force=args.force)
        if result is True:
            ok += 1
        elif result is False:
            fail += 1

    print(f"\n===== DONE =====")
    print(f"  OK: {ok}  |  Failed/skipped: {fail}")

    # Write manifest
    import json
    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "fpocket_binary": fpocket_bin,
        "fpocket_version": fpocket_version,
        "stride": args.stride,
        "project_root": str(project_root),
        "jobs": [f"{s}/{r}" for s, r in jobs],
        "ok": ok,
        "failed": fail,
    }
    manifest_path = pocket_root / "fpocket_frames_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"  Manifest → {manifest_path}")

    print(f"\nNext: python 09_scripts/05_pocket_detection/05_parse_fpocket.py")


if __name__ == "__main__":
    main()
