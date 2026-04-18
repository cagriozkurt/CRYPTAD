"""
analyze_persistence_gate.py  —  CRYPTAD
Unbiased persistence gate analysis (Step 3.4b) for cryptic pocket candidates.

For each persistence gate replica:
  1. Loads pocket_data.csv produced by run_fpocket_frames.py --run persist
  2. Counts frames where a pocket centroid is within SITE_RADIUS of the
     candidate's centroid → "open frames"
  3. Reports persistence per replica and pass/fail (≥20% in ≥1 replica)

Gate PASS criterion:
  pocket_persistence ≥ 0.20 in at least 1 of 3 replicas → confirmed cryptic
  pocket for virtual screening receptor preparation.

Usage:
    # S1 BIN1 BAR (default):
    python3 09_scripts/06_persistence_gate/04_analyze_persistence_gate.py

    # S3 PICALM ANTH:
    python3 09_scripts/06_persistence_gate/04_analyze_persistence_gate.py --sys S3_PICALM_ANTH
"""

import argparse
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
_script = Path(__file__).resolve()
_default_root = _script.parents[2]

def parse_args():
    parser = argparse.ArgumentParser(description="Unbiased persistence gate analysis")
    parser.add_argument("--sys", default="S1_BIN1_BAR",
                        help="System name from cryptic_candidates.csv (default: S1_BIN1_BAR)")
    parser.add_argument("--radius", type=float, default=8.0,
                        help="Pocket search radius in Å (default: 8.0)")
    parser.add_argument("--threshold", type=float, default=0.20,
                        help="Persistence threshold for PASS (default: 0.20)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    return parser.parse_args()

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; CYAN  = "\033[96m"; RED    = "\033[91m"
def log(m, c=RESET): print(f"{c}{m}{RESET}")

N_REPS = 3

# ---------------------------------------------------------------------------
def load_pocket_csv(gatedir, rep):
    """Load pocket_data.csv for persistence gate replica rep (1-indexed).

    Returns (df, total_frames) where total_frames is counted BEFORE filtering
    zero-centroid rows, so frames where fpocket found no pockets are included
    in the denominator and persistence values are not inflated.
    Returns (None, 0) if the file does not exist.
    """
    p = gatedir / f"rep{rep}" / "pocket_data.csv"
    if not p.exists():
        return None, 0
    df = pd.read_csv(p)
    total_frames = df["frame"].nunique()
    df = df[(df["cx"].abs() + df["cy"].abs() + df["cz"].abs()) > 0]
    return df, total_frames




def persistence_near_centroid(df, cx, cy, cz, total_frames, radius=8.0):
    """
    Fraction of frames that have a pocket centroid within `radius` Å of (cx,cy,cz).
    A frame counts if ANY pocket in that frame is within the radius.
    """
    dists = np.sqrt((df["cx"] - cx)**2 + (df["cy"] - cy)**2 + (df["cz"] - cz)**2)
    near  = df[dists <= radius]
    n_open = near["frame"].nunique()
    return n_open / total_frames if total_frames > 0 else 0.0


# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    sys_name     = args.sys
    site_radius  = args.radius
    pass_persist = args.threshold

    project_root = args.project_root.resolve()
    cands_csv = project_root / "03_pocket_analysis" / "metadynamics" / "cryptic_candidates.csv"
    gatedir   = project_root / "03_pocket_analysis" / "fpocket_results" / sys_name / "persist"
    outdir    = project_root / "03_pocket_analysis" / "persistence_gate"
    outdir.mkdir(parents=True, exist_ok=True)

    if not cands_csv.exists():
        log(f"[ERROR] {cands_csv} not found. "
            f"Run 09_scripts/07_conformational/03_analyze_metad_pockets.py first.", RED)
        return

    cands = pd.read_csv(cands_csv)
    cands = cands[(cands["system"] == sys_name) & (cands["tier"] != "Noise")]
    if cands.empty:
        log(f"No cryptic candidates found for {sys_name} in cryptic_candidates.csv.", YELLOW)
        return

    log(f"\n{BOLD}Unbiased Persistence Gate — {sys_name}{RESET}")
    log(f"Candidates: {len(cands)}  |  Replicas: {N_REPS}  |  Radius: {site_radius} Å  |  Threshold: {pass_persist:.0%}\n")

    # Load all replica data
    rep_dfs = {}
    rep_frames = {}
    for rep in range(1, N_REPS + 1):
        df, total = load_pocket_csv(gatedir, rep)
        if df is None:
            log(f"  [SKIP] rep{rep}: pocket_data.csv not found — run fpocket first", RED)
        else:
            rep_dfs[rep] = df
            rep_frames[rep] = total
            log(f"  rep{rep}: {total} frames, {len(df)} pocket observations")

    if not rep_dfs:
        log("\nNo replica data found. Run fpocket on persistence gate trajectories first.", RED)
        log(f"Expected: {gatedir}/rep{{1,2,3}}/pocket_data.csv")
        return

    # ── Per-candidate analysis ────────────────────────────────────────────────
    results = []
    log(f"\n{'Site':>6}  {'Persist R1':>10}  {'Persist R2':>10}  {'Persist R3':>10}  {'Max':>6}  {'Gate'}")
    log(f"{'─'*6}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*6}  {'─'*6}")

    for _, row in cands.iterrows():
        site_id = int(row["site"])
        cx, cy, cz = row["cx"], row["cy"], row["cz"]

        rep_persists = {}
        for rep, df in rep_dfs.items():
            p = persistence_near_centroid(df, cx, cy, cz, rep_frames[rep], radius=site_radius)
            rep_persists[rep] = p

        max_p = max(rep_persists.values()) if rep_persists else 0.0
        gate  = "PASS" if max_p >= pass_persist else "FAIL"

        r_cols = "  ".join(
            f"{rep_persists.get(r, float('nan')):>10.3f}" for r in range(1, N_REPS + 1)
        )
        color = GREEN if gate == "PASS" else YELLOW
        log(f"{site_id:>6}  {r_cols}  {max_p:>6.3f}  {gate}", color)

        results.append({
            "site":    site_id,
            **{f"persist_rep{r}": rep_persists.get(r, np.nan) for r in range(1, N_REPS + 1)},
            "max_persist": round(max_p, 4),
            "gate":   gate,
            "cx": cx, "cy": cy, "cz": cz,
        })

    results_df = pd.DataFrame(results)
    out_csv = outdir / f"persistence_gate_{sys_name}.csv"
    results_df.to_csv(out_csv, index=False)
    log(f"\n→ Results saved: {out_csv}")

    passed  = results_df[results_df["gate"] == "PASS"]
    failed  = results_df[results_df["gate"] == "FAIL"]
    log(f"\n{'═'*60}", BOLD)
    log(f"  PASS: {len(passed)}  |  FAIL: {len(failed)}", BOLD)
    if not passed.empty:
        log(f"\n  Sites passing unbiased persistence gate:", GREEN)
        for _, r in passed.iterrows():
            log(f"    ★ site {int(r['site']):>3}  max_persist={r['max_persist']:.3f}  "
                f"pos=({r['cx']:.1f},{r['cy']:.1f},{r['cz']:.1f})", GREEN)
        log(f"\n  → Proceed to receptor preparation (Step 4.1) for passed sites", GREEN)
    else:
        log(f"\n  No sites passed. Consider:", YELLOW)
        log(f"    1. Inspect trajectories in PyMOL — pocket may be present but fpocket misses it", YELLOW)
        log(f"    2. Try MDpocket/POVME3 as orthogonal cavity metric (§3.4 fallback)", YELLOW)
        log(f"    3. Report as 'transiently sampled' — still VS-eligible if max_persist ≥ 0.10", YELLOW)

    # ── Persistence timeline figure ───────────────────────────────────────────
    fig, axes = plt.subplots(1, len(cands), figsize=(4 * len(cands), 4),
                             constrained_layout=True)
    if len(cands) == 1:
        axes = [axes]
    fig.suptitle(f"{sys_name.replace('_', ' ')} — Unbiased Persistence Gate",
                 fontsize=11, fontweight="bold")

    for ax, (_, row) in zip(axes, cands.iterrows()):
        site_id = int(row["site"])
        cx, cy, cz = row["cx"], row["cy"], row["cz"]
        gate_row = results_df[results_df["site"] == site_id].iloc[0]

        for rep, df in rep_dfs.items():
            all_frames  = sorted(df["frame"].unique())
            dists       = np.sqrt((df["cx"] - cx)**2 + (df["cy"] - cy)**2 + (df["cz"] - cz)**2)
            open_frames = set(df[dists <= site_radius]["frame"].unique())
            times_ns, open_flag = [], []
            for fr in all_frames:
                t_ns = df[df["frame"] == fr]["time_ns"].iloc[0] if "time_ns" in df.columns else fr * 0.1
                times_ns.append(t_ns)
                open_flag.append(1 if fr in open_frames else 0)
            ax.plot(times_ns, open_flag, alpha=0.6, lw=0.8, label=f"rep{rep}")

        ax.axhline(pass_persist, color="crimson", ls="--", lw=0.8,
                   label=f"threshold ({pass_persist:.0%})")
        ax.set_title(f"site {site_id}  [max={gate_row['max_persist']:.2f}  {gate_row['gate']}]",
                     fontsize=9, fontweight="bold" if gate_row["gate"] == "PASS" else "normal")
        ax.set_xlabel("Time (ns)", fontsize=8)
        ax.set_ylabel("Pocket open (0/1)", fontsize=8)
        ax.set_ylim(-0.05, 1.15)
        ax.legend(fontsize=7)

    out_png = outdir / f"persistence_gate_{sys_name}.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"  Plot → {out_png}", GREEN)

    # Write manifest — record all parameters that affect gate outcome
    manifest = {
        "generated_at":   time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":         "09_scripts/06_persistence_gate/04_analyze_persistence_gate.py",
        "project_root":   str(project_root),
        "system":         sys_name,
        "site_radius_A":  site_radius,
        "pass_threshold": pass_persist,
        "n_candidates":   len(cands),
        "n_reps_loaded":  len(rep_dfs),
        "n_pass":         int(len(results_df[results_df["gate"] == "PASS"])),
        "n_fail":         int(len(results_df[results_df["gate"] == "FAIL"])),
        "output_csv":     str(out_csv),
        "output_png":     str(out_png),
    }
    manifest_path = outdir / f"persistence_gate_{sys_name}_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"  Manifest → {manifest_path}")


if __name__ == "__main__":
    main()
