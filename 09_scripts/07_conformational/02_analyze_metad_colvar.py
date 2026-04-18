"""
analyze_metad_colvar.py  —  CRYPTAD
Post-simulation metadynamics COLVAR analysis for all 5 systems.

Checks:
  1. CV statistics (min/max/mean/std) — gross unfolding detection
  2. Basin revisitation count (research plan §9 Step 3.3 criterion: >10 crossings)
  3. Bias accumulation (metad.bias growth — confirms Gaussians were deposited)
  4. Generates CV trajectory plots → 03_pocket_analysis/metadynamics/

Usage:
    source .venv/bin/activate
    python3 09_scripts/07_conformational/02_analyze_metad_colvar.py
    python3 09_scripts/07_conformational/02_analyze_metad_colvar.py --sys S3_PICALM_ANTH
"""

import argparse
import json
import os
import time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_script = Path(__file__).resolve()
_default_root = _script.parents[2]

SYSTEMS = {
    "S1_BIN1_BAR": {
        "path": "BIN1/metadynamics/S1_BAR",
        "cvs":  ["cv1 (jaw, nm)", "cv2 (hinge, rad)"],
        "ref":  [4.000, 2.693],
    },
    "S2_BIN1_SH3": {
        "path": "BIN1/metadynamics/S2_SH3",
        "cvs":  ["cv_site4 (nm)", "cv_site5 (nm)"],
        "ref":  [1.258, 1.293],
    },
    "S3_PICALM_ANTH": {
        "path": "PICALM/metadynamics/S3_ANTH",
        "cvs":  ["cv_site1 (nm)", "cv_site4 (nm)"],
        "ref":  [1.233, 1.545],
    },
    "S4_CD2AP_SH3-2": {
        "path": "CD2AP/metadynamics/S4_SH3-2",
        "cvs":  ["cv_site1 (nm)", "cv_site2 (nm)"],
        "ref":  [1.346, 1.267],
    },
    "S5_CD2AP_SH3-1": {
        "path": "CD2AP/metadynamics/S5_SH3-1",
        "cvs":  ["cv_site1 (nm)", "cv_site4 (nm)"],
        "ref":  [1.224, 1.496],
    },
}

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; CYAN  = "\033[96m"; RED    = "\033[91m"
def log(m, c=RESET): print(f"{c}{m}{RESET}")

def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD metadynamics COLVAR analysis")
    parser.add_argument("--sys", default=None,
                        help="Single system name (default: all 5 systems)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    rundir = project_root / "02_md_simulations"
    outdir = project_root / "03_pocket_analysis" / "metadynamics"
    outdir.mkdir(parents=True, exist_ok=True)

    if args.sys and args.sys not in SYSTEMS:
        import sys as _sys
        _sys.exit(f"[ERROR] Unknown system '{args.sys}'. "
                  f"Choose from: {list(SYSTEMS.keys())}")
    systems = {args.sys: SYSTEMS[args.sys]} if args.sys else SYSTEMS

    results = {}

    for sname, cfg in systems.items():
        colvar_path = rundir / cfg["path"] / "COLVAR"
        if not colvar_path.exists():
            log(f"[SKIP] {sname}: COLVAR not found", RED)
            continue

        log(f"\n{CYAN}▶ {sname}{RESET}")

        data    = np.loadtxt(colvar_path, comments="#")
        time_ns = data[:, 0] / 1000.0
        cv1     = data[:, 1]
        cv2     = data[:, 2]
        bias    = data[:, 3]

        # ── CV statistics ──────────────────────────────────────────────────
        for cv, label, ref in [(cv1, cfg["cvs"][0], cfg["ref"][0]),
                               (cv2, cfg["cvs"][1], cfg["ref"][1])]:
            log(f"  {label}:")
            log(f"    min={cv.min():.3f}  max={cv.max():.3f}  "
                f"mean={cv.mean():.3f}  std={cv.std():.3f}")
            if cv.max() > ref * 3.0:
                log(f"    ⚠  MAX exceeds 3× ref ({ref:.3f}) — inspect for unfolding", YELLOW)
            else:
                log(f"    ✓  No extreme excursion", GREEN)

        # ── Basin revisitation (CV1) — research plan criterion: >10 crossings
        thresh    = cfg["ref"][0] * 1.2
        crossings = int(np.sum(np.abs(np.diff((cv1 > thresh).astype(int)))))
        ok        = crossings > 10
        log(f"  CV1 basin revisitations: {crossings}  "
            f"{'✓ PASS (>10)' if ok else '⚠ FAIL (≤10)'}",
            GREEN if ok else YELLOW)

        # ── Bias growth ────────────────────────────────────────────────────
        log(f"  Final accumulated bias: {bias[-1]:.1f} kJ/mol")

        results[sname] = {
            "total_ns":  float(time_ns[-1]),
            "cv1_min": float(cv1.min()), "cv1_max": float(cv1.max()),
            "cv1_mean": float(cv1.mean()), "cv1_std": float(cv1.std()),
            "cv2_min": float(cv2.min()), "cv2_max": float(cv2.max()),
            "cv2_mean": float(cv2.mean()), "cv2_std": float(cv2.std()),
            "crossings": crossings, "bias_final": float(bias[-1]),
        }

        # ── Plot ───────────────────────────────────────────────────────────
        fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
        fig.suptitle(f"{sname}  —  WTMetaD COLVAR  ({time_ns[-1]:.1f} ns)",
                     fontsize=13, fontweight="bold")

        axes[0].plot(time_ns, cv1, lw=0.3, color="steelblue", alpha=0.8)
        axes[0].axhline(cfg["ref"][0], color="gray",   ls="--", lw=0.8, label="t=0")
        axes[0].axhline(thresh,        color="orange",  ls=":",  lw=0.8, label="×1.2 threshold")
        axes[0].set_ylabel(cfg["cvs"][0], fontsize=9)
        axes[0].legend(fontsize=7, loc="upper right")

        axes[1].plot(time_ns, cv2, lw=0.3, color="firebrick", alpha=0.8)
        axes[1].axhline(cfg["ref"][1], color="gray", ls="--", lw=0.8)
        axes[1].set_ylabel(cfg["cvs"][1], fontsize=9)

        axes[2].plot(time_ns, bias, lw=0.5, color="darkgreen", alpha=0.8)
        axes[2].set_ylabel("metad.bias (kJ/mol)", fontsize=9)
        axes[2].set_xlabel("Time (ns)", fontsize=10)

        plt.tight_layout()
        outfile = outdir / f"colvar_{sname}.png"
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        plt.close()
        log(f"  Plot → {outfile}", GREEN)

    # ── Summary table ──────────────────────────────────────────────────────
    if results:
        log(f"\n{'═'*72}", BOLD)
        log("  CRYPTAD — Metadynamics COLVAR Summary", BOLD)
        log(f"{'═'*72}", BOLD)
        log(f"  {'System':<22} {'ns':>6}  {'CV1 range (nm/rad)':>18}  "
            f"{'CV2 range':>14}  {'Revisit':>8}  {'Bias':>10}", BOLD)
        log(f"  {'-'*22} {'-'*6}  {'-'*18}  {'-'*14}  {'-'*8}  {'-'*10}")
        for sname, r in results.items():
            cv1r = f"{r['cv1_min']:.2f}–{r['cv1_max']:.2f}"
            cv2r = f"{r['cv2_min']:.2f}–{r['cv2_max']:.2f}"
            flag = "✓" if r["crossings"] > 10 else "⚠"
            log(f"  {sname:<22} {r['total_ns']:>6.1f}  {cv1r:>18}  {cv2r:>14}  "
                f"  {flag}{r['crossings']:>5}  {r['bias_final']:>8.0f} kJ/mol")

        log(f"\n  Plots saved to: {outdir}")
        log(f"  Next: python3 09_scripts/07_conformational/01_analyze_fes_convergence.py\n")

        # Write manifest
        manifest = {
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "script":       "09_scripts/07_conformational/02_analyze_metad_colvar.py",
            "project_root": str(project_root),
            "systems":      list(systems.keys()),
            "results":      results,
        }
        manifest_path = outdir / "colvar_manifest.json"
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        log(f"  Manifest → {manifest_path}")
    else:
        log(f"\n[!] No COLVAR data found in {rundir}.", YELLOW)
        log(f"    Ensure metadynamics runs are complete and synced from TRUBA.", YELLOW)


if __name__ == "__main__":
    main()
