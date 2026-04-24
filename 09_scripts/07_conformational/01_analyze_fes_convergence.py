"""
analyze_fes_convergence.py  —  CRYPTAD
FES block convergence analysis after plumed sum_hills.

Reads fes_67ns.dat, fes_133ns.dat, fes_200ns.dat per system and:
  1. Overlays 1D projections (CV1 and CV2 marginals) for visual inspection
  2. Computes ΔΔF between t=133 ns and t=200 ns for the global minimum basin
     Criterion: ΔΔF < 1 kJ/mol → converged
  3. Reports the lowest free-energy basin coordinates (CV1*, CV2*)
  4. Saves convergence overlay plots → 03_pocket_analysis/metadynamics/

Usage:
    source .venv/bin/activate
    python3 09_scripts/07_conformational/01_analyze_fes_convergence.py
    python3 09_scripts/07_conformational/01_analyze_fes_convergence.py --sys S3_PICALM_ANTH

Prerequisites:
    - plumed sum_hills run on TRUBA (09_scripts/03_metadynamics/02_run_sum_hills.sh)
    - fes_67ns.dat, fes_133ns.dat, fes_200ns.dat in each run directory
    - rsync'd to local (02_md_simulations/.../metadynamics/{SYS}/)
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
from scipy.ndimage import gaussian_filter

_script = Path(__file__).resolve()
_default_root = _script.parents[2]

SYSTEMS = {
    "S1_BIN1_BAR": {
        "path": "BIN1/metadynamics/S1_BAR",
        "cvs":  ["cv1 jaw (nm)", "cv2 hinge (rad)"],
        "ref":  [4.000, 2.693],
        "threshold": {"axis": "v", "value": 2.8, "label": "CV1 > 2.8 nm"},
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
        "threshold": {"axis": "h", "value": 3.0, "label": "CV2 > 3.0 nm"},
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


def load_fes(filepath):
    """
    Load a plumed sum_hills 2D FES file.
    Returns (cv1_vals, cv2_vals, fes_2d) where fes_2d is shaped (n1, n2).
    Skips header lines (starting with #) and blank lines.
    PLUMED writes row-major: outer loop = CV1, inner loop = CV2.
    Returns None if file missing.
    """
    if not Path(filepath).exists():
        return None
    data = np.loadtxt(filepath, comments="#")
    if data.ndim == 1 or data.shape[0] == 0:
        return None

    cv1_all = data[:, 0]
    cv2_all = data[:, 1]
    fes_all = data[:, 2]

    cv1_vals = np.unique(cv1_all)
    cv2_vals = np.unique(cv2_all)
    n1, n2 = len(cv1_vals), len(cv2_vals)

    # Reshape: PLUMED writes outer=CV1, inner=CV2
    try:
        fes_2d = fes_all.reshape(n1, n2)
    except ValueError:
        # Fallback: truncate to fit
        fes_2d = fes_all[: n1 * n2].reshape(n1, n2)

    return cv1_vals, cv2_vals, fes_2d


def marginal_1d(cv_vals, fes_2d, axis):
    """
    Boltzmann-weighted 1D projection along one CV axis.
    axis=0 → project over CV2 (gives F as function of CV1)
    axis=1 → project over CV1 (gives F as function of CV2)
    Returns free-energy profile shifted to min=0.
    """
    kT = 2.577  # kJ/mol at 310 K (kB=8.314e-3 kJ/mol/K * 310 K)
    # Convert FES to probability, marginalize, convert back
    prob_2d = np.exp(-fes_2d / kT)
    prob_2d = np.nan_to_num(prob_2d, nan=0.0, posinf=0.0, neginf=0.0)
    prob_1d = prob_2d.sum(axis=1 - axis)  # sum over the other axis
    prob_1d = np.clip(prob_1d, 1e-300, None)
    fes_1d  = -kT * np.log(prob_1d)
    fes_1d -= fes_1d.min()
    return fes_1d


def delta_delta_F(fes_133, fes_200):
    """
    Compute ΔΔF between t=133 and t=200 ns over the low-energy basin.
    Basin defined as regions where fes_200 < 5 kJ/mol (above global minimum).
    ΔΔF = mean |F_200 - F_133| over basin points.
    """
    if fes_133 is None or fes_200 is None:
        return None
    _, _, fes2d_133 = fes_133
    _, _, fes2d_200 = fes_200
    basin_mask = fes2d_200 < 5.0
    if basin_mask.sum() == 0:
        basin_mask = fes2d_200 < 10.0
    diff = np.abs(fes2d_200[basin_mask] - fes2d_133[basin_mask])
    return float(diff.mean())


def main():
    parser = argparse.ArgumentParser(description="CRYPTAD FES block convergence analysis")
    parser.add_argument("--sys", default=None,
                        help="Single system name (default: all 5 systems)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    rundir = project_root / "02_md_simulations"
    outdir = project_root / "03_pocket_analysis" / "metadynamics"
    outdir.mkdir(parents=True, exist_ok=True)

    systems = {args.sys: SYSTEMS[args.sys]} if args.sys else SYSTEMS
    if args.sys and args.sys not in SYSTEMS:
        import sys as _sys
        _sys.exit(f"[ERROR] Unknown system '{args.sys}'. "
                  f"Choose from: {list(SYSTEMS.keys())}")

    results = {}

    for sname, cfg in systems.items():
        run_path = rundir / cfg["path"]
        f67  = load_fes(run_path / "fes_67ns.dat")
        f133 = load_fes(run_path / "fes_133ns.dat")
        f200 = load_fes(run_path / "fes_200ns.dat")

        missing = [t for t, f in [("67ns", f67), ("133ns", f133), ("200ns", f200)] if f is None]
        if missing:
            log(f"[SKIP] {sname}: FES files missing — {missing}", RED)
            log(f"         Run 09_scripts/03_metadynamics/02_run_sum_hills.sh on TRUBA first, "
                f"then rsync.", YELLOW)
            continue

        log(f"\n{CYAN}▶ {sname}{RESET}")

        cv1_vals, cv2_vals, fes200 = f200
        cv1_vals_133, _, fes133 = f133

        # ── Global minimum basin coordinates (from 200 ns FES) ──────────────
        idx_min = np.unravel_index(np.argmin(fes200), fes200.shape)
        cv1_star = cv1_vals[idx_min[0]]
        cv2_star = cv2_vals[idx_min[1]]
        log(f"  Global FES minimum (200 ns): CV1*={cv1_star:.3f}, CV2*={cv2_star:.3f}")
        log(f"    vs. starting: CV1_ref={cfg['ref'][0]:.3f}, CV2_ref={cfg['ref'][1]:.3f}")

        # ── Cryptic pocket indicator ─────────────────────────────────────────
        d_cv1 = abs(cv1_star - cfg["ref"][0])
        d_cv2 = abs(cv2_star - cfg["ref"][1])
        if d_cv1 > 0.3 * cfg["ref"][0] or d_cv2 > 0.3 * cfg["ref"][1]:
            log(f"  ★ Basin displaced from start (ΔCV1={d_cv1:.3f}, ΔCV2={d_cv2:.3f}) "
                f"— potential alternate conformation sampled", GREEN)
        else:
            log(f"  Basin near starting point (ΔCV1={d_cv1:.3f}, ΔCV2={d_cv2:.3f})", YELLOW)

        # ── ΔΔF convergence criterion ────────────────────────────────────────
        ddF = delta_delta_F(f133, f200)
        converged = ddF is not None and ddF < 1.0
        log(f"  ΔΔF (133→200 ns, basin): {ddF:.2f} kJ/mol  "
            f"{'✓ CONVERGED (<1 kJ/mol)' if converged else '⚠ NOT CONVERGED (≥1 kJ/mol)'}",
            GREEN if converged else YELLOW)
        if not converged:
            log("    → Interpret as enhanced sampling only; no quantitative ΔG claims.", YELLOW)

        results[sname] = {
            "cv1_star": cv1_star, "cv2_star": cv2_star,
            "ddF": ddF, "converged": converged,
            "d_cv1": d_cv1, "d_cv2": d_cv2,
        }

        # ── Figure: 1D marginal overlays + 2D FES (200 ns) ──────────────────
        fig = plt.figure(figsize=(14, 13))
        fig.suptitle(f"{sname}  —  FES Block Convergence  (67 / 133 / 200 ns)",
                     fontsize=13, fontweight="bold")

        gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35,
                              height_ratios=[1, 1.8])
        ax_cv1 = fig.add_subplot(gs[0, 0])
        ax_cv2 = fig.add_subplot(gs[0, 1])
        ax_2d  = fig.add_subplot(gs[1, :])

        for fes, label, color in [(f67, "67 ns", "gray"), (f133, "133 ns", "steelblue"),
                                   (f200, "200 ns", "crimson")]:
            cv_v, _, fes2d = fes
            fes1d = marginal_1d(cv_v, fes2d, axis=0)
            ax_cv1.plot(cv_v, fes1d, lw=1.5, color=color, label=label)
        ax_cv1.axvline(cfg["ref"][0], color="black", ls="--", lw=0.8, alpha=0.6, label="t=0 ref")
        ax_cv1.set_xlabel(cfg["cvs"][0], fontsize=9)
        ax_cv1.set_ylabel("F (kJ/mol)", fontsize=9)
        ax_cv1.set_title("CV1 marginal FES", fontsize=10)
        ax_cv1.legend(fontsize=7)
        ax_cv1.set_ylim(bottom=0)

        for fes, label, color in [(f67, "67 ns", "gray"), (f133, "133 ns", "steelblue"),
                                   (f200, "200 ns", "crimson")]:
            cv_v, cv2_v, fes2d = fes
            fes1d = marginal_1d(cv2_v, fes2d, axis=1)
            ax_cv2.plot(cv2_v, fes1d, lw=1.5, color=color, label=label)
        ax_cv2.axvline(cfg["ref"][1], color="black", ls="--", lw=0.8, alpha=0.6, label="t=0 ref")
        ax_cv2.set_xlabel(cfg["cvs"][1], fontsize=9)
        ax_cv2.set_ylabel("F (kJ/mol)", fontsize=9)
        ax_cv2.set_title("CV2 marginal FES", fontsize=10)
        ax_cv2.legend(fontsize=7)
        ax_cv2.set_ylim(bottom=0)

        fes_smooth = gaussian_filter(np.clip(fes200, 0, 40), sigma=0.8)
        cf = ax_2d.contourf(cv1_vals, cv2_vals, fes_smooth.T,
                            levels=20, cmap="RdYlBu_r", vmin=0, vmax=40)
        plt.colorbar(cf, ax=ax_2d, label="F (kJ/mol)", pad=0.02)
        ax_2d.contour(cv1_vals, cv2_vals, fes_smooth.T,
                      levels=[2, 5, 10, 15, 20, 30], colors="k", linewidths=0.4, alpha=0.5)
        ax_2d.scatter([cfg["ref"][0]], [cfg["ref"][1]], c="white", s=80, zorder=5,
                      marker="x", linewidths=2, label="t=0 ref")
        ax_2d.scatter([cv1_star], [cv2_star], c="yellow", s=100, zorder=5,
                      marker="*", edgecolors="k", linewidths=0.5, label="global min")

        # Cyan dashed threshold line marking the cryptic pocket open region
        thr = cfg.get("threshold")
        if thr:
            if thr["axis"] == "v":
                ax_2d.axvline(thr["value"], color="cyan", lw=1.2, ls="--",
                              alpha=0.85, label=thr["label"])
            else:
                ax_2d.axhline(thr["value"], color="cyan", lw=1.2, ls="--",
                              alpha=0.85, label=thr["label"])
        ax_2d.set_xlabel(cfg["cvs"][0], fontsize=10)
        ax_2d.set_ylabel(cfg["cvs"][1], fontsize=10)
        ax_2d.set_title(f"2D FES at 200 ns  |  ΔΔF={ddF:.2f} kJ/mol  "
                        f"({'CONVERGED' if converged else 'NOT CONVERGED'})", fontsize=10)
        ax_2d.legend(fontsize=8, loc="upper right")

        outfile = outdir / f"fes_convergence_{sname}.png"
        plt.savefig(outfile, dpi=600, bbox_inches="tight")
        plt.close()
        log(f"  Plot → {outfile}", GREEN)

    # ── Summary table ─────────────────────────────────────────────────────────
    if results:
        log(f"\n{'═'*72}", BOLD)
        log("  CRYPTAD — FES Block Convergence Summary", BOLD)
        log(f"{'═'*72}", BOLD)
        log(f"  {'System':<22} {'CV1* (nm/rad)':>14}  {'CV2* (nm/rad)':>14}  "
            f"{'ΔΔF':>8}  {'Status':>12}", BOLD)
        log(f"  {'-'*22} {'-'*14}  {'-'*14}  {'-'*8}  {'-'*12}")
        for sname, r in results.items():
            status = "CONVERGED" if r["converged"] else "not conv."
            log(f"  {sname:<22} {r['cv1_star']:>14.3f}  {r['cv2_star']:>14.3f}  "
                f"{r['ddF']:>7.2f}  {status:>12}",
                GREEN if r["converged"] else YELLOW)

        log(f"\n  Interpretation:")
        log(f"  • CONVERGED (ΔΔF < 1 kJ/mol): quantitative ΔG claims allowed")
        log(f"  • NOT CONVERGED: report as 'enhanced sampling'; use FES for")
        log(f"    basin identification only, not quantitative binding affinities")
        log(f"\n  Plots → {outdir}")
        log(f"  Next: python3 09_scripts/07_conformational/02_analyze_fes_2d.py\n")

        # Write manifest
        manifest = {
            "generated_at":   time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "script":         "09_scripts/07_conformational/01_analyze_fes_convergence.py",
            "project_root":   str(project_root),
            "systems":        list(systems.keys()),
            "kT_kJ_mol":      2.577,
            "temperature_K":  310,
            "basin_threshold_kJ_mol": 5.0,
            "convergence_criterion_kJ_mol": 1.0,
            "results": {
                sname: {k: (float(v) if isinstance(v, (float, np.floating)) else v)
                        for k, v in r.items()}
                for sname, r in results.items()
            },
        }
        manifest_path = outdir / "fes_convergence_manifest.json"
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        log(f"  Manifest → {manifest_path}")
    else:
        log(f"\n[!] No FES data found. Run on TRUBA:", YELLOW)
        log(f"    bash 09_scripts/03_metadynamics/02_run_sum_hills.sh", YELLOW)
        log(f"    rsync output back to local", YELLOW)


if __name__ == "__main__":
    main()
