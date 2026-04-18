#!/usr/bin/env python3
"""
CRYPTAD — Trajectory QC plots
Reads .xvg files from 03_trajectory_qc/ and produces:
  - qc_rmsd.png       : Backbone RMSD for all 15 runs
  - qc_rg.png         : Protein Rg for all 15 runs
  - qc_energy.png     : Temperature, Pressure, Potential for all 15 runs
  - qc_summary.tsv    : Numerical summary (mean ± std for production phase)

Usage:
    python 09_scripts/04_trajectory_qc/05_plot_qc.py [--qc-dir PATH] [--project-root PATH]

Options:
    --qc-dir        Explicit path to 03_trajectory_qc/ output directory.
                    Default: <project-root>/03_trajectory_qc
    --project-root  Explicit path to CRYPTAD project root.
                    Default: inferred as three directories above this script.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
COLORS = ["#2196F3", "#FF5722", "#4CAF50"]   # blue, orange, green per replica

# Acceptance thresholds
RMSD_WARN_ANG = 3.0    # Å — flag if mean RMSD exceeds this
RG_DRIFT_ANG  = 1.5    # Å — flag if inter-replica Rg range exceeds this
TEMP_TARGET   = 310.0  # K
TEMP_TOL      = 5.0    # K
PRESS_TOL     = 200.0  # bar — mean |P| flagged if above this


# ---------------------------------------------------------------------------
# XVG parser
# ---------------------------------------------------------------------------
def parse_xvg(path: Path):
    """Return (time_array, values_array) skipping # and @ header lines."""
    times, values = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            times.append(float(cols[0]))
            values.append([float(c) for c in cols[1:]])
    t = np.array(times)
    v = np.array(values)
    if v.ndim == 2 and v.shape[1] == 1:
        v = v[:, 0]
    return t, v


def equilibrated_stats(y):
    """Mean and std of the second half of the array."""
    half = len(y) // 2
    return float(np.mean(y[half:])), float(np.std(y[half:]))


# ---------------------------------------------------------------------------
# Per-system panel helper
# ---------------------------------------------------------------------------
def fill_panel(ax, qc_root, sys_name, xvg_name, col_idx=0, scale=1.0):
    """Plot all 3 replicas. Returns list of (mean, std) or None per run."""
    stats = []
    found_any = False
    for i, run in enumerate(RUNS):
        fpath = qc_root / sys_name / run / xvg_name
        if not fpath.exists():
            stats.append(None)
            continue
        found_any = True
        t, v = parse_xvg(fpath)
        y = (v[:, col_idx] if v.ndim == 2 else v) * scale
        ax.plot(t, y, color=COLORS[i], alpha=0.85, linewidth=0.6, label=run)
        stats.append(equilibrated_stats(y))
    if not found_any:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", color="gray", fontsize=8)
    return stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    # Script lives at <root>/09_scripts/04_trajectory_qc/05_plot_qc.py
    # parents[0] = 04_trajectory_qc/  parents[1] = 09_scripts/  parents[2] = root
    _script = Path(__file__).resolve()
    _default_root = _script.parents[2]

    parser = argparse.ArgumentParser(description="CRYPTAD trajectory QC plots.")
    parser.add_argument("--qc-dir", default=None,
                        help="Path to 03_trajectory_qc/ directory. "
                             "Default: <project-root>/03_trajectory_qc")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    qc_root = Path(args.qc_dir).resolve() if args.qc_dir else project_root / "03_trajectory_qc"
    out_dir = qc_root

    if not qc_root.exists():
        sys.exit(f"[ERROR] Not found: {qc_root}\n"
                 "Run 09_scripts/04_trajectory_qc/03_qc_trajectories.sh on TRUBA first, "
                 "then rsync 03_trajectory_qc/ back.")

    print(f"Reading from: {qc_root}")
    summary_rows = []

    # -----------------------------------------------------------------------
    # Figure 1: RMSD
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 5, figsize=(18, 3.5))
    fig.suptitle("Backbone RMSD vs Frame 0", fontsize=11, fontweight="bold")

    for ax, sys_name in zip(axes, SYSTEMS):
        stats = fill_panel(ax, qc_root, sys_name, "rmsd_backbone.xvg", scale=10.0)  # nm→Å
        ax.set_title(sys_name.replace("_", "\n"), fontsize=8, fontweight="bold")
        ax.set_xlabel("Time (ns)", fontsize=8)
        ax.set_ylabel("RMSD (Å)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="upper left")
        for i, s in enumerate(stats):
            if s is None:
                continue
            flag = "HIGH" if s[0] > RMSD_WARN_ANG else "OK"
            if flag == "HIGH":
                ax.axhline(RMSD_WARN_ANG, color="red", linestyle="--",
                           linewidth=0.8, alpha=0.5)
                ax.text(0.98, 0.95, "⚠ HIGH", transform=ax.transAxes,
                        ha="right", va="top", fontsize=7, color="red")
            summary_rows.append(dict(system=sys_name, run=RUNS[i],
                                     quantity="RMSD_backbone_A",
                                     mean=round(s[0], 3), std=round(s[1], 3), flag=flag))

    fig.tight_layout()
    p = out_dir / "qc_rmsd.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"  Saved: {p}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Figure 2: Radius of gyration
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 5, figsize=(18, 3.5))
    fig.suptitle("Protein Radius of Gyration", fontsize=11, fontweight="bold")

    for ax, sys_name in zip(axes, SYSTEMS):
        stats = fill_panel(ax, qc_root, sys_name, "rg.xvg", col_idx=0, scale=10.0)  # nm→Å
        ax.set_title(sys_name.replace("_", "\n"), fontsize=8, fontweight="bold")
        ax.set_xlabel("Time (ns)", fontsize=8)
        ax.set_ylabel("Rg (Å)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=6, loc="upper left")
        valid = [s for s in stats if s]
        if valid:
            rng = max(s[0] for s in valid) - min(s[0] for s in valid)
            if rng > RG_DRIFT_ANG:
                ax.text(0.98, 0.95, f"⚠ {rng:.1f}Å spread",
                        transform=ax.transAxes, ha="right", va="top",
                        fontsize=7, color="red")
        for i, s in enumerate(stats):
            if s is None:
                continue
            summary_rows.append(dict(system=sys_name, run=RUNS[i], quantity="Rg_A",
                                     mean=round(s[0], 3), std=round(s[1], 3), flag="OK"))

    fig.tight_layout()
    p = out_dir / "qc_rg.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"  Saved: {p}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # Figure 3: Energy — 3 rows (T, P, Pot) × 5 systems
    # -----------------------------------------------------------------------
    fig, axes_e = plt.subplots(3, 5, figsize=(18, 9))
    fig.suptitle("Thermodynamic Properties (Production MD)",
                 fontsize=11, fontweight="bold")

    energy_cfg = [
        ("temperature.xvg", "Temperature (K)",   1.0, "Temperature", TEMP_TARGET, TEMP_TOL),
        ("pressure.xvg",    "Pressure (bar)",     1.0, "Pressure",    0.0,         PRESS_TOL),
        ("potential.xvg",   "Potential (kJ/mol)", 1.0, "Potential",   None,        None),
    ]

    for row, (xvg, ylabel, scale, qty, target, tol) in enumerate(energy_cfg):
        for col, sys_name in enumerate(SYSTEMS):
            ax = axes_e[row, col]
            stats = fill_panel(ax, qc_root, sys_name, xvg, scale=scale)
            if row == 0:
                ax.set_title(sys_name.replace("_", "\n"), fontsize=8, fontweight="bold")
            if col == 0:
                ax.set_ylabel(ylabel, fontsize=8)
            ax.set_xlabel("Time (ps)", fontsize=7)
            ax.tick_params(labelsize=6)
            if row < 2:
                ax.legend(fontsize=5, loc="upper right")
            for i, s in enumerate(stats):
                if s is None:
                    continue
                flag = "OK"
                if target is not None and tol is not None:
                    if abs(s[0] - target) > tol:
                        flag = "FLAG"
                        ax.axhline(target, color="green", linestyle="--",
                                   linewidth=0.8, alpha=0.6)
                        ax.text(0.98, 0.95, "⚠", transform=ax.transAxes,
                                ha="right", va="top", fontsize=9, color="red")
                summary_rows.append(dict(system=sys_name, run=RUNS[i], quantity=qty,
                                         mean=round(s[0], 3), std=round(s[1], 3), flag=flag))

    fig.tight_layout()
    p = out_dir / "qc_energy.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    print(f"  Saved: {p}")
    plt.close(fig)

    # -----------------------------------------------------------------------
    # TSV summary
    # -----------------------------------------------------------------------
    tsv = out_dir / "qc_summary.tsv"
    with open(tsv, "w") as fh:
        fh.write("system\trun\tquantity\tmean\tstd\tflag\n")
        for r in summary_rows:
            fh.write(f"{r['system']}\t{r['run']}\t{r['quantity']}\t"
                     f"{r['mean']}\t{r['std']}\t{r['flag']}\n")
    print(f"  Saved: {tsv}")

    # -----------------------------------------------------------------------
    # Console flags
    # -----------------------------------------------------------------------
    flags = [r for r in summary_rows if r["flag"] != "OK"]
    print("\n===== QC FLAGS =====")
    if not flags:
        print("  All checks PASSED.")
    else:
        for r in flags:
            print(f"  [{r['flag']}]  {r['system']:20s} {r['run']}  "
                  f"{r['quantity']:20s}  mean={r['mean']:.2f}  std={r['std']:.2f}")
    print("\nIf all OK → proceed to pocket detection: 09_scripts/05_pocket_detection/")


if __name__ == "__main__":
    main()
