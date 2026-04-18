#!/usr/bin/env python3
"""
CRYPTAD — Per-chain RMSD for S1 BIN1 BAR dimer.
Reads rmsd_chain_A.xvg and rmsd_chain_B.xvg from
03_trajectory_qc/S1_BIN1_BAR/{run}/ and produces qc_perchain.png.

Usage:
    python 09_scripts/04_trajectory_qc/06_plot_qc_perchain.py [--qc-dir PATH] [--project-root PATH]

Options:
    --qc-dir        Path to 03_trajectory_qc/S1_BIN1_BAR/ directory.
                    Default: <project-root>/03_trajectory_qc/S1_BIN1_BAR
    --project-root  Path to CRYPTAD project root.
                    Default: inferred as three directories above this script.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS   = ["run1", "run2", "run3"]
COLORS = ["#2196F3", "#FF5722", "#4CAF50"]
THRESH_WARN = 5.0   # Å — per-chain RMSD "elevated" threshold
THRESH_FAIL = 8.0   # Å — per-chain RMSD "high / possible unfolding"


def parse_xvg(path: Path):
    t, v = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            cols = line.split()
            if len(cols) < 2:
                continue
            t.append(float(cols[0]))
            v.append(float(cols[1]))
    return np.array(t), np.array(v)


def equilibrated_stats(v: np.ndarray):
    h = len(v) // 2
    return float(np.mean(v[h:])), float(np.std(v[h:]))


def main():
    # Script lives at <root>/09_scripts/04_trajectory_qc/06_plot_qc_perchain.py
    # parents[0] = 04_trajectory_qc/  parents[1] = 09_scripts/  parents[2] = root
    _script = Path(__file__).resolve()
    _default_root = _script.parents[2]

    parser = argparse.ArgumentParser(description="CRYPTAD per-chain RMSD plot for S1 BIN1 BAR.")
    parser.add_argument("--qc-dir", default=None,
                        help="Path to 03_trajectory_qc/S1_BIN1_BAR/. "
                             "Default: <project-root>/03_trajectory_qc/S1_BIN1_BAR")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    qc_s1 = (Path(args.qc_dir).resolve() if args.qc_dir
              else project_root / "03_trajectory_qc" / "S1_BIN1_BAR")

    if not qc_s1.exists():
        sys.exit(f"[ERROR] Not found: {qc_s1}\n"
                 "Run 09_scripts/04_trajectory_qc/04_qc_perchain_s1.sh on TRUBA first, "
                 "then rsync 03_trajectory_qc/S1_BIN1_BAR/ back.")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
    fig.suptitle("S1 BIN1 BAR — Per-Chain Backbone RMSD", fontsize=11, fontweight="bold")

    print("\n===== S1 BIN1 BAR — Per-chain RMSD =====")
    print(f"{'Run':<8} {'Chain':<8} {'Mean (Å)':<12} {'Std (Å)':<10} {'Status'}")
    print("-" * 50)

    all_ok = True
    for ax, chain in zip(axes, ["A", "B"]):
        ax.set_title(f"Chain {chain}", fontsize=10, fontweight="bold")
        ax.set_xlabel("Time (ns)", fontsize=9)
        ax.set_ylabel("RMSD (Å)", fontsize=9)
        ax.axhline(THRESH_WARN, color="orange", linestyle="--", linewidth=0.8,
                   alpha=0.7, label=f"{THRESH_WARN} Å warn")
        ax.axhline(THRESH_FAIL, color="red",    linestyle="--", linewidth=0.8,
                   alpha=0.7, label=f"{THRESH_FAIL} Å fail")

        for i, run in enumerate(RUNS):
            fpath = qc_s1 / run / f"rmsd_chain_{chain}.xvg"
            if not fpath.exists():
                print(f"{run:<8} {'Chain '+chain:<8} {'NO DATA':>12}")
                continue
            t, v = parse_xvg(fpath)
            y = v * 10.0   # nm → Å
            ax.plot(t, y, color=COLORS[i], alpha=0.85, linewidth=0.7, label=run)
            m, s = equilibrated_stats(y)
            if   m >= THRESH_FAIL: status = "HIGH — possible unfolding"
            elif m >= THRESH_WARN: status = "ELEVATED — check structure"
            else:                  status = "OK"
            if status != "OK":
                all_ok = False
            print(f"{run:<8} {'Chain '+chain:<8} {m:>10.2f}   {s:>8.2f}   {status}")

        ax.legend(fontsize=7, loc="upper left")

    fig.tight_layout()
    out = qc_s1 / "qc_perchain.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {out}")

    print("\n===== VERDICT =====")
    if all_ok:
        print("  Both chains individually stable (<5 Å RMSD).")
        print("  High whole-dimer RMSD = inter-chain relative motion (BAR dynamics).")
        print("  S1 BIN1 BAR: GO for pocket analysis.")
        print("  The inter-chain motion may expose interface pockets — flag for inspection.")
    else:
        print("  One or both chains show elevated per-chain RMSD.")
        print("  Inspect structure visually (VMD/PyMOL) before proceeding.")
        print("  S1 BIN1 BAR: CONDITIONAL — visual check required.")


if __name__ == "__main__":
    main()
