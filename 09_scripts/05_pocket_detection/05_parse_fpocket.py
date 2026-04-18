#!/usr/bin/env python3
"""
CRYPTAD — fpocket per-frame output parser and pocket persistence analysis

For each system/run, reads pocket_data.csv and:
  1. Clusters pocket centroids spatially (complete-linkage hierarchical, eps=8 Å)
     to identify recurring pocket sites across all frames.
     Complete linkage avoids DBSCAN chain-linking: every point in a cluster
     is guaranteed to be within eps of every other point in that cluster.
  2. Computes per-site persistence, mean score, druggability, volume.
  3. Performs cross-replica consistency check (≥2/3 runs).
  4. Classifies pockets: Tier 1 (persistent + druggable) / Tier 2 / Noise.
  5. Outputs summary CSVs and per-system figures.

Output in 03_pocket_analysis/pocket_sites/:
  {SYS}_{RUN}_sites.csv   — per-site metrics for each run
  {SYS}_summary.csv       — cross-replica summary per system
  {SYS}_analysis.png      — persistence timeline + druggability figure
  selected_pockets.csv    — Tier-1/2 pockets with ≥2 replicas for VS

Usage:
  python 09_scripts/05_pocket_detection/05_parse_fpocket.py
  python 09_scripts/05_pocket_detection/05_parse_fpocket.py --sys S4_CD2AP_SH3-2
  python 09_scripts/05_pocket_detection/05_parse_fpocket.py --eps 6 --tier1 0.15
  python 09_scripts/05_pocket_detection/05_parse_fpocket.py --runs metad
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import cdist

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

CLUSTER_EPS   = 8.0   # Å — complete-linkage diameter threshold per site
CROSS_DIST    = 10.0  # Å — cross-replica site matching (also complete linkage)
MAX_SAMPLE    = 5000  # max rows for direct hierarchical clustering;
                      # larger datasets are subsampled then all points assigned
TIER1_PERSIST = 0.20  # fraction of frames for Tier 1
TIER2_PERSIST = 0.05  # fraction of frames for Tier 2
TIER1_DRUGG   = 0.10  # mean druggability score for Tier 1
MIN_VOLUME    = 100.0 # Å³ — discard sites with mean volume below this

RUN_COLORS  = ["#2196F3", "#FF5722", "#4CAF50"]
TIER_COLORS = {"Tier1": "#E53935", "Tier2": "#FB8C00", "Noise": "#BDBDBD"}


# ---------------------------------------------------------------------------
# Per-run site identification
# ---------------------------------------------------------------------------
def load_run(csv_path):
    """Load pocket_data.csv. Returns DataFrame and total frame count.

    total_frames is counted BEFORE filtering zero-centroid rows so that
    frames where fpocket found no pockets (cx=cy=cz=0) are included in
    the denominator; otherwise persistence values would be inflated.
    """
    df = pd.read_csv(csv_path)
    total_frames = df["frame"].nunique()
    # Drop rows where centroid extraction failed (all zeros)
    df = df[df["cx"].abs() + df["cy"].abs() + df["cz"].abs() > 0]
    return df, total_frames


def _cluster_coords(coords, eps):
    """
    Complete-linkage hierarchical clustering.
    Returns 0-indexed integer label array.
    Within each cluster, max pairwise distance ≤ eps (no chain-linking).
    """
    if len(coords) == 1:
        return np.array([0])
    labels = fclusterdata(coords, t=eps, metric="euclidean",
                          criterion="distance", method="complete")
    return labels - 1  # fclusterdata is 1-indexed


def cluster_run(df, eps=CLUSTER_EPS, max_sample=MAX_SAMPLE):
    """
    Cluster pocket centroids across all frames using complete-linkage HC.
    For datasets larger than max_sample rows: cluster a random subsample,
    then assign every remaining point to the nearest cluster centroid.
    Returns df with 'site' column added.
    """
    df = df.copy()
    coords = df[["cx", "cy", "cz"]].values

    if len(coords) <= max_sample:
        df["site"] = _cluster_coords(coords, eps)
    else:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(coords), max_sample, replace=False)
        sample = coords[idx]
        sample_labels = _cluster_coords(sample, eps)

        # Compute centroid of each cluster from the sample
        n_clusters = int(sample_labels.max()) + 1
        centroids = np.vstack([
            sample[sample_labels == k].mean(axis=0)
            for k in range(n_clusters)
        ])

        # Assign all points to nearest centroid
        dists = cdist(coords, centroids)
        df["site"] = dists.argmin(axis=1)

    return df


def compute_site_metrics(df, total_frames):
    """
    Per-site aggregation. Returns DataFrame with one row per site.
    persistence = frames where ≥1 pocket belongs to site / total_frames.
    """
    rows = []
    for site_id, grp in df.groupby("site"):
        n_frames = grp["frame"].nunique()
        rows.append({
            "site":        int(site_id),
            "persistence": round(n_frames / total_frames, 4),
            "n_frames":    n_frames,
            "mean_score":  round(grp["score"].mean(), 4),
            "mean_drugg":  round(grp["drugg_score"].mean(), 4),
            "mean_volume": round(grp["volume"].mean(), 1),
            "cx":          round(grp["cx"].mean(), 2),
            "cy":          round(grp["cy"].mean(), 2),
            "cz":          round(grp["cz"].mean(), 2),
            "n_obs":       len(grp),
        })
    sites = pd.DataFrame(rows).sort_values("persistence", ascending=False)
    return sites.reset_index(drop=True)


def classify_tiers(sites, tier1_persist=TIER1_PERSIST,
                   tier2_persist=TIER2_PERSIST,
                   tier1_drugg=TIER1_DRUGG,
                   min_volume=MIN_VOLUME):
    """Add 'tier' column to site metrics DataFrame."""
    def _tier(row):
        if row["mean_volume"] < min_volume:
            return "Noise"
        if row["persistence"] >= tier1_persist and row["mean_drugg"] >= tier1_drugg:
            return "Tier1"
        if row["persistence"] >= tier2_persist:
            return "Tier2"
        return "Noise"
    sites = sites.copy()
    sites["tier"] = sites.apply(_tier, axis=1)
    return sites


# ---------------------------------------------------------------------------
# Cross-replica consistency
# ---------------------------------------------------------------------------
def match_cross_replica(site_list, max_dist=CROSS_DIST):
    """
    Match pocket sites across the 3 replicas using complete-linkage
    hierarchical clustering on site centroids (avoids chain-linking).
    Returns summary DataFrame with n_replicas column per group.
    """
    combined = []
    for run_idx, sites in enumerate(site_list):
        for _, row in sites[sites["tier"] != "Noise"].iterrows():
            combined.append({
                "run_idx":     run_idx,
                "site":        row["site"],
                "persistence": row["persistence"],
                "mean_score":  row["mean_score"],
                "mean_drugg":  row["mean_drugg"],
                "mean_volume": row["mean_volume"],
                "cx": row["cx"], "cy": row["cy"], "cz": row["cz"],
                "tier": row["tier"],
            })

    if not combined:
        return pd.DataFrame()

    comb = pd.DataFrame(combined)
    coords = comb[["cx", "cy", "cz"]].values

    if len(coords) == 1:
        comb["group"] = 0
    else:
        labels = fclusterdata(coords, t=max_dist, metric="euclidean",
                              criterion="distance", method="complete")
        comb["group"] = labels - 1

    rows = []
    for gid, grp in comb.groupby("group"):
        n_replicas = int(grp["run_idx"].nunique())
        rows.append({
            "group_id":    int(gid),
            "n_replicas":  n_replicas,
            "best_tier":   "Tier1" if "Tier1" in grp["tier"].values else "Tier2",
            "persistence": round(grp["persistence"].mean(), 4),
            "mean_score":  round(grp["mean_score"].mean(), 4),
            "mean_drugg":  round(grp["mean_drugg"].mean(), 4),
            "mean_volume": round(grp["mean_volume"].mean(), 1),
            "cx":          round(grp["cx"].mean(), 2),
            "cy":          round(grp["cy"].mean(), 2),
            "cz":          round(grp["cz"].mean(), 2),
            "runs":        ",".join(str(r+1) for r in sorted(grp["run_idx"].unique())),
        })

    return pd.DataFrame(rows).sort_values(
        ["n_replicas", "persistence"], ascending=[False, False]
    ).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_system(sys_name, df_list, sites_list, run_names, outpath):
    """
    Three-panel figure per system:
      Left   — n_pockets per frame over time (one line per run)
      Centre — persistence bar chart (top 20 non-Noise sites, coloured by tier)
      Right  — druggability vs persistence scatter (bubble area ∝ volume)
    """
    fig = plt.figure(figsize=(16, 5), constrained_layout=True)
    fig.suptitle(sys_name.replace("_", " "), fontsize=12, fontweight="bold")
    gs = gridspec.GridSpec(1, 3, figure=fig)

    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1])
    ax2 = fig.add_subplot(gs[2])

    # Panel 0 — n_pockets over time
    ax0.set_title("Pockets per frame", fontsize=9)
    ax0.set_xlabel("Time (ns)", fontsize=8)
    ax0.set_ylabel("n_pockets", fontsize=8)
    for i, (df, rname) in enumerate(zip(df_list, run_names)):
        npk = df.groupby("frame")[["time_ns", "n_pockets"]].first()
        ax0.plot(npk["time_ns"], npk["n_pockets"],
                 color=RUN_COLORS[i], alpha=0.55, linewidth=0.5, label=rname)
    ax0.legend(fontsize=7)

    # Panel 1 — persistence bar chart
    all_sites = pd.concat(
        [s.assign(run=r) for s, r in zip(sites_list, run_names)],
        ignore_index=True
    )
    non_noise = all_sites[all_sites["tier"] != "Noise"]
    top = non_noise.nlargest(20, "persistence")

    ax1.set_title("Site persistence (top 20)", fontsize=9)
    ax1.set_xlabel("Site (run:id)", fontsize=8)
    ax1.set_ylabel("Persistence", fontsize=8)
    if not top.empty:
        labels = [f"{r}:s{int(s)}" for r, s in zip(top["run"], top["site"])]
        colors = [TIER_COLORS.get(t, "#BDBDBD") for t in top["tier"]]
        ax1.bar(range(len(top)), top["persistence"], color=colors)
        ax1.set_xticks(range(len(top)))
        ax1.set_xticklabels(labels, rotation=90, fontsize=6)
        ax1.axhline(TIER1_PERSIST, color="#E53935", linestyle="--",
                    linewidth=0.8, label=f"Tier1 ({TIER1_PERSIST})")
        ax1.axhline(TIER2_PERSIST, color="#FB8C00", linestyle="--",
                    linewidth=0.8, label=f"Tier2 ({TIER2_PERSIST})")
        ax1.legend(fontsize=7)

    # Panel 2 — druggability vs persistence scatter
    ax2.set_title("Druggability vs persistence", fontsize=9)
    ax2.set_xlabel("Persistence", fontsize=8)
    ax2.set_ylabel("Mean druggability score", fontsize=8)
    if not non_noise.empty:
        sizes = np.clip(non_noise["mean_volume"] / 5, 5, 200)
        scatter_colors = [TIER_COLORS.get(t, "#BDBDBD") for t in non_noise["tier"]]
        ax2.scatter(non_noise["persistence"], non_noise["mean_drugg"],
                    s=sizes, c=scatter_colors, alpha=0.7, edgecolors="none")
        ax2.axvline(TIER1_PERSIST, color="#E53935", linestyle="--", linewidth=0.8)
        ax2.axhline(TIER1_DRUGG,   color="#E53935", linestyle=":",  linewidth=0.8,
                    label="Tier1 thresholds")
        ax2.legend(fontsize=7)

    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Per-system processor
# ---------------------------------------------------------------------------
def process_system(sys_name, pocket_root, sites_dir, args, runs):
    print(f"\n{'='*60}")
    print(f"  {sys_name}")
    print(f"{'='*60}")

    df_list    = []
    sites_list = []
    run_names  = []

    for run in runs:
        csv_path = pocket_root / sys_name / run / "pocket_data.csv"
        if not csv_path.exists():
            print(f"  [SKIP] {run}: pocket_data.csv not found")
            continue

        df, total_frames = load_run(csv_path)
        print(f"  {run}: {total_frames} frames, {len(df)} pocket observations")

        df_clust = cluster_run(df, eps=args.eps, max_sample=MAX_SAMPLE)
        n_sites  = df_clust["site"].nunique()
        sites    = compute_site_metrics(df_clust, total_frames)
        sites    = classify_tiers(sites,
                                  tier1_persist=args.tier1,
                                  tier2_persist=args.tier2,
                                  tier1_drugg=args.drugg,
                                  min_volume=args.min_volume)

        t1 = (sites["tier"] == "Tier1").sum()
        t2 = (sites["tier"] == "Tier2").sum()
        print(f"         {n_sites} clusters → Tier1: {t1}  Tier2: {t2}  "
              f"(top persist: {sites['persistence'].iloc[0]:.3f})")

        out_csv = sites_dir / f"{sys_name}_{run}_sites.csv"
        sites.to_csv(out_csv, index=False)

        df_list.append(df_clust)
        sites_list.append(sites)
        run_names.append(run)

    if not sites_list:
        print(f"  [SKIP] No data for {sys_name}")
        return None

    # Cross-replica consistency
    summary = match_cross_replica(sites_list, max_dist=args.cross_dist)
    if not summary.empty:
        summary.insert(0, "system", sys_name)
        summary.to_csv(sites_dir / f"{sys_name}_summary.csv", index=False)
        t1_cross = (summary["best_tier"] == "Tier1").sum()
        multi    = (summary["n_replicas"] >= 2).sum()
        print(f"  Cross-replica groups: {len(summary)} total  |  "
              f"Tier1: {t1_cross}  |  ≥2 replicas: {multi}")
    else:
        print(f"  No cross-replica sites found")
        summary = pd.DataFrame()

    plot_path = sites_dir / f"{sys_name}_analysis.png"
    plot_system(sys_name, df_list, sites_list, run_names, plot_path)
    print(f"  Figure → {plot_path.name}")

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    _script = Path(__file__).resolve()
    _default_root = _script.parents[2]

    parser = argparse.ArgumentParser(
        description="CRYPTAD fpocket persistence analysis")
    parser.add_argument("--sys",        default=None,
                        help="Single system (e.g. S4_CD2AP_SH3-2); default: all")
    parser.add_argument("--runs",       nargs="+", default=list(RUNS),
                        help="Run names to process (default: run1 run2 run3). "
                             "Use --runs metad to process metadynamics data.")
    parser.add_argument("--eps",        type=float, default=CLUSTER_EPS,
                        help=f"Complete-linkage cluster diameter in Å "
                             f"(default: {CLUSTER_EPS})")
    parser.add_argument("--cross-dist", type=float, default=CROSS_DIST,
                        help=f"Cross-replica matching diameter in Å "
                             f"(default: {CROSS_DIST})")
    parser.add_argument("--tier1",      type=float, default=TIER1_PERSIST,
                        help=f"Tier-1 persistence threshold (default: {TIER1_PERSIST})")
    parser.add_argument("--tier2",      type=float, default=TIER2_PERSIST,
                        help=f"Tier-2 persistence threshold (default: {TIER2_PERSIST})")
    parser.add_argument("--drugg",      type=float, default=TIER1_DRUGG,
                        help=f"Tier-1 druggability threshold "
                             f"(default: {TIER1_DRUGG})")
    parser.add_argument("--min-volume", type=float, default=MIN_VOLUME,
                        help=f"Minimum mean volume in Å³ (default: {MIN_VOLUME})")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    pocket_root  = project_root / "03_pocket_analysis" / "fpocket_results"
    sites_dir    = project_root / "03_pocket_analysis" / "pocket_sites"
    sites_dir.mkdir(parents=True, exist_ok=True)

    systems = [args.sys] if args.sys else SYSTEMS

    print("CRYPTAD fpocket persistence analysis")
    print(f"  Clustering: complete-linkage HC, eps={args.eps} Å  "
          f"(max cluster diameter — no chain-linking)")
    print(f"  Cross-replica match radius: {args.cross_dist} Å")
    print(f"  Tiers: Tier1≥{args.tier1} persist + drugg≥{args.drugg}  |  "
          f"Tier2≥{args.tier2} persist  |  min_vol≥{args.min_volume} Å³")
    print(f"  Runs: {args.runs}")

    all_summaries = []
    selected = pd.DataFrame()
    for sys_name in systems:
        summary = process_system(sys_name, pocket_root, sites_dir, args, args.runs)
        if summary is not None and not summary.empty:
            all_summaries.append(summary)

    # Write selected_pockets.csv — Tier1 and Tier2 with ≥2 replicas
    if all_summaries:
        combined = pd.concat(all_summaries, ignore_index=True)
        selected = combined[
            combined["best_tier"].isin(["Tier1", "Tier2"]) &
            (combined["n_replicas"] >= 2)
        ].sort_values(
            ["system", "best_tier", "n_replicas", "persistence"],
            ascending=[True, True, False, False]
        )
        selected.to_csv(sites_dir / "selected_pockets.csv", index=False)

        print(f"\n{'='*60}")
        print(f"Selected pockets (Tier1+Tier2, ≥2 replicas): {len(selected)}")
        print(f"\n{'System':<22} {'Tier':<7} {'Rep':<5} "
              f"{'Persist':<9} {'Drugg':<7} {'Vol':>7}  Coords")
        print("-" * 75)
        for _, row in selected.iterrows():
            print(f"{row['system']:<22} {row['best_tier']:<7} "
                  f"{row['n_replicas']:<5} {row['persistence']:<9.3f} "
                  f"{row['mean_drugg']:<7.3f} {row['mean_volume']:>7.1f}  "
                  f"({row['cx']:.1f}, {row['cy']:.1f}, {row['cz']:.1f})")

    print(f"\n===== DONE =====")
    print(f"Output: {sites_dir}")

    # Write manifest — record all parameters that affect output
    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "project_root": str(project_root),
        "systems": systems,
        "runs": args.runs,
        "cluster_eps_A": args.eps,
        "cross_dist_A": args.cross_dist,
        "tier1_persist": args.tier1,
        "tier2_persist": args.tier2,
        "tier1_drugg": args.drugg,
        "min_volume_A3": args.min_volume,
        "max_sample": MAX_SAMPLE,
        "n_selected": int(len(selected)),
    }
    manifest_path = sites_dir / "parse_fpocket_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"  Manifest → {manifest_path}")

    print(f"\nNext: python 09_scripts/06_persistence_gate/01_apply_persistence_gate.py")


if __name__ == "__main__":
    main()
