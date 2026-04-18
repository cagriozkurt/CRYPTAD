"""
analyze_metad_pockets.py  —  CRYPTAD
Pocket persistence analysis on WTMetaD trajectories and cryptic pocket
identification by comparison against production MD pocket sites.

For each system:
  1. Loads fpocket_results/{SYS}/metad/pocket_data.csv
  2. Clusters centroids (complete-linkage HC, eps=8 Å) → metad pocket sites
  3. Computes persistence, druggability, volume per site
  4. Cross-references against production MD pocket_sites/{SYS}_summary.csv
       - Metad site with nearest production site > 10 Å → CRYPTIC CANDIDATE
       - Metad site with nearest production site ≤ 10 Å → CANONICAL (already seen)
  5. Outputs:
       03_pocket_analysis/metadynamics/metad_pockets_{SYS}.csv
       03_pocket_analysis/metadynamics/cryptic_candidates.csv
       03_pocket_analysis/metadynamics/metad_pockets_{SYS}.png

Usage:
    source .venv/bin/activate
    python3 09_scripts/07_conformational/03_analyze_metad_pockets.py
    python3 09_scripts/07_conformational/03_analyze_metad_pockets.py --sys S3_PICALM_ANTH
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fclusterdata
from scipy.spatial.distance import cdist

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

CLUSTER_EPS    = 8.0    # Å — complete-linkage diameter
CRYPTIC_DIST   = 10.0   # Å — metad site farther than this from all prod sites = cryptic
TIER1_PERSIST  = 0.20   # ≥20% frames
TIER2_PERSIST  = 0.05   # ≥5% frames
TIER1_DRUGG    = 0.10
MIN_VOLUME     = 100.0  # Å³

SYSTEMS = [
    "S1_BIN1_BAR",
    "S2_BIN1_SH3",
    "S3_PICALM_ANTH",
    "S4_CD2AP_SH3-2",
    "S5_CD2AP_SH3-1",
]

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; CYAN  = "\033[96m"; RED    = "\033[91m"
def log(m, c=RESET): print(f"{c}{m}{RESET}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def cluster_pockets(df, eps=CLUSTER_EPS):
    """Complete-linkage HC clustering on centroids. Returns df with 'site' col."""
    df = df.copy()
    coords = df[["cx", "cy", "cz"]].values
    if len(coords) == 0:
        df["site"] = pd.Series(dtype=int)
        return df
    if len(coords) == 1:
        df["site"] = 0
        return df
    labels = fclusterdata(coords, t=eps, metric="euclidean",
                          criterion="distance", method="complete")
    df["site"] = labels - 1
    return df


def site_metrics(df, total_frames):
    """Aggregate per-site statistics."""
    rows = []
    for sid, grp in df.groupby("site"):
        n_fr = grp["frame"].nunique()
        rows.append({
            "site":        int(sid),
            "persistence": round(n_fr / total_frames, 4),
            "n_frames":    n_fr,
            "mean_score":  round(grp["score"].mean(), 4),
            "mean_drugg":  round(grp["drugg_score"].mean(), 4),
            "mean_volume": round(grp["volume"].mean(), 1),
            "cx":          round(grp["cx"].mean(), 2),
            "cy":          round(grp["cy"].mean(), 2),
            "cz":          round(grp["cz"].mean(), 2),
        })
    sites = pd.DataFrame(rows).sort_values("persistence", ascending=False)
    return sites.reset_index(drop=True)


def classify(row):
    if row["mean_volume"] < MIN_VOLUME:
        return "Noise"
    if row["persistence"] >= TIER1_PERSIST and row["mean_drugg"] >= TIER1_DRUGG:
        return "Tier1"
    if row["persistence"] >= TIER2_PERSIST:
        return "Tier2"
    return "Noise"


def load_production_sites(sys_name, sites_dir):
    """Load production MD cross-replica summary for comparison."""
    path = sites_dir / f"{sys_name}_summary.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def find_cryptic(metad_sites, prod_sites, threshold=CRYPTIC_DIST):
    """
    For each metad site, find the distance to the nearest production site.
    Sites farther than `threshold` Å from all production sites = CRYPTIC CANDIDATE.
    Returns metad_sites DataFrame with 'nearest_prod_dist' and 'classification' cols.
    """
    metad_sites = metad_sites.copy()
    if prod_sites is None or prod_sites.empty:
        metad_sites["nearest_prod_dist"] = np.nan
        metad_sites["classification"] = "unknown (no prod data)"
        return metad_sites

    prod_coords  = prod_sites[["cx", "cy", "cz"]].values
    metad_coords = metad_sites[["cx", "cy", "cz"]].values

    if len(prod_coords) == 0:
        metad_sites["nearest_prod_dist"] = np.nan
        metad_sites["classification"] = "unknown (no prod sites)"
        return metad_sites

    dists = cdist(metad_coords, prod_coords).min(axis=1)
    metad_sites["nearest_prod_dist"] = dists.round(1)
    metad_sites["classification"] = [
        "CRYPTIC" if d > threshold else "canonical"
        for d in dists
    ]
    return metad_sites


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD metadynamics pocket persistence and cryptic candidate analysis")
    parser.add_argument("--sys", default=None,
                        help="Single system name (default: all 5 systems)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    pocket_root  = project_root / "03_pocket_analysis" / "fpocket_results"
    sites_dir    = project_root / "03_pocket_analysis" / "pocket_sites"
    outdir       = project_root / "03_pocket_analysis" / "metadynamics"
    outdir.mkdir(parents=True, exist_ok=True)

    if args.sys and args.sys not in SYSTEMS:
        sys.exit(f"[ERROR] Unknown system '{args.sys}'. "
                 f"Choose from: {SYSTEMS}")
    systems = [args.sys] if args.sys else SYSTEMS

    all_cryptic = []

    for sname in systems:
        csv_path = pocket_root / sname / "metad" / "pocket_data.csv"
        if not csv_path.exists():
            log(f"[SKIP] {sname}: metad/pocket_data.csv not found", RED)
            continue

        log(f"\n{CYAN}▶ {sname}{RESET}")

        df_raw = pd.read_csv(csv_path)
        # Count total frames BEFORE filtering zero-centroid rows so that frames
        # where fpocket found no pockets are included in the denominator;
        # otherwise persistence values would be inflated.
        total_frames = df_raw["frame"].nunique()
        df_raw = df_raw[(df_raw["cx"].abs() + df_raw["cy"].abs() + df_raw["cz"].abs()) > 0]
        log(f"  {total_frames} frames, {len(df_raw)} pocket observations")

        # Cluster
        df = cluster_pockets(df_raw)
        n_sites = df["site"].nunique()

        # Site metrics
        sites = site_metrics(df, total_frames)
        sites["tier"] = sites.apply(classify, axis=1)

        t1 = (sites["tier"] == "Tier1").sum()
        t2 = (sites["tier"] == "Tier2").sum()
        log(f"  {n_sites} clusters → Tier1: {t1}  Tier2: {t2}  Noise: {n_sites-t1-t2}")
        if len(sites) > 0:
            log(f"  Top site: persistence={sites['persistence'].iloc[0]:.3f}  "
                f"drugg={sites['mean_drugg'].iloc[0]:.3f}  "
                f"vol={sites['mean_volume'].iloc[0]:.0f} Å³")

        # Cross-reference against production MD
        prod = load_production_sites(sname, sites_dir)
        sites = find_cryptic(sites, prod)

        # Report
        non_noise = sites[sites["tier"] != "Noise"]
        cryptic   = non_noise[non_noise["classification"] == "CRYPTIC"]
        canonical = non_noise[non_noise["classification"] == "canonical"]
        log(f"  Cryptic candidates (>{CRYPTIC_DIST} Å from prod): {len(cryptic)}  "
            f"Canonical: {len(canonical)}",
            GREEN if len(cryptic) > 0 else YELLOW)

        for _, row in cryptic.iterrows():
            log(f"    ★ site {int(row['site']):>3}  persist={row['persistence']:.3f}  "
                f"drugg={row['mean_drugg']:.3f}  vol={row['mean_volume']:.0f} Å³  "
                f"Δprod={row['nearest_prod_dist']:.1f} Å  "
                f"pos=({row['cx']:.1f},{row['cy']:.1f},{row['cz']:.1f})  [{row['tier']}]",
                GREEN)

        # Save per-system CSV
        sites.insert(0, "system", sname)
        out_csv = outdir / f"metad_pockets_{sname}.csv"
        sites.to_csv(out_csv, index=False)
        log(f"  → {out_csv.name}")

        # Collect cryptic for global summary
        if not cryptic.empty:
            cdf = cryptic.copy()
            cdf.insert(0, "system", sname)
            all_cryptic.append(cdf)

        # ── Figure ──────────────────────────────────────────────────────────
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
        fig.suptitle(f"{sname}  —  Metadynamics Pocket Analysis",
                     fontsize=12, fontweight="bold")

        npk = df_raw.groupby("frame")[["time_ns", "n_pockets"]].first()
        axes[0].plot(npk["time_ns"], npk["n_pockets"], lw=0.4, color="steelblue", alpha=0.8)
        axes[0].set_xlabel("Time (ns)", fontsize=9)
        axes[0].set_ylabel("n_pockets", fontsize=9)
        axes[0].set_title("Pockets per frame", fontsize=10)

        top = non_noise.nlargest(25, "persistence") if not non_noise.empty else pd.DataFrame()
        if not top.empty:
            bar_colors = ["#E53935" if t == "Tier1" else
                          ("#FB8C00" if t == "Tier2" else "#BDBDBD") for t in top["tier"]]
            hatches = ["//" if c == "CRYPTIC" else "" for c in top["classification"]]
            bars = axes[1].bar(range(len(top)), top["persistence"],
                               color=bar_colors, edgecolor="k", linewidth=0.3)
            for bar, h in zip(bars, hatches):
                bar.set_hatch(h)
            axes[1].axhline(TIER1_PERSIST, color="#E53935", ls="--", lw=0.8,
                            label=f"Tier1 ({TIER1_PERSIST})")
            axes[1].axhline(TIER2_PERSIST, color="#FB8C00", ls="--", lw=0.8,
                            label=f"Tier2 ({TIER2_PERSIST})")
            axes[1].set_xticks(range(len(top)))
            axes[1].set_xticklabels([f"s{int(s)}" for s in top["site"]],
                                    rotation=90, fontsize=7)
            axes[1].legend(fontsize=7)
            axes[1].set_title("Site persistence (// = cryptic candidate)", fontsize=10)
            axes[1].set_ylabel("Persistence", fontsize=9)

        if not non_noise.empty:
            for cls, color, marker in [("canonical", "steelblue", "o"),
                                        ("CRYPTIC",   "crimson",   "*")]:
                sub = non_noise[non_noise["classification"] == cls]
                if sub.empty:
                    continue
                sz = np.clip(sub["mean_volume"] / 5, 10, 300)
                axes[2].scatter(sub["persistence"], sub["mean_drugg"],
                                s=sz, c=color, marker=marker,
                                alpha=0.75, edgecolors="k", linewidths=0.3,
                                label=cls, zorder=3 if cls == "CRYPTIC" else 2)
            axes[2].axvline(TIER1_PERSIST, color="#E53935", ls="--", lw=0.8)
            axes[2].axhline(TIER1_DRUGG,   color="#E53935", ls=":",  lw=0.8)
            axes[2].set_xlabel("Persistence", fontsize=9)
            axes[2].set_ylabel("Mean druggability", fontsize=9)
            axes[2].set_title("Druggability vs persistence", fontsize=10)
            axes[2].legend(fontsize=8)

        outpng = outdir / f"metad_pockets_{sname}.png"
        fig.savefig(outpng, dpi=150, bbox_inches="tight")
        plt.close(fig)
        log(f"  Plot → {outpng.name}", GREEN)

    # ── Global cryptic candidates summary ────────────────────────────────────
    log(f"\n{'═'*72}", BOLD)
    log("  CRYPTAD — Metadynamics Cryptic Pocket Candidates", BOLD)
    log(f"{'═'*72}", BOLD)

    if all_cryptic:
        combined = pd.concat(all_cryptic, ignore_index=True)
        combined = combined[combined["tier"] != "Noise"].sort_values(
            ["system", "persistence"], ascending=[True, False]
        )
        out_cands = outdir / "cryptic_candidates.csv"
        combined.to_csv(out_cands, index=False)

        log(f"\n  {'System':<22} {'Site':>5}  {'Tier':<7}  {'Persist':>8}  "
            f"{'Drugg':>6}  {'Vol (Å³)':>9}  {'Δprod (Å)':>10}  {'Position'}", BOLD)
        log(f"  {'-'*22} {'-'*5}  {'-'*7}  {'-'*8}  {'-'*6}  {'-'*9}  {'-'*10}  {'-'*20}")
        for _, row in combined.iterrows():
            log(f"  {row['system']:<22} {int(row['site']):>5}  {row['tier']:<7}  "
                f"{row['persistence']:>8.3f}  {row['mean_drugg']:>6.3f}  "
                f"{row['mean_volume']:>9.0f}  {row['nearest_prod_dist']:>10.1f}  "
                f"({row['cx']:.1f}, {row['cy']:.1f}, {row['cz']:.1f})",
                GREEN if row["tier"] == "Tier1" else YELLOW)

        log(f"\n  Total cryptic candidates: {len(combined)}  "
            f"(Tier1: {(combined['tier']=='Tier1').sum()}  "
            f"Tier2: {(combined['tier']=='Tier2').sum()})")
        log(f"  Saved → {out_cands}")

        # Write manifest
        manifest = {
            "generated_at":        time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "script":              "09_scripts/07_conformational/03_analyze_metad_pockets.py",
            "project_root":        str(project_root),
            "systems":             systems,
            "cluster_eps_A":       CLUSTER_EPS,
            "cryptic_dist_A":      CRYPTIC_DIST,
            "tier1_persist":       TIER1_PERSIST,
            "tier2_persist":       TIER2_PERSIST,
            "tier1_drugg":         TIER1_DRUGG,
            "min_volume_A3":       MIN_VOLUME,
            "n_cryptic_total":     int(len(combined)),
            "n_cryptic_tier1":     int((combined["tier"] == "Tier1").sum()),
            "n_cryptic_tier2":     int((combined["tier"] == "Tier2").sum()),
            "output_candidates":   str(out_cands),
        }
        manifest_path = outdir / "metad_pockets_manifest.json"
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        log(f"  Manifest → {manifest_path}")
    else:
        log("  No cryptic candidates found above thresholds.", YELLOW)
        log("  Consider lowering TIER2_PERSIST or CRYPTIC_DIST constants.")

    log(f"\n  Next: python3 09_scripts/06_persistence_gate/04_analyze_persistence_gate.py\n")


if __name__ == "__main__":
    main()
