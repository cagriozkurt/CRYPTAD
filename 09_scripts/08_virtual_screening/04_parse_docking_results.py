"""
parse_docking_results.py  —  CRYPTAD
Parse AutoDock Vina output and build cross-pocket hit table (Step 4.3 post-processing).

For each ligand × pocket combination:
  - Extracts best Vina score from PDBQT REMARK VINA RESULT line
  - Joins with ligand manifest (name, MW, LogP, TPSA, etc.)
  - Computes ligand efficiency (LE = score / heavy atom count)
  - Flags consensus hits (strong score in ≥2 pockets)
  - Deduplicates salt forms (keeps best-scoring representative)

Output:
  04_virtual_screening/docking_results/
    ├── docking_scores_all.csv      full ligand × pocket score table
    ├── top50_S1_site688.csv        top 50 for lead BIN1 pocket
    ├── top50_S1_site680.csv
    ├── top50_S1_site1813.csv
    ├── top50_S3_site473.csv
    ├── top_hits_consensus.csv      hits scoring ≤ HIT_THRESHOLD in ≥2 pockets
    └── score_distributions.png    violin/hist plot per pocket

Hit threshold: ≤ -7.0 kcal/mol (standard VS criterion)
Consensus:     ≤ -7.0 in ≥ 2 pockets

Usage:
  python3 09_scripts/08_virtual_screening/04_parse_docking_results.py
  python3 09_scripts/08_virtual_screening/04_parse_docking_results.py --threshold -7.5
  python3 09_scripts/08_virtual_screening/04_parse_docking_results.py --project-root /path/to/CRYPTAD
"""

import argparse
import json
import re
import time
from pathlib import Path
from typing import Optional

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem

# ---------------------------------------------------------------------------
_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; RED   = "\033[91m"; CYAN   = "\033[96m"
def log(m, c=RESET): print(f"{c}{m}{RESET}", flush=True)

HIT_THRESHOLD  = -7.0   # kcal/mol
CONSENSUS_N    = 2       # must hit in this many pockets
TOP_N          = 50      # compounds per pocket for MM-GBSA


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def parse_best_score(pdbqt_path: Path) -> Optional[float]:
    """Return best Vina score (most negative) from a PDBQT output file."""
    try:
        for line in pdbqt_path.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT"):
                return float(line.split()[3])
    except Exception:
        pass
    return None


def heavy_atom_count(smiles: str) -> Optional[int]:
    """Return number of heavy atoms from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol.GetNumHeavyAtoms() if mol else None
    except Exception:
        return None


def normalise_name(name: str) -> str:
    """Strip common salt suffixes for deduplication grouping."""
    suffixes = [
        r"\s+hydrochloride$", r"\s+hcl$", r"\s+hydrobromide$",
        r"\s+tosylate$", r"\s+mesylate$", r"\s+acetate$",
        r"\s+sodium$", r"\s+potassium$", r"\s+maleate$",
        r"\s+fumarate$", r"\s+tartrate$", r"\s+phosphate$",
        r"\s+sulfate$", r"\s+citrate$", r"\s+succinate$",
    ]
    n = name.lower().strip()
    for s in suffixes:
        n = re.sub(s, "", n, flags=re.IGNORECASE)
    return n.strip()


# ---------------------------------------------------------------------------
# Load manifest
# ---------------------------------------------------------------------------
def load_manifest(manifest_path: Path) -> pd.DataFrame:
    df = pd.read_csv(manifest_path)
    df = df[df["status"] == "PASS"].copy()
    df["hac"] = df["smiles"].apply(heavy_atom_count)
    df["name_norm"] = df["name"].apply(normalise_name)
    return df.set_index("id")


# ---------------------------------------------------------------------------
# Load scores for all pockets
# ---------------------------------------------------------------------------
def load_scores(pockets: list, dockdir: Path) -> pd.DataFrame:
    """Return DataFrame with one row per ligand, one column per pocket score."""
    scores = {}
    for pocket in pockets:
        pocket_dir = dockdir / pocket
        pocket_scores = {}
        for pdbqt in pocket_dir.glob("*_out.pdbqt"):
            lig_id = pdbqt.name.replace("_out.pdbqt", "")
            score  = parse_best_score(pdbqt)
            if score is not None:
                pocket_scores[lig_id] = score
        scores[pocket] = pocket_scores
        log(f"  {pocket}: {len(pocket_scores)} scores loaded", CYAN)

    df = pd.DataFrame(scores)
    df.index.name = "id"
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="CRYPTAD docking results parser")
    parser.add_argument("--threshold", type=float, default=HIT_THRESHOLD,
                        help=f"Hit score threshold (default: {HIT_THRESHOLD})")
    parser.add_argument("--consensus-n", type=int, default=CONSENSUS_N,
                        help=f"Min pockets for consensus hit (default: {CONSENSUS_N})")
    parser.add_argument("--top-n", type=int, default=TOP_N,
                        help=f"Top N per pocket for MM-GBSA list (default: {TOP_N})")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    dockdir  = project_root / "04_virtual_screening" / "docking_results"
    manifest_path = project_root / "04_virtual_screening" / "ligands" / "library_manifest.csv"

    log(f"\n{BOLD}CRYPTAD — Docking Results Parser (Step 4.3){RESET}")

    # Discover pockets
    pockets = sorted(p.name for p in dockdir.iterdir()
                     if p.is_dir() and p.name != "logs")
    log(f"Pockets: {pockets}\n")

    # Load data
    manifest = load_manifest(manifest_path)
    scores_df = load_scores(pockets, dockdir)

    # Merge with manifest
    df = scores_df.join(manifest[["name", "name_norm", "smiles", "mw",
                                   "logp", "tpsa", "hbd", "hba",
                                   "rotbonds", "hac"]], how="left")

    # Cross-pocket summary columns
    score_cols = pockets
    df["best_score"]   = df[score_cols].min(axis=1)
    df["mean_score"]   = df[score_cols].mean(axis=1).round(3)
    df["n_pockets_hit"] = (df[score_cols] <= args.threshold).sum(axis=1)

    # Ligand efficiency on best pocket score
    df["LE"] = (df["best_score"] / df["hac"]).round(4)

    # Per-pocket rank (1 = best)
    for p in pockets:
        df[f"rank_{p}"] = df[p].rank(method="min").astype("Int64")

    # ── Deduplication: keep best-scoring representative per normalised name ──
    df_dedup = (df.reset_index()
                  .sort_values("best_score")
                  .groupby("name_norm", as_index=False)
                  .first()
                  .set_index("id"))

    # ── Save full table ───────────────────────────────────────────────────────
    out_all = dockdir / "docking_scores_all.csv"
    col_order = (["name", "mw", "logp", "tpsa", "hbd", "hba", "rotbonds", "hac"]
                 + score_cols
                 + ["best_score", "mean_score", "LE", "n_pockets_hit"]
                 + [f"rank_{p}" for p in pockets])
    df[col_order].sort_values("best_score").to_csv(out_all)
    log(f"\n  Full table ({len(df)} compounds) → {out_all.relative_to(project_root)}", GREEN)

    # ── Top N per pocket ──────────────────────────────────────────────────────
    for pocket in pockets:
        top = (df_dedup[[pocket, "name", "mw", "logp", "tpsa", "LE",
                         "best_score", "mean_score", "n_pockets_hit"]]
               .dropna(subset=[pocket])
               .sort_values(pocket)
               .head(args.top_n))
        out = dockdir / f"top{args.top_n}_{pocket}.csv"
        top.to_csv(out)
        log(f"  Top {args.top_n} {pocket}: best={top[pocket].iloc[0]:.3f}  "
            f"worst={top[pocket].iloc[-1]:.3f}  → {out.name}", GREEN)

    # ── Consensus hits ────────────────────────────────────────────────────────
    consensus = (df_dedup[df_dedup["n_pockets_hit"] >= args.consensus_n]
                 .sort_values("mean_score"))
    out_cons = dockdir / "top_hits_consensus.csv"
    consensus[col_order].to_csv(out_cons)
    log(f"\n  Consensus hits (≤{args.threshold} in ≥{args.consensus_n} pockets): "
        f"{len(consensus)} compounds → {out_cons.name}", GREEN)

    # ── Print consensus table ─────────────────────────────────────────────────
    log(f"\n{'─'*90}", BOLD)
    log(f"  {'Name':<28} {'MW':>6} {'LogP':>5}  "
        + "  ".join(f"{p[-6:]:>7}" for p in pockets)
        + f"  {'Mean':>7}  {'LE':>6}  Pockets", BOLD)
    log(f"{'─'*90}", BOLD)
    for cid, row in consensus.head(30).iterrows():
        scores_str = "  ".join(
            f"{row[p]:>7.3f}" if pd.notna(row[p]) else f"{'NA':>7}"
            for p in pockets
        )
        hits = int(row["n_pockets_hit"])
        color = GREEN if hits >= 3 else (CYAN if hits == 2 else RESET)
        log(f"  {str(row['name']):<28} {row['mw']:>6.1f} {row['logp']:>5.2f}  "
            f"{scores_str}  {row['mean_score']:>7.3f}  {row['LE']:>6.4f}  {hits}/4",
            color)

    # ── Score distribution plot ───────────────────────────────────────────────
    fig, axes = plt.subplots(1, len(pockets), figsize=(4 * len(pockets), 5),
                             sharey=True, constrained_layout=True)
    fig.suptitle("CRYPTAD — Vina Score Distributions", fontsize=11, fontweight="bold")

    for ax, pocket in zip(axes, pockets):
        vals = df[pocket].dropna().values
        ax.hist(vals, bins=60, color="#4C72B0", edgecolor="none", alpha=0.85)
        ax.axvline(args.threshold, color="crimson", lw=1.2, ls="--",
                   label=f"threshold ({args.threshold})")
        n_hits = (vals <= args.threshold).sum()
        ax.set_title(f"{pocket}\n(n={len(vals)}, hits={n_hits})", fontsize=9)
        ax.set_xlabel("Vina score (kcal/mol)", fontsize=8)
        ax.legend(fontsize=7)

    axes[0].set_ylabel("Count", fontsize=8)
    out_plot = dockdir / "score_distributions.png"
    fig.savefig(out_plot, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log(f"\n  Plot → {out_plot.relative_to(project_root)}", GREEN)

    # ── Final summary ─────────────────────────────────────────────────────────
    log(f"\n{'═'*60}", BOLD)
    log(f"  Total ligands scored:  {len(df)}", BOLD)
    log(f"  Hits (≤{args.threshold} kcal/mol):", BOLD)
    for pocket in pockets:
        n = (df[pocket] <= args.threshold).sum()
        log(f"    {pocket}: {n}", CYAN)
    log(f"  Consensus (≥{args.consensus_n} pockets): {len(consensus)}", BOLD)
    log(f"\n  Next: Step 4.4 — MM-GBSA rescoring of top {args.top_n} hits/pocket", GREEN)
    log(f"  python3 09_scripts/08_virtual_screening/05_setup_mmgbsa.py", GREEN)

    # ── Manifest ──────────────────────────────────────────────────────────────
    manifest_out = {
        "generated_at":    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":          "09_scripts/08_virtual_screening/04_parse_docking_results.py",
        "project_root":    str(project_root),
        "hit_threshold":   args.threshold,
        "consensus_n":     args.consensus_n,
        "top_n":           args.top_n,
        "pockets":         pockets,
        "n_ligands_scored": int(len(df)),
        "n_consensus_hits": int(len(consensus)),
        "hits_per_pocket": {p: int((df[p] <= args.threshold).sum()) for p in pockets},
        "outputs": {
            "all_scores":  str(out_all),
            "consensus":   str(out_cons),
            "plot":        str(out_plot),
        },
    }
    manifest_path_out = dockdir / "parse_docking_manifest.json"
    with open(manifest_path_out, "w") as fh:
        json.dump(manifest_out, fh, indent=2)
    log(f"Manifest → {manifest_path_out}")


if __name__ == "__main__":
    main()
