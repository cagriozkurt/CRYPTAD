"""
CRYPTAD — Step 5.4: Selectivity Assessment
==========================================
Dock the 56 unique final hits against 4 off-target / canonical-site receptors:

  1. BIN1_SH3_canonical  — BIN1 SH3 domain peptide-binding groove (1MUZ)
  2. PICALM_PIP2         — PICALM ANTH domain PIP2-binding basic patch (7JXV)
  3. BACE1               — BACE-1 aspartyl protease active site (2QMG)
  4. AChE                — acetylcholinesterase catalytic gorge (4EY7, human)

Goal: confirm hits are selective for the cryptic pockets they were found in,
NOT general binders of these canonical / off-target sites.

Selection criterion (plan §Step 5.4):
  For each cryptic-site hit, selectivity is claimed if off-target Vina score is
  ≥ 2 kcal/mol less favourable (closer to 0) than the cryptic-site MM-GBSA ΔG.
  This is a qualitative screen; report ΔΔG and flag borderline cases.

Outputs (04_virtual_screening/selectivity/):
  receptors/{receptor_id}/receptor.pdbqt + vina_config.txt
  results/scores_all.csv
  results/selectivity_summary.csv
  06_figures/fig_selectivity.png
  results/selectivity_manifest.json

Usage:
  conda run -n vina_env python 09_scripts/08_virtual_screening/08_selectivity_docking.py
  conda run -n vina_env python 09_scripts/08_virtual_screening/08_selectivity_docking.py --project-root /path/to/CRYPTAD
"""

import argparse
import glob as _glob
import json
import logging
import subprocess
import shutil
import sys
import time
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch

warnings.filterwarnings("ignore", category=UserWarning)
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

RESET = "\033[0m"; GREEN = "\033[92m"; YELLOW = "\033[93m"
BOLD  = "\033[1m"; RED   = "\033[91m"; CYAN   = "\033[96m"
def log(m, c=RESET): print(f"{c}{m}{RESET}", flush=True)


# ── Helper: download PDB ──────────────────────────────────────────────────────
def download_pdb(pdb_id: str, dest: Path):
    if dest.exists():
        log(f"  [SKIP] {pdb_id}.pdb already downloaded", CYAN)
        return
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    log(f"  Downloading {pdb_id} from RCSB...")
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        urllib.request.urlretrieve(url, dest)
        log(f"  ✓ {dest.name} ({dest.stat().st_size // 1024} KB)", GREEN)
    except Exception as e:
        raise RuntimeError(f"Failed to download {pdb_id}: {e}")


# ── Helper: compute centroid from PDB ────────────────────────────────────────
def get_grid_centre(pdb_path: Path, chain: str, ligand_resname,
                    pocket_residues: list) -> tuple:
    """
    Grid box centroid calculation.
    Priority 1: if ligand_resname given, use all HETATM heavy-atom coords (no chain filter).
    Priority 2: Cα of pocket_residues in the specified chain.
                If pocket_residues is empty [], use ALL Cα in the chain (whole domain centre).
    Returns (cx, cy, cz).
    """
    lig_coords = []
    res_coords = []
    use_all_ca = not pocket_residues
    res_set    = set(pocket_residues) if pocket_residues else set()

    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "HETATM" and ligand_resname:
                rn = line[17:20].strip()
                if rn == ligand_resname:
                    try:
                        lig_coords.append([float(line[30:38]),
                                           float(line[38:46]),
                                           float(line[46:54])])
                    except ValueError:
                        pass
            elif rec == "ATOM":
                ch = line[21].strip()
                if chain and ch != chain:
                    continue
                if line[12:16].strip() != "CA":
                    continue
                try:
                    resid = int(line[22:26].strip())
                    if use_all_ca or resid in res_set:
                        res_coords.append([float(line[30:38]),
                                           float(line[38:46]),
                                           float(line[46:54])])
                except (ValueError, IndexError):
                    pass

    coords = lig_coords if lig_coords else res_coords
    if not coords:
        raise RuntimeError(
            f"No atoms found for centroid computation in {pdb_path.name} "
            f"(chain={chain}, ligand={ligand_resname}, "
            f"residues={'all' if use_all_ca else res_set})")
    c = np.mean(coords, axis=0)
    return float(c[0]), float(c[1]), float(c[2])


# ── Helper: clean PDB ─────────────────────────────────────────────────────────
def clean_pdb(src: Path, dst: Path):
    """Keep only ATOM/TER from the first MODEL (handles NMR multi-model ensembles)."""
    lines          = []
    in_first_model = False
    seen_model     = False

    with open(src) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "MODEL":
                if seen_model:
                    break
                seen_model = in_first_model = True
                continue
            if rec == "ENDMDL":
                in_first_model = False
                continue
            if rec not in ("ATOM", "TER"):
                continue
            if rec == "ATOM":
                alt = line[16]
                if alt not in (" ", "A"):
                    continue
                line = line[:16] + " " + line[17:]
                if line[17:20] in ("HSD", "HSE", "HSP"):
                    line = line[:17] + "HIS" + line[20:]
            lines.append(line)

    with open(dst, "w") as f:
        f.writelines(lines)


# ── Helper: prepare receptor (pdb2pqr + meeko) ───────────────────────────────
def prepare_receptor(rec: dict, rec_dir: Path) -> Path:
    rid    = rec["id"]
    outdir = rec_dir / rid
    outdir.mkdir(parents=True, exist_ok=True)
    pdbqt  = outdir / "receptor.pdbqt"
    cfg    = outdir / "vina_config.txt"

    if pdbqt.exists() and cfg.exists():
        log(f"  [SKIP] {rid}: already prepared", CYAN)
        return pdbqt

    log(f"\n{'─'*60}", BOLD)
    log(f"  Preparing: {rid}  ({rec['label']})", BOLD)

    src_pdb = rec["pdb_src"]

    if rec["download"]:
        download_pdb(rec["pdb_id"], src_pdb)

    if not src_pdb.exists():
        raise FileNotFoundError(f"PDB not found: {src_pdb}")

    cx, cy, cz = get_grid_centre(
        src_pdb, rec["chain"],
        rec.get("ligand_resname"), rec["pocket_residues"]
    )
    log(f"  Grid centre: ({cx:.2f}, {cy:.2f}, {cz:.2f})")

    clean = outdir / "receptor_clean.pdb"
    clean_pdb(src_pdb, clean)

    pqr   = outdir / "receptor.pqr"
    pdb_H = outdir / "receptor_H.pdb"
    pdb2pqr_cmd = shutil.which("pdb2pqr") or shutil.which("pdb2pqr30")
    if not pdb2pqr_cmd:
        raise RuntimeError("pdb2pqr not found in PATH")
    r = subprocess.run(
        [pdb2pqr_cmd, "--ff", "AMBER", "--pdb-output", str(pdb_H),
         "--drop-water", "--nodebump", str(clean), str(pqr)],
        capture_output=True, text=True
    )
    if not pdb_H.exists():
        raise RuntimeError(f"pdb2pqr failed for {rid}: {r.stderr[-300:]}")
    log(f"  ✓ pdb2pqr")

    mk = shutil.which("mk_prepare_receptor.py")
    if not mk:
        raise RuntimeError("mk_prepare_receptor.py not found")
    r = subprocess.run(
        [mk, "-i", str(pdb_H), "-o", str(pdbqt.with_suffix("")), "-p"],
        capture_output=True, text=True
    )
    if not pdbqt.exists():
        raise RuntimeError(f"meeko failed for {rid}: {r.stderr[-300:]}")
    log(f"  ✓ Meeko PDBQT ({pdbqt.stat().st_size//1024} KB)", GREEN)

    # Write vina config with filename-only receptor to avoid space-in-path issues
    cfg.write_text(
        f"receptor = {pdbqt.name}\n\n"
        f"center_x = {cx:.3f}\n"
        f"center_y = {cy:.3f}\n"
        f"center_z = {cz:.3f}\n\n"
        f"size_x = {rec['box_x']}\n"
        f"size_y = {rec['box_y']}\n"
        f"size_z = {rec['box_z']}\n\n"
        f"exhaustiveness = 8\n"
        f"num_modes       = 9\n"
        f"energy_range    = 3\n"
        f"cpu             = 8\n"
    )
    log(f"  ✓ Vina config written", GREEN)
    return pdbqt


# ── Helper: run vina for one ligand ──────────────────────────────────────────
def _parse_vina_score(pdbqt_out: Path) -> float:
    """Extract best (lowest) score from Vina output PDBQT."""
    try:
        with open(pdbqt_out) as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT"):
                    return float(line.split()[3])
    except Exception:
        pass
    return float("nan")


def dock_ligand(ligand_pdbqt: Path, rec_id: str, rec_dir: Path,
                res_dir: Path, vina_bin: str) -> float:
    """
    Run vina CLI in a subprocess (crash-isolated — avoids SIGABRT from Python API
    when a ligand atom falls outside the grid).
    Returns best docking score (kcal/mol). NaN on failure.
    """
    outdir    = res_dir / rec_id
    outdir.mkdir(exist_ok=True)
    out_pdbqt = outdir / f"{ligand_pdbqt.stem}_out.pdbqt"

    if out_pdbqt.exists():
        return _parse_vina_score(out_pdbqt)

    cfg = rec_dir / rec_id / "vina_config.txt"
    receptor_pdbqt = rec_dir / rec_id / "receptor.pdbqt"

    # Parse config — DO NOT use --config flag: vina's config parser fails on
    # paths containing spaces. Pass all parameters as CLI flags instead.
    cfg_data = {}
    for line in cfg.read_text().splitlines():
        if "=" in line and not line.strip().startswith("#"):
            k, _, v = line.partition("=")
            cfg_data[k.strip()] = v.strip()

    cmd = [
        vina_bin,
        "--receptor",       str(receptor_pdbqt),
        "--ligand",         str(ligand_pdbqt),
        "--out",            str(out_pdbqt),
        "--center_x",       cfg_data["center_x"],
        "--center_y",       cfg_data["center_y"],
        "--center_z",       cfg_data["center_z"],
        "--size_x",         cfg_data["size_x"],
        "--size_y",         cfg_data["size_y"],
        "--size_z",         cfg_data["size_z"],
        "--exhaustiveness", cfg_data.get("exhaustiveness", "8"),
        "--num_modes",      cfg_data.get("num_modes", "9"),
        "--cpu",            "8",
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode != 0 or not out_pdbqt.exists():
            return float("nan")
        return _parse_vina_score(out_pdbqt)
    except Exception as e:
        log(f"    [WARN] dock failed for {ligand_pdbqt.stem} / {rec_id}: {e}", YELLOW)
        return float("nan")


# ── Main ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD Selectivity Assessment (Step 5.4)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    sel_dir  = project_root / "04_virtual_screening" / "selectivity"
    rec_dir  = sel_dir / "receptors"
    res_dir  = sel_dir / "results"
    fig_dir  = project_root / "06_figures"
    ligpdbqt = project_root / "04_virtual_screening" / "ligands" / "pdbqt"
    mmgbsa   = project_root / "04_virtual_screening" / "mmgbsa_results" / "final_hits.csv"

    for d in [rec_dir, res_dir, fig_dir]:
        d.mkdir(parents=True, exist_ok=True)

    vina_bin = shutil.which("vina")
    if not vina_bin:
        log("[ERROR] 'vina' not found in PATH. "
            "Activate the vina_env conda environment before running.", RED)
        sys.exit(1)

    # ── Receptor definitions ──────────────────────────────────────────────────
    receptors = [
        {
            "id":       "BIN1_SH3_canonical",
            "label":    "BIN1 SH3 c-Myc binding groove",
            "pdb_id":   "1MUZ",
            "pdb_src":  project_root / "01_structures" / "BIN1" / "pdb" / "1MUZ.pdb",
            "download": False,
            "pocket_residues": [],
            "ligand_resname":  None,
            "chain": "A",
            "box_x": 25, "box_y": 25, "box_z": 25,
        },
        {
            "id":       "PICALM_ubiquitin_site",
            "label":    "PICALM ANTH ubiquitin-binding site",
            "pdb_id":   "7JXV",
            "pdb_src":  project_root / "01_structures" / "PICALM" / "pdb" / "7JXV.pdb",
            "download": False,
            "pocket_residues": [],
            "ligand_resname":  None,
            "chain": "B",
            "box_x": 25, "box_y": 25, "box_z": 25,
        },
        {
            "id":       "BACE1",
            "label":    "BACE-1 aspartyl protease active site",
            "pdb_id":   "2QMG",
            "pdb_src":  rec_dir / "BACE1" / "2QMG.pdb",
            "download": True,
            "pocket_residues": [93, 289, 71, 72, 169],
            "ligand_resname":  "SC6",
            "chain": "A",
            "box_x": 25, "box_y": 25, "box_z": 25,
        },
        {
            "id":       "AChE",
            "label":    "AChE catalytic gorge",
            "pdb_id":   "4EY7",
            "pdb_src":  rec_dir / "AChE" / "4EY7.pdb",
            "download": True,
            "pocket_residues": [203, 447, 334, 86, 286],
            "ligand_resname":  "E20",
            "chain": "A",
            "box_x": 25, "box_y": 25, "box_z": 25,
        },
    ]

    # ── Load final hits (unique compounds) ────────────────────────────────────
    df_hits    = pd.read_csv(mmgbsa)
    unique_ids = df_hits["chembl_id"].unique()
    log(f"\n{BOLD}CRYPTAD — Selectivity Assessment (Step 5.4){RESET}")
    log(f"  {len(df_hits)} hit rows  |  {len(unique_ids)} unique compounds  |"
        f"  {len(receptors)} off-target receptors")

    # ── Step 1: Receptor Preparation ─────────────────────────────────────────
    log(f"\n{BOLD}── Step 1: Receptor Preparation ──{RESET}")
    for rec in receptors:
        try:
            prepare_receptor(rec, rec_dir)
        except Exception as e:
            log(f"  [ERROR] {rec['id']}: {e}", RED)
            sys.exit(1)

    # ── Step 2: Docking ───────────────────────────────────────────────────────
    log(f"\n{BOLD}── Step 2: Docking {len(unique_ids)} compounds × {len(receptors)} receptors ──{RESET}")
    total = len(unique_ids) * len(receptors)
    done  = 0

    records = []
    for rec in receptors:
        for chembl_id in unique_ids:
            lig = ligpdbqt / f"{chembl_id}.pdbqt"
            if not lig.exists():
                log(f"  [SKIP] ligand not found: {chembl_id}", YELLOW)
                score = float("nan")
            else:
                score = dock_ligand(lig, rec["id"], rec_dir, res_dir, vina_bin)
            records.append({
                "chembl_id":      chembl_id,
                "receptor_id":    rec["id"],
                "receptor_label": rec["label"],
                "vina_score":     score,
            })
            done += 1
            if done % 20 == 0:
                log(f"  {done}/{total} docking runs complete...", CYAN)

    scores_df  = pd.DataFrame(records)
    scores_csv = res_dir / "scores_all.csv"
    scores_df.to_csv(scores_csv, index=False)
    log(f"\n  Scores saved: {scores_csv.relative_to(project_root)}", GREEN)

    # ── Step 3: Selectivity Summary ───────────────────────────────────────────
    log(f"\n{BOLD}── Step 3: Selectivity Summary ──{RESET}")

    cryptic_vina = {}
    top50_pattern = str(project_root / "04_virtual_screening" / "docking_results" / "top50_*.csv")
    for fpath in _glob.glob(top50_pattern):
        try:
            tdf = pd.read_csv(fpath)
            id_col    = "id" if "id" in tdf.columns else tdf.columns[0]
            score_col = "best_score" if "best_score" in tdf.columns else None
            if score_col:
                for _, row in tdf.iterrows():
                    cid = str(row[id_col])
                    sc  = float(row[score_col]) if pd.notna(row[score_col]) else float("nan")
                    if cid not in cryptic_vina or (not np.isnan(sc) and sc < cryptic_vina.get(cid, 0)):
                        cryptic_vina[cid] = sc
        except Exception:
            pass

    cryptic_mmgbsa = (
        df_hits.groupby("chembl_id")["dG_mean"].min().reset_index()
        .rename(columns={"dG_mean": "cryptic_mmgbsa_best"})
    )

    pivot = scores_df.pivot_table(
        index="chembl_id", columns="receptor_id", values="vina_score", aggfunc="min"
    )
    pivot.columns.name = None
    pivot = pivot.reset_index()
    for rid in [r["id"] for r in receptors]:
        if rid not in pivot.columns:
            pivot[rid] = float("nan")

    summary = cryptic_mmgbsa.merge(pivot, on="chembl_id", how="left")
    summary["cryptic_vina_best"] = summary["chembl_id"].map(cryptic_vina)

    rec_ids = [r["id"] for r in receptors]
    for rid in rec_ids:
        summary[f"ddG_{rid}"] = summary[rid] - summary["cryptic_vina_best"]

    summary["n_selective"]   = sum(
        (summary[f"ddG_{rid}"] >= 2.0).astype(int) for rid in rec_ids
    )
    summary["fully_selective"] = summary["n_selective"] == len(rec_ids)

    name_map = df_hits.set_index("chembl_id")["name"].to_dict()
    summary["name"] = summary["chembl_id"].map(name_map)

    sel_csv = res_dir / "selectivity_summary.csv"
    summary.to_csv(sel_csv, index=False)
    log(f"  Summary saved: {sel_csv.relative_to(project_root)}", GREEN)

    n_sel = int(summary["fully_selective"].sum())
    log(f"\n  Compounds fully selective (ΔΔG ≥ 2 kcal/mol for all 4 targets): "
        f"{n_sel} / {len(summary)}", GREEN if n_sel > 0 else YELLOW)
    log(f"  Compounds with ≥1 borderline off-target (ΔΔG < 2):               "
        f"{len(summary) - n_sel} / {len(summary)}")

    # ── Step 4: Publication Figure ────────────────────────────────────────────
    log(f"\n{BOLD}── Step 4: Selectivity Figure ──{RESET}")

    ddg_cols   = [f"ddG_{r['id']}" for r in receptors]
    rec_labels = [r["label"] for r in receptors]

    summary["mean_ddG"] = summary[ddg_cols].mean(axis=1)
    plot_df = summary.sort_values("mean_ddG", ascending=False).reset_index(drop=True)
    mat = plot_df[ddg_cols].values

    fig, axes = plt.subplots(1, 2, figsize=(15, max(7, len(plot_df) * 0.28 + 2)),
                             gridspec_kw={"width_ratios": [3, 1]})

    ax   = axes[0]
    norm = TwoSlopeNorm(vmin=-5, vcenter=2, vmax=12)
    im   = ax.imshow(mat, cmap="RdYlGn", norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(rec_labels)))
    ax.set_xticklabels(rec_labels, rotation=25, ha="right", fontsize=9)
    ax.set_yticks(range(len(plot_df)))
    compound_labels = [f"{row['name'][:22]} ({row['chembl_id']})"
                       for _, row in plot_df.iterrows()]
    ax.set_yticklabels(compound_labels, fontsize=7)

    for i in range(len(plot_df)):
        for j in range(len(rec_labels)):
            v = mat[i, j]
            if np.isnan(v):
                ax.text(j, i, "N/A", ha="center", va="center", fontsize=6, color="gray")
            else:
                color = "white" if abs(v) > 6 else "#222222"
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6.5, color=color)

    plt.colorbar(im, ax=ax, label="ΔΔG (kcal/mol)\n(off-target Vina − cryptic-site Vina)",
                 fraction=0.03, pad=0.02)
    ax.set_title(
        "Selectivity Assessment — ΔΔG per off-target receptor (Vina vs Vina)\n"
        "Green (≥2 kcal/mol) = selective for cryptic site; "
        "Red (<2) = borderline/not selective at Vina level",
        fontsize=9, fontweight="bold")
    ax.axvline(-0.5, color="#777777", lw=0.5)
    for i, row in plot_df.iterrows():
        if not row["fully_selective"]:
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(rec_labels), 1,
                                       fill=False, edgecolor="#E74C3C", lw=1.0))

    ax2 = axes[1]
    colors = ["#2ECC71" if v else "#E74C3C" for v in plot_df["fully_selective"]]
    bars = ax2.barh(range(len(plot_df)), plot_df["mean_ddG"], color=colors,
                    edgecolor="white", height=0.7)
    ax2.axvline(2.0, color="black", lw=1.0, ls="--", alpha=0.6, label="2 kcal/mol threshold")
    ax2.set_xlabel("Mean ΔΔG (kcal/mol)", fontsize=9)

    # Y-axis: dashed gridlines at every row to aid cross-referencing with heatmap
    ax2.set_yticks(range(len(plot_df)))
    ax2.set_yticklabels([])
    ax2.yaxis.grid(True, linestyle="--", linewidth=0.4, color="#bbbbbb", alpha=0.7)
    ax2.set_axisbelow(True)

    # Value annotation at end of each bar
    x_range = ax2.get_xlim()
    offset  = (x_range[1] - x_range[0]) * 0.02
    for bar, val in zip(bars, plot_df["mean_ddG"]):
        if not np.isnan(val):
            ax2.text(
                bar.get_width() + offset,
                bar.get_y() + bar.get_height() / 2,
                f"{val:.1f}",
                va="center", ha="left", fontsize=6.5, color="#333333",
            )
    # Expand xlim to prevent annotation clipping
    ax2.set_xlim(left=ax2.get_xlim()[0],
                 right=ax2.get_xlim()[1] + (ax2.get_xlim()[1] - ax2.get_xlim()[0]) * 0.12)

    # Match y-axis extent and direction to the heatmap (imshow inverts y)
    ax2.set_ylim(ax.get_ylim())
    ax2.invert_yaxis()

    ax2.set_title("Mean ΔΔG\nacross 4 targets", fontsize=9)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(fontsize=8, loc="upper right")

    legend_elements = [
        Patch(color="#2ECC71", label=f"Fully selective ({n_sel} compounds)"),
        Patch(color="#E74C3C", label=f"Borderline ({len(summary)-n_sel} compounds)"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=2,
               fontsize=9, frameon=False, bbox_to_anchor=(0.35, -0.02))

    plt.suptitle(
        f"CRYPTAD Selectivity Screen — {len(unique_ids)} final hits vs. "
        f"4 off-target/canonical receptors\n"
        f"ΔΔG = off-target Vina − cryptic-site Vina  "
        f"(Vina uncertainty ±1–2 kcal/mol; BACE-1/AChE selectivity requires MM-GBSA confirmation)",
        fontsize=10, fontweight="bold", y=1.01)
    plt.tight_layout()
    out_fig = fig_dir / "fig_selectivity.png"
    plt.savefig(out_fig, dpi=600, bbox_inches="tight")
    plt.close()
    log(f"  Figure saved: {out_fig.relative_to(project_root)}", GREEN)

    # Publication copies
    from PIL import Image as _Image
    pub_dir = fig_dir / "publication"
    pub_dir.mkdir(parents=True, exist_ok=True)
    pub_png = pub_dir / "Figure_10_selectivity.png"
    pub_tif = pub_dir / "Figure_10_selectivity.tif"
    import shutil as _shutil
    _shutil.copy(out_fig, pub_png)
    img = _Image.open(pub_png).convert("RGB")
    img.save(pub_tif, dpi=(600, 600), compression="tiff_lzw")
    log(f"  Publication PNG: {pub_png.relative_to(project_root)}", GREEN)
    log(f"  Publication TIF: {pub_tif.relative_to(project_root)}  "
        f"({img.width}×{img.height} px @ 600 dpi)", GREEN)

    log(f"\n{BOLD}Done. Selectivity assessment complete.{RESET}")
    log(f"  Scores:  {scores_csv}")
    log(f"  Summary: {sel_csv}")
    log(f"  Figure:  {out_fig}")

    # ── Manifest ──────────────────────────────────────────────────────────────
    manifest = {
        "generated_at":    time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script":          "09_scripts/08_virtual_screening/08_selectivity_docking.py",
        "project_root":    str(project_root),
        "vina_bin":        vina_bin,
        "n_compounds":     int(len(unique_ids)),
        "n_receptors":     len(receptors),
        "n_fully_selective": n_sel,
        "receptors":       [{"id": r["id"], "pdb_id": r["pdb_id"]} for r in receptors],
        "outputs": {
            "scores":   str(scores_csv),
            "summary":  str(sel_csv),
            "figure":   str(out_fig),
        },
    }
    manifest_path = res_dir / "selectivity_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"  Manifest: {manifest_path.relative_to(project_root)}")


if __name__ == "__main__":
    main()
