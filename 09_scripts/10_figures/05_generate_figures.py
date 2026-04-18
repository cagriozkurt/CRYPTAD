"""
CRYPTAD — Publication figures for Phase 3 virtual screening results.

Figure 1: Screening funnel (3,127 → 7 stable hits)
Figure 2: Hit properties scatter (LE vs dG, coloured by pocket)
Figure 3: Contact occupancy heatmap — site473 (PICALM)
Figure 4: Contact occupancy heatmap — site688 + site680 (BIN1)
Figure 5: Final hits summary table

Outputs: 06_figures/fig{1-5}_*.png  (300 dpi)

Usage:
  python3 09_scripts/10_figures/05_generate_figures.py
  python3 09_scripts/10_figures/05_generate_figures.py --project-root /path/to/CRYPTAD
"""

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings("ignore")

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

STABLE_HITS = {
    ("S3_site473", "CHEMBL175"),
    ("S3_site473", "CHEMBL341812"),
    ("S3_site473", "CHEMBL797"),
    ("S1_site688", "CHEMBL1162"),
    ("S1_site688", "CHEMBL3707377"),
    ("S1_site688", "CHEMBL4105630"),
    ("S1_site688", "CHEMBL1200624"),
    ("S1_site680", "CHEMBL4105630"),
}

SHORT_NAMES = {
    "CHEMBL175":    "Dexibuprofen",
    "CHEMBL341812": "Ibufenac",
    "CHEMBL797":    "Phensuximide",
    "CHEMBL1162":   "Norethindrone",
    "CHEMBL3707377":"Segesterone Ac.",
    "CHEMBL4105630":"Zuranolone",
    "CHEMBL1200624":"Ethynodiol DA",
    "CHEMBL1201151":"Mestranol",
    "CHEMBL976":    "Praziquantel",
    "CHEMBL691":    "Ethinyl Estr.",
}

POCKET_COLORS = {
    "S1_site688":  "#1f77b4",
    "S1_site680":  "#aec7e8",
    "S1_site1813": "#d62728",
    "S3_site473":  "#2ca02c",
}

POCKET_LABELS = {
    "S1_site688":  "BIN1 site688 (primary)",
    "S1_site680":  "BIN1 site680",
    "S1_site1813": "BIN1 site1813 [dissociated]",
    "S3_site473":  "PICALM site473",
}


def make_contact_heatmap(ax, contacts, pocket, itype, title, cmap_color):
    sub = contacts[
        (contacts["pocket"] == pocket) &
        (contacts["itype"]  == itype)
    ].copy()
    if sub.empty:
        ax.text(0.5, 0.5, "No contacts", ha="center", va="center",
                transform=ax.transAxes, color="gray")
        ax.set_title(title)
        return None

    sub["compound"]  = sub["chembl_id"].map(SHORT_NAMES).fillna(sub["chembl_id"])
    sub["res_label"] = sub["resname"] + sub["resid"].astype(str)
    pivot = sub.pivot_table(index="res_label", columns="compound",
                            values="occupancy", fill_value=0)
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

    cmap = LinearSegmentedColormap.from_list("occ", ["white", cmap_color])
    im   = ax.imshow(pivot.values, cmap=cmap, vmin=0, vmax=1,
                     aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right", fontsize=9)
    ax.set_yticks(range(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=9)

    for i in range(len(pivot.index)):
        for j in range(len(pivot.columns)):
            v = pivot.values[i, j]
            if v >= 0.20:
                ax.text(j, i, f"{v*100:.0f}%", ha="center", va="center",
                        fontsize=7.5,
                        color="white" if v > 0.6 else "#333333")
    ax.set_title(title, fontweight="bold", fontsize=10)
    return im


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD publication figures (virtual screening)")
    parser.add_argument("--project-root", type=Path, default=_default_root,
                        help=f"Path to CRYPTAD project root (default: {_default_root})")
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    outdir       = project_root / "06_figures"
    outdir.mkdir(parents=True, exist_ok=True)

    mmgbsa_csv   = project_root / "04_virtual_screening" / "mmgbsa_results" / "mmgbsa_all.csv"
    admet_csv    = project_root / "04_virtual_screening" / "mmgbsa_results" / "final_hits.csv"
    contacts_csv = project_root / "03_pocket_analysis" / "binding_mode" / "contacts_all.csv"

    plt.rcParams.update({
        "font.family":    "DejaVu Sans",
        "font.size":      10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "figure.dpi":     150,
    })

    # ── Load data ─────────────────────────────────────────────────────────────
    mmgbsa   = pd.read_csv(mmgbsa_csv)
    admet    = pd.read_csv(admet_csv)
    contacts = pd.read_csv(contacts_csv)

    mmgbsa["stable"] = mmgbsa.apply(
        lambda r: (r["pocket"], r["chembl_id"]) in STABLE_HITS, axis=1)
    admet["stable"] = admet.apply(
        lambda r: (r["pocket"], r["chembl_id"]) in STABLE_HITS, axis=1)

    # ── Figure 1 — Screening Funnel ───────────────────────────────────────────
    stages = [
        ("ChEMBL approved\nsmall molecules",        3127, "#9ecae1"),
        ("CNS/BBB filter\n+ PAINS removal",         1211, "#6baed6"),
        ("Docking hits\n(≤−7.0 in ≥2 pockets)",      37, "#4292c6"),
        ("MM-GBSA complete\n(LE ≥ 0.3)",            152, "#2171b5"),
        ("ADMET pass\n(all 6 CNS criteria)",          10, "#08519c"),
        ("MD stable\nbinders",                         7, "#08306b"),
    ]

    fig, ax = plt.subplots(figsize=(7, 5))
    bar_h   = 0.55
    max_n   = stages[0][1]
    for i, (label, n, color) in enumerate(stages):
        width = n / max_n
        ax.barh(i, width, height=bar_h, color=color,
                edgecolor="white", linewidth=0.8)
        ax.text(width + 0.01, i, f"{n:,}", va="center", ha="left",
                fontsize=10, fontweight="bold", color="#222222")
        ax.text(-0.01, i, label, va="center", ha="right",
                fontsize=9, color="#333333")

    ax.set_xlim(-0.55, 1.18)
    ax.set_ylim(-0.6, len(stages) - 0.4)
    ax.set_yticks([]); ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(len(stages) - 1):
        ax.annotate("", xy=(0.0, i + 1 - bar_h / 2 - 0.05),
                    xytext=(0.0, i + bar_h / 2 + 0.05),
                    arrowprops=dict(arrowstyle="-|>", color="#666666", lw=1.2))

    ax.set_title("CRYPTAD — Virtual Screening Funnel", fontweight="bold", pad=12)
    plt.tight_layout()
    plt.savefig(outdir / "fig1_screening_funnel.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Fig 1 done")

    # ── Figure 2 — Hit properties scatter ────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for pocket, grp in mmgbsa.groupby("pocket"):
        ax.scatter(-grp["dG_mean"], grp["LE"],
                   color=POCKET_COLORS.get(pocket, "#cccccc"),
                   alpha=0.15, s=18, zorder=1)

    stable_mmgbsa = mmgbsa[mmgbsa["stable"]]
    for pocket, grp in stable_mmgbsa.groupby("pocket"):
        ax.scatter(-grp["dG_mean"], grp["LE"],
                   color=POCKET_COLORS[pocket], s=90,
                   edgecolors="black", linewidths=0.7,
                   zorder=3, label=POCKET_LABELS[pocket])
        for _, row in grp.iterrows():
            label = SHORT_NAMES.get(row["chembl_id"], row["chembl_id"])
            ax.annotate(label, (-row["dG_mean"], row["LE"]),
                        textcoords="offset points", xytext=(5, 3),
                        fontsize=7.5, color="#222222", zorder=4)

    ax.axhline(0.3, color="#aaaaaa", linestyle="--", lw=1,
               label="LE = 0.3 threshold")
    ax.set_xlabel("−ΔG$_{MM-GBSA}$ (kcal/mol)", fontsize=10)
    ax.set_ylabel("Ligand Efficiency (LE, kcal/mol/heavy atom)", fontsize=10)
    ax.set_title("MM-GBSA Hit Properties — All Pockets", fontweight="bold")
    ax.legend(fontsize=8, framealpha=0.9, loc="upper left")
    ax.set_xlim(left=0); ax.set_ylim(bottom=0)
    plt.tight_layout()
    plt.savefig(outdir / "fig2_hit_properties.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Fig 2 done")

    # ── Figure 3 — Contact heatmap: PICALM site473 ──────────────────────────���
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    fig.suptitle("PICALM site473 — Contact Occupancy (≥20%, 101 frames)",
                 fontweight="bold", fontsize=11)

    im1 = make_contact_heatmap(axes[0], contacts, "S3_site473", "Hydrophobic",
                               "Hydrophobic contacts", "#2ca02c")
    im2 = make_contact_heatmap(axes[1], contacts, "S3_site473", "HBond",
                               "H-bond contacts", "#ff7f0e")

    for im, ax in [(im1, axes[0]), (im2, axes[1])]:
        if im is not None:
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                         label="Occupancy")
    plt.tight_layout()
    plt.savefig(outdir / "fig3_contacts_site473.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Fig 3 done")

    # ── Figure 4 — Contact heatmap: BIN1 site688 + site680 ───────────────────
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    fig.suptitle("BIN1 BAR — Contact Occupancy (≥20%, 101 frames)",
                 fontweight="bold", fontsize=11)

    make_contact_heatmap(axes[0], contacts, "S1_site688", "Hydrophobic",
                         "site688 — Hydrophobic", "#1f77b4")

    sub680 = contacts[contacts["pocket"] == "S1_site680"].copy()
    sub680["compound"]  = sub680["chembl_id"].map(SHORT_NAMES).fillna(sub680["chembl_id"])
    sub680["res_label"] = sub680["resname"] + sub680["resid"].astype(str)
    if not sub680.empty:
        pivot680 = sub680.pivot_table(index="res_label", columns="itype",
                                      values="occupancy", fill_value=0)
        pivot680 = pivot680.loc[
            pivot680.max(axis=1).sort_values(ascending=False).index]
        cmap680 = LinearSegmentedColormap.from_list("occ", ["white", "#aec7e8"])
        im680   = axes[1].imshow(pivot680.values, cmap=cmap680, vmin=0, vmax=1,
                                 aspect="auto", interpolation="nearest")
        axes[1].set_xticks(range(len(pivot680.columns)))
        axes[1].set_xticklabels(pivot680.columns, rotation=20,
                                ha="right", fontsize=9)
        axes[1].set_yticks(range(len(pivot680.index)))
        axes[1].set_yticklabels(pivot680.index, fontsize=9)
        for i in range(len(pivot680.index)):
            for j in range(len(pivot680.columns)):
                v = pivot680.values[i, j]
                if v >= 0.20:
                    axes[1].text(j, i, f"{v*100:.0f}%", ha="center", va="center",
                                 fontsize=8,
                                 color="white" if v > 0.6 else "#333333")
        axes[1].set_title("site680 — Zuranolone contacts",
                          fontweight="bold", fontsize=10)
        plt.colorbar(im680, ax=axes[1], fraction=0.046, pad=0.04,
                     label="Occupancy")

    plt.tight_layout()
    plt.savefig(outdir / "fig4_contacts_bin1.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Fig 4 done")

    # ── Figure 5 — Final hits summary table ──────────────────────────────────
    admet_dict  = {row["chembl_id"]: row for _, row in admet.iterrows()}
    stable_rows = []

    for pocket, chembl_id in sorted(STABLE_HITS, key=lambda x: (x[0], x[1])):
        name = SHORT_NAMES.get(chembl_id, chembl_id)
        mg   = mmgbsa[(mmgbsa["pocket"] == pocket) &
                      (mmgbsa["chembl_id"] == chembl_id)]
        if mg.empty:
            continue
        mg = mg.iloc[0]

        hb = contacts[(contacts["pocket"]    == pocket) &
                      (contacts["chembl_id"] == chembl_id) &
                      (contacts["itype"]     == "HBond")]
        hb_str = "; ".join(
            f"{r['resname']}{r['resid']} ({r['occupancy']*100:.0f}%)"
            for _, r in hb.sort_values("occupancy", ascending=False).head(2).iterrows()
        ) or "none"

        admet_row = admet_dict.get(chembl_id, {})
        stable_rows.append({
            "Pocket":         pocket.replace("S1_", "BIN1 ").replace("S3_", "PICALM "),
            "Compound":       name,
            "dG\n(kcal/mol)": f"{mg['dG_mean']:.1f}",
            "LE":             f"{mg['LE']:.3f}",
            "CNS-MPO":        f"{mg['CNS_MPO_6']:.1f}",
            "MW":             f"{mg['MW']:.0f}",
            "BBB\n(logBB)":   f"{float(admet_row.get('logBBB_pred', 0)):.2f}",
            "Key H-bonds":    hb_str,
        })

    df_table = pd.DataFrame(stable_rows)
    fig, ax  = plt.subplots(figsize=(13, 4.2))
    ax.axis("off")

    col_widths = [0.10, 0.12, 0.08, 0.06, 0.08, 0.06, 0.09, 0.41]
    table = ax.table(
        cellText=df_table.values,
        colLabels=df_table.columns,
        cellLoc="center",
        loc="center",
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.7)

    for j in range(len(df_table.columns)):
        table[0, j].set_facecolor("#2c3e50")
        table[0, j].set_text_props(color="white", fontweight="bold")

    pocket_row_colors = {
        "PICALM site473": "#d5f5e3",
        "BIN1 site688":   "#d6eaf8",
        "BIN1 site680":   "#eaf4fb",
    }
    for i, row in enumerate(stable_rows):
        bg = pocket_row_colors.get(row["Pocket"], "#f9f9f9")
        for j in range(len(df_table.columns)):
            table[i + 1, j].set_facecolor(bg)
            table[i + 1, j].set_edgecolor("#cccccc")

    for i in range(len(stable_rows)):
        table[i + 1, len(df_table.columns) - 1].set_text_props(ha="left")

    ax.set_title("CRYPTAD — Final 7 Stable Hits: Key Properties",
                 fontweight="bold", fontsize=11, pad=14)

    legend_patches = [
        mpatches.Patch(color="#d5f5e3", label="PICALM site473"),
        mpatches.Patch(color="#d6eaf8", label="BIN1 site688"),
        mpatches.Patch(color="#eaf4fb", label="BIN1 site680"),
    ]
    ax.legend(handles=legend_patches, loc="lower right",
              fontsize=8, framealpha=0.9, bbox_to_anchor=(1, -0.02))

    plt.tight_layout()
    plt.savefig(outdir / "fig5_final_hits_table.png", dpi=300, bbox_inches="tight")
    plt.close()
    print("Fig 5 done")

    print(f"\nAll figures saved to {outdir.relative_to(project_root)}/")


if __name__ == "__main__":
    main()
