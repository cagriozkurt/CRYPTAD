"""
CRYPTAD — Assemble three-panel pLDDT Figure S1.

Reads AF3 model_0 CIF files for BIN1, PICALM, and CD2AP, extracts per-residue
pLDDT from Cα B-factors, and produces a single three-panel figure with panel
labels (a), (b), (c) positioned inside the top-left corner of each panel to
avoid overlap with y-axis tick labels.

Outputs:
  06_figures/figS2_plddt.png              (300 dpi, for local reference)
  06_figures/publication/Figure_S1_plddt.png   (300 dpi)
  06_figures/publication/Figure_S1_plddt.tif   (300 dpi)

Usage:
  python3 09_scripts/10_figures/08_assemble_plddt_panel.py [--project-root PATH]
"""

import argparse
from pathlib import Path

import gemmi
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

_SCRIPT_PATH = Path(__file__).resolve()
DEFAULT_ROOT  = _SCRIPT_PATH.parents[2]

TARGETS = [
    ("BIN1",   "fold_cryptad_bin1_isoform1_model_0.cif"),
    ("PICALM", "fold_cryptad_picalm_full_model_0.cif"),
    ("CD2AP",  "fold_cryptad_cd2ap_sh3tandem_model_0.cif"),
]

PANEL_LABELS = ["(a)", "(b)", "(c)"]


def extract_plddt(cif_path: Path):
    """Return (residue_numbers, plddt_values) from Cα B-factors in a CIF."""
    st = gemmi.read_structure(str(cif_path))
    residues, plddt = [], []
    for model in st:
        for chain in model:
            for res in chain:
                if res.is_water():
                    continue
                for atom in res:
                    if atom.name == "CA":
                        residues.append(res.seqid.num)
                        plddt.append(atom.b_iso)
                        break
    return np.array(residues), np.array(plddt)


def draw_plddt_panel(ax, target, residues, plddt, panel_label):
    """Draw a single pLDDT panel on ax with the given panel label."""
    # Confidence colour bands (matching AlphaFold Database viewer)
    ax.axhspan(90, 100, alpha=0.15, color="#0053D6")
    ax.axhspan(70,  90, alpha=0.15, color="#65CBF3")
    ax.axhspan(50,  70, alpha=0.15, color="#FFDB13")
    ax.axhspan( 0,  50, alpha=0.15, color="#FF7D45")

    ax.plot(residues, plddt, color="#333333", lw=1.2, zorder=3)
    ax.axhline(70, color="steelblue", lw=0.8, ls="--", alpha=0.6)
    ax.axhline(50, color="orange",    lw=0.8, ls="--", alpha=0.6)

    ax.set_xlim(residues[0], residues[-1])
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue number", fontsize=11)
    ax.set_ylabel("pLDDT", fontsize=11)
    ax.set_title(f"{target} — AlphaFold3 per-residue confidence (pLDDT)",
                 fontsize=11, fontweight="bold")

    # Shade disordered stretches (pLDDT < 50)
    low_mask = plddt < 50
    if low_mask.any():
        runs, in_run, start = [], False, 0
        for i, v in enumerate(low_mask):
            if v and not in_run:
                in_run, start = True, residues[i]
            elif not v and in_run:
                in_run = False
                runs.append((start, residues[i - 1]))
        if in_run:
            runs.append((start, residues[-1]))
        for s, e in runs:
            ax.axvspan(s, e, alpha=0.25, color="red", zorder=2)
            ax.text((s + e) / 2, 5, f"{s}–{e}", ha="center",
                    fontsize=7, color="darkred")

    # Panel label — placed inside the top-left corner of the axes to avoid
    # overlapping y-axis tick labels (axes-coordinate positioning)
    ax.text(0.01, 0.97, panel_label,
            transform=ax.transAxes,
            fontsize=13, fontweight="bold",
            va="top", ha="left",
            color="black")

    return ax


def main():
    parser = argparse.ArgumentParser(
        description="Assemble three-panel pLDDT Figure S1 for CRYPTAD.")
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()

    root      = args.project_root.resolve()
    pub_dir   = root / "06_figures" / "publication"
    figs_dir  = root / "06_figures"
    pub_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9))
    fig.subplots_adjust(hspace=0.45)

    for ax, (target, cif_name), label in zip(axes, TARGETS, PANEL_LABELS):
        cif_path = root / "01_structures" / target / "alphafold" / cif_name
        if not cif_path.exists():
            raise FileNotFoundError(f"CIF not found: {cif_path}")
        residues, plddt = extract_plddt(cif_path)
        draw_plddt_panel(ax, target, residues, plddt, label)

    # Shared legend — bottom-right of the figure
    legend_patches = [
        mpatches.Patch(color="#0053D6", alpha=0.4, label="Very high >90"),
        mpatches.Patch(color="#65CBF3", alpha=0.4, label="Confident 70–90"),
        mpatches.Patch(color="#FFDB13", alpha=0.4, label="Low 50–70"),
        mpatches.Patch(color="#FF7D45", alpha=0.4, label="Very low <50"),
    ]
    axes[-1].legend(handles=legend_patches, loc="lower right",
                    fontsize=8, framealpha=0.9)

    # Save
    for out_path in [
        figs_dir  / "figS2_plddt.png",
        pub_dir   / "Figure_S1_plddt.png",
        pub_dir   / "Figure_S1_plddt.tif",
    ]:
        fig.savefig(out_path, dpi=600, bbox_inches="tight")
        print(f"Saved: {out_path.relative_to(root)}")

    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
