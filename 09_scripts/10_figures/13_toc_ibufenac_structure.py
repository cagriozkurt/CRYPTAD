"""
CRYPTAD — TOC Graphic: Ibufenac 2D chemical structure

Generates a clean vector-style 2D depiction of Ibufenac
(CHEMBL341812; p-isobutylphenylacetic acid; CC(C)Cc1ccc(CC(=O)O)cc1)
for use as a TOC graphic component.

Color scheme (matches TOC FES and scatter panels):
  Background : dark navy   #0a1628
  C/C bonds  : white
  O atoms    : amber/gold  (matches FES minimum colour)
  No explicit H (except implicit on OH shown via label)

Outputs (06_figures/publication/):
  TOC_ibufenac.svg           — vector, dark navy background
  TOC_ibufenac.png           — rasterised, dark navy background, 300 dpi
  TOC_ibufenac_light.svg     — vector, white background (for reference)
  TOC_ibufenac_light.png

Usage:
  python3 09_scripts/10_figures/13_toc_ibufenac_structure.py
  python3 09_scripts/10_figures/13_toc_ibufenac_structure.py --project-root /path
"""

import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

_script       = Path(__file__).resolve()
_default_root = _script.parents[2]

SMILES = "CC(C)Cc1ccc(CC(=O)O)cc1"

# Colours as (R, G, B, A) floats 0–1
NAVY    = (0.039, 0.086, 0.157, 1.0)
WHITE   = (1.000, 1.000, 1.000, 1.0)
AMBER   = (0.950, 0.720, 0.150, 1.0)   # matches FES global-minimum colour
TRANSPARENT = (0.0, 0.0, 0.0, 0.0)

# Canvas size in pixels (vector SVG — pixel size only affects viewport)
W, H = 560, 380


def make_mol():
    mol = Chem.MolFromSmiles(SMILES)
    rdDepictor.Compute2DCoords(mol)
    rdDepictor.NormalizeDepiction(mol)
    # Flip horizontally so isobutyl group faces left (aesthetic preference)
    rdDepictor.StraightenDepiction(mol)
    return mol


def draw(mol, bg_colour: tuple, bond_colour: tuple,
         o_colour: tuple, width=W, height=H) -> tuple[str, bytes]:
    """
    Returns (svg_text, png_bytes) for the given colour scheme.
    """
    # ── SVG ──────────────────────────────────────────────────────────────────
    svg_drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = svg_drawer.drawOptions()
    opts.backgroundColour         = bg_colour
    opts.bondLineWidth             = 2.8
    opts.baseFontSize              = 0.70     # relative to bond length
    opts.padding                   = 0.12
    opts.addAtomIndices            = False
    opts.addStereoAnnotation       = False
    opts.updateAtomPalette({
        6:  bond_colour[:3],   # C → bond_colour (white or black)
        8:  o_colour[:3],      # O → amber or standard red
        1:  bond_colour[:3],   # H (if shown)
    })
    svg_drawer.DrawMolecule(mol)
    svg_drawer.FinishDrawing()
    svg_text = svg_drawer.GetDrawingText()

    # ── PNG (Cairo) ───────────────────────────────────────────────────────────
    png_drawer = rdMolDraw2D.MolDraw2DCairo(width * 2, height * 2)  # 2× for 300 dpi
    opts2 = png_drawer.drawOptions()
    opts2.backgroundColour         = bg_colour
    opts2.bondLineWidth             = 2.8
    opts2.baseFontSize              = 0.70
    opts2.padding                   = 0.12
    opts2.addAtomIndices            = False
    opts2.addStereoAnnotation       = False
    opts2.updateAtomPalette({
        6:  bond_colour[:3],
        8:  o_colour[:3],
        1:  bond_colour[:3],
    })
    png_drawer.DrawMolecule(mol)
    png_drawer.FinishDrawing()
    png_bytes = png_drawer.GetDrawingText()

    return svg_text, png_bytes


def main():
    parser = argparse.ArgumentParser(
        description="CRYPTAD TOC — Ibufenac 2D structure")
    parser.add_argument("--project-root", type=Path, default=_default_root)
    args   = parser.parse_args()
    root   = args.project_root.resolve()
    outdir = root / "06_figures/publication"
    outdir.mkdir(parents=True, exist_ok=True)

    mol = make_mol()
    print(f"Molecule: {Chem.MolToSmiles(mol)}")
    print(f"Heavy atoms: {mol.GetNumAtoms()}")

    # Dark version (TOC)
    svg, png = draw(mol, bg_colour=NAVY, bond_colour=WHITE, o_colour=AMBER)
    (outdir / "TOC_ibufenac.svg").write_text(svg)
    (outdir / "TOC_ibufenac.png").write_bytes(png)
    print("Saved: TOC_ibufenac.svg / .png")

    # Light reference version
    RED = (0.80, 0.10, 0.10, 1.0)
    svg_l, png_l = draw(mol, bg_colour=(1,1,1,1),
                        bond_colour=(0,0,0,1), o_colour=RED)
    (outdir / "TOC_ibufenac_light.svg").write_text(svg_l)
    (outdir / "TOC_ibufenac_light.png").write_bytes(png_l)
    print("Saved: TOC_ibufenac_light.svg / .png")

    print("Done.")


if __name__ == "__main__":
    main()
