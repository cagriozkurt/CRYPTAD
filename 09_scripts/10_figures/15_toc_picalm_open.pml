# CRYPTAD — TOC Graphic: PICALM ANTH open-pocket conformation
#
# Uses a metadynamics frame where cv_site4 > 3.5 nm (t ≈ 112340 ps),
# showing the site473 cryptic pocket physically open.
# Cavity region (pocket-lining residues) coloured amber.
#
# Usage (command-line, no GUI):
#   cd "/Volumes/PortableSSD/untitled folder/CRYPTAD"
#   pymol -cq 09_scripts/10_figures/15_toc_picalm_open.pml
#
# Output: 06_figures/publication/TOC_picalm_open.png
#
# To fine-tune view interactively:
#   pymol 09_scripts/10_figures/15_toc_picalm_open.pml
# then rotate, run: ray 1200, 900  and  png TOC_picalm_open.png

# ── Load open-state structure ─────────────────────────────────────────────────
load 01_structures/PICALM/prepared/7JXV_ANTH_open_state.pdb, anth_open

# Fix CHARMM histidine naming so PyMOL colours/selects correctly
alter anth_open, resn = 'HIS' if resn in ('HSD', 'HSE', 'HSP') else resn
rebuild

# Remove solvent (none expected, but safety net)
remove solvent

# ── Colors ────────────────────────────────────────────────────────────────────
set bg_rgb, [0.039, 0.086, 0.157]          # dark navy #0a1628

set_color muted_grey,   [0.42, 0.46, 0.54]
set_color pip2_teal,    [0.12, 0.62, 0.70]
set_color pocket_amber, [0.95, 0.72, 0.15]

# ── Representation ────────────────────────────────────────────────────────────
hide everything, anth_open

show surface, anth_open
color muted_grey, anth_open

# Teal PIP2-helix accent (same face as ground-state render)
select pip2_helix, anth_open and resi 19-40
color pip2_teal, pip2_helix
deselect

# Amber cavity: site473 pocket-lining residues
# (holo resids: LEU219, TYR237, PHE240, LEU267, SER268, GLN269, ALA270, LEU278)
select site473, anth_open and resi 219+237+240+267+268+269+270+278
color pocket_amber, site473
deselect

# ── Surface quality ───────────────────────────────────────────────────────────
set surface_quality,    1
set surface_solvent,    off
set transparency,       0.0      # fully opaque: makes the cavity edge crisp

# ── Lighting & rendering ──────────────────────────────────────────────────────
set ambient,            0.35
set direct,             0.65
set reflect,            0.25
set spec_reflect,       0.30
set spec_power,         180
set two_sided_lighting, off
set ray_shadow,         on
set ray_shadow_decay_factor, 0.1
set antialias,          2
set ray_trace_mode,     1        # normal with outlines
set ray_trace_color,    black
set ray_opaque_background, on

# ── View: match ground-state orientation ─────────────────────────────────────
orient anth_open

# Same rotation as ground-state script — keeps panels comparable in TOC
# Adjust these if the pocket cavity is not facing the viewer
rotate x, 15
rotate y, -20
rotate z, 5

# ── Render ────────────────────────────────────────────────────────────────────
ray 1200, 900
png 06_figures/publication/TOC_picalm_open.png, dpi=300

quit
