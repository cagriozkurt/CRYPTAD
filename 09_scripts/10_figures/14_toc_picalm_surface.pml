# CRYPTAD — TOC Graphic: PICALM ANTH ground-state surface
#
# Produces a muted grey molecular surface on dark navy background
# representing the "undrugged" ground state (no open cryptic pocket).
#
# Usage (command-line, no GUI):
#   pymol -cq 09_scripts/10_figures/14_toc_picalm_surface.pml
#
# Output: 06_figures/publication/TOC_picalm_ground.png
#         06_figures/publication/TOC_picalm_ground_open.png  (with pocket)
#
# After generating, open PyMOL GUI interactively to fine-tune the view:
#   pymol 09_scripts/10_figures/14_toc_picalm_surface.pml
# then rotate to taste, run: ray 1200, 900   and   png TOC_picalm_ground.png

# ── Load structure ────────────────────────────────────────────────────────────
load 01_structures/PICALM/prepared/7JXV_ANTH_prepared.pdb, anth

# Remove waters and heteroatoms
remove solvent
remove hetatm

# ── Colors ───────────────────────────────────────────────────────────────────
# Background: dark navy #0a1628
set bg_rgb, [0.039, 0.086, 0.157]

# Muted steel-grey for ground state surface
set_color muted_grey,   [0.42, 0.46, 0.54]

# Teal accent for the KKK/PIP2-helix (helix α1, resids 19–40)
# — subtle, gives the membrane-interaction face a hint of identity
set_color pip2_teal,    [0.12, 0.62, 0.70]

# ── Representation ───────────────────────────────────────────────────────────
hide everything, anth

# Surface
show surface, anth
color muted_grey, anth

# Subtle teal tint on PIP2-helix
select pip2_helix, anth and resi 19-40
color pip2_teal, pip2_helix
deselect

# ── Surface quality ───────────────────────────────────────────────────────────
set surface_quality,    1
set surface_solvent,    off
set transparency,       0.08     # very slight transparency for depth

# ── Lighting & rendering ─────────────────────────────────────────────────────
set ambient,            0.35
set direct,             0.65
set reflect,            0.25
set spec_reflect,       0.30
set spec_power,         180
set two_sided_lighting, off
set ray_shadow,         on
set ray_shadow_decay_factor, 0.1
set antialias,          2
set ray_trace_mode,     1        # 1 = normal with outlines (clean for TOC)
set ray_trace_color,    black
set ray_opaque_background, on

# ── View: orient to show membrane-interaction face ───────────────────────────
orient anth

# Tilt slightly: rotate ~15° around x to bring PIP2 helix into view
rotate x, 15
rotate y, -20
rotate z, 5

# ── Render and save ───────────────────────────────────────────────────────────
ray 1200, 900
png 06_figures/publication/TOC_picalm_ground.png, dpi=300

# ── Second render: highlight site473 pocket region (open-state hint) ─────────
# For the "open state + Ibufenac" TOC panel:
# Color site473 pocket residues amber to mark the binding site
set_color pocket_amber, [0.95, 0.72, 0.15]

select site473_resid, anth and resi 219+237+240+267+268+269+270+278
color pocket_amber, site473_resid
deselect

ray 1200, 900
png 06_figures/publication/TOC_picalm_pocket.png, dpi=300

quit
