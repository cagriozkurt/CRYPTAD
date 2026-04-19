# CRYPTAD — TOC Graphic: Ibufenac docked inside PICALM site473
#
# Source: ref_nw.pdb from complex MD rep1 (100 % engagement)
#   Chain A = PICALM ANTH (resids 19-286)
#   Chain B / MOL 287 = Ibufenac (CHEMBL341812)
#
# View 1 (full protein): semi-transparent grey surface + amber pocket patch +
#   ligand visible as bold sticks through the surface
# View 2 (close-up): pocket residues as sticks + ligand sticks, surface clipped
#
# Usage:
#   cd "/Volumes/PortableSSD/untitled folder/CRYPTAD"
#   pymol -cq 09_scripts/10_figures/16_toc_picalm_docked.pml
#
# Output: 06_figures/publication/TOC_picalm_docked.png
#         06_figures/publication/TOC_picalm_docked_closeup.png
#
# Fine-tune interactively:
#   pymol 09_scripts/10_figures/16_toc_picalm_docked.pml

# ── Load ──────────────────────────────────────────────────────────────────────
load 02_md_simulations/complex_md/S3_site473_CHEMBL341812/rep1/ref_nw.pdb, complex

alter complex, resn = 'HIS' if resn in ('HSD', 'HSE', 'HSP') else resn
rebuild
remove solvent

# Strip protein hydrogens (surface doesn't need them; removes spurious spheres)
remove complex and elem H and not resn MOL

# ── Selections ────────────────────────────────────────────────────────────────
select prot,       complex and chain A
select ibufenac,   complex and resn MOL
select pip2_helix, complex and chain A and resi 19-40
select site473,    complex and chain A and resi 219+237+240+267+268+269+270+278
select lig_heavy,  ibufenac and not elem H

# ── Colors ────────────────────────────────────────────────────────────────────
set bg_rgb, [0.039, 0.086, 0.157]

set_color muted_grey,   [0.42, 0.46, 0.54]
set_color pip2_teal,    [0.12, 0.62, 0.70]
set_color pocket_amber, [0.95, 0.72, 0.15]
set_color lig_white,    [0.95, 0.95, 0.95]

# ── Representations ───────────────────────────────────────────────────────────
hide everything, complex

# Protein surface
show surface, prot
color muted_grey,   prot
color pip2_teal,    pip2_helix
color pocket_amber, site473

set transparency,       0.10     # mostly opaque — pocket opening reveals ligand
set surface_quality,    1
set surface_solvent,    off

# Open a window into the pocket by removing surface from site473 residues
# The ligand then shows directly through the cavity opening
hide surface, site473

# Ligand: CPK spheres — fully exposed through the pocket window
show spheres, lig_heavy
set sphere_scale,  0.55

# Bright warm yellow C (contrasts grey surface + amber patch surroundings)
set_color lig_yellow, [1.00, 0.95, 0.40]
color lig_yellow,   lig_heavy
color pocket_amber, ibufenac and elem O

# ── Lighting ──────────────────────────────────────────────────────────────────
set ambient,            0.35
set direct,             0.65
set reflect,            0.25
set spec_reflect,       0.35
set spec_power,         180
set two_sided_lighting, on
set ray_shadow,         on
set ray_shadow_decay_factor, 0.1
set antialias,          2
set ray_trace_mode,     1
set ray_trace_color,    black
set ray_opaque_background, on

# ── View 1: full protein ──────────────────────────────────────────────────────
orient prot
rotate x, 15
rotate y, -20
rotate z, 5

ray 1200, 900
png 06_figures/publication/TOC_picalm_docked.png, dpi=300

# ── View 2: close-up — cartoon protein + stick pocket residues + ligand ───────
# Switch to cartoon for depth; show pocket residues and ligand as sticks
hide surface, prot
show cartoon, prot
color muted_grey, prot
color pip2_teal,  pip2_helix

# Pocket residues as sticks for atomic context
show sticks, site473
color muted_grey, site473
util.cnc site473                 # colour-by-element (grey C, red O, blue N)

# Ligand sticks — bold yellow/amber on top
set stick_radius, 0.18
color lig_yellow,   lig_heavy
color pocket_amber, ibufenac and elem O

# Cartoon quality
set cartoon_fancy_helices, 1
set cartoon_smooth_loops,  1

# Zoom centred on ligand, tight crop
zoom ibufenac, 10
rotate y, 20    # tilt to open pocket face toward viewer

ray 1200, 900
png 06_figures/publication/TOC_picalm_docked_closeup.png, dpi=300

quit
