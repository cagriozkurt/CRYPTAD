# ============================================================
# CRYPTAD — Apo vs Cryptic: BIN1 BAR domain
# Apo:     2FIC (crystal, tau-axis closed state)
# Cryptic: site688_frame302.pdb (t=151 ns, vol=286 Å³)
#          Primary lead — BIN1/tau-binding interface
# Pocket centroid: (78.03, 92.58, 103.52)
# ============================================================

reinitialize
bg_color white

# ── Load structures ──────────────────────────────────────────
load /Volumes/PortableSSD/untitled folder/CRYPTAD/01_structures/BIN1/pdb/2FIC.pdb, apo
load /Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis/cryptic_inspection/S1_BIN1_BAR/site688_frame302.pdb, cryptic

# ── Align apo onto cryptic (apo moves, cryptic is reference) ─
cealign cryptic and name CA, apo and name CA

# ── Global display settings ──────────────────────────────────
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_tube_radius, 0.3
set cartoon_oval_length, 1.2
set depth_cue, 0
set fog, 0
set ray_shadows, 0
set ambient, 0.6
set direct, 0.8
set specular, 0.15
set shininess, 5
set surface_quality, 1
set mesh_quality, 2

# ── Apo structure (closed state, neutral) ────────────────────
hide everything, apo
show cartoon, apo
color gray70, apo
set cartoon_transparency, 0.6, apo

# ── Cryptic structure (open state, warm palette) ─────────────
hide everything, cryptic

# Remove N/C-terminal capping groups (CAY/CAY added by pdb2gmx)
select cryptic_protein, cryptic and (name CA+C+N+O+CB)
show cartoon, cryptic

# Chain A: warm coral  (the tau-interaction helix face)
color salmon, cryptic and chain A
# Chain B: complementary slate
color lightblue, cryptic and chain B

# ── Pocket cavity ────────────────────────────────────────────
# Semi-transparent sphere at centroid = approximate cavity volume (286 Å³, r≈4.1 Å)
pseudoatom pocket_cen, pos=[78.03, 92.58, 103.52]
hide everything, pocket_cen
show spheres, pocket_cen
set sphere_scale, 4.1, pocket_cen
color tv_orange, pocket_cen
set sphere_transparency, 0.45, pocket_cen

# ── Pocket-lining residues (within 5 Å of centroid) ──────────
select pocket_res, (cryptic) and (byres (all within 5.0 of pocket_cen))
show sticks, pocket_res
show sticks, pocket_res and (not name C+N+O+CA)   # sidechains cleaner
color tv_orange, pocket_res
set stick_radius, 0.15

# Backbone of pocket region: slightly highlighted
show cartoon, pocket_res
set cartoon_transparency, 0.0, pocket_res
color brightorange, pocket_res and chain A
color limon, pocket_res and chain B

# ── Apo: show the SAME residues as thin sticks for contrast ──
# (same sequence positions — shows that pocket is absent in apo)
select apo_same_res, apo and (byres (all within 5.0 of pocket_cen))
show sticks, apo_same_res and (not name C+N+O+CA)
color gray50, apo_same_res
set stick_radius, 0.10

# ── Surface — cryptic only, pocket region ────────────────────
# Solvent-accessible surface for just the pocket shell residues
create pocket_shell, cryptic and (byres (all within 7.0 of pocket_cen))
show surface, pocket_shell
set surface_color, tv_orange, pocket_shell
set transparency, 0.55, pocket_shell
set surface_mode, 3    # by_object atom mode for cleaner look

# ── Labels ───────────────────────────────────────────────────
# Label the pocket centroid
pseudoatom lbl_pocket, pos=[78.03, 99.0, 103.52]
label lbl_pocket, "site688 (286 Å³)"
set label_size, 14
set label_color, tv_orange
hide everything, lbl_pocket
show labels, lbl_pocket

# ── Cleanup selections ───────────────────────────────────────
deselect

# ── View ─────────────────────────────────────────────────────
# Zoom on pocket region with context
zoom pocket_cen, 25
# Fine-tune: rotate to show pocket opening clearly
turn y, 15
turn x, -10

# Save this as a scene
scene apo_vs_cryptic, store

# ── Ray-trace and save ───────────────────────────────────────
# Uncomment to auto-render at 300 dpi:
# set ray_opaque_background, 1
# ray 2400, 1800
# png /Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis/cryptic_inspection/S1_BIN1_BAR/apo_vs_cryptic_BIN1_BAR.png, dpi=300

# ── Legend (as comments for figure caption) ──────────────────
# Gray transparent cartoon    = 2FIC apo (crystal structure, closed)
# Salmon/lightblue cartoon    = site688_frame302 cryptic (open, t=151 ns)
# Orange sphere (transparent) = site688 pocket cavity (286 Å³)
# Orange sticks               = pocket-lining residues (cryptic)
# Gray sticks                 = same residues in apo (no pocket)
