# ============================================================
# CRYPTAD — Apo vs Cryptic: PICALM ANTH domain
# Apo:     7JXV chain A (crystal, Aβ-axis closed state)
#          [chain B = ubiquitin, excluded from overlay]
# Cryptic: site473_frame384.pdb (t=192 ns, vol=271 Å³)
#          Confirmed lead — PICALM/Aβ-endocytosis interface
# Pocket centroid: (65.4, 64.1, 59.5)
# ============================================================

reinitialize
bg_color white

# ── Load structures ──────────────────────────────────────────
load /Volumes/PortableSSD/untitled folder/CRYPTAD/01_structures/PICALM/pdb/7JXV.pdb, apo_7jxv
load /Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis/cryptic_inspection/S3_PICALM_ANTH/site473_frame384.pdb, cryptic

# Isolate ANTH domain from 7JXV (chain A only; chain B = ubiquitin)
create apo, apo_7jxv and chain A
delete apo_7jxv

# ── Align apo onto cryptic ───────────────────────────────────
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

# ── Cryptic structure (open state, cool palette) ─────────────
hide everything, cryptic
show cartoon, cryptic

# PICALM ANTH: cool teal palette  (Aβ axis)
color teal, cryptic

# ── Pocket cavity ────────────────────────────────────────────
# Semi-transparent sphere at centroid (vol=271 Å³, r≈4.0 Å)
pseudoatom pocket_cen, pos=[65.4, 64.1, 59.5]
hide everything, pocket_cen
show spheres, pocket_cen
set sphere_scale, 4.0, pocket_cen
color cyan, pocket_cen
set sphere_transparency, 0.45, pocket_cen

# ── Pocket-lining residues (within 5 Å of centroid) ──────────
select pocket_res, (cryptic) and (byres (all within 5.0 of pocket_cen))
show sticks, pocket_res
color cyan, pocket_res
set stick_radius, 0.15

# Backbone of pocket region: highlighted (no transparency)
show cartoon, pocket_res
set cartoon_transparency, 0.0, pocket_res
color aquamarine, pocket_res

# ── Apo: same residues as thin gray sticks for contrast ──────
select apo_same_res, apo and (byres (all within 5.0 of pocket_cen))
show sticks, apo_same_res and (not name C+N+O+CA)
color gray50, apo_same_res
set stick_radius, 0.10

# ── Surface — cryptic pocket shell ───────────────────────────
create pocket_shell, cryptic and (byres (all within 7.0 of pocket_cen))
show surface, pocket_shell
set surface_color, cyan, pocket_shell
set transparency, 0.55, pocket_shell
set surface_mode, 3

# ── Labels ───────────────────────────────────────────────────
pseudoatom lbl_pocket, pos=[65.4, 70.5, 59.5]
label lbl_pocket, "site473 (271 Å³)"
set label_size, 14
set label_color, teal
hide everything, lbl_pocket
show labels, lbl_pocket

# ── Cleanup ──────────────────────────────────────────────────
deselect

# ── View ─────────────────────────────────────────────────────
zoom pocket_cen, 22
turn y, -15
turn x, 10

scene apo_vs_cryptic, store

# ── Ray-trace and save ───────────────────────────────────────
# set ray_opaque_background, 1
# ray 2400, 1800
# png /Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis/cryptic_inspection/S3_PICALM_ANTH/apo_vs_cryptic_PICALM_ANTH.png, dpi=300

# ── Legend ───────────────────────────────────────────────────
# Gray transparent cartoon   = 7JXV apo chain A (crystal, closed)
# Teal cartoon               = site473_frame384 cryptic (open, t=192 ns)
# Cyan sphere (transparent)  = site473 pocket cavity (271 Å³)
# Cyan/aquamarine sticks     = pocket-lining residues (cryptic)
# Gray sticks                = same residues in apo (no pocket)
