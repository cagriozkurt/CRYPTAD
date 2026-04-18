# CRYPTAD — Pocket Inspection: S3_PICALM_ANTH
# 5 top cross-replica Tier-2 sites (sorted: replicas desc, persistence desc)
#
# Biological context:
#   Look for sites near the PI(4,5)P2-binding loops (helix α0, β1-α1 linker) and the clathrin-box region.
#
# Color guide:
#   firebrick    — Site 1 | persist=8.1% | drugg=0.058 | vol=274 A3 | run1, run2, run3
#   marine       — Site 2 | persist=8.0% | drugg=0.049 | vol=290 A3 | run1, run2, run3
#   forest       — Site 3 | persist=7.9% | drugg=0.083 | vol=333 A3 | run1, run2, run3
#   gold         — Site 4 | persist=7.4% | drugg=0.040 | vol=277 A3 | run1, run2, run3
#   purple       — Site 5 | persist=7.4% | drugg=0.039 | vol=296 A3 | run1, run2, run3
#
# Usage: File > Run Script... > select this .pml
#
# Quick exploration commands:
#   hide all, open_frames      show only the reference + centroid markers
#   show sticks, pocket1_env   show stick repr for site-1 pocket residues
#   set sphere_scale, 1.0, centroids   smaller centroid spheres
#   color white, reference     make protein background white for figures

reinitialize

# Reference structure (protein_ref.pdb, frame 0 of production)
load ref.pdb, reference
show_as cartoon, reference
color gray70, reference
show surface, reference
set surface_color, gray85, reference
set transparency, 0.55, reference
set cartoon_smooth_loops, 1, reference

# --- Site 1 ---  persist=8.1%  drugg=0.058  vol=274 A3  (run1, run2, run3)
load site1_best.pdb, site1_open
align site1_open, reference
show_as cartoon, site1_open
color firebrick, site1_open
set cartoon_transparency, 0.40, site1_open
pseudoatom centroid_1, pos=[41.260, 42.150, 42.590]
show spheres, centroid_1
set sphere_scale, 2.5, centroid_1
color firebrick, centroid_1
label centroid_1, "S1 (8.1%)"
select pocket1_env, (byres reference within 8 of centroid_1)
show surface, pocket1_env
color firebrick, pocket1_env
set transparency, 0.20, pocket1_env

# --- Site 2 ---  persist=8.0%  drugg=0.049  vol=290 A3  (run1, run2, run3)
load site2_best.pdb, site2_open
align site2_open, reference
show_as cartoon, site2_open
color marine, site2_open
set cartoon_transparency, 0.40, site2_open
pseudoatom centroid_2, pos=[47.060, 37.550, 43.970]
show spheres, centroid_2
set sphere_scale, 2.5, centroid_2
color marine, centroid_2
label centroid_2, "S2 (8.0%)"
select pocket2_env, (byres reference within 8 of centroid_2)
show surface, pocket2_env
color marine, pocket2_env
set transparency, 0.20, pocket2_env

# --- Site 3 ---  persist=7.9%  drugg=0.083  vol=333 A3  (run1, run2, run3)
load site3_best.pdb, site3_open
align site3_open, reference
show_as cartoon, site3_open
color forest, site3_open
set cartoon_transparency, 0.40, site3_open
pseudoatom centroid_3, pos=[38.580, 54.420, 49.900]
show spheres, centroid_3
set sphere_scale, 2.5, centroid_3
color forest, centroid_3
label centroid_3, "S3 (7.9%)"
select pocket3_env, (byres reference within 8 of centroid_3)
show surface, pocket3_env
color forest, pocket3_env
set transparency, 0.20, pocket3_env

# --- Site 4 ---  persist=7.4%  drugg=0.040  vol=277 A3  (run1, run2, run3)
load site4_best.pdb, site4_open
align site4_open, reference
show_as cartoon, site4_open
color gold, site4_open
set cartoon_transparency, 0.40, site4_open
pseudoatom centroid_4, pos=[49.250, 55.700, 53.120]
show spheres, centroid_4
set sphere_scale, 2.5, centroid_4
color gold, centroid_4
label centroid_4, "S4 (7.4%)"
select pocket4_env, (byres reference within 8 of centroid_4)
show surface, pocket4_env
color gold, pocket4_env
set transparency, 0.20, pocket4_env

# --- Site 5 ---  persist=7.4%  drugg=0.039  vol=296 A3  (run1, run2, run3)
load site5_best.pdb, site5_open
align site5_open, reference
show_as cartoon, site5_open
color purple, site5_open
set cartoon_transparency, 0.40, site5_open
pseudoatom centroid_5, pos=[39.750, 47.820, 54.400]
show spheres, centroid_5
set sphere_scale, 2.5, centroid_5
color purple, centroid_5
label centroid_5, "S5 (7.4%)"
select pocket5_env, (byres reference within 8 of centroid_5)
show surface, pocket5_env
color purple, pocket5_env
set transparency, 0.20, pocket5_env

# Display settings
set ray_shadow, 0
set bg_rgb, white
set depth_cue, 0
set specular, 0.3
set cartoon_fancy_helices, 1
set label_size, 14
set label_color, black
zoom all, buffer=5

# Group objects for easy show/hide in the PyMOL object panel
group centroids, centroid_*
group pocket_envs, pocket*_env
group open_frames, site*_open

# System-specific orientation hint:
# Look for sites near the PI(4,5)P2-binding loops (helix α0, β1-α1 linker) and the clathrin-box region.
