# CRYPTAD — Pocket Inspection: S5_CD2AP_SH3-1
# 4 top cross-replica Tier-2 sites (sorted: replicas desc, persistence desc)
#
# Biological context:
#   Look for sites in or near the RT-loop groove (canonical SH3 binding cleft, Pro-rich peptide face).
#
# Color guide:
#   firebrick    — Site 1 | persist=5.5% | drugg=0.066 | vol=299 A3 | run1, run2, run3
#   marine       — Site 2 | persist=5.4% | drugg=0.077 | vol=317 A3 | run1, run2, run3
#   forest       — Site 3 | persist=5.4% | drugg=0.071 | vol=315 A3 | run1, run3
#   gold         — Site 4 | persist=5.1% | drugg=0.071 | vol=304 A3 | run1, run2
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

# --- Site 1 ---  persist=5.5%  drugg=0.066  vol=299 A3  (run1, run2, run3)
load site1_best.pdb, site1_open
align site1_open, reference
show_as cartoon, site1_open
color firebrick, site1_open
set cartoon_transparency, 0.40, site1_open
pseudoatom centroid_1, pos=[29.040, 28.650, 19.880]
show spheres, centroid_1
set sphere_scale, 2.5, centroid_1
color firebrick, centroid_1
label centroid_1, "S1 (5.5%)"
select pocket1_env, (byres reference within 8 of centroid_1)
show surface, pocket1_env
color firebrick, pocket1_env
set transparency, 0.20, pocket1_env

# --- Site 2 ---  persist=5.4%  drugg=0.077  vol=317 A3  (run1, run2, run3)
load site2_best.pdb, site2_open
align site2_open, reference
show_as cartoon, site2_open
color marine, site2_open
set cartoon_transparency, 0.40, site2_open
pseudoatom centroid_2, pos=[27.030, 26.150, 36.600]
show spheres, centroid_2
set sphere_scale, 2.5, centroid_2
color marine, centroid_2
label centroid_2, "S2 (5.4%)"
select pocket2_env, (byres reference within 8 of centroid_2)
show surface, pocket2_env
color marine, pocket2_env
set transparency, 0.20, pocket2_env

# --- Site 3 ---  persist=5.4%  drugg=0.071  vol=315 A3  (run1, run3)
load site3_best.pdb, site3_open
align site3_open, reference
show_as cartoon, site3_open
color forest, site3_open
set cartoon_transparency, 0.40, site3_open
pseudoatom centroid_3, pos=[26.340, 20.780, 28.710]
show spheres, centroid_3
set sphere_scale, 2.5, centroid_3
color forest, centroid_3
label centroid_3, "S3 (5.4%)"
select pocket3_env, (byres reference within 8 of centroid_3)
show surface, pocket3_env
color forest, pocket3_env
set transparency, 0.20, pocket3_env

# --- Site 4 ---  persist=5.1%  drugg=0.071  vol=304 A3  (run1, run2)
load site4_best.pdb, site4_open
align site4_open, reference
show_as cartoon, site4_open
color gold, site4_open
set cartoon_transparency, 0.40, site4_open
pseudoatom centroid_4, pos=[32.740, 36.650, 29.050]
show spheres, centroid_4
set sphere_scale, 2.5, centroid_4
color gold, centroid_4
label centroid_4, "S4 (5.1%)"
select pocket4_env, (byres reference within 8 of centroid_4)
show surface, pocket4_env
color gold, pocket4_env
set transparency, 0.20, pocket4_env

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
# Look for sites in or near the RT-loop groove (canonical SH3 binding cleft, Pro-rich peptide face).
