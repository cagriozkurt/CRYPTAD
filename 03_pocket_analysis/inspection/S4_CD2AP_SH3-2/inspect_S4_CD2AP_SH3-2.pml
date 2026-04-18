# CRYPTAD — Pocket Inspection: S4_CD2AP_SH3-2
# 2 top cross-replica Tier-2 sites (sorted: replicas desc, persistence desc)
#
# Biological context:
#   Look for sites in or near the RT-loop groove (canonical SH3 binding cleft, Pro-rich peptide face).
#
# Color guide:
#   firebrick    — Site 1 | persist=6.2% | drugg=0.103 | vol=315 A3 | run1, run2, run3
#   marine       — Site 2 | persist=5.2% | drugg=0.072 | vol=295 A3 | run2, run3
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

# --- Site 1 ---  persist=6.2%  drugg=0.103  vol=315 A3  (run1, run2, run3)
load site1_best.pdb, site1_open
align site1_open, reference
show_as cartoon, site1_open
color firebrick, site1_open
set cartoon_transparency, 0.40, site1_open
pseudoatom centroid_1, pos=[37.000, 27.260, 30.560]
show spheres, centroid_1
set sphere_scale, 2.5, centroid_1
color firebrick, centroid_1
label centroid_1, "S1 (6.2%)"
select pocket1_env, (byres reference within 8 of centroid_1)
show surface, pocket1_env
color firebrick, pocket1_env
set transparency, 0.20, pocket1_env

# --- Site 2 ---  persist=5.2%  drugg=0.072  vol=295 A3  (run2, run3)
load site2_best.pdb, site2_open
align site2_open, reference
show_as cartoon, site2_open
color marine, site2_open
set cartoon_transparency, 0.40, site2_open
pseudoatom centroid_2, pos=[28.400, 20.880, 26.860]
show spheres, centroid_2
set sphere_scale, 2.5, centroid_2
color marine, centroid_2
label centroid_2, "S2 (5.2%)"
select pocket2_env, (byres reference within 8 of centroid_2)
show surface, pocket2_env
color marine, pocket2_env
set transparency, 0.20, pocket2_env

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
