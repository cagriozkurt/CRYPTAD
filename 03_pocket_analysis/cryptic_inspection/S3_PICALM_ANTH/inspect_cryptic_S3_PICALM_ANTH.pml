# CRYPTAD — Cryptic pocket inspection
# System: S3_PICALM_ANTH
# Colour legend: magenta sphere = cryptic site centroid
reinitialize

load /Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis/cryptic_inspection/S3_PICALM_ANTH/site473_frame384.pdb, site473
hide everything, site473
show cartoon, site473
color slate, site473
pseudoatom crypt_site473, pos=[65.41,64.13,59.47]
show spheres, crypt_site473
set sphere_scale, 3.0, crypt_site473
color magenta, crypt_site473

bg_color white
set ray_shadows, 0
orient
