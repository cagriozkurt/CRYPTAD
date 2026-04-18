# CRYPTAD binding mode — LAROTRECTINIB SULFATE @ S1_site680
# Rep frame 67 | 101 frames | contacts ≥ 20% occupancy

load /Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis/binding_mode/S1_site680_CHEMBL3989939_rep.pdb, complex

bg_color white
set cartoon_fancy_helices, 1
hide everything
show cartoon, complex and polymer
show sticks, complex and resname MOL
color slate, complex and polymer
color yellow, complex and resname MOL
set stick_radius, 0.15

select contacts, complex and polymer and resi 192+193+194+195+196+198+199+202+295+341+345+349
show sticks, contacts
color cyan, contacts
set stick_radius, 0.12, contacts

select hbond_res, complex and polymer and resi 194+195+345
color tv_orange, hbond_res

distance hbonds, complex and resname MOL, hbond_res, 3.5, mode=2
set dash_color, firebrick
set dash_width, 2.5

label hbond_res and name CA, "%s%s" % (resn, resi)
set label_size, 11
set label_color, black

center complex and resname MOL
zoom complex and resname MOL, 8
set ray_shadows, 0
set ray_opaque_background, 1
bg_color white
ray 1200, 1200
png /Volumes/PortableSSD/untitled folder/CRYPTAD/03_pocket_analysis/binding_mode/S1_site680_CHEMBL3989939.png, dpi=300
quit
