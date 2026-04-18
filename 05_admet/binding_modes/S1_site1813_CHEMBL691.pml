# CRYPTAD binding mode — ETHINYL ESTRADIOL @ S1_site1813
# Rep frame 73 | 101 frames | contacts ≥ 20% occupancy

load /Volumes/PortableSSD/untitled folder/CRYPTAD/05_admet/binding_modes/S1_site1813_CHEMBL691_rep.pdb, complex

bg_color white
set cartoon_fancy_helices, 1
hide everything
show cartoon, complex and polymer
show sticks, complex and resname MOL LIG
color slate, complex and polymer
color yellow, complex and resname MOL LIG
set stick_radius, 0.15

select contacts, complex and polymer and resi none
show sticks, contacts
color cyan, contacts
set stick_radius, 0.12, contacts

select hbond_res, complex and polymer and resi 172
color tv_orange, hbond_res

distance hbonds, complex and resname MOL LIG, hbond_res, 3.5, mode=2
set dash_color, firebrick
set dash_width, 2.5

label hbond_res and name CA, "%s%s" % (resn, resi)
set label_size, 11
set label_color, black

center complex and resname MOL LIG
zoom complex and resname MOL LIG, 8
set ray_shadows, 0
