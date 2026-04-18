from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/site473_frame384.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/site473_frame384.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [1008,1077,685,690,1197,763,1381,1379,708,710,712,698,700,702,657,704,707,714,716,1247,1271,1313,650,652,1308,1249,1847,1852,1078,1363,1838,1377,1906,1910] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.302,0.278,0.702]
select surf_pocket2, protein and id [1894,2342,2345,1923,2007,2009,2014,1976,1920,2207,1979,1984,1872,2538,1890,1899,2037,2043,2059,2063,2034,2057,2113,2471,2420,2283,2422,2464] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.631,0.361,0.902]
select surf_pocket3, protein and id [3111,4047,3073,3130,3131,3902,3925,3853,4025,3856,3858,2500,2503,2507,3113,3118,2569,3608,3610,3619,3688,3600,3602,3844,3582,3588,3599,2548,2555,3077,2512,2498] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.678,0.278,0.702]
select surf_pocket4, protein and id [1048,1088,1094,985,973,2906,966,980,1012,1013,1090,2971,2981,2984,2986,2985,2977,2914,2916,2918,2915,2904,1755,1760,1829,2073,2075] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.682]
select surf_pocket5, protein and id [2555,2596,2998,3065,3068,3075,3077,2943,2945,2872,2587,2590,2600,2602,2604,2592,2601,2673,2870,2560,2993,2923,2935,2937,2997,2611,3902] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.341]
select surf_pocket6, protein and id [3616,3848,4056,4068,4085,3801,3806,3815,3547,4089,3629,3693,4163,4094,3708,3712,4159] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.522,0.361]
select surf_pocket7, protein and id [2657,2726,2661,2733,2637,2653,1650,1729,1734,2666,2750,1660,1664,1666,2641,1613,1646,1656,1658,2820,1604,1606,1609,1610,2868] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.596,0.278]
select surf_pocket8, protein and id [4298,4302,3325,4258,4246,3308,4255,3437,3305,3354,3324] 
set surface_color,  pcol8, surf_pocket8 
   

deselect

orient
