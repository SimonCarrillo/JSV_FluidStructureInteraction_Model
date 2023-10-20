Main programs" to compute the natural frequencies and mode shapes in displacement and stress. Both in-vacuum and air-backed vibrations are computed in each code:

FSI_main.m
For positive r and s

FSI_main_aluminumplate_pureinplane.m
For r=0 or s=0 cases



Sub-functions called by the main program: 

For air-backed vibrations:
Fluid_Vacumm_plate_ij_zero.m (Never used; same as Vacumm_plate_ij_zero.m)
Fluid_Vacumm_plate_nonzero.m

For in-vacuum vibrations: 
Vacumm_plate_ij_zero.m
Vacumm_plate_nonzero.m


Definition of characteristic functions, F and G, for displacement and stress shape functions U and B
new_pagano_non.m

Computation of vibration frequencies, comparison between COMSOL and theory
comsol_freq.m