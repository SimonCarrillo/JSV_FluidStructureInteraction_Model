Code for reproducing results in Simon Carrillo Segura, Peng Zhang, Maurizio Porfiri,
Three-dimensional exact solution of free vibrations of a simply supported rectangular plate in contact with a fluid,
Journal of Sound and Vibration,
Volume 534,
2022

Manuscript can be found here: https://doi.org/10.1016/j.jsv.2022.117007

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
