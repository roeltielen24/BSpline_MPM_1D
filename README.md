These files contain a working code for simulating very simple soil-deformation problems using Powell-Sabin spline MPM on triangulations. 
The main file to be run are Oedometer_Bspline_2D or Vib_bar_Bspline_2D, for a soil column under self weight or a vibrating bar with fixed ends respectively.
The time stepping is performed in MPM_2D_B_spline. From this script, all the other scripts are called. 
Finally, the file Basis_functions_Bezier_ordinates creates the PS-splines from a triangulation, 
and value_Bpline2D can evaluate these basis functions in the positions of material points.   
