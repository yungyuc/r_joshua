
How To Run "CAVITY_2D_STAGGERED.F90":
     Run "PISO Code" using visual studio 2010 & intel fortran compiler. 


Files Includes With This "cavity2d":

CAVITY_2D_STAGGERED.F90 Fortran source code 
u-y_re100		         Comparison the Ghia et al.(1982) and PISO Code of u-velocity along the vertical centreline through 
                        geometric center for Re=100
x-v_re100		         Comparison the Ghia et al.(1982) and PISO Code of v-velocity along the horizontal centreline through 
                        geometric center for Re=100
schematic_diagram	    The schematic diagram of two dimension cavity flow


Design Decisions & Project Issues:
     We introduced a two-dimensional Pressure-Implicit with Splitting Operators (PISO) scheme to simulate the laminar incompressible flow in a unit square cavity whose top wall moves with a uniform velocity serves as a model problem to 
test and evaluate the unified wall-boundary condition.

Analysis Results
     The results "u-y_re100" & "x-v_re100" compare well with the numerical results obtained by Ghia et al. (1982).
