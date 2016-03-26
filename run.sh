#! /bin/bash
#
CELLREF=1000
LIMITER=min-bee
# kevin-1D (reference solution 1)
~/solid-1D/kevinref/Elastic1DUnigrid barton2 $CELLREF

# # cns-1D (reference solution 2)
~/projects/1D_Riemann/1D_RiemannSingle.sh $CELLREF $LIMITER

# solid-1D (my) code
~/solid-1D/src/Elastic1D ~/solid-1D/settings/barton1D.cfg 

# plotting
~/solid-1D/plot_solid-1D.sh $CELLREF

# 10000 cells takes ~ 273 seconds on 4 cores
