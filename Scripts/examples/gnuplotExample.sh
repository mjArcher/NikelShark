#! /bin/bash
#param 1 = test number, param 2 = limiter
#g++ -o SlicData/runSlic SlicData/SlicSolve.cpp SlicData/originalTest.cpp
#gfortran -o ExactData/exact ExactData/exact.f

#for loop around here - for plotting 5 different tests

SLICFILE="Slic/Limiters/t_$1_$2.dat"
EXACTFILE="../Exact/exact_$1.edat"

#./ExactData/exact ExactData/exact_1.in ExactData/exact_1.edat
#./SlicData/runSlic 
#parameters can be passed into gnuplot code - see laptop code

gnuplot -persist << PLOT
set term postscript landscape enhanced colour size 10,7
set output '| ps2pdf - test_$1_$2.pdf'

set xrange [ 0.0 : 1.0 ]
set nokey
set lmargin 10
set rmargin 2
set pointsize 0.5

set multiplot
set ylabel "density"
set size 0.5,0.45
set origin 0.0,0.5
set bmargin 1
set tmargin 0
set format x ""
set xlabel ""
plot "$SLICFILE" u 1:2 w p pt 19, "$EXACTFILE" u 1:2 w l 

set ylabel "velocity"
set size 0.5,0.45
set origin 0.5,0.5
set bmargin 1
set tmargin 0
set format x ""
set xlabel ""
plot "$SLICFILE" u 1:3 w p pt 19, "$EXACTFILE" u 1:3 w l

set ylabel "pressure"
set size 0.5,0.45
set xlabel "position"
set origin 0.00,0.05
set bmargin 1
set tmargin 1
plot "$SLICFILE" u 1:4 w p pt 19, "$EXACTFILE" u 1:4 w l

set ylabel "internal energy"
set size 0.5,0.45
set xlabel "position"
set origin 0.5,0.05
set bmargin 1
set tmargin 1
plot "$SLICFILE" u 1:5 w p pt 19, "$EXACTFILE" u 1:5 w l

set nomultiplot

quit 
PLOT
