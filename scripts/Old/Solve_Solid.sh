#! /bin/bash

if [ $1 -eq 1 ];
then 
# make ./kevinElastic/Elastic1DUnigrid
# rm -rf ./kevinElastic/out/Barton1_200/*
# ./Elastic1DUnigrid barton1 200
g++ -O3 -I ~/Libraries/eigen/ -I ~/Libraries/ UtilityFunctions.cpp SlicSolve.cpp -o ss
./ss
fi 
#gnuplot output


# plot "$RES1" u 1:2 w p pt 19 ps 0.4 lw 0.2 lc rgb "blue", "$RES2" u 1:2 w p pt 19 ps 0.4 lw 0.2 linecolor rgb "green", "$RES3" u 1:2 w p pt 19 ps 0.4 lw 0.2 linecolor rgb "red", "$EXACT" u 1:2 w l lw 2 lt 1 linecolor rgb "black"

fileBase="/lsc/zeushome/ma595/Dropbox/2013-2014/Code/Solid/"
output=$fileBase"ElasticTest.pdf"
#initial
out1D200_ini=$fileBase"Output/out1D200.dat"
out1D200_k_ini=$fileBase"KevinElastic/out/Barton1_200/out000.tec"
#final
out1D200=$fileBase"Output/out1D200.dat"
out1D200_k=$fileBase"KevinElastic/out/Barton1_200/out001.tec"

echo "looking in "$out1D200" for output"
gnuplot -persist << PLOT
set term postscript portrait enhanced colour size 7,10
set output '| ps2pdf - $output'
set style line 2 lc rgb '#5e9c36' pt 6 ps 1.5 lt 1 lw 5 # --- green

set xrange [ 0.0 : 0.01 ]
set nokey
set lmargin 10
set rmargin 2
set pointsize 0.5

set style line 1 lt 2 lw 2 ps 0.5 pt 19 lc rgb "red"
# set style line 1 lt 2 lw 2 pt 3 ps 0.5 pt 19 lc rgb "black"
# set style line 3 lt 1 lw 5 lc rgb "red"
# set style line 4 lt 1 lc rgb "black"
# set style line 3 
# set style line 4 lt 1 lc rgb "black"

set multiplot
set ylabel "density"
set xlabel "position"
set size 0.5,0.30
set origin 0.0,0.66
plot "$out1D200_ini" u 1:2 w p pt 19 lc rgb "red", "$out1D200_k_ini" u 1:2 w l ls 2 

set size 0.5,0.30
set ylabel "u1"
set size 0.5,0.3
set origin 0.5,0.66
plot "$out1D200_ini" u 1:3 w p pt 19 lc rgb "red", "$out1D200_k_ini" u 1:3 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "u2"
set origin 0.00,0.33
plot "$out1D200_ini" u 1:4 w p pt 19 lc rgb "red", "$out1D200_k_ini" u 1:4 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "u3"
set origin 0.5,0.33
plot "$out1D200_ini" u 1:5 w p pt 19 lc rgb "red", "$out1D200_k_ini" u 1:5 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "sigma11"
set origin 0.00,0.05
plot "$out1D200_ini" u 1:6 w p pt 19 lc rgb "red", "$out1D200_k_ini" u 1:6 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "sigma12"
set origin 0.5,0.05
plot "$out1D200_ini" u 1:7 w p pt 19 lc rgb "red", "$out1D200_k_ini" u 1:7 w p pt 19 lc rgb "black"

set nomultiplot
set multiplot
# set output '| ps2pdf - $output'

set ylabel "density"
set xlabel "position"
set size 0.5,0.30
set origin 0.0,0.66
plot "$out1D200" u 1:2 w p pt 19 lc rgb "red", "$out1D200_k" u 1:2 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "u1"
set size 0.5,0.3
set origin 0.5,0.66
plot "$out1D200" u 1:3 w p pt 19 lc rgb "red", "$out1D200_k" u 1:3 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "u2"
set origin 0.00,0.33
plot "$out1D200" u 1:4 w p pt 19 lc rgb "red", "$out1D200_k" u 1:4 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "sigma11"
set origin 0.5,0.33
plot "$out1D200" u 1:5 w p pt 19 lc rgb "red", "$out1D200_k" u 1:5 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "sigma12"
set origin 0.00,0.05
plot "$out1D200" u 1:6 w p pt 19 lc rgb "red", "$out1D200_k" u 1:6 w p pt 19 lc rgb "black"

set size 0.5,0.30
set ylabel "sigma13"
set origin 0.5,0.05
plot "$out1D200" u 1:7 w p pt 19 lc rgb "red", "$out1D200_k" u 1:7 w p pt 19 lc rgb "black"

set nomultiplot
quit 
PLOT

evince $output
 
#for comparison plot Kevin's results for same test
#param 1 : compile % run
# ./KevinElastic/scriptTests.sh 0 


#----------------------------------------------------------------------------------------


# # gnuplot -persist << PLOT
# set multiplot layout 2,2
# set nokey
# set yrange[0:1]
# set xrange[0:1]
# #set autoscale
# #set size square
# set pointsize 3

# set xlabel "position" font "Helvetica,5"
# set ylabel "density" font "Helvetica,11"
# plot 'test_01_0.dat' using 1:2 with points pointtype 7 pointsize 0.5

# set ylabel "velocity"
# #set yrange[-0.1:1.6]
# #set xrange[0:1]
# plot 'test_01_0.dat' using 1:3 with points pointtype 7 pointsize 0.5

# set ylabel "pressure"
# #set yrange[0:1]
# #Set xrange[0:1]
# plot 'test_01_0.dat' using 1:4 with points pointtype 7 pointsize 0.5

# set ylabel "internal energy" 
# #set yrange[1.8:3.8]
# #set xrange[0:1]
# plot 'test_01_0.dat' using 1:5 with points pointtype 7 pointsize 0.5

# quit

# PLOT


