
gnuplot -persist << PLOT
set xrange [ 0.0 : 0.01 ]
set nokey
set lmargin 10
set rmargin 2
set pointsize 0.5

set multiplot
set ylabel "density"
set xlabel "position"
set size 0.5,0.30
set origin 0.0,0.66
plot "$out1D200" u 1:2 w p pt 19

set size 0.5,0.30
set ylabel "u1"
set size 0.5,0.3
set origin 0.5,0.66
plot "$out1D200" u 1:3 w p pt 19

set size 0.5,0.30
set ylabel "u2"
set origin 0.00,0.33
plot "$out1D200" u 1:4 w p pt 19

set size 0.5,0.30
set ylabel "sigma11"
set origin 0.5,0.33
plot "$out1D200" u 1:5 w p pt 19

set size 0.5,0.30
set ylabel "sigma12"
set origin 0.00,0.05
plot "$out1D200" u 1:6 w p pt 19

set size 0.5,0.30
set ylabel "sigma13"
set origin 0.5,0.05
plot "$out1D200" u 1:7 w p pt 19


PLOT
