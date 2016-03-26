#! /usr/bin/env gnuplot 
set term pdfcairo enhanced font 'Helvetica, 12' size 21cm, 10cm dashed
set style line 11 lc rgb '#808080' lt 1 
set style line 12 lc rgb'#808080' lt 0 lw 1 # dotted line 
set style line 1 lw 5 pt 19 ps 0.2 lt 1
set style line 2 lw 5 pt 1 ps 0.2 linetype 2
set style data points
set tics nomirror  
set grid
set rmarg 5  
set output "gnuplottest.pdf"
set samples 1000
set xrange [-2*pi:2*pi]
set border linewidth 1.5
set style function lines
set style data lines 
plot cos(x), sin(x) lc rgb "blue"
