#! /usr/bin/env gnuplot 
set term pdfcairo noenhanced font 'Helvetica, 12' size 25cm, 12cm dashed
set style line 11 lc rgb '#808080' lt 1 
set style line 12 lc rgb'#808080' lt 0 lw 1 # dotted line 

set style line 1 lw 5 pt 19 ps 0.8 lt 1
set style line 2 lw 5 pt 1 ps 0.2 linetype 2

set output "~/solid-1D/out.pdf"
dataOut='/local/data/public/ma595/output/solid-1D/barton1D'
plot dataOut u 1:2 w l ls 11

plot dataOut index 1 u 1:2 w l ls 1 
