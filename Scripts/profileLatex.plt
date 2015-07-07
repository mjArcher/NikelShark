#!/usr/bin/gnuplot

reset
set style line 1 lt 2 lw 2 ps 0.3 pt 19 lc rgb "red" 
set style line 2 ps 0.5 pt 19 lc rgb "black"
set style line 3 lc rgb '#5e9c36' pt 19 ps 0.5 lt 1 lw 0.2 # --- green
set style line 12 lc rgb'#808080' lt 0 lw 1

set nokey
set lmargin 1
set rmargin 1
set border linewidth 1.5 1+2 back
# set ticscale 0.5 0.25
set xtics 0.002 nomirror 
set mxtics 0.001
set ytics nomirror
set ylabel font "Arial, 12"
#set ylabel font "helvetica, 10" 
set grid back ls 12


# set style line 1 lt 2 lw 2 pt 3 ps 0.5 pt 19 lc rgb "black"
# set style line 3 lt 1 lw 5 lc rgb "red"
# set style line 4 lt 1 lc rgb "black"
# set style line 3 
# set style line 4 lt 1 lc rgb "black"

# set xlabel "position" font "Times-Roman, 12"

set format y "%.3g"
