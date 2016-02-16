reset
# epslatex
set terminal epslatex size 16.5cm,12cm color colortext 9 standalone  \
header '\definecolor{tics}{rgb}{0.5,0.5,0.5}'

#21cm, 14.8

# define axis
# remove border on top and right and set color to gray
set pointsize 1
#set size 1,1
#set origin 0.0,0.0
#set bmargin 1
#set tmargin 0
set size 1, 0.75

set style line 11 lc rgb '#808080' lt 1
set border 3 back ls 11
set tics nomirror
# define grid
set style line 12 lc rgb'#808080' lt 0 lw 1
set grid back ls 12

# color definitions
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1.5 lt 1 lw 5 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1.5 lt 1 lw 5 # --- green

#set key bottom right
set key top left
set format '\color{black}$%g$'
