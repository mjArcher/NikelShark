#! /usr/bin/env gnuplot 
set term pdfcairo enhanced font 'Helvetica, 12' size 21cm, 29cm dashed
set style line 11 lc rgb '#808080' lt 1 
set style line 12 lc rgb'#808080' lt 0 lw 1 # dotted line 

load '~/cns_work/work/scripts/output/1D/Profile.plt' 
set style line 1 lw 5 pt 19 ps 0.2 lt 1
set style line 2 lw 5 pt 1 ps 0.2 linetype 2
set style data points

set output "`echo $pdfname`"
set tics nomirror  
set rmarg 5  
# plot "`echo $elasticout`" u 1:2 w l ls 11

# plot "`echo $elasticout`" index 10 u 1:2 w l ls 1 

lis=1
lis2=2
i=10 # index 

set multiplot layout 3,2
set lmargin 10
# set key center top 
# unset key 
set ylabel "{/Symbol r} (g/cm^3) " offset 0
plot "`echo $elasticout`" index i u 1:2 w p ls lis title "my", "`echo $refout`" u ($1/100):2 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:2 w p ls lis2 lc rgb "blue" title "kev"

set ylabel "u_1 (km/s)" offset 0
plot "`echo $elasticout`" index i u 1:3 w p ls lis title "my", "`echo $refout`" u ($1/100):4 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:3 w p ls lis2 lc rgb "blue" title "kev"

set ylabel "u_2 (km/s)" offset 0
plot "`echo $elasticout`" index i u 1:4 w p ls lis title "my", "`echo $refout`" u ($1/100):5 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:4 w p ls lis2 lc rgb "blue" title "kev"

set ylabel "u_3 (km/s)" offset 0 
plot "`echo $elasticout`" index i u 1:5 w p ls lis title "my", "`echo $refout`" u ($1/100):6 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:5 w p ls lis2 lc rgb "blue" title "kev"

set bmargin 4
set ylabel " {/Symbol s_{11}} (GPa)" offset 0
plot "`echo $elasticout`" index i u 1:6 w p ls lis title "my", "`echo $refout`" u ($1/100):7 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:6 w p ls lis2 lc rgb "blue" title "kev" 

set ylabel " {/Symbol s_{12}} (GPa) " offset 0
plot "`echo $elasticout`" index i u 1:7 w p ls lis title "my", "`echo $refout`" u ($1/100):8 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:7 w p ls lis2 lc rgb "blue" title "kev" 


set multiplot layout 3,2
# unset key

set ylabel " {/Symbol s_{13}} (GPa) " offset 0
plot "`echo $elasticout`" index i  u 1:8 w p ls lis title "my", "`echo $refout`" u ($1/100):9 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:8 w p ls lis2 lc rgb "blue" title "kev" 

set ylabel " {/Symbol s_{22}} (GPa) " offset 0
plot "`echo $elasticout`" index i u 1:10 w p ls lis title "my", "`echo $refout`" u ($1/100):11 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:10 w p ls lis2 lc rgb "blue" title "kev" 

set ylabel " {/Symbol s_{23}} (GPa) " offset 0
plot "`echo $elasticout`" index i u 1:11 w p ls lis title "my", "`echo $refout`" u ($1/100):12 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:11 w p ls lis2 lc rgb "blue" title "kev" 

set ylabel " {/Symbol s_{33}} (GPa) " offset 0
plot "`echo $elasticout`" index i u 1:14 w p ls lis title "my", "`echo $refout`" u ($1/100):15 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:14 w p ls lis2 lc rgb "blue" title "kev" 

set ylabel "S" 
plot "`echo $elasticout`" index i u 1:15 w p ls lis title "my", "`echo $refout`" u ($1/100):16 w p ls lis2 lc rgb "green" title "cns", "`echo $refkevout`" u 1:15 w p ls lis2 lc rgb "blue" title "kev" 

