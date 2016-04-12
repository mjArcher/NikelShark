#! /usr/bin/env gnuplot 
set term pdfcairo enhanced font 'Helvetica, 12' size 21cm, 29cm dashed
set style line 11 lc rgb '#808080' lt 1 
set style line 12 lc rgb'#808080' lt 0 lw 1 # dotted line 

load '~/cns_work/work/scripts/output/1D/Profile.plt' 
set style line 1 lw 5 pt 19 ps 0.2 lt 1
set style line 2 lw 5 pt 1 ps 0.2 linetype 2

set output "~/solid-1D/"."`echo $pdfout`".".pdf"

set tics nomirror  
set rmarg 5  

lis=1
i=1 # index 

set multiplot layout 3,2
set lmargin 10
# set key center top 
unset key 
set ylabel "{/Symbol r} (g/cm^3) " offset 0
plot "`echo $elasticout`" index 10 u 1:2 w p ls lis 

set ylabel "u_1 (km/s)" offset 0
plot "`echo $elasticout`" index 10 u 1:3 w p ls lis

set ylabel "u_2 (km/s)" offset 0
plot "`echo $elasticout`" index 10 u 1:4 w p ls lis

set ylabel "u_3 (km/s)" offset 0 
plot "`echo $elasticout`" index 10 u 1:5 w p ls lis

set bmargin 4
set ylabel " {/Symbol s_{11}} (GPa)" offset 0
plot "`echo $elasticout`"index 10  u 1:6 w p ls lis

set ylabel " {/Symbol s_{12}} (GPa) " offset 0
plot "`echo $elasticout`" index 10 u 1:7 w p ls lis


set multiplot layout 3,2
unset key

set ylabel " {/Symbol s_{13}} (GPa) " offset 0
plot "`echo $elasticout`" index 10 u 1:8 w p ls lis

set ylabel " {/Symbol s_{22}} (GPa) " offset 0
plot "`echo $elasticout`" index 10 u 1:10 w p ls lis

set ylabel " {/Symbol s_{23}} (GPa) " offset 0
plot "`echo $elasticout`" index 10 u 1:11 w p ls lis

set ylabel " {/Symbol s_{33}} (GPa) " offset 0
plot "`echo $elasticout`" index 10 u 1:14 w p ls lis

set ylabel "S" 
plot "`echo $elasticout`" index 10 u 1:15 w p ls lis

