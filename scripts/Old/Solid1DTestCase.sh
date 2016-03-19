#! /bin/bash
#Format is ./Scripts/SlicSolve.sh 1(run) 200(test) 2(test Case)
res=$3 #resolution
testCase=$4 #test case
kevDir="./../../KevinElastic"
output="./Barton1D"$res"_"$testCase".pdf"

#Latex output
output1="./Barton1D_"$testCase"Latex_1"
output2="./Barton1D_"$testCase"Latex_2"
output3="./Barton1D_"$testCase"Latex_3"
#initial
out1D_ini="./Output/out1D"$res"_"$testCase"ini.dat"
out1D_k_ini=$kevDir"/out/Barton"$testCase"_"$res"/out000.tec"
#final
out1D="./Output/out1D"$res"_"$testCase".dat"
out1D_k=$kevDir"/out/Barton"$testCase"_"$res"/out001.tec"

if [ $1 -eq 1 ];
then 
# g++ -O3 -Wall -I ~/Libraries/eigen/ -I ~/Libraries/ UtilityFunctions.cpp SlicSolve.cpp -o ss
make Elastic1D
./Elastic1D barton$testCase $res $out1D_ini $out1D

fi 

# Compute Kevin's results for comparison
if [ $2 -eq 1 ];
then 
cd ./../../KevinElastic
rm -rf ./out/Barton$testCase_$res/*
make ./Elastic1DUnigrid
./Elastic1DUnigrid barton$testCase $res
cd ../myCode/SimpleElastic/
fi
#gnuplot output

echo "looking in "$out1D_ini" for initial output"
echo "looking in "$out1D" for final output"
echo "looking in "$out1D_k_ini" for initial exact output"
echo "looking in "$out1D_k" for final exact output"

lis=1
ltp="p"
lisE=2
ltpE="l"

gnuplot -persist << PLOT
set term postscript portrait enhanced colour size 21cm,29cm
set output '| ps2pdf - $output'

set xrange [ 0.0 : 0.01 ]
load './Scripts/profile.plt'
set multiplot layout 3,2
set size 0.4, 0.3
# yLabel = "density u1 u2 sigma11 sigma12 sigma13"
set ylabel "{/Symbol r} (g/cm^3) " offset 2 
plot "$out1D_ini" u 1:2 w $ltp ls $lis, "$out1D_k_ini" u 1:2 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "u_1 (km/s)" offset 2 
plot "$out1D_ini" u 1:3 w $ltp ls $lis, "$out1D_k_ini" u 1:3 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "u_2 (km/s)" 
plot "$out1D_ini" u 1:4 w $ltp ls $lis, "$out1D_k_ini" u 1:4 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "u_3 (km/s)" offset 2 
plot "$out1D_ini" u 1:5 w $ltp ls $lis, "$out1D_k_ini" u 1:5 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel " {/Symbol s_{11}} (GPa)" offset 2 
plot "$out1D_ini" u 1:6 w $ltp ls $lis, "$out1D_k_ini" u 1:6 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel " {/Symbol s_{12}} (GPa) " offset 2 
plot "$out1D_ini" u 1:7 w $ltp ls $lis, "$out1D_k_ini" u 1:7 w $ltpE ls $lisE


set nomultiplot
set multiplot
set format y "%.2g"
set size 0.4, 0.3
set origin 0.0, 0.65
set pointsize 0.3
set ylabel  "{/Symbol s_{13}} (GPa)" offset 2 
plot "$out1D_ini" u 1:8 w $ltp ls $lis, "$out1D_k_ini" u 1:8 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.5, 0.65
set ylabel "{/Symbol s_{22}} (GPa)" offset 2 
plot "$out1D_ini" u 1:9 w $ltp ls $lis, "$out1D_k_ini" u 1:9 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.0, 0.35
set ylabel "{/Symbol s_{23}} (GPa)" offset 2 
plot "$out1D_ini" u 1:10 w $ltp ls $lis, "$out1D_k_ini" u 1:10 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.5, 0.35
set ylabel "{/Symbol s_{33}} (GPa)" offset 2 
plot "$out1D_ini" u 1:11 w $ltp ls $lis, "$out1D_k_ini" u 1:11 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.25, 0.05
set ylabel "S (kJ/gK)"  
plot "$out1D_ini" u 1:12 w $ltp ls $lis, "$out1D_k_ini" u 1:12 w $ltpE ls $lisE

set nomultiplot 
set size 0.4,0.30
set multiplot layout 3,2

set ylabel "{/Symbol r} (g/cm^3)"
set size 0.4, 0.3
plot "$out1D" u 1:2 w $ltp ls $lis, "$out1D_k" u 1:2 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "u_1 (km/s)"
plot "$out1D" u 1:3 w $ltp ls $lis, "$out1D_k" u 1:3 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "u_2 (km/s)"
plot "$out1D" u 1:4 w $ltp ls $lis, "$out1D_k" u 1:4 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "u_3 (km/s)"
plot "$out1D" u 1:5 w $ltp ls $lis, "$out1D_k" u 1:5 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "{/Symbol s_{11}} (GPa)" 
plot "$out1D" u 1:6 w $ltp ls $lis, "$out1D_k" u 1:6 w $ltpE ls $lisE

set size 0.4, 0.3
set ylabel "{/Symbol s_{12}} (GPa)"
plot "$out1D" u 1:7 w $ltp ls $lis, "$out1D_k" u 1:7 w $ltpE ls $lisE

set nomultiplot
# set size 1, 1
set multiplot
# set multiplot layout 2,2
set format y "%.2g"
set size 0.4, 0.3
set origin 0.0, 0.65
set pointsize 0.3
yLabel = "density u1 u2 sigma11 sigma12 sigma13"
set ylabel  "{/Symbol s_{13}} (GPa)" offset 2 
plot "$out1D" u 1:8 w $ltp ls $lis, "$out1D_k" u 1:8 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.5, 0.65
set ylabel "{/Symbol s_{22}} (GPa)" offset 2 
plot "$out1D" u 1:9 w $ltp ls $lis, "$out1D_k" u 1:9 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.0, 0.35
set ylabel "{/Symbol s_{23}} (GPa)" offset 2 
plot "$out1D" u 1:10 w $ltp ls $lis, "$out1D_k" u 1:10 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.5, 0.35
set ylabel "{/Symbol s_{33}} (GPa)" offset 2 
plot "$out1D" u 1:11 w $ltp ls $lis, "$out1D_k" u 1:11 w $ltpE ls $lisE

set size 0.4, 0.3
set origin 0.25, 0.05
set ylabel "S (kJ/gK)"  
plot "$out1D" u 1:12 w $ltp ls $lis, "$out1D_k" u 1:12 w $ltpE ls $lisE

set nomultiplot
quit 
PLOT

ox1=0.075
ox2=0.55
oy1=0.65
oy2=0.35
oy3=0.05

gnuplot -persist << PLOT
set terminal epslatex size 21cm,29cm color colortext 9 standalone
load './Scripts/profileLatex.plt'
set output '$output1.tex'
set xrange [ 0.0 : 0.01 ]

set multiplot
set size 0.4, 0.3
set origin $ox1, $oy1 
yLabel = "density u1 u2 sigma11 sigma12 sigma13"
set ylabel ' $ \rho$ (g/cm$^3$) density' offset 1
plot "$out1D_ini" u 1:2 w $ltp ls $lis, "$out1D_k_ini" u 1:2 w $ltp ls $lis

set size 0.4, 0.3
set origin $ox2, $oy1 
set ylabel ' $ u_1 $ ' offset 1
plot "$out1D_ini" u 1:3 w $ltp ls $lis, "$out1D_k_ini" u 1:3 w $ltp ls $lis

set size 0.4, 0.3
set origin $ox1, $oy2
set ylabel ' $ u_2$ ' offset 1
plot "$out1D_ini" u 1:4 w $ltp ls $lis, "$out1D_k_ini" u 1:4 w $ltp ls $lis

set size 0.4, 0.3
set origin $ox2, $oy2
set ylabel ' $ u_3 $' offset 1
plot "$out1D_ini" u 1:5 w $ltp ls $lis, "$out1D_k_ini" u 1:5 w $ltp ls $lis

set size 0.4, 0.3
set origin $ox1, $oy3
set ylabel ' $\sigma_{11} $ ' offset 1
plot "$out1D_ini" u 1:6 w $ltp ls $lis, "$out1D_k_ini" u 1:6 w $ltp ls $lis

set size 0.4, 0.3
set origin $ox2, $oy3
set ylabel ' $\sigma_{12}$ ' offset 1
plot "$out1D_ini" u 1:7 w $ltp ls $lis, "$out1D_k_ini" u 1:7 w $ltp ls $lis

set nomultiplot
quit 
PLOT

latex $output1.tex
dvipdfm -o $output1.pdf $output1.dvi

rm  ./*.aux
rm  ./*.eps 
rm  ./*.log 
rm  ./*.dvi 
rm  ./*.tex 

evince $output
# evince $output1.pdf
# gnuplot -persist << PLOT
# set terminal epslatex size 21cm,29cm color colortext 9 standalone
# # load './Scripts/latexplot.plt'
# set output '$output1.tex'

# set style line 1 ps 0.5 pt 19 lc rgb "black"
# set style line 2 lc rgb '#5e9c36' pt 19 ps 0.5 lt 1 lw 0.2 # --- green

# set xrange [ 0.0 : 0.01 ]
# set nokey

# set lmargin 1
# set rmargin 1
# set border linewidth 1.5
# set xtics nomirror font "arial, 10"
# set ytics nomirror font "arial, 10"

# set style line 1 lt 2 lw 2 ps 0.5 pt 19 lc rgb "red"
# set style line 12 lc rgb'#808080' lt 0 lw 1
# set grid back ls 12

# # set xlabel "position" font "Times-Roman, 12"

# set multiplot layout 3,2
# set format y "%.3g"
# set size 0.4, 0.3
# set pointsize 0.3
# yLabel = "density u1 u2 sigma11 sigma12 sigma13"
# set ylabel ' $ \rho$ (g/cm$^3$) ' offset 2 
# plot "$out1D_ini" u 1:2 w p pt 19 lc rgb "red", "$out1D_k_ini" u 1:2 w $ltp ls $lis

# set size 0.4, 0.3
# set ylabel ' $ u_1 $ ' offset 2 
# plot "$out1D_ini" u 1:3 w p pt 19 lc rgb "red", "$out1D_k_ini" u 1:3 w $ltp ls $lis

# set size 0.4, 0.3
# set ylabel ' $ u_2$ ' offset 2 
# plot "$out1D_ini" u 1:4 w p pt 19 lc rgb "red", "$out1D_k_ini" u 1:4 w $ltp ls $lis

# set size 0.4, 0.3
# set ylabel ' $ u_3 $' offset 2 
# plot "$out1D_ini" u 1:5 w p pt 19 lc rgb "red", "$out1D_k_ini" u 1:5 w $ltp ls $lis

# set size 0.4, 0.3
# set ylabel ' $\sigma_{11} $ ' offset 2 
# plot "$out1D_ini" u 1:6 w p pt 19 lc rgb "red", "$out1D_k_ini" u 1:6 w $ltp ls $lis

# set size 0.4, 0.3
# set ylabel ' $\sigma_{12}$ ' offset 2 
# plot "$out1D_ini" u 1:7 w p pt 19 lc rgb "red", "$out1D_k_ini" u 1:7 w $ltp ls $lis


# # set size 0.4,0.30
# # set multiplot layout 3,2

# # set ylabel "density"
# # set size 0.4, 0.3
# # plot "$out1D" u 1:2 w p pt 19 lc rgb "red", "$out1D_k" u 1:2 w $ltp ls $lis

# # set size 0.4, 0.3
# # set ylabel "u1"
# # plot "$out1D" u 1:3 w p pt 19 lc rgb "red", "$out1D_k" u 1:3 w $ltp ls $lis

# # set size 0.4, 0.3
# # set ylabel "u2"
# # plot "$out1D" u 1:4 w p pt 19 lc rgb "red", "$out1D_k" u 1:4 w $ltp ls $lis

# # set size 0.4, 0.3
# # set ylabel "u3"
# # plot "$out1D" u 1:5 w p pt 19 lc rgb "red", "$out1D_k" u 1:5 w $ltp ls $lis

# # set size 0.4, 0.3
# # set ylabel "sigma11"
# # plot "$out1D" u 1:6 w p pt 19 lc rgb "red", "$out1D_k" u 1:6 w $ltp ls $lis

# # set size 0.4, 0.3
# # set ylabel "sigma12"
# # plot "$out1D" u 1:7 w p pt 19 lc rgb "red", "$out1D_k" u 1:7 w $ltp ls $lis

# set nomultiplot
# set size 1, 1
# set multiplot
# # set lmargin 5

# set format y "%.2g"
# set size 0.4, 0.3
# set pointsize 0.3
# set border 1+2 back 


# set origin $ox1, $oy1 
# yLabel = "density u1 u2 sigma11 sigma12 sigma13"
# set ylabel ' $\sigma_{13}$ (GPa) ' offset 1 
# plot "$out1D" u 1:8 w p pt 19 lc rgb "red", "$out1D_k" u 1:8 w $ltp ls $lis

# set size 0.4, 0.3
# set origin $ox2, $oy1 
# set ylabel ' $\sigma_{22}$ (GPa) ' offset 1 
# plot "$out1D" u 1:9 w p pt 19 lc rgb "red", "$out1D_k" u 1:9 w $ltp ls $lis

# set size 0.4, 0.3
# set origin $ox1, $oy2
# set ylabel ' $\sigma_{23}$ (GPa) ' offset 1 
# plot "$out1D" u 1:10 w p pt 19 lc rgb "red", "$out1D_k" u 1:10 w $ltp ls $lis

# set size 0.4, 0.3
# set origin $ox2, $oy2
# set ylabel ' $\sigma_{33}$ (GPa) ' offset 1 
# plot "$out1D" u 1:11 w p pt 19 lc rgb "red", "$out1D_k" u 1:11 w $ltp ls $lis

# set size 0.4, 0.3
# set origin 0.3, 0.05
# set ylabel ' S (kJ/gK) '  
# plot "$out1D" u 1:12 w p pt 19 lc rgb "red", "$out1D_k" u 1:12 w $ltp ls $lis

# set nomultiplot


# quit 
# PLOT

 
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


