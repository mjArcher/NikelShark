#! /bin/bash
base="~/Dropbox/Project/Code/"
LOutDir="~/Dropbox/Project/ReportPlots/"


plotD_R=$LOutDir"data/displacementProfileRalf.dat" #ralf displacement profile
plotD_GN=$LOutDir"data/displacementProfile.dat" #displacement profile from Giordano
plotD_GE=$LOutDir"data/displacementProfileExp.dat" #experimental displacement profile
plotD_GP=$LOutDir"data/pressureProfile.dat" #giordano experimental pressure profile
plotD_MP="/local/data/public/ma595/lsc_amr/LSC_AMR_MultiMaterial/output/New/elasticBeamPlatform/sensorWOAMR.dat"


plotD_M="data/BeambeamAgain_profile.dat"
plotD_GRCP=$LOutDir"data/beamPlatformWOAMRGPressure.dat" #without amr pressure
plotD_GRCD=$LOutDir"data/beamPlatformWOAMRGDefl.dat" #without amr group deflection+velocity
plotD_GRCPAMR=$LOutDir"data/beamPlatform3MHGPressure.dat" #with amr group pressure
plotD_GRCDAMR=$LOutDir"data/beamPlatform3MHGDefl.dat" #with amr group deflection+velocity


outName0a="myDeflection"
outName0b="myVelocity"
outName1="compareDefl"
outName2="myplotCompare"
outName3="pressureBeamLoad"
outName4="pressureBeamLoadAMR"
outName5="pressurePt"

linewidth=2
pointsize=1
pointtype1=19
pointtype2=1
pointtype3=2
time=0.005
gnuplot -persist << PLOT
load '~/Dropbox/Project/Code/ver5/Beam/latexplot.plt'
set key at screen 0.14,0.95


set ylabel ' $\textbf{deflection}$ [m]' offset +5.0
set style line 3 lt 1 lw 5 pt 3
set style line 4 lt 3 lw 5 pt 3

set terminal epslatex size 16.5cm,18cm color colortext 9 standalone
set output '$outName0a.tex'
set size 1,1
set multiplot
set xtics format " " 
set xtics 0.001
set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.55
set tmargin at screen 0.95
plot "$plotD_M" u 1:2 w l ls 3 lc rgb "red" title 'Deflection (SLIC)',\
"$plotD_GN" u 1:2 w l ls 3 lc rgb "black" title 'Giordano (CARBUR)',\
"$plotD_R" u 1:2 w l ls 4 lc rgb "black" title 'Deiterding (MUSCL)'
set xtics format '\color{black}$%g$ '
set xlabel ' $\textbf{time}$ [s] '
set ylabel ' $\textbf{velocity}$ [m/s]' offset +5.0
set xrange [ 0.0 : $time ]
set lmargin at screen 0.15
set rmargin at screen 0.90
set bmargin at screen 0.10
set tmargin at screen 0.5
set key at screen 0.13,0.5

#set output '$outName0b.tex'
plot "$plotD_M" u 1:3 w l ls 3 lc rgb "red" title 'Velocity profile'

unset multiplot


load '~/Dropbox/Project/Code/ver5/Beam/latexplot.plt'
set key at screen 0.15,0.7
set xtics format '\color{black}$%g$ '
set xlabel ' $\textbf{time}$ [s] '
set ylabel ' $\textbf{deflection}$ [m]' offset +5.0
set style line 3 lt 1 lw 3 pt 3
set style line 4 lt 3 lw 3 pt 3

set output '$outName1.tex'
plot "$plotD_R" u 1:2 w l ls 3 lc rgb "black" title 'Deiterding (MUSCL)',\
 "$plotD_GN" u 1:2 w l ls 4 lc rgb "black" title 'Giordano (CARBUR) ',\
 "$plotD_GE" u 1:2 w p pt 5 ps $pointsize lw 0.2 lc rgb "black" title 'Giordano (Experiment)' 

set output '$outName2.tex'
plot "$plotD_M" u 1:2 w l ls 3 lc rgb "red" title 'Deflection SLIC',\
"$plotD_GRCD" u 1:2 w l ls 3 lc rgb "green" title 'WAF',\
"$plotD_GRCDAMR" u 1:2 w l ls 3 lc rgb "blue" title 'MUSCL + AMR, r $\times$ 3',\
"$plotD_GN" u 1:2 w l ls 3 lc rgb "black" title 'Giordano (CARBUR)',\
"$plotD_R" u 1:2 w l ls 4 lc rgb "black" title 'Deiterding (MUSCL + AMR)'

set yrange [ -100000 : 200000 ]
set output '$outName3.tex'
set key at screen 0.18,0.7
set ylabel ' $\textbf{Pressure load}$ [Pa] ' offset +5.0
plot for [ii=1:1] "$plotD_GRCP" u 1:100 w l ls 3 lc rgb "red" title '(100)',\
for [ii=1:1] "$plotD_GRCP" u 1:50 w l ls 3 lc rgb "green" title '(50)',\
for [ii=1:1] "$plotD_GRCP" u 1:2 w l ls 3 lc rgb "blue" title '(1)'

#plot "$plotD_GRCP u 1:2 w l ls 3 lc rgb "black" notitle

set output '$outName4.tex'
plot for [ii=1:99] "$plotD_GRCPAMR" u 1:ii w l ls 3 lc rgb "black" notitle

unset yrange
set yrange [ 100000 : 280000 ]
set key at screen 0.2,0.95
set output '$outName5.tex'
set xrange [ 0.0 : 0.0032 ]
set ylabel ' $\textbf{Pressure}$ [Pa] ' offset +5.0
plot "$plotD_MP" u 1:2 index 3 w l ls 3 lc rgb "red" title 'Pressure (SLIC)',\
"$plotD_GP" u 1:2 w l ls 3 lc "black" title 'Giordano (CARBUR)'


PLOT

#ps2pdf -dEPSCrop  $outName1-inc.eps $outName1-inc.pdf

latex $outName0a.tex
dvipdfm -o $outName0a.pdf $outName0a.dvi

#latex $outName0b.tex
#dvipdfm -o $outName0b.pdf $outName0b.dvi

latex $outName1.tex
dvipdfm -o $outName1.pdf $outName1.dvi

latex $outName2.tex
dvipdfm -o $outName2.pdf $outName2.dvi

latex $outName3.tex
dvipdfm -o $outName3.pdf $outName3.dvi

latex $outName4.tex
dvipdfm -o $outName4.pdf $outName4.dvi

latex $outName5.tex
dvipdfm -o $outName5.pdf $outName5.dvi


rm  ./*.aux
rm  ./*.eps 
rm  ./*.log 
rm  ./*.dvi 
rm  ./*.tex 
