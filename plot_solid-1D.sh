#! /bin/bash
export elasticout="/local/data/public/ma595/output/solid-1D/barton1D"
export refout="~/cns/output/entropy_1D/RiemannA_1.dat"
export pdfname="solid-1D.pdf"
gnuplot-4.6 ./scripts/gnuplot/plot-riemann_multi.gnu


