#! /bin/bash

# ensure that openmp is disabled 
cd ~/solid-1D
./make-clean.sh
./make.sh
HOME="/home/raid/ma595"

/lsc/opt/bin/valgrind-3.10 --tool=callgrind --log-file="$HOME/solid-1D/profiling/logdat/mylog.out" \
  --callgrind-out-file="$HOME/solid-1D/profiling/logdat/mycall.out" \
  ~/solid-1D/src/Elastic1D ~/solid-1D/settings/barton1D.cfg 

cd ~/solid-1D/kevinref/
make clean
make Elastic1DUnigrid
/lsc/opt/bin/valgrind-3.10 --tool=callgrind --log-file="$HOME/solid-1D/profiling/logdat/kevlog.out" \
  --callgrind-out-file="$HOME/solid-1D/profiling/logdat/kevcall.out" \
  ~/solid-1D/kevinref/Elastic1DUnigrid barton1 100
  

#run code
