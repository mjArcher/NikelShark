#! /bin/bash

#make all unit tests 
cd ~/solid-1D/src
make SystemUnitTests

cd ~/solid-1D/kevinref
make KevinTests 

cd ~/solid-1D/

echo "Kevin's tests" 
./kevinref/KevinTests
echo "Tests" 
./src/SystemUnitTests


